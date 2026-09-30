#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <algorithm>
#include <unordered_map>
#include <chrono>

struct Clause {
    std::vector<int> literals;
};

class ADLSATSolverCorrected {
private:
    int num_vars;
    int num_clauses;
    std::vector<Clause> clauses;
    std::vector<double> beliefs;
    std::unordered_map<int, bool> fixed_vars;

    // Adjacency structure: var -> list of (clause_idx, literal_pos)
    std::vector<std::vector<std::pair<int, size_t>>> var_to_clauses;
    std::vector<std::vector<double>> eta_msg;
    std::vector<std::unordered_map<int, double>> mu_msg;
    std::mt19937 rng;

    int compute_2_adic_valuation(double delta_mu, int k_precision = 16) {
        if (std::abs(delta_mu) < 1e-12) return 20;
        int quantized_m = static_cast<int>(std::round(delta_mu * (1 << k_precision)));
        if (quantized_m == 0) return 20;

        int v2 = 0, abs_m = std::abs(quantized_m);
        while ((abs_m & 1) == 0 && abs_m != 0) {
            v2++;
            abs_m >>= 1;
        }
        return v2 - k_precision;
    }

    double apply_2_adic_forcing(double mu_curr, double mu_prev, double base_gamma = 0.5) {
        double delta_mu = mu_curr - mu_prev;
        int v2 = compute_2_adic_valuation(delta_mu);
        double forcing_power = 1.0 + (base_gamma * std::pow(2.0, v2 / 5.0));

        double p1 = std::pow(std::max(mu_curr, 1e-10), forcing_power);
        double p0 = std::pow(std::max(1.0 - mu_curr, 1e-10), forcing_power);
        return p1 / (p1 + p0 + 1e-15);
    }

public:
    ADLSATSolverCorrected(int vars, const std::vector<Clause>& cl)
        : num_vars(vars), num_clauses(cl.size()), clauses(cl),
          beliefs(vars + 1, 0.5), rng(42) {
        eta_msg.resize(num_clauses);
        var_to_clauses.resize(num_vars + 1);

        for (size_t c = 0; c < clauses.size(); ++c) {
            eta_msg[c].resize(clauses[c].literals.size(), 0.5);
            for (size_t l = 0; l < clauses[c].literals.size(); ++l) {
                int var = std::abs(clauses[c].literals[l]);
                var_to_clauses[var].push_back({static_cast<int>(c), l});
            }
        }

        mu_msg.resize(num_vars + 1);
        for (size_t c = 0; c < clauses.size(); ++c)
            for (int lit : clauses[c].literals)
                mu_msg[std::abs(lit)][c] = 0.5;
    }

    std::vector<Clause> phase1_decimation(int max_iters = 100, double threshold = 0.90) {
        std::vector<double> mu_prev = beliefs;

        for (int iter = 0; iter < max_iters; ++iter) {
            mu_prev = beliefs;

            for (size_t a = 0; a < clauses.size(); ++a) {
                for (size_t j_idx = 0; j_idx < clauses[a].literals.size(); ++j_idx) {
                    double prod_term = 1.0;
                    for (size_t k_idx = 0; k_idx < clauses[a].literals.size(); ++k_idx) {
                        if (j_idx == k_idx) continue;
                        int lit_k = clauses[a].literals[k_idx];
                        int var_k = std::abs(lit_k);
                        double mu_val = mu_msg[var_k][a];
                        double term = (lit_k > 0) ? mu_val : (1.0 - mu_val);
                        prod_term *= (1.0 - term);
                    }
                    eta_msg[a][j_idx] = 1.0 - prod_term;
                }
            }

            for (int i = 1; i <= num_vars; ++i) {
                if (fixed_vars.count(i)) continue;

                double total_prod_pos = 1.0, total_prod_neg = 1.0;
                for (const auto& entry : var_to_clauses[i]) {
                    int c = entry.first;
                    size_t l_idx = entry.second;
                    double eta = eta_msg[c][l_idx];
                    total_prod_pos *= eta;
                    total_prod_neg *= (1.0 - eta);
                }

                double raw_mu = total_prod_pos / (total_prod_pos + total_prod_neg + 1e-15);
                beliefs[i] = apply_2_adic_forcing(raw_mu, mu_prev[i]);

                for (auto& pair : mu_msg[i]) pair.second = beliefs[i];

                if (beliefs[i] > threshold) fixed_vars[i] = true;
                else if (beliefs[i] < (1.0 - threshold)) fixed_vars[i] = false;
            }

            if (fixed_vars.size() >= static_cast<size_t>(num_vars * 0.85)) break;
        }

        std::vector<Clause> residual;
        for (const auto& clause : clauses) {
            bool satisfied = false;
            Clause reduced;
            for (int lit : clause.literals) {
                int var = std::abs(lit);
                if (fixed_vars.count(var)) {
                    if ((lit > 0 && fixed_vars[var]) || (lit < 0 && !fixed_vars[var])) {
                        satisfied = true;
                        break;
                    }
                } else {
                    reduced.literals.push_back(lit);
                }
            }
            if (!satisfied && !reduced.literals.empty()) residual.push_back(reduced);
        }
        return residual;
    }

    bool phase2_walksat(const std::vector<Clause>& residual, int max_flips = 200000, double p_noise = 0.4) {
        if (residual.empty()) return true;

        std::unordered_map<int, bool> assignment;
        for (const auto& cl : residual)
            for (int lit : cl.literals)
                if (!assignment.count(std::abs(lit)))
                    assignment[std::abs(lit)] = (rng() % 2 == 0);

        auto eval_clause = [&](const Clause& cl) {
            for (int lit : cl.literals) {
                bool val = assignment[std::abs(lit)];
                if ((lit > 0 && val) || (lit < 0 && !val)) return true;
            }
            return false;
        };

        std::uniform_real_distribution<double> dist(0.0, 1.0);
        std::vector<int> unsat_indices;
        unsat_indices.reserve(residual.size());

        for (int flip = 0; flip < max_flips; ++flip) {
            unsat_indices.clear();
            for (size_t i = 0; i < residual.size(); ++i)
                if (!eval_clause(residual[i])) unsat_indices.push_back(i);

            if (unsat_indices.empty()) return true;

            int target_idx = unsat_indices[rng() % unsat_indices.size()];
            const Clause& target_clause = residual[target_idx];

            if (dist(rng) > p_noise) {
                int var_to_flip = std::abs(target_clause.literals[rng() % target_clause.literals.size()]);
                int min_broken = 1e9;
                for (int lit : target_clause.literals) {
                    int v = std::abs(lit);
                    assignment[v] = !assignment[v];

                    int broken = 0;
                    for (const auto& cl : residual)
                        if (!eval_clause(cl)) broken++;
                    assignment[v] = !assignment[v];

                    if (broken < min_broken) {
                        min_broken = broken;
                        var_to_flip = v;
                    }
                }
                assignment[var_to_flip] = !assignment[var_to_flip];
            } else {
                int var_to_flip = std::abs(target_clause.literals[rng() % target_clause.literals.size()]);
                assignment[var_to_flip] = !assignment[var_to_flip];
            }
        }
        return false;
    }
};
