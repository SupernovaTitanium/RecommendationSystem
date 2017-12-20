//
// Created by Xun Zheng on 2/23/17.
//

#ifndef CBGF_BERNOULLIMETROPOLIS_H
#define CBGF_BERNOULLIMETROPOLIS_H


#include "Model.h"

/**
 * Simulated annealing with Metropolis-Gibbs, using Bernoulli proposal.
 */
class BernoulliMetropolis {
public:
    /**
     * Set the model pointer.
     *
     * @param model pointer to Model
     */
    void init_from_model(Model *model) {
        model_ = model;
    }

    /**
     * Run the main algorithm, using provided configuration file.
     *
     * @param config file containing runtime configs
     */
    void run(const Config &config) {
        // Print pretty header
        Metrics metrics;
        metrics.print_header();

        // Evaluate before starting
        metrics.epoch = 0;
        metrics.epoch_time = 0;
        metrics.elapsed_time = 0;
        model_->sparsity_metrics(&metrics);
        model_->prediction_metrics_fast(&metrics);
        metrics.print_full();

        // Start
        for (int epoch = 1; epoch <= config.max_epoch; ++epoch) {
            // Update parameters
            T_ = config.T0 / (1.0 + epoch);
            delta_ = config.delta0 / log(1.0 + epoch);
            accept_count_ = 0;
            // Metropolis Gibbs
            auto t1 = std::chrono::steady_clock::now();
            for (int i = 0; i < model_->graph.num_nodes(); ++i) {
                mh_step(i);
            }
            auto t2 = std::chrono::steady_clock::now();
            // Report some statistics
            metrics.epoch = epoch;
            metrics.epoch_time = std::chrono::duration<double>(t2 - t1).count();
            metrics.elapsed_time += metrics.epoch_time;
            metrics.T = T_;
            metrics.delta = delta_;
            metrics.mh_rate = 1.0 * accept_count_ / model_->graph.num_nodes();
            if (epoch % config.eval_every == 0) {
                model_->sparsity_metrics(&metrics);
                model_->prediction_metrics_fast(&metrics);
                metrics.print_full();
            } else {
                model_->sparsity_metrics(&metrics);
                metrics.print_simple();
            }
        }
    }

    /**
     * Perform Metropolis step for a single node.
     * Complexity: O(nnz/n * nnz/k + k)
     *
     * @param node_id node id
     */
    void mh_step(int node_id) {
        // Propose
        SArray<double> mu = neighborhood_average(node_id);
        SArray<int> proposal = draw_bernoulli(mu);
        // Accept
        double log_ratio = log_bernoulli_ratio(mu, model_->Z.row(node_id), proposal);
        double h_new = model_->single_objective_fast(node_id, proposal);
        double h_old = model_->single_objective_fast(node_id, model_->Z.row(node_id));
        double log_accept_prob = (h_new - h_old) / T_ + log_ratio;
        if (log(rand_.next_double()) < log_accept_prob) {
            model_->Z.update_col_from_row(node_id, model_->Z.row(node_id), proposal);
            model_->Z.update_row(node_id, proposal);
            ++accept_count_;
        }
    }

    /**
     * Average features of the given node's neighbors, in dense form (i.e. length K).
     * Since it will be used in log(mu/(1-mu)), 0 or 1's should be prevented, so add a reducing slack value.
     * Complexity: O(d * nnz/n + k)
     *
     * @param node_id node id
     * @return average of its neighbors
     */
    SArray<double> neighborhood_average(int node_id) {
        SArray<double> mu(model_->Z.col_size());
        for (int j : model_->graph.neighbors(node_id)) {
            for (int k : model_->Z.row(j)) {
                mu[k] += 1.0;
            }
        }
        int d = model_->graph.degree(node_id);
        for (int k = 0; k < mu.size(); ++k) {
            if (mu[k] == 0.0) {
                mu[k] = delta_; // avoid numerical problem
            } else if (mu[k] == d) {
                mu[k] = 1.0 - delta_; // avoid numerical problem
            } else {
                mu[k] /= d;
            }
        }
        return mu;
    }

    /**
     * Draw a binary vector, where k-th entry follows Bernoulli(mu_k).
     * Complexity: O(k)
     *
     * @param mu vector of mean parameter of Bernoulli
     * @return sampled binary vector
     */
    SArray<int> draw_bernoulli(const SArray<double> &mu) {
        SArray<int> sample;
        for (int k = 0; k < mu.size(); ++k) {
            if (rand_.next_double() < mu[k]) {
                sample.push_back(k);
            }
        }
        return sample;
    }

    /**
     * Log of likelihood ratio of two Bernoulli's: log(q(z_old; mu) / q(z_new; mu)) = theta' * (z_old - z_new),
     * where theta = log(mu / (1-mu)) is the natural parameter of Bernoulli.
     * Since mu is dense, z_old and z_new are sparse, use 2-way merge.
     * Complexity: O(nnz/n)
     *
     * @param mu mean parameter of Bernoulli
     * @param z_old old binary vector
     * @param z_new new binary vector
     * @return log of likelihood ratio
     */
    double log_bernoulli_ratio(const SArray<double> &mu, const SArray<int> &z_old, const SArray<int> &z_new) {
        double log_ratio = 0.0;
        auto it_old = z_old.begin();
        auto it_new = z_new.begin();
        while (it_old != z_old.end() && it_new != z_new.end()) {
            if (*it_old < *it_new) {
                // theta_{k_old} * 1
                log_ratio += log(mu[*it_old] / (1.0 - mu[*it_old]));
                ++it_old;
            } else if (*it_old > *it_new) {
                // theta_{k_new} * -1
                log_ratio -= log(mu[*it_new] / (1.0 - mu[*it_new]));
                ++it_new;
            } else {
                // theta_k * 0, do nothing
                ++it_old;
                ++it_new;
            }
        }
        while (it_old != z_old.end()) {
            // theta_{k_old} * 1
            log_ratio += log(mu[*it_old] / (1.0 - mu[*it_old]));
            ++it_old;
        }
        while (it_new != z_new.end()) {
            // theta_{k_new} * -1
            log_ratio -= log(mu[*it_new] / (1.0 - mu[*it_new]));
            ++it_new;
        }
        return log_ratio;
    }

    // Print short info
    void print_info() {
        LOG("Algorithm: BernoulliMetropolis");
    }

private:
    Model *model_; // pointer since algorithm modifies the model
    Random rand_; // random engine
    double T_; // temperature for simulated annealing
    double delta_; // slack for mu, avoid numerical problems when mu = 0 or mu = 1
    int accept_count_; // count MH proposal acceptance
};


#endif //CBGF_BERNOULLIMETROPOLIS_H
