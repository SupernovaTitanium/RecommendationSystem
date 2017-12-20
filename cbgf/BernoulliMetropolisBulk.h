//
// Created by Xun Zheng on 2/24/17.
//

#ifndef CBGF_BERNOULLIMETROPOLISBULK_H
#define CBGF_BERNOULLIMETROPOLISBULK_H


#include "Model.h"

/**
 * Simulated annealing with Metropolis-Gibbs, using Bernoulli proposal, IN PARALLEL.
 * Defer updates to the end of each epoch and perform bulk update together.
 */
class BernoulliMetropolisBulk {
public:
    /**
     * Set the model pointer, init structures for bulk update.
     *
     * @param model pointer to Model
     */
    void init_from_model(Model *model) {
        model_ = model;
        is_accepted_.resize(model_->graph.num_nodes()); // accept = 1, reject = 0
        accepted_proposals_.resize(model_->graph.num_nodes());
    }

    /**
     * Run the main algorithm IN PARALLEL, using provided configuration file.
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
        model_->prediction_metrics_fast_parallel(config.cpu, &metrics);
        metrics.print_full();

        // Start
        for (int epoch = 1; epoch <= config.max_epoch; ++epoch) {
            // Update parameters
            T_ = config.T0 / (1.0 + epoch);
            delta_ = config.delta0 / log(1.0 + epoch);
            // Metropolis Gibbs
            auto t1 = std::chrono::steady_clock::now();
            std::atomic<int> atomic_id(0);
            std::vector<std::thread> thrs(config.cpu);
            for (auto &thr : thrs) {
                thr = std::thread([this, &atomic_id]() {
                    Random local_rand;
                    while (true) {
                        int node_id = atomic_id++;
                        if (node_id >= model_->graph.num_nodes()) {
                            break;
                        }
                        mh_step(node_id, &local_rand);
                    }
                });
            }
            for (auto &thr : thrs) {
                thr.join();
            }
            // Bulk update
            for (int i = 0; i < model_->graph.num_nodes(); ++i) {
                if (is_accepted_[i] == 1) {
                    model_->Z.update_row(i, accepted_proposals_[i]);
                }
            }
            model_->Z.reset_cols_from_rows();  // Important!!
            auto t2 = std::chrono::steady_clock::now();
            // Report some statistics
            metrics.epoch = epoch;
            metrics.epoch_time = std::chrono::duration<double>(t2 - t1).count();
            metrics.elapsed_time += metrics.epoch_time;
            metrics.T = T_;
            metrics.delta = delta_;
            metrics.mh_rate = std::accumulate(is_accepted_.begin(), is_accepted_.end(), 0.0) / is_accepted_.size();
            if (epoch % config.eval_every == 0) {
                model_->sparsity_metrics(&metrics);
                model_->prediction_metrics_fast_parallel(config.cpu, &metrics);
                metrics.print_full();
            } else {
                model_->sparsity_metrics(&metrics);
                metrics.print_simple();
            }
        }
    }

    /**
     * Perform Metropolis step for a single node.
     *
     * @param node_id node id
     */
    void mh_step(int node_id, Random *local_rand) {
        // Propose
        SArray<double> mu = neighborhood_average(node_id);
        SArray<int> proposal = draw_bernoulli(mu, local_rand);
        // Accept
        double log_ratio = log_bernoulli_ratio(mu, model_->Z.row(node_id), proposal);
        double h_new = model_->single_objective_fast(node_id, proposal);
        double h_old = model_->single_objective_fast(node_id, model_->Z.row(node_id));
        double log_accept_prob = (h_new - h_old) / T_ + log_ratio;
        if (log(local_rand->next_double()) < log_accept_prob) {
            is_accepted_[node_id] = 1;
            accepted_proposals_[node_id] = proposal; // zero copy
        } else {
            is_accepted_[node_id] = 0;
            accepted_proposals_[node_id].clear();
        }
    }

    /**
     * Average features of the given node's neighbors, in dense form (i.e. length K).
     * Since it will be used in log(mu/(1-mu)), 0 or 1's should be prevented, so add a reducing slack value.
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
     *
     * @param mu vector of mean parameter of Bernoulli
     * @return sampled binary vector
     */
    SArray<int> draw_bernoulli(const SArray<double> &mu, Random *local_rand) {
        SArray<int> sample;
        for (int k = 0; k < mu.size(); ++k) {
            if (local_rand->next_double() < mu[k]) {
                sample.push_back(k);
            }
        }
        return sample;
    }

    /**
     * Log of likelihood ratio of two Bernoulli's: log(q(z_old; mu) / q(z_new; mu)) = theta' * (z_old - z_new),
     * where theta = log(mu / (1-mu)) is the natural parameter of Bernoulli.
     * Since mu is dense, z_old and z_new are sparse, use 2-way merge.
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
        LOG("Algorithm: BernoulliMetropolisBulk");
    }

private:
    Model *model_; // pointer since algorithm modifies the model
    double T_; // temperature for simulated annealing
    double delta_; // slack for mu, avoid numerical problems when mu = 0 or mu = 1
    SArray<int> is_accepted_; // record acceptance decision, accept = 1, reject = 0
    std::vector<SArray<int>> accepted_proposals_; // store accepted proposals for bulk update
};


#endif //CBGF_BERNOULLIMETROPOLISBULK_H
