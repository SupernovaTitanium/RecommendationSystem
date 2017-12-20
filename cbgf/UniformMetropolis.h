//
// Created by Xun Zheng on 2/25/17.
//

#ifndef CBGF_UNIFORMMETROPOLIS_H
#define CBGF_UNIFORMMETROPOLIS_H


#include "Model.h"

/**
 * Simulated annealing with Metropolis-Gibbs, using Uniform proposal.
 */
class UniformMetropolis {
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
     * Uniform proposal typically leads to dense Z, so use non-fast versions of
     * model functions that does not use column view.
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
        model_->prediction_metrics(&metrics);
        metrics.print_full();

        // Start
        for (int epoch = 1; epoch <= config.max_epoch; ++epoch) {
            // Update parameters
            T_ = config.T0 / (1.0 + epoch);
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
            metrics.mh_rate = 1.0 * accept_count_ / model_->graph.num_nodes();
            if (epoch % config.eval_every == 0) {
                model_->sparsity_metrics(&metrics);
                model_->prediction_metrics(&metrics);
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
    void mh_step(int node_id) {
        // Propose
        SArray<int> proposal = draw_uniform();
        // Accept
        double h_new = model_->single_objective(node_id, proposal);
        double h_old = model_->single_objective(node_id, model_->Z.row(node_id));
        double log_accept_prob = (h_new - h_old) / T_;
        if (log(rand_.next_double()) < log_accept_prob) {
            model_->Z.update_row(node_id, proposal); // does not need col update
            ++accept_count_;
        }
    }

    /**
     * Draw a random binary vector of length K
     *
     * @return random sparse binary vector
     */
    SArray<int> draw_uniform() {
        SArray<int> sample;
        for (int k = 0; k < model_->Z.col_size(); ++k) {
            if (rand_.next_double() < 0.5) {
                sample.push_back(k);
            }
        }
        return sample;
    }

    // Print short info
    void print_info() {
        LOG("Algorithm: UniformMetropolis");
    }

private:
    Model *model_; // pointer since algorithm modifies the model
    double T_; // temperature for simulated annealing
    int accept_count_; // count MH proposal acceptance
    Random rand_; // random engine
};


#endif //CBGF_UNIFORMMETROPOLIS_H
