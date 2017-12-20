//
// Created by Xun Zheng on 10/10/16.
//

#ifndef CBGF_MODEL_H
#define CBGF_MODEL_H


#include <atomic>
#include <mutex>
#include <thread>
#include "Config.h"
#include "Graph.h"
#include "SparseBinaryMatrix.h"
#include "util.h"
#include "Metrics.h"

/**
 * Binary graph factorization model.
 */
struct Model {
    // Make members public for convenience
    Graph graph;
    SparseBinaryMatrix Z;
    double weight;

    /**
     * Initialize the model from config.
     * Load graph and randomly initialize Z.
     *
     * @param config config
     */
    void init_from_config(const Config &config) {
        // Load graph
        LOG("Loading graph");
        if (!config.edgelist_file.empty()) {
            graph.init_from_edgelist(config.edgelist_file);
        } else {
            CHECK(false, "No input graph specified");
        }
        LOG("Loading graph: done");

        // Initialize Z
        Z.resize(graph.num_nodes(), config.K);
        if (!config.init_from_rows_file.empty()) {
            LOG("Initializing Z from rows");
            Z.init_from_rows(config.init_from_rows_file);
            LOG("Initializing Z from rows: done");
        } else if (!config.init_from_cols_file.empty()) {
            LOG("Initializing Z from columns");
            Z.init_from_cols(config.init_from_cols_file, graph.name_id_map());
            LOG("Initializing Z from columns: done");
        } else {
            LOG("Initializing Z randomly");
            Z.init_randomly();
            LOG("Initializing Z randomly: done");
        }

        // Determine weight on true positives
        if (config.weight == 0) {
            weight = graph.num_nonedges() / graph.num_edges();
        } else {
            weight = config.weight;
        }
        LOG("Setting weight = %g", weight);
    }

    /**
     * Evaluate objective for a single term.
     * Enumerate all nodes, thus complexity is O(n * nnz/row).
     *
     * @param node_id node id
     * @param feature factor of the vertex
     * @return objective value for node_id
     */
    double single_objective(int node_id, const SArray<int> &feature) {
        double num_true_positive = 0.0, num_true_negative = 0.0;
        for (int i = 0; i < graph.num_nodes(); ++i) {
            int a = (graph.is_neighbors(i, node_id) || i == node_id);
            int a_hat = has_common_element(Z.row(i), feature);
            num_true_positive += (a == 1 && a_hat == 1);
            num_true_negative += (a == 0 && a_hat == 0);
        }
        return weight * num_true_positive + num_true_negative;
    }

    /**
     * Evaluate full objective and some predictive metrics.
     * Enumerate all nodes, thus complexity is O(n^2 * nnz/row).
     *
     * @param mt metrics to be filled
     */
    void prediction_metrics(Metrics *mt) {
        auto t1 = std::chrono::steady_clock::now();
        double num_true_positive = 0.0, num_true_negative = 0.0, num_positive = 0.0;
        for (int i = 0; i < graph.num_nodes(); ++i) {
            for (int j = i; j < graph.num_nodes(); ++j) { // Upper triangle
                int a = (graph.is_neighbors(i, j) || i == j);
                int a_hat = has_common_element(Z.row(i), Z.row(j));
                num_true_positive += (a == 1 && a_hat == 1);
                num_true_negative += (a == 0 && a_hat == 0);
                num_positive += (a_hat == 1);
            }
        }
        mt->obj = 0.5 * (weight * num_true_positive + num_true_negative) / graph.num_nonedges();
        mt->norm_acc = 0.5 * (num_true_positive / graph.num_edges() + num_true_negative / graph.num_nonedges());
        mt->accuracy = 2 * (num_true_positive + num_true_negative) / graph.num_nodes() / (graph.num_nodes() + 1.0);
        mt->precision = num_true_positive / num_positive;
        mt->recall = num_true_positive / graph.num_edges();
        mt->f1 = 2 * mt->precision * mt->recall / (mt->precision + mt->recall);
        auto t2 = std::chrono::steady_clock::now();
        mt->eval_time = std::chrono::duration<double>(t2 - t1).count();
    }

    /**
     * Evaluate objective for a single term.
     * Utilize sparsity of Z. Complexity: O(nnz/n * (d + nnz/k)).
     *
     * @param node_id node id
     * @param feature factor of the vertex
     * @return objective value for node_id
     */
    double single_objective_fast(int node_id, const SArray<int> &feature) {
        double num_true_positive = 0.0;
        for (int j : graph.neighbors(node_id)) {
            num_true_positive += has_common_element(feature, Z.row(j));
        }
        num_true_positive += has_common_element(feature, Z.row(node_id));
        double num_positive = Z.count_positive_prediction(feature);
        double num_negative = graph.num_nodes() - num_positive;
        double num_false_negative = graph.degree(node_id) + 1 - num_true_positive;
        double num_true_negative = num_negative - num_false_negative;
        return weight * num_true_positive + num_true_negative;
    }

    /**
     * Evaluate full objective and some predictive metrics.
     * Equivalent to calling single_objective() for each node, i.e. every entry in the adjacency matrix.
     * However, we want results on upper-triangular only.
     * Therefore we add additional diagonal prediction and then divide the results by 2.
     *
     * @param mt metrics to be filled
     */
    void prediction_metrics_fast(Metrics *mt) {
        auto t1 = std::chrono::steady_clock::now();
        double tp = 0.0, tn = 0.0, pos = 0.0;
        for (int i = 0; i < graph.num_nodes(); ++i) {
            double num_true_positive = 0.0;
            for (int j : graph.neighbors(i)) {
                num_true_positive += has_common_element(Z.row(i), Z.row(j));
            }
            int self_positive = has_common_element(Z.row(i), Z.row(i));
            num_true_positive += self_positive;
            double num_positive = Z.count_positive_prediction(Z.row(i));
            double num_negative = graph.num_nodes() - num_positive;
            double num_false_negative = graph.degree(i) + 1 - num_true_positive;
            double num_true_negative = num_negative - num_false_negative;
            tp += num_true_positive + self_positive; // add additional diagonal
            tn += num_true_negative; // additional diagonal never appears since condition = 0
            pos += num_positive + self_positive; // add additional diagonal
        }
        tp /= 2; tn /= 2; pos /= 2;
        mt->obj = 0.5 * (weight * tp + tn) / graph.num_nonedges();
        mt->norm_acc = 0.5 * (tp / graph.num_edges() + tn / graph.num_nonedges());
        mt->accuracy = 2 * (tp + tn) / graph.num_nodes() / (graph.num_nodes() + 1.0);
        mt->precision = tp / pos;
        mt->recall = tp / graph.num_edges();
        mt->f1 = 2 * mt->precision * mt->recall / (mt->precision + mt->recall);
        auto t2 = std::chrono::steady_clock::now();
        mt->eval_time = std::chrono::duration<double>(t2 - t1).count();
    }

    /**
     * Evaluate full objective and some predictive metrics, in parallel.
     *
     * @param cpu num of threads to use
     * @param mt metrics to be filled
     */
    void prediction_metrics_fast_parallel(int cpu, Metrics *mt) {
        auto t1 = std::chrono::steady_clock::now();
        double tp = 0.0, tn = 0.0, pos = 0.0;
        std::atomic<int> atomic_id(0);
        std::mutex mut;
        std::vector<std::thread> thrs(cpu);
        for (auto &thr : thrs) {
            thr = std::thread([&]() {
                double local_tp = 0.0, local_tn = 0.0, local_pos = 0.0;
                while (true) {
                    int node_id = atomic_id++;
                    if (node_id >= graph.num_nodes()) {
                        break;
                    }
                    double num_true_positive = 0.0;
                    for (int j : graph.neighbors(node_id)) {
                        num_true_positive += has_common_element(Z.row(node_id), Z.row(j));
                    }
                    int self_positive = has_common_element(Z.row(node_id), Z.row(node_id));
                    num_true_positive += self_positive;
                    double num_positive = Z.count_positive_prediction(Z.row(node_id));
                    double num_negative = graph.num_nodes() - num_positive;
                    double num_false_negative = graph.degree(node_id) + 1 - num_true_positive;
                    double num_true_negative = num_negative - num_false_negative;
                    local_tp += num_true_positive + self_positive; // add additional diagonal
                    local_tn += num_true_negative; // additional diagonal never appears since condition = 0
                    local_pos += num_positive + self_positive; // add additional diagonal
                }
                std::lock_guard<std::mutex> guard(mut);
                tp += local_tp;
                tn += local_tn;
                pos += local_pos;
            });
        }
        for (auto &thr : thrs) {
            thr.join();
        }
        tp /= 2; tn /= 2; pos /= 2;
        mt->obj = 0.5 * (weight * tp + tn) / graph.num_nonedges();
        mt->norm_acc = 0.5 * (tp / graph.num_edges() + tn / graph.num_nonedges());
        mt->accuracy = 2 * (tp + tn) / graph.num_nodes() / (graph.num_nodes() + 1.0);
        mt->precision = tp / pos;
        mt->recall = tp / graph.num_edges();
        mt->f1 = 2 * mt->precision * mt->recall / (mt->precision + mt->recall);
        auto t2 = std::chrono::steady_clock::now();
        mt->eval_time = std::chrono::duration<double>(t2 - t1).count();
    }

    /**
     * Report some sparsity metrics of Z.
     *
     * @param mt metrics to be filled
     */
    void sparsity_metrics(Metrics *mt) {
        mt->nnz_per_row = Z.nnz_per_row();
        mt->nnz_per_col = Z.nnz_per_col();
        mt->empty_row_ratio = Z.empty_row_ratio();
    }

    /**
     * Write Z in the output file.
     * If no output filename is specified, do nothing.
     *
     * @param config config file that specifies the output filename
     */
    void dump_matrix(const Config &config) {
        if (!config.output_by_rows_file.empty()) {
            LOG("Writing rows of Z to %s", config.output_by_rows_file.c_str());
            Z.output_by_rows(config.output_by_rows_file);
            LOG("Writing rows of Z: done");
        }
        if (!config.output_by_cols_file.empty()) {
            LOG("Writing cols of Z to %s", config.output_by_cols_file.c_str());
            Z.output_by_cols(config.output_by_cols_file, graph.id_name_map());
            LOG("Writing cols of Z: done");
        }
    }

    // Print short info about graph and Z
    void print_info() {
        graph.print_info();
        Z.print_info();
    }
};


#endif //CBGF_MODEL_H
