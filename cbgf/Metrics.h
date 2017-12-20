//
// Created by Xun Zheng on 10/10/16.
//

#ifndef CBGF_METRICS_H
#define CBGF_METRICS_H


#include <string>
#include "logging.h"

/**
 * Some metrics of the algorithm, model, and prediction.
 */
struct Metrics {
    // Epoch and time
    int epoch = 0;
    double epoch_time = 0;
    double elapsed_time = 0;

    // Algorithm statistics
    double T = 0;
    double delta = 0;
    double mh_rate = 0;

    // Sparsity patterns
    double nnz_per_row = 0;
    double nnz_per_col = 0;
    double empty_row_ratio = 0;

    // objectives
    double obj = 0;
    double norm_acc = 0;

    // Prediction statistics
    double accuracy = 0;
    double precision = 0;
    double recall = 0;
    double f1 = 0;
    double eval_time = 0;

    void print_header() {
        LOG("%s", std::string(6 + (1 + 12) * 15, '-').c_str());
        LOG("%6s %12s %12s  | %9s %12s %12s  | %9s %12s %12s  | %9s %12s  | %9s %12s %12s %12s %12s",
            "epoch", "seconds", "elapsed",
            "T", "delta", "mh_rate",
            "nnz/row", "nnz/col", "empty_row",
            "obj", "norm_acc",
            "accuracy", "precision", "recall", "F1", "eval_time");
    }

    void print_simple() {
        LOG("%6d %12.4g %12.6g  | %9.4g %12.4g %12.4g  | %9.4g %12.4g %12.4g  | %9s %12s  | %9s %12s %12s %12s %12s",
            epoch, epoch_time, elapsed_time,
            T, delta, mh_rate,
            nnz_per_row, nnz_per_col, empty_row_ratio,
            "-", "-",
            "-", "-", "-", "-", "-");
    }

    void print_full() {
        LOG("%6d %12.4g %12.6g  | %9.4g %12.4g %12.4g  | %9.4g %12.4g %12.4g  | %9.4g %12.4g  | %9.4g %12.4g %12.4g %12.4g %12.4g",
            epoch, epoch_time, elapsed_time,
            T, delta, mh_rate,
            nnz_per_row, nnz_per_col, empty_row_ratio,
            obj, norm_acc,
            accuracy, precision, recall, f1, eval_time);
    }
};


#endif //CBGF_METRICS_H
