//
// Created by Xun Zheng on 10/10/16.
//

#ifndef CBGF_CONFIG_H
#define CBGF_CONFIG_H


#include <string>
#include <fstream>
#include <thread>
#include "logging.h"
#include "FileParser.h"

/**
 * Stores all the configurations required by the algorithm.
 */
struct Config {
    int cpu = 1; // num of threads to use
    int K = 100; // number of communities
    int max_epoch = 500; // max num of sweep over the vertices
    double T0 = 1.0; // temperature schedule: T = T0 / (1 + epoch)
    double delta0 = 1e-8; // delta schedule: delta = delta0 / log(1 + epoch)
    double weight = 0; // weight on true positives, default is #nonedge/#edge
    int eval_every = 1; // evaluate metrics every X epochs
    std::string algorithm = ""; // algorithm type: {uniform, bernoulli, bulk}
    std::string edgelist_file = ""; // graph stored as edge list
    std::string init_from_rows_file = ""; // file where each line is a row of Z
    std::string init_from_cols_file = ""; // file where each line is a column of Z
    std::string output_by_rows_file = ""; // one line per vertex, indicating all communities this vertex is a part of
    std::string output_by_cols_file = ""; // one line per community, indicating all members of that community

    /**
     * Reads config from a file, where each line is a key-value pair separated by space.
     *
     * @param filename config filename
     */
    void init_from_file(const std::string &filename) {
        FileParser fp(filename);
        while (fp.next_line()) {
            auto fields = fp.fields();
            // Skip empty lines and comments
            if (fields.empty() || *fields[0] == '%') {
                continue;
            }
            CHECK(fields.size() == 2, "Got %d fields (expect 2) at line %d", fields.size(), fp.line());
            std::string key = fields[0];
            if (key == "cpu") { cpu = str2int(fields[1]); }
            else if (key == "K") { K = str2int(fields[1]); }
            else if (key == "maxEpoch") { max_epoch = str2int(fields[1]); }
            else if (key == "T0") { T0 = str2double(fields[1]); }
            else if (key == "delta0") { delta0 = str2double(fields[1]); }
            else if (key == "weight") { weight = str2double(fields[1]); }
            else if (key == "evalEvery") { eval_every = str2int(fields[1]); }
            else if (key == "algorithm") { algorithm = fields[1]; }
            else if (key == "edgeListFile") { edgelist_file = fields[1]; }
            else if (key == "initFromRowsFile") { init_from_rows_file = fields[1]; }
            else if (key == "initFromColsFile") { init_from_cols_file = fields[1]; }
            else if (key == "outputByRowsFile") { output_by_rows_file = fields[1]; }
            else if (key == "outputByColsFile") { output_by_cols_file = fields[1]; }
            else { CHECK(false, "Unknown config name: %s", fields[0]); }
        }
    }

    /**
     * Log current configs.
     */
    void print_info() {
        LOG("cpu = %d", cpu);
        LOG("K = %d", K);
        LOG("maxEpoch = %d", max_epoch);
        LOG("T0 = %g", T0);
        LOG("delta0 = %g", delta0);
        LOG("weight = %g", weight);
        LOG("evalEvery = %d", eval_every);
        LOG("algorithm = %s", algorithm.c_str());
        LOG("edgeListFile = %s", edgelist_file.c_str());
        LOG("initFromRowsFile = %s", init_from_rows_file.c_str());
        LOG("initFromColsFile = %s", init_from_cols_file.c_str());
        LOG("outputByRowsFile = %s", output_by_rows_file.c_str());
        LOG("outputByColsFile = %s", output_by_cols_file.c_str());
    }
};


#endif //CBGF_CONFIG_H
