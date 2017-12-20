#include "logging.h"
#include "Config.h"
#include "Model.h"
#include "UniformMetropolis.h"
#include "BernoulliMetropolis.h"
#include "BernoulliMetropolisBulk.h"


int main(int argc, char *argv[]) {
    CHECK(argc == 2, "Usage: ./program config_file");

    Config config;
    config.init_from_file(argv[1]);
    config.print_info();

    Model model;
    model.init_from_config(config);
    model.print_info();

    if (config.algorithm == "uniform") {
        UniformMetropolis alg;
        alg.init_from_model(&model);
        alg.print_info();
        alg.run(config);
    } else if (config.algorithm == "bernoulli") {
        BernoulliMetropolis alg;
        alg.init_from_model(&model);
        alg.print_info();
        alg.run(config);
    } else if (config.algorithm == "bulk") {
        BernoulliMetropolisBulk alg;
        alg.init_from_model(&model);
        alg.print_info();
        alg.run(config);
    } else {
        LOG("Unkown algorithm: %s", config.algorithm.c_str());
    }

    model.dump_matrix(config);

    return 0;
}
