cmake_minimum_required(VERSION 3.6)
project(cbgf)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=c++11")
find_package(Threads REQUIRED)

set(SOURCE_FILES main.cpp Config.h FileParser.h Graph.h Metrics.h Model.h SArray.h SparseBinaryMatrix.h Random.h logging.h util.h BernoulliMetropolis.h BernoulliMetropolisBulk.h UniformMetropolis.h)
add_executable(cbgf ${SOURCE_FILES})
target_link_libraries(cbgf Threads::Threads)
