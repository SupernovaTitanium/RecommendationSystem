//
// Created by Xun Zheng on 10/10/16.
//

#ifndef CBGF_GRAPH_H
#define CBGF_GRAPH_H


#include <unordered_map>
#include <sstream>
#include "SArray.h"
#include "FileParser.h"

/**
 * Undirected graph using adjacency list representation.
 * Node id's are 0-indexed.
 */
class Graph {
public:
    /**
     * Clear everything.
     */
    void clear() {
        N_ = 0;
        E_ = 0;
        num_edges_ = 0;
        num_nonedges_ = 0;
        adjacency_list_.clear();
        num_loaded_nodes_ = 0;
        name_id_map_.clear();
        id_name_map_.clear();
    }

    /**
     * Initialize graph from edge list file (e.g. https://snap.stanford.edu/data/ca-AstroPh.txt.gz).
     * Nodes can be represented by any string names, not just integers.
     * Lines starting with '#' are treated as comments.
     * Empty lines are skipped.
     *
     * @param filename edge list file
     */
    void init_from_edgelist(const std::string &filename) {
        clear();
        FileParser fp(filename);
        while (fp.next_line()) {
            auto fields = fp.fields();
            // Ignore empty lines and comments starting with #
            if (fields.empty() || *fields[0] == '#') {
                continue;
            }
            int node_id0 = add_node(fields[0]);
            int node_id1 = add_node(fields[1]);
            if (node_id0 == node_id1) { continue; } // skip self loop
            adjacency_list_[node_id0].push_back(node_id1);
            adjacency_list_[node_id1].push_back(node_id0);
        }
        for (auto &nbrs : adjacency_list_) {
            std::sort(nbrs.begin(), nbrs.end());
            auto last = std::unique(nbrs.begin(), nbrs.end());
            int compact_size = (int) (last - nbrs.begin());
            nbrs.resize(compact_size);
            E_ += compact_size;
            nbrs.shrink_to_fit();
        }
        N_ = num_loaded_nodes_;
        E_ = E_ / 2; // excluding self-loops
        num_edges_ = E_ + N_; // including self-loops
        num_nonedges_ = N_ * ((N_ + 1.0) / 2 - num_edges_ / N_); // avoid int overflow
    }

    /**
     * Try adding a node with a name (string) and map it to a 0-based contiguous id.
     * Return new id if the node name is new, otherwise return existing id.
     *
     * @param name the name (string) of the node
     * @return unique id of the node, 0-index
     */
    int add_node(const std::string &name) {
        auto it = name_id_map_.find(name);
        if (it == name_id_map_.end()) { // new name
            int node_id = num_loaded_nodes_++;
            name_id_map_[name] = node_id;
            id_name_map_[node_id] = name;
            adjacency_list_.emplace_back(); // push an empty SArray
            if (num_loaded_nodes_ % 100000 == 0) {
                LOG("Nodes loaded: %d", num_loaded_nodes_);
            }
            return node_id;
        } else { // existing name
            return it->second;
        }
    }

    // Getter
    int num_nodes() {
        return N_;
    }

    // Getter, includes self-loops
    double num_edges() {
        return num_edges_;
    }

    // Getter
    double num_nonedges() {
        return num_nonedges_;
    }

    // Getter
    std::unordered_map<int, std::string> id_name_map() {
        return id_name_map_;
    };

    // Getter
    std::unordered_map<std::string, int> name_id_map() {
        return name_id_map_;
    };

    /**
     * Get degree of a node.
     *
     * @param i node id
     * @return degree of the node (number of neighbors)
     */
    int degree(int i) {
        return adjacency_list_[i].size();
    }

    /**
     * Get neighbors of a node.
     *
     * @param i node id
     * @return neighbors of the node
     */
    SArray<int> neighbors(int i) {
        return adjacency_list_[i];
    }

    /**
     * Check if two nodes are neighbors by doing binary search on the smaller adjacency list.
     * If i == j, return false too.
     *
     * @param i node id
     * @param j node id
     * @return true iff two nodes are neighbors
     */
    bool is_neighbors(int i, int j) {
        return degree(i) < degree(j) ? adjacency_list_[i].bsearch(j) : adjacency_list_[j].bsearch(i);
    }

    // Print short info
    void print_info() {
        LOG("Graph: N = %d, E = %d, one = %g, zero = %g, one/(one+zero) = %g, zero/one = %g",
            N_, E_, num_edges_, num_nonedges_, num_edges_ / N_ * 2.0 / (N_ + 1.0), num_nonedges_ / num_edges_);
    }

    // Debug
    void debug() {
        print_info();
        int max_line = 10;
        int max_fields = 10;
        for (int i = 0; i < std::min(max_line, N_); ++i) {
            const auto &nbrs = adjacency_list_[i];
            std::stringstream ss;
            ss << "node " << i << "(" << id_name_map_[i] << "):";
            ss << " [" << nbrs.size() << "/" << nbrs.capacity() << "]";
            if (nbrs.size() <= 2 * max_fields) {
                for (int j = 0; j < nbrs.size(); ++j) {
                    ss << " " << nbrs[j] << "(" << id_name_map_[nbrs[j]] << ")";
                }
            } else {
                for (int j = 0; j < max_fields; ++j) {
                    ss << " " << nbrs[j] << "(" << id_name_map_[nbrs[j]] << ")";
                }
                ss << " ...";
                for (int j = nbrs.size() - max_fields; j < nbrs.size(); ++j) {
                    ss << " " << nbrs[j] << "(" << id_name_map_[nbrs[j]] << ")";
                }
            }
            LOG("%s", ss.str().c_str());
        }
    }

private:
    int N_; // number of nodes
    int E_; // number of edges, excluding self-loops
    double num_edges_; // number of edges, including self-loops
    double num_nonedges_; // number of zeros in the upper triangular of adjacency matrix
    std::vector<SArray<int>> adjacency_list_; // adjacency list, each in ascending order
    int num_loaded_nodes_ = 0; // used in graph loading
    std::unordered_map<std::string, int> name_id_map_; // node name to unique id
    std::unordered_map<int, std::string> id_name_map_; // unique id to node name
};


#endif //CBGF_GRAPH_H
