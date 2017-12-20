//
// Created by Xun Zheng on 10/10/16.
//

#ifndef CBGF_SPARSEBINARYMATRIX_H
#define CBGF_SPARSEBINARYMATRIX_H


#include <unordered_map>
#include <vector>
#include <fstream>
#include <queue>
#include "SArray.h"
#include "Random.h"
#include "FileParser.h"

/**
 * Sparse N x K binary matrix, represented as N sparse vectors.
 * Some algorithms may also use K sparse vectors as column view in addition, but some don't.
 * So all member functions should rely solely on rows, unless specified otherwise.
 */
class SparseBinaryMatrix {
public:
    /**
     * Empty constructor.
     */
    SparseBinaryMatrix() {}

    /**
     * Constructor with size.
     * @param N row size
     * @param K column size
     */
    SparseBinaryMatrix(int N, int K) {
        resize(N, K);
    }

    /**
     * Initialize a matrix of a given size.
     *
     * @param N row size
     * @param K column size
     */
    void resize(int N, int K) {
        N_ = N;
        K_ = K;
        rows_.resize(N_);
        cols_.resize(K_);
    }

    /**
     * Set one nonzero randomly for each row.
     */
    void init_randomly() {
        Random r;
        for (int n = 0; n < N_; ++n) {
            int k = r.next_int(K_);
            rows_[n] = {k};
        }
        reset_cols_from_rows();
    }

    /**
     * Init from the file that contains row view of the matrix.
     * Each line contains nonzero column id's of a row, separated by spaces.
     * Assume file content is 1-indexed.
     *
     * @param filename row file
     */
    void init_from_rows(const std::string &filename) {
        FileParser fp(filename);
        int n = 0;
        while (fp.next_line()) {
            auto fields = fp.fields();
            SArray<int> row(fields.size());
            for (int i = 0; i < fields.size(); ++i) {
                int k = str2int(fields[i]) - 1; // convert to 0-index
                CHECK(k >= 0 && k < K_, "Invalid k = %d (K = %d), at line %d", k, K_, fp.line());
                row[i] = k;
            }
            std::sort(row.begin(), row.end());
            rows_[n] = row;
            ++n;
        }
        CHECK(n == N_, "Invalid num rows = %d (N = %d)", n, N_);
        reset_cols_from_rows();
    }

    /**
     * Init from the file that contains column view of the matrix.
     * Each line contains nonzero node names of a column, separated by spaces.
     * Translate node names to unique id's using the dict obtained from graph.
     *
     * @param filename column file
     */
    void init_from_cols(const std::string &filename, const std::unordered_map<std::string, int> &name_id_map) {
        FileParser fp(filename);
        int k = 0;
        while (fp.next_line()) {
            auto fields = fp.fields();
            for (int i = 0; i < fields.size(); ++i) {
                int n = name_id_map.at(fields[i]);
                CHECK(n >= 0 && n < N_, "Invalid n = %d (N_ = %d), at line %d", n, N_, fp.line());
                rows_[n].push_back(k);
            }
            ++k;
        }
        CHECK(k == K_, "Invalid num cols = %d (K_ = %d)", k, K_);
        for (auto &row : rows_) {
            std::sort(row.begin(), row.end());
            row.shrink_to_fit();
        }
        reset_cols_from_rows();
    }

    /**
     * Adapt columns to rows.
     * Erase all elements, reconstruct from rows, then sort.
     */
    void reset_cols_from_rows() {
        for (auto &col : cols_) {
            col.clear();
        }
        for (int n = 0; n < N_; ++n) {
            for (int k : rows_[n]) {
                cols_[k].push_back(n);
            }
        }
        for (auto &col : cols_) {
            std::sort(col.begin(), col.end());
            col.shrink_to_fit();
        }
    }

    /**
     * Write nonzeros of each row as a line in a file.
     * Nonzeros are separated by spaces, sorted in ascending order, 1-indexed.
     *
     * @param filename row file
     */
    void output_by_rows(const std::string &filename) {
        // Copy the matrix, make entries 1-indexed
        std::vector<std::vector<int>> rows(N_);
        for (int n = 0; n < N_; ++n) {
            for (int k : rows_[n]) {
                rows[n].push_back(k + 1);
            }
        }
        output_nested_container(filename, rows_);
    }

    /**
     * Write nonzeros of each column as a line in a file.
     * Nonzeros are shown in their original names, separated by spaces, not in a particular order.
     *
     * @param filename col file
     */
    void output_by_cols(const std::string &filename, const std::unordered_map<int, std::string> &id_name_map) {
        // Build column view, translate node id to node name
        std::vector<std::vector<std::string>> cols(K_);
        for (int n = 0; n < N_; ++n) {
            for (int k : rows_[n]) {
                cols[k].push_back(id_name_map.at(n));
            }
        }
        output_nested_container(filename, cols);
    }

    /**
     * Get n-th row.
     *
     * @param n row id
     * @return n-th row
     */
    SArray<int> row(int n) {
        return rows_[n];
    }

    /**
     * Average number of nonzeros per row.
     *
     * @return average number of nonzeros per row
     */
    double nnz_per_row() {
        double count = 0.0;
        for (const auto &row : rows_) {
            count += row.size();
        }
        return count / N_;
    }

    /**
     * Average number of nonzeros per column.
     *
     * @return average number of nonzeros per column
     */
    double nnz_per_col() {
        double count = 0.0;
        for (const auto &row : rows_) {
            count += row.size();
        }
        return count / K_;
    }

    /**
     * Ratio of empty rows among all rows.
     * Empty row means no assignment has been made.
     *
     * @return proportion of empty rows
     */
    double empty_row_ratio() {
        double count = 0;
        for (const auto &row : rows_) {
            if (row.empty()) {
                ++count;
            }
        }
        return count / N_;
    }

    /**
     * Sum of all rows, i.e. vector of column size
     *
     * @return SArray of column size
     */
    SArray<int> rowsum() {
        SArray<int> count(K_);
        for (const auto &row : rows_) {
            for (const auto &k : row) {
                count[k] += 1;
            }
        }
        return count;
    }

    /**
     * Set a row to a new value.
     * Never use ```Z.row(42) = new_row;``` !!!
     *
     * @param n row id
     * @param new_row value of the new row
     */
    void update_row(int n, const SArray<int> &new_row) {
        rows_[n] = new_row;
    }

    /**
     * Given old and new values of the n-th row, update columns accordingly.
     * Remove n from column id's contained in old row, and add n to column id's in new row.
     * Since rows are sorted, use 2-way merge.
     * Also since columns are sorted, pay attention to the order after add/removal of n.
     * Each add/removal in a column has complexity O(nnz/col).
     * Note, only columns are updated. New row should be updated separately.
     *
     * @param n row id
     * @param old_row old row
     * @param new_row new row
     */
    void update_col_from_row(int n, const SArray<int> &old_row, const SArray<int> &new_row) {
        auto it_old = old_row.begin();
        auto it_new = new_row.begin();
        while (it_old != old_row.end() && it_new != new_row.end()) {
            if (*it_old < *it_new) {
                // remove n from col[k_old]
                int k = *it_old;
                auto before = cols_[k];
                std::remove(cols_[k].begin(), cols_[k].end(), n);
                auto after1 = cols_[k];
                cols_[k].resize(cols_[k].size() - 1);
                auto after2 = cols_[k];
                ++it_old;
            } else if (*it_old > *it_new) {
                // add n to col[k_new]
                int k = *it_new;
                auto before = cols_[k];
                cols_[k].push_back(n);
                auto after1 = cols_[k];
                std::inplace_merge(cols_[k].begin(), cols_[k].end() - 1, cols_[k].end());
                auto after2 = cols_[k];
                ++it_new;
            } else {
                // do nothing
                ++it_old;
                ++it_new;
            }
        }
        while (it_old != old_row.end()) {
            // remove n from col[k_old]
            int k = *it_old;
            std::remove(cols_[k].begin(), cols_[k].end(), n);
            cols_[k].resize(cols_[k].size() - 1);
            ++it_old;
        }
        while (it_new != new_row.end()) {
            // add n to col[k_new]
            int k = *it_new;
            cols_[k].push_back(n);
            std::inplace_merge(cols_[k].begin(), cols_[k].end() - 1, cols_[k].end());
            ++it_new;
        }
    }

    /**
     * Compute nnz(y), where y = Z * vec, in complexity O(nnz/row * nnz/col) by exploiting sparsity.
     * Let c = nonzeros(vec), then y = sum_{i in c} Z.col(i).
     * Since columns are sorted, summation can be done by k-way merge.
     * Also since only nnz(y) is needed, no need to keep all sums, memorizing the most recent one suffices.
     *
     * @param vec sparse vector with true size K
     * @return nnz(Z * vec)
     */
    int count_positive_prediction(const SArray<int> &vec) {
        using IterPair = std::pair<int *, int *>; // (iterator,end)
        auto gt = [](IterPair &a, IterPair &b) { return *(a.first) > *(b.first); };
        std::priority_queue<IterPair, std::vector<IterPair>, decltype(gt)> pq(gt);
        for (int k : vec) {
            if (!cols_[k].empty()) {
                pq.emplace(cols_[k].begin(), cols_[k].end());
            }
        }
        int prev_n = -1, num_positive = 0;
        while (!pq.empty()) {
            IterPair item = pq.top();
            int n = *(item.first);
            if (n != prev_n) {
                prev_n = n;
                ++num_positive;
            }
            pq.pop();
            if (item.first + 1 != item.second) {
                pq.emplace(item.first + 1, item.second);
            }
        }
        return num_positive;
    }

    // Getter
    int row_size() {
        return N_;
    }

    // Getter
    int col_size() {
        return K_;
    }

    // Print short info
    void print_info() {
        LOG("SparseBinaryMatrix: N = %d, K = %d", N_, K_);
    }

    // Debug
    void debug() {
        print_info();
        for (int n = 0; n < N_; ++n) {
            LOG("row %d: %s", n, rows_[n].str().c_str());
        }
    }

private:
    /**
     * Write the content of each inner container as a line in a file.
     * Items are separated by spaces.
     *
     * @param filename filename
     * @param nested_container a container containing containers
     */
    template<typename T>
    void output_nested_container(const std::string &filename, const T &nested_container) {
        std::ofstream ofs(filename);
        CHECK(ofs.is_open(), "Failed to open file: %s", filename.c_str());
        for (const auto &container : nested_container) {
            for (auto it = container.begin(); it != container.end(); ++it) {
                if (it != container.begin()) {
                    ofs << " ";
                }
                ofs << *it;
            }
            ofs << std::endl;
        }
        ofs.close();
    }

private:
    int N_; // row size
    int K_; // column size
    std::vector<SArray<int>> rows_; // rows, sorted in ascending order
    std::vector<SArray<int>> cols_; // cols, sorted in ascending order
};


#endif //CBGF_SPARSEBINARYMATRIX_H
