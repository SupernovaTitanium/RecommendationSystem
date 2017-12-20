//
// Created by Xun Zheng on 10/10/16.
//

#ifndef CBGF_SARRAY_H
#define CBGF_SARRAY_H


#include <memory>
#include <sstream>
#include <vector>
#include <algorithm>
#include <cstring>
#include "logging.h"

/**
 * Shared array using smart pointers.
 * Safe to use as return value since copy constructor is zero-copy.
 * Only use with primitive types as it relies on memcpy() and memset().
 * Thanks to dmlc/ps-lite which served as a reference to this implementation.
 */
template<typename T>
class SArray {
public:
    /**
     * Constructs an empty array.
     */
    SArray() {}

    /**
     * Empty destructor. Actual release of memory is handled by ptr_.
     * If the reference count becomes zero, deleter of ptr_ is called to release the memory.
     */
    ~SArray() {}

    /**
     * Construct array with certain size, initialized by zeros.
     *
     * @param size initial size
     */
    SArray(int size) {
        resize(size);
    }

    /**
     * List constructor.
     *
     * @param li initializer list
     */
    SArray(const std::initializer_list<T> &li) {
        int new_size = (int) li.size();
        resize(new_size);
        T *a = begin();
        for (const T &value : li) {
            *a = value;
            ++a;
        }
    }

//    /**
//     * Initialize from std::vector<int>.
//     *
//     * @param vec STL vector
//     */
//    SArray(const std::vector<T>& vec) {
//        int vec_size = (int) vec.size();
//        resize(vec_size);
//        memcpy(array(), vec.data(), vec_size * sizeof(T));
//    }

    /**
     * Copy construct, with zero copy.
     *
     * @param other object to copy from
     */
    SArray(const SArray<T> &other) {
        *this = other;
    }

    /**
     * Copy other to myself, with zero copy.
     *
     * @param other object to copy from
     */
    void operator=(const SArray<T> &other) {
        // Avoid self assignment
        if (this == &other) {
            return;
        }
        size_ = other.size();
        capacity_ = other.capacity();
        ptr_ = other.ptr();
    }

    /**
     * Reset myself to a new array, modifies both size_ and capacity_ to the length of new array.
     *
     * @param new_array new array to point to
     * @param new_size length of the new array
     */
    void reset(T *new_array, int new_size) {
        size_ = new_size;
        capacity_ = new_size;
        ptr_.reset(new_array, [](T *p) { delete[] p; });
    }

    /**
     * Resize and fill new entries to zero if needed.
     *
     * @param new_size new size
     */
    void resize(int new_size) {
        int old_size = size_;
        // Update size, allocate more memory if needed
        if (new_size <= capacity_) {
            size_ = new_size;
        } else {
            T *new_array = new T[new_size];
            memcpy(new_array, array(), size_ * sizeof(T)); // copy old content
            reset(new_array, new_size);
        }
        // If there are new values to fill
        if (old_size - new_size < 0) {
            memset(array() + old_size, 0, (new_size - old_size) * sizeof(T));
        }
    }

    /**
     * Increase capacity if needed.
     *
     * @param new_capacity new capacity
     */
    void reserve(int new_capacity) {
        if (capacity_ >= new_capacity) {
            return;
        }
        int old_size = size_;
        resize(new_capacity);
        size_ = old_size;
    }

    /**
     * Shrink capacity if needed by copying to a new compact array.
     */
    void shrink_to_fit() {
        if (size_ == capacity_) {
            return;
        }
        T *new_array = new T[size_];
        memcpy(new_array, array(), size_ * sizeof(T)); // copy old content
        reset(new_array, size_);
    }

    /**
     * Push an item to the back of the array, increase size_ by 1.
     * Doubles memory if needed.
     *
     * @param value value
     */
    inline void push_back(const T &value) {
        if (size_ == capacity_) {
            reserve(size_ == 0 ? 4 : size_ * 2);
        }
        array()[size_++] = value;
    }

    /**
     * Binary search for a value, assuming array is sorted in ascending order.
     *
     * @param value value to search for
     * @return true iff array contains value
     */
    bool bsearch(const T &value) {
        return std::binary_search(begin(), end(), value);
    }

    /**
     * Detach underlying array, set both size_ and capacity_ to zero.
     */
    void clear() {
        size_ = 0;
        capacity_ = 0;
        ptr_.reset();
    }

    /**
     * Test if underlying array is empty.
     *
     * @return true iff size_ is 0.
     */
    bool empty() const {
        return size_ == 0;
    }

    // Index operator
    inline T &operator[](int i) {
//        CHECK(i < size_, "Invalid index = %d (size_ = %d)", i, size_);
        return array()[i];
    }

    // Index operator (const)
    inline const T &operator[](int i) const {
//        CHECK(i < size_, "Invalid index = %d (size_ = %d)", i, size_);
        return array()[i];
    }

    // Iterator
    T *begin() { return array(); }

    // Iterator (const)
    const T *begin() const { return array(); }

    // Iterator
    T *end() { return array() + size_; }

    // Iterator (const)
    const T *end() const { return array() + size_; }

    // Getter
    int size() const { return size_; }

    // Getter
    int capacity() const { return capacity_; }

    // Getter
    std::shared_ptr<T> &ptr() { return ptr_; }

    // Getter
    const std::shared_ptr<T> &ptr() const { return ptr_; }

    // Getter
    T *array() const { return ptr_.get(); } // actual underlying array

    // Debug str, only contain first and last 5 elements if the size exceeds 10
    std::string str() {
        int m = 10;
        std::stringstream ss;
        ss << "[" << size_ << "/" << capacity_ << "]";
        if (size_ <= 2 * m) {
            for (int i = 0; i < size_; ++i) { ss << " " << array()[i]; }
        } else {
            for (int i = 0; i < m; ++i) { ss << " " << array()[i]; }
            ss << " ...";
            for (int i = size_ - m; i < size_; ++i) { ss << " " << array()[i]; }
        }
        ss << " (use_count: " << ptr_.use_count() << ")";
        return ss.str();
    }

private:
    int size_ = 0;
    int capacity_ = 0;
    std::shared_ptr<T> ptr_;
};


#endif //CBGF_SARRAY_H
