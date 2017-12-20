//
// Created by Xun Zheng on 10/10/16.
//

#ifndef CBGF_UTIL_H
#define CBGF_UTIL_H


#include "SArray.h"

/**
 * Check if two sparse arrays have at least one common element.
 *
 * @param a1 a sorted array
 * @param a2 a sorted array
 * @return true iff two arrays have at least one common element
 */
bool has_common_element(const SArray<int> &a1, const SArray<int> &a2) {
    auto head1 = a1.begin();
    auto head2 = a2.begin();
    while (head1 != a1.end() && head2 != a2.end()) {
        if (*head1 < *head2) { ++head1; }
        else if (*head1 > *head2) { ++head2; }
        else { return true; }
    }
    return false;
}


#endif //CBGF_UTIL_H
