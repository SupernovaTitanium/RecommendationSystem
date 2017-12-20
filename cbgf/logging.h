//
// Created by Xun Zheng on 10/10/16.
//

#ifndef CBGF_LOGGING_H
#define CBGF_LOGGING_H


#include <cstdio>
#include <cstdlib>
#include <ctime>

/**
 * A logging macro with prefix of the form "mm/dd hh:mm:ss ". Also appends newline.
 */
#define LOG(fmt, args...) do {                                     \
    time_t time_since_epoch = time(NULL);                          \
    struct tm* tm_info = localtime(&time_since_epoch);             \
    fprintf(stdout, "%02d/%02d %02d:%02d:%02d " fmt "\n",          \
            1+tm_info->tm_mon, tm_info->tm_mday, tm_info->tm_hour, \
            tm_info->tm_min, tm_info->tm_sec, ##args);             \
} while (0)

/**
 * Assertion followed by an explanation.
 */
#define CHECK(clause, fmt, args...) do {                                \
    if (!(clause)) {                                                    \
        LOG("Check failed at %s:%d: " fmt, __FILE__, __LINE__, ##args); \
        exit(EXIT_FAILURE);                                             \
    }                                                                   \
} while (0)


#endif //CBGF_LOGGING_H
