//
// Created by Xun Zheng on 10/10/16.
//

#ifndef CBGF_FILEPARSER_H
#define CBGF_FILEPARSER_H


#include <cstdio>
#include <string>
#include <cstdlib>
#include "logging.h"
#include "SArray.h"

/**
 * Some common space characters to serve as delimeters for Tokenizer.
 */
static const char SPACES[] = " \f\n\r\t\v";

/**
 * Buffered file reader that parses space-separaed fields for each line.
 * Typical usage:
 *
 * FileParser fp(filename);
 * while (fp.next_line()) {
 *     auto fields = fp.fields();
 *     // do something
 * }
 */
class FileParser {
public:
    /**
     * Construct parser by opening a file.
     *
     * @param filename filename
     */
    FileParser(const std::string &filename) {
        file_ = fopen(filename.c_str(), "r");
        CHECK(file_ != nullptr, "Open failed: %s", filename.c_str());
    }

    /**
     * Destructor. Close the file and release memory.
     */
    ~FileParser() {
        fclose(file_);
        free(line_str_);
    }

    /**
     * Check and read the next line.
     * If there is next line, split it to prepare for space-separated fields.
     *
     * @return true iff there is next line
     */
    bool next_line() {
        auto r = getline(&line_str_, &num_bytes_, file_);
        if (r == -1) {
            return false;
        } else {
            ++line_;
            fields_.clear();
            fields_.reserve(r); // line length is an upper bound of num fields
            char *token = strtok(line_str_, SPACES);
            while (token != nullptr) {
                fields_.push_back(token);
                token = strtok(nullptr, SPACES);
            }
            return true;
        }
    }

    /**
     * Get the current line number of the file. Useful in debugging.
     *
     * @return current line number
     */
    int line() {
        return line_;
    }

    /**
     * Get space-separated fields of the current line.
     * Each entry is a C-style string that lies in line_str_, so no need to release memory.
     *
     * @return space-separated fields of the current line
     */
    SArray<char *> fields() {
        return fields_;
    }

private:
    FILE *file_ = nullptr; // file pointer
    char *line_str_ = nullptr; // buffer for line content
    size_t num_bytes_ = 0; // number of bytes read in getline()
    int line_ = 0; // current line number
    SArray<char *> fields_; // space-separated fields for the current line
};

/**
 * Convert C-style string to int. Just a shorthand for strtol().
 *
 * @param input C-style string
 * @return integer
 */
int str2int(char *input) {
    return (int) strtol(input, nullptr, 10);
}

/**
 * Convert C-style string to double. Just a shorthand for strtod().
 * @param input C-style string
 * @return double
 */
double str2double(char *input) {
    return strtod(input, nullptr);
}


#endif //CBGF_FILEPARSER_H
