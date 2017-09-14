/*
 * utils.h
 *
 *  Created on: Jan 10, 2014
 *      Author: randomizer
 */

#ifndef UTILS_H_
#define UTILS_H_

#include <string>

namespace utils {

/**
 * Checks if a folder exists
 * @param foldername path to the folder to check.
 * @return true if the folder exists, false otherwise.
 */
bool folder_exists(std::string foldername);

/**
 * Recursive wrapper for mkdir.
 * @param[in] path the full path of the directory to create.
 * @return zero on success, otherwise -1.
 */
int mkdir_recursive(const char *path);

extern const char SEPARATOR;

} /* namespace utils */

#endif /* UTILS_H_ */
