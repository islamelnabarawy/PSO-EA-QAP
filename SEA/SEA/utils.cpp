/*
 * utils.cpp
 *
 *  Created on: Jan 10, 2014
 *      Author: randomizer
 */

#include "utils.h"

#include <sstream>
#include <iostream>
#include <cstdlib>

#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>

#if defined(_MSC_VER) || defined(WIN32)
#include <direct.h> // for _mkdir() on Windows
#endif

#ifdef _MSC_VER
	// VC++ uses _mkdir() instead of mkdir()
	#define mkdir(x) _mkdir(x)
#elif !defined(WIN32)
	// unified mkdir() definition with a single argument
	int mkdir(const char* path) { return mkdir(path, 0755); }
#endif

namespace utils {

#if defined(_MSC_VER) || defined(WIN32)
const char SEPARATOR = '\\'; // Windows directory separator
#else
const char SEPARATOR = '/'; // Linux directory separator
#endif

/**
 * Checks if a folder exists
 * @param foldername path to the folder to check.
 * @return true if the folder exists, false otherwise.
 */
bool folder_exists(std::string foldername)
{
	struct stat st;
	if (stat(foldername.c_str(), &st) == -1) return false;
	return st.st_mode & S_IFDIR;
}

/**
 * Recursive wrapper for mkdir.
 * @param[in] path the full path of the directory to create.
 * @return zero on success, otherwise -1.
 */
int mkdir_recursive(const char *path)
{
	std::string current_level = "";
    std::string level;
    std::stringstream ss(path);

    // split path using slash as a separator
    while (std::getline(ss, level, '/'))
    {
        current_level += level; // append folder to the current level

        // create current level
        if (!folder_exists(current_level) && mkdir(current_level.c_str()) != 0) {
			std::cerr << "Error creating folder: " << current_level.c_str() << std::endl;
            return -1;
		}

        current_level += SEPARATOR; // don't forget to append a separator
    }

    return 0;
}

} /* namespace utils */
