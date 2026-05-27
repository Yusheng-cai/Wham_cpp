#pragma once

#include <iostream>
#include <stdio.h>
#include <string>
#include <unistd.h>

namespace FileSystem {
std::string getCurrentPath();

std::string joinPath(const std::string& path1, const std::string& path2);
}; // namespace FileSystem