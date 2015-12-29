/*************************************************************************
eKATS1.0 (source code generated 2015-12-15)
Copyright (c) Pierre-François Marteau (eKATS project).

>>> SOURCE LICENSE >>>
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation (www.fsf.org); either version 2 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

A copy of the GNU General Public License is available at
http://www.fsf.org/licensing/licenses
>>> END OF LICENSE >>>
*************************************************************************/
#ifndef UTILS_H
#define UTILS_H
#include <memory>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <dirent.h>
#include <sys/stat.h>
#include "UniformTimeseries.h"

std::vector<int> shuffleIndex(int N);
void saveFileMatrix(std::string fileN, std::vector<std::vector<double>> mat);

#endif
