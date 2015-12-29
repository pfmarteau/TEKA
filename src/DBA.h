/**************************************************************************************************
Translated from java to C++ from source code available at https://github.com/fpetitjean/DBA
For more details, c.f. F. Petitjean & P. Gançarski, “Summarizing a Set of Time Series by Averaging:
from Steiner Sequence to Compact Multiple Alignment,” Theoretical Computer Science; 2012.
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
**************************************************************************************************/


#ifndef DBA_H
#define DBA_H

#include <iostream>
#include <cfloat>
#include <stdlib.h>
#include <stdio.h>
#include "UniformTimeseries.h"
#include "DTW.h"
#include "DistanceL2.h"

#define PI 3.14159265359
#define NIL -1
#define DIAGONAL 0
#define LEFT 1
#define UP 2

class DBA {
public:
	static const long serialVersionUID = 1L;

private:

	/**
	 * This attribute is used in order to initialize only once the matrixes
	 */
	static const int MAX_SEQ_LENGTH = 2000;

	/**
	 * store the cost of the alignment
	 */
	std::vector<std::vector<double>> costMatrix;
	
	/**
	 * store the warping path
	 */
	std::vector<std::vector<int>>  pathMatrix; //[MAX_SEQ_LENGTH][MAX_SEQ_LENGTH];

	/**
	 * Store the length of the optimal path in each cell
	 */
	std::vector<std::vector<int>> optimalPathLength; //[MAX_SEQ_LENGTH][MAX_SEQ_LENGTH];
	
public:
	DBA();
	UniformTimeseries  computeAverage(UniformTimeseries C, std::vector<UniformTimeseries> sequences);
	static void dba(std::vector<double> C, std::vector<std::vector<double>> sequences);
	static void test();

};

#endif
