/*************************************************************************
eKATS0.1 (source code generated 2015-12-15)
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
#ifndef KCENTROID_H
#define KCENTROID_H

#include <iostream>
#include <cfloat>
#include <stdlib.h>
#include <stdio.h>
#include "UniformTimeseries.h"
#include "DTW.h"
#include "KDTW.h"
#include "GramMatrix.h"
#include "DBA.h"
#include "kmedoids.h"
#include "utils.h"

std::vector<UniformTimeseries>  kcentroid_idba(std::vector<UniformTimeseries> dsTrain, int nclass, int nclust);
std::vector<UniformTimeseries> kadd2by2div2(std::vector<UniformTimeseries> ds, double sigma, int lab);
UniformTimeseries getKCentroid(std::vector<UniformTimeseries> ds, double sigma, int lab, int nShuf);
UniformTimeseries getKDTWCentroid(UniformTimeseries C, std::vector<UniformTimeseries> ds, double sigma);
UniformTimeseries getDTWCentroid(UniformTimeseries averageSequence, std::vector<UniformTimeseries> ds);
std::vector<UniformTimeseries>  kcentroid_kdtw_pwa(std::vector<UniformTimeseries> dsTrain, int nclass, int nclust, double sigma);
std::vector<UniformTimeseries>  kcentroid_ikdba(std::vector<UniformTimeseries> dsTrain, int nclass, int nclust, double sigma);

#endif
