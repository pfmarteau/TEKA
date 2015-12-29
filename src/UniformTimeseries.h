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
#ifndef UNIFORMTIMESERIES_H
#define UNIFORMTIMESERIES_H
#include <iostream>
#include <fstream>
#include <sstream>
#include <cassert>
#include <math.h>
#include <cfloat>
#include <map>
#include "Point.h"

enum {_DTW, _KDTW};
enum {_UCR, _MTDS };

class UniformTimeseries {
public:
    std::vector<Point> data; // the sample array
    double label;			 // label (category)
    int subject;			 // subject/source identifier
    UniformTimeseries();
    UniformTimeseries(const std::vector<Point> &data);
    UniformTimeseries(const std::vector<Point> &data, double lab);
    UniformTimeseries(const std::vector<Point> &data, int subject, double lab);
    UniformTimeseries(int L, int D);
    UniformTimeseries reverse();
    std::vector<Point>::size_type size() const {
        return data.size();
    }
    std::vector<double>  sumSamples();
    double norm();
    double norml1();
    double maxOnDim(int k);
    double minOnDim(int k);
    void scalMult(double scal);
    static std::vector<UniformTimeseries> parse_ucr_file(std::string &filename, int limit=0);
    static std::vector<UniformTimeseries> parse_mdts_file(std::string &filename, int dim, int limit, const bool withSubject=true);
    void saveFile(std::string file);
    static void saveUniformDatasetSVM(std::string fname, std::vector<UniformTimeseries> ds);
    static void saveUniformDatasetUCR(std::string fname, std::vector<UniformTimeseries> ds);
    static void saveUniformDatasetInColumns(std::string fname, std::vector<UniformTimeseries> ds);
    static std::vector<UniformTimeseries> extractTSfromLabelDS(std::vector<UniformTimeseries> ds, int lab);
    static std::vector<UniformTimeseries> relabelDS(std::vector<UniformTimeseries>& ds, std::map<int,int>& tab);
    void print();
    };

#endif
