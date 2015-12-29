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

#ifndef KDTW_H
#define KDTW_H

#include <limits>
#include <cmath>
#include <cfloat>
#include "UniformTimeseries.h"
#include "matrix.h"


/*
 * KDTW corresponds to the implementation of the time elastic p.d. kernel close to the DTW elastic measure defined in [1].
 * The .cpp file just contains the main implementation of the algorithm.
 * [1] Marteau, P.F., Gibet, S., On Recursive Edit Distance Kernels with Application to Time Series Classification,
 *     Accepted for publication (June 2014) in IEEE Trans. on Neural Networks and Learning Systems
 */
class KDTW {
public:

	static double kernel(const UniformTimeseries &ts1,
                                 const UniformTimeseries &ts2,
                                 double sigma);
	static double kernel_t(const UniformTimeseries &ts1,
                                 const UniformTimeseries &ts2,
                                 double sigma, double nu);
	static matrix<double> kernel_mat(const UniformTimeseries &ts1,
                                const UniformTimeseries &ts2,
                                double sigma);
	static matrix<double> kernel_AMA(UniformTimeseries ts1, UniformTimeseries ts2, double sigma);
	static matrix<double>  kdistance_mat(const UniformTimeseries &ts1,
                                const UniformTimeseries &ts2,
                                double sigma);
	static double kdistance(const UniformTimeseries &ts1,
                                const UniformTimeseries &ts2,
                                double sigma);
    static UniformTimeseries add_kernelDBA(const UniformTimeseries C, std::vector<UniformTimeseries> sequences,
    		double sigma, double lab);
    static UniformTimeseries add_kernelT(const UniformTimeseries &ts1, const UniformTimeseries &ts2,
    		matrix<double> mat, double gscal, double scal2);
    static UniformTimeseries add_kernelTbis(const UniformTimeseries &ts1, const UniformTimeseries &ts2,
    		matrix<double> mat, double gscal, double scal2);
    static UniformTimeseries averagePairT(UniformTimeseries &ts1, UniformTimeseries &ts2, double sigma);
};




#endif
