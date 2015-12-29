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
#include <limits>
#include <cassert>
#include <stdexcept>
#include <algorithm>
#include "DTW.h"
#include "DistanceL2.h"

using namespace std;

static const double INFTY = numeric_limits<double>::infinity();

// Singleton for use throughout the whole program.
std::shared_ptr<DTW> DTW_DIST = std::shared_ptr<DTW>(new DTW());

// Main implementation of DTW.
double DTW::dist_hint(
        const UniformTimeseries &ts1,
        const UniformTimeseries &ts2,
        double early_abandon)
{
    const bool enable_LB = true;

    if (ts1.size() == 0 || ts2.size() == 0) {
        throw std::runtime_error("cannot calculate DTW on empty timeseries");
    }

    const size_t l1 = ts1.size();
    const size_t l2 = ts2.size();

    // LB Kim -- WARNING: PROBABLY INCORRECT FOR PHI_STAIRS
    double lb_kim = DistanceL2().dist_squared(ts1.data[0], ts2.data[0]);
    if (l1 > 1 || l2 > 1)
        lb_kim += DistanceL2().dist_squared(ts1.data[l1-1], ts2.data[l2-1]);

    if (enable_LB && lb_kim > early_abandon) {
        //DEBUG_NOENDL("\\");
        return INFTY;
    }

    const size_t size = (l1+1)*(l2+1);

    auto mat = vector<double>(size, INFTY);
    const size_t k = l2+1;

    auto minsrow = vector<double>();
    mat[0*k + 0] = 0.0;
    for (size_t i1 = 0; i1 < l1; i1++) {
        double min_of_row = INFTY;
        double prev_min_of_row = INFTY;

        for (size_t i2 = 0; i2 < l2; i2++) {
            // Current index in timeseries is at [i1] (resp [i2]).
            // Current index in the matrix is at [i1+1, i2+1].

            // Cost = distance between the two points.
            const double cost = DistanceL2().dist_squared(ts1.data[i1], ts2.data[i2]);

            // Compute all three paths coming to (i,j).
            const double choice1 = mat[ i1   *k + i2+1];
            const double choice2 = mat[(i1+1)*k + i2  ];
            const double choice3 = mat[ i1   *k + i2  ];

            const double minchoice = min(min(choice1, choice2), choice3);
            //DEBUG("  minchoice[" << i1 << "," << i2 << "] = " << minchoice);

            double tmp_result = cost + minchoice;
            mat[(i1+1)*k + i2+1] = tmp_result;

            if (tmp_result < min_of_row)
                min_of_row = tmp_result;
        }

        // Early abandoning DTW: the minimum of all values we calculated is
        // a lower bound.
        //DEBUG("min DTW on row " << i1 << " is = \t " << min_of_row);
        if (enable_LB && min_of_row > early_abandon) {
            //DEBUG_NOENDL("/");
            return INFTY;
        }
        assert(min_of_row <= prev_min_of_row); // assert decreasing
        minsrow.push_back(min_of_row);
        prev_min_of_row = min_of_row;
    }

    assert(l1*k + l2 + 1 == size);
    double result = mat[l1*k + l2];
    assert(lb_kim <= result);
    assert(std::all_of(minsrow.cbegin(), minsrow.cend(), [=](double x) { return x <= result; }));
    //DEBUG("final DTW ==> \t " << result);
    return result;
}
