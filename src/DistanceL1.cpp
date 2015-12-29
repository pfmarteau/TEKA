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
#include <cassert>
#include <cmath>
#include "DistanceL1.h"

/*
 * L1 norm.
 */
double DistanceL1::dist(const Point &p1, const Point &p2) {
    assert(p1.size() == p2.size());
    int n = p1.size();

    double sum = 0.0;
    for (int i = 0; i < n; i++) {
        double dp = p2[i] - p1[i];
        sum += fabs(dp);
    }

    return sum;
}


