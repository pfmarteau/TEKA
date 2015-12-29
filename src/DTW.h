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
#ifndef DTW_H
#define DTW_H

#include <limits>
#include <cmath>
#include <cfloat>
#include "UniformTimeseries.h"


/*
 * DTW is an elastic distance.
 * The .cpp file just contains the main implementation of the algorithm. Other
 * functions (squared, hint...) are just derived versions of this.
 *
 * TODO: "hint" will be useful in the future, for early abandoned computations.
 * But it is not yet implemented.
 */
class DTW {
public:
    double dist(const UniformTimeseries &ts1,
                const UniformTimeseries &ts2) {
        return dist_hint(ts1, ts2, std::numeric_limits<double>::infinity());
    }

    double dist_hint(const UniformTimeseries &ts1,
                     const UniformTimeseries &ts2,
                     double early_abandon);

    double dist_squared(const UniformTimeseries &ts1,
                        const UniformTimeseries &ts2)
    { return dist_squared_hint(ts1, ts2, std::numeric_limits<double>::infinity()); }

    double dist_squared_hint(const UniformTimeseries &ts1,
                             const UniformTimeseries &ts2,
                             double early_abandon)
    {
        auto d = dist_hint(ts1, ts2, std::sqrt(early_abandon));
        return d*d;
    }

};




#endif
