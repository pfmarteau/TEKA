/****************************************************************************************************
Translated from java to C++ from source code available at https://github.com/fpetitjean/DBA
pfm 15/12/2015
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
******************************************************************************************************/
#include "DBA.h"
using namespace std;

	double distanceTo(double a, double b) {
		return (a - b) * (a - b);
	}
	double dist_squared(Point p1, Point p2) {
	    int n = p1.size();
	    double sum = 0.0;
	    for (int i = 0; i < n; i++) {
	        double dp = p2[i] - p1[i];
	        sum += dp*dp;
	    }
	    return sum;
	}

	double Min3(double a, double b, double c) {
		if (a < b) {
			if (a < c) {
				return a;
			} else {
				return c;
			}
		} else {
			if (b < c) {
				return b;
			} else {
				return c;
			}
		}
	}

	int ArgMin3(double a, double b, double c) {
		if (a < b) {
			if (a < c) {
				return 0;
			} else {
				return 2;
			}
		} else {
			if (b < c) {
				return 1;
			} else {
				return 2;
			}
		}
	}

	double barycenter(vector<double> tab) {
		if (tab.size() < 1) {
			cout << "RuntimeException('empty double tab');" << endl;
		}
		double sum = 0.0;
		sum = 0.0;
		for (unsigned int i=0; i< tab.size(); i++) {
			sum += (tab[i]);
		}
		return sum / tab.size();
	}
	
	std::vector<double> barycenterP(vector<Point> tab) {
		if (tab.size() < 1) {
			cout << "RuntimeException('empty double tab');" << endl;
		}
		vector<double> sum(tab[0].size());
		for(unsigned int i=0; i<tab[0].size(); i++ )
		    sum[i] = 0.0;
		for (unsigned int i=0; i< tab.size(); i++) {
			for(unsigned int j=0; j<tab[0].size(); j++)
			    sum[j] += (tab[i][j]);
		}
		for(unsigned int i=0; i<tab[0].size(); i++ )
			sum[i]=sum[i] / tab.size();
		return sum;
	}

DBA::DBA(){
	costMatrix.resize(MAX_SEQ_LENGTH, std::vector<double>(MAX_SEQ_LENGTH));
	pathMatrix.resize(MAX_SEQ_LENGTH, std::vector<int>(MAX_SEQ_LENGTH));
	optimalPathLength.resize(MAX_SEQ_LENGTH, std::vector<int>(MAX_SEQ_LENGTH));
}
/*
 * DTW Barycenter Averaging (DBA)
 * C: average sequence to update
 * sequence set of time series to average
 */
UniformTimeseries  DBA::computeAverage(UniformTimeseries C, std::vector<UniformTimeseries> sequences) {
	costMatrix.resize(MAX_SEQ_LENGTH);
	for(unsigned int i=0; i<costMatrix.size(); i++)
		costMatrix[i].resize(MAX_SEQ_LENGTH);

	std::vector<vector<Point>> tupleAssociation;
	tupleAssociation.resize(C.size());

	unsigned int nbTuplesAverageSeq, i, j, indiceRes;
	double res = 0.0;
	unsigned int centerLength = C.size();
	unsigned int seqLength;

	UniformTimeseries T;	
	for(unsigned int n=0; n<sequences.size(); n++){
		T=sequences[n];
		seqLength = T.size();
		costMatrix[0][0] = distanceTo(C.data[0][0], T.data[0][0]);
		pathMatrix[0][0] = NIL;
		optimalPathLength[0][0] = 0;
		for (unsigned int i = 1; i < centerLength; i++) {
			costMatrix[i][0] = costMatrix[i - 1][0] + dist_squared(C.data[i], T.data[0]);
			pathMatrix[i][0] = UP;
			optimalPathLength[i][0] = i;
		}
		for (unsigned int j = 1; j < seqLength; j++) {
			costMatrix[0][j] = costMatrix[0][j - 1] + dist_squared(T.data[j], C.data[0]);
			pathMatrix[0][j] = LEFT;
			optimalPathLength[0][j] = j;
		}

			for (unsigned int i = 1; i < centerLength; i++) {
				for (j = 1; j < seqLength; j++) {
					indiceRes = ArgMin3(costMatrix[i - 1][j - 1], costMatrix[i][j - 1], costMatrix[i - 1][j]);
					pathMatrix[i][j] = indiceRes;
					switch (indiceRes) {
						case DIAGONAL:
							res = costMatrix[i - 1][j - 1];
							optimalPathLength[i][j] = optimalPathLength[i - 1][j - 1] + 1;
							break;
						case LEFT:
							res = costMatrix[i][j - 1];
							optimalPathLength[i][j] = optimalPathLength[i][j - 1] + 1;
							break;
						case UP:
							res = costMatrix[i - 1][j];
							optimalPathLength[i][j] = optimalPathLength[i - 1][j] + 1;
							break;
					}
					costMatrix[i][j] = res + dist_squared(C.data[i], T.data[j]);
				}
			}
			nbTuplesAverageSeq = optimalPathLength[centerLength - 1][seqLength - 1] + 1;

			i = centerLength - 1;
			j = seqLength - 1;

			for (int t = nbTuplesAverageSeq - 1; t >= 0; t--) {
				tupleAssociation[i].push_back(T.data[j]);
				switch (pathMatrix[i][j]) {
					case DIAGONAL:
						i = i - 1;
						j = j - 1;
						//cout << "diag";
						break;
					case LEFT:
						j = j - 1;
						//cout << "left";
						break;
					case UP:
						i = i - 1;
						//cout << "up";
						break;
				}

			}
		}
		for (unsigned int t = 0; t < centerLength; t++) {
			C.data[t] = barycenterP(tupleAssociation[t]);
		}
	return C;
	}

