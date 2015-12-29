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
#include "GramMatrix.h"

using namespace std;

/*
 * return the min of a matrix (vector of vectors of doubles).
 */
double matmin(std::vector<std::vector<double>> mat){
double min = 1e300;
  for ( unsigned int  i = 0 ; i < mat.size(); i++ ){
     for(unsigned int  j=0; j < mat[i].size(); j++){
	if(min>mat[i][j])
		min=mat[i][j];
     }
  }
return min;
}

/*
 * return the max of a matrix (vector of vectors of doubles).
 */
double matmax(std::vector<std::vector<double>> mat){
double max = -1e300;
  for ( unsigned int  i = 0 ; i < mat.size(); i++ ){
     for(unsigned int  j=0; j < mat[i].size(); j++){
	if(max<mat[i][j])
		max=mat[i][j];
     }
  }
return max;
}

/*
 * transform a matrix (vector of vectors of doubles) in a pointer of pointer of double.
 */
double **vecmat2doubleptr(matrix<double> M){
	double **out=(double **)calloc(M.nr(), sizeof(double *));
	for ( unsigned int i = 0 ; i < M.nr(); i++ ){
		out[i]=(double *)calloc(M.nc(), sizeof(double ));
		for(unsigned int j=0; j < M.nc(); j++){
			out[i][j]=M(i,j);
		}
	}
return out;
}

/*
 * compute the DTW gram matrix from a vector of uniformly sampled time series.
 */
matrix<double> ComputeDTWGramMatrix(std::vector<UniformTimeseries> ds){
	matrix<double> distMat(ds.size(), ds.size());
	for ( unsigned int i = 0 ; i < ds.size(); i++ )
		for(unsigned int j=0; j < ds.size(); j++)
			distMat(i,j)=0;
	for ( unsigned int i = 0 ; i < ds.size(); i++ ){
			for(unsigned int j=0; j < ds.size(); j++){
						distMat(i,j)=DTW().dist_squared(ds[i], ds[j]);
			}
	}
	return (distMat);
}

/*
 * compute the DTW gram matrix, for the class label <nclass>, from a vector of uniformly sampled time series.
 */
matrix<double> ComputeDTWGramMatrix(std::vector<UniformTimeseries> ds, int nclass){
	matrix<double>  distMat(ds.size(), ds.size());
	for ( unsigned int i = 0 ; i < ds.size(); i++ )
		for(unsigned int j=0; j < ds.size(); j++)
			distMat(i,j)=0;
	for (int n=1; n<= nclass; n++){
		//cout << " class= " << n << endl;
		for ( unsigned int i = 0 ; i < ds.size(); i++ ){
			if(ds[i].label==n){
				for(unsigned int j=0; j < ds.size(); j++){
					if(ds[j].label==n){
						distMat(i,j)=DTW().dist_squared(ds[i], ds[j]);
					}
				}
			}
		}
	}
	return (distMat);
}

/*
 * compute the KDTW gram matrix from a vector of uniformly sampled time series with <sigma> value.
 */
matrix<double> ComputeKDTWGramMatrix(std::vector<UniformTimeseries> ds, double sigma){
	matrix<double> GM(ds.size(), ds.size());
for (unsigned int i = 0 ; i < ds.size(); i++ ){
   for(unsigned int j=i; j < ds.size(); j++){
	   GM(i,j)=KDTW().kernel(ds[i], ds[j], sigma);
	   GM(j,i)=GM(i,j);
   }
}
return (GM);
}

/*
 * compute the KDTW gram matrix, for the class label <nclass>, from a vector of uniformly sampled time series with <sigma> value.
 */
matrix<double>  ComputeKDTWGramMatrix(std::vector<UniformTimeseries> ds, double sigma, int nclass){
	matrix<double>  GM(ds.size(), ds.size());
	for (unsigned int i = 0 ; i < ds.size(); i++ )
		for(unsigned int j=0; j < ds.size(); j++)
			GM(i,j)=-DBL_MAX;
	for (int n=1; n<= nclass; n++){
		for (unsigned int i = 0 ; i < ds.size(); i++ ){
			if(ds[i].label==n){
				for(unsigned int j=i; j < ds.size(); j++){
					if(ds[j].label==n){
						GM(i,j)=KDTW().kernel(ds[i], ds[j], sigma);
						GM(j,i)=GM(i,j);
					}
				}
			}
		}
	}
	return (GM);
}

/*
 * compute the KDTW distance gram matrix, for the class label <nclass>, from a vector of uniformly sampled time series with <sigma> value.
 */
matrix<double> ComputeKDistanceMatrix(std::vector<UniformTimeseries> ds, double sigma, std::vector<double> label,
		int nclass){
	matrix<double> distMat(ds.size(), ds.size());
	for (unsigned int i = 0 ; i < ds.size(); i++)
		for(unsigned int j=0; j < ds.size(); j++)
			distMat(i,j)=DBL_MAX;
	for (int n=1; n<= nclass; n++){
		for (unsigned  int i = 0 ; i < ds.size(); i++ ){
			if(label[i]==n){
				for(unsigned int  j=i; j < ds.size(); j++){
					if(label[j]==n){
						distMat(i,j)=KDTW().kdistance(ds[i], ds[j], sigma);
						distMat(j,i)=distMat(i,j);
					}
				}
			}
		}
	}
	return (distMat);
}

/*
 * replace each element of the input matrix  by its logarithm (neperian).
 */
std::vector<std::vector<double>> ComputeLogMatrix(std::vector<std::vector<double>> GM){
unsigned int sz=GM.size();
for (unsigned int  i = 0 ; i < sz; i++ ){
   for(unsigned int  j=i; j < sz; j++){
	   GM[i][j]=log(GM[i][j]+1e-300);
      GM[j][i]=GM[i][j];
   }
}
return (GM);
}

/*
 * replace each element of the input matrix  by its 'normalized' logarithm (neperian).
 */
std::vector<std::vector<double>> ComputeLogNormalizedMatrix(std::vector<std::vector<double>> GM){
unsigned int  sz=GM.size();
for (unsigned int  i = 0 ; i < sz; i++ ){
   for(unsigned int  j=i; j < sz; j++){
	   GM[i][j]=log(GM[i][j]+1e-300);
	   GM[j][i]=GM[i][j];
   }
}
double logmax=matmax(GM);
double logmin=matmin(GM);
for (unsigned int  i = 0 ; i < sz; i++ ){
   for(unsigned int  j=i; j < sz; j++){
	   GM[i][j]=exp((GM[i][j]-logmin)/(logmax-logmin));
	   GM[j][i]=GM[i][j];
   }
}
return (GM);
}

/*
 * replace each element of the input matrix M by its 't' power, with t=1.0/(log(max(M)-log(min(M)).
 */
std::vector<std::vector<double>> ComputePowerNormalizedMatrix(std::vector<std::vector<double>> GM, double *t){
unsigned int  sz = GM.size();
double logmin = 1e300;
double logmax = -1e300, x;
for ( unsigned int  i = 0 ; i < sz; i++ ){
   for(unsigned int  j=i; j < sz; j++){
      x=log(GM[i][j]+1e-300);
      if (logmin>x)
      	  logmin=x;
      if (logmax<x)
      	  logmax=x;
   }
}
*t=1.0/(logmax-logmin);
for ( unsigned int  i = 0 ; i < sz; i++ ){
   for(unsigned int  j=i; j < sz; j++){
	   GM[i][j]=pow(GM[i][j],*t);
	   GM[j][i]=GM[i][j];
   }
}
return (GM);
}





