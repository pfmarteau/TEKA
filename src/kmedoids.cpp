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
#include "kmedoids.h"
using namespace std;

/*
 * Returns the index of the Medoid (the time series that is at a minimum cumulated distance of the other elements)
 * matrix<double>  distMat: matrix of pairwise distances
 * int *clusterid: array of labelled cluster ids
 * int nc: number of distinct labels (categories)
 *
 */
int getMedoid_min(matrix<double>  distMat, int *clusterid, int nc){
	int ii=0;
	double x;
	double dtwv[distMat.nr()];
	for(unsigned int i=0; i<distMat.nr(); i++){
		dtwv[i]=DBL_MAX;
		if(clusterid[i]==nc){
		   x=0; ii=0;
		   for(unsigned int j=0; j<distMat.nc(); j++){
			   if(clusterid[j]==nc){
				   x+=distMat(i,j);
				   ii++;
			   }
		   }
		   if(ii>0)
		      dtwv[i]=x/(double)ii;
		   else
			   dtwv[i]=DBL_MAX;
		}
	}
	double min=DBL_MAX;
	int imin=0;
	for(unsigned int i=0; i<distMat.nr(); i++){
		if(min>dtwv[i] && clusterid[i]==nc){
			imin=i;
			min=dtwv[i];
		}
	}
return imin;
}

/*
 * Returns the index of the Medoid (the time series that is at a maximum cumulated kernel value of the other elements)
 * matrix<double>  GramMat: matrix of pairwise kernel values (Gram matrix)
 * int *clusterid: array of labelled cluster ids
 * int nc: number of distinct labels (categories)
 */
int getMedoid_max(matrix<double> GramMat, int *clusterid, int nc){
	int ii=0;
	double x;
	double kdtwv[GramMat.nr()];
	for(unsigned int i=0; i<GramMat.nr(); i++){
		kdtwv[i]=-DBL_MAX;
		if(clusterid[i]==nc){
		   x=0; ii=0;
		   for(unsigned int j=0; j<GramMat.nc(); j++){
			   if(clusterid[j]==nc){
				   x+=GramMat(i,j);
				   ii++;
			   }
		   }
		   if(ii>0)
		      kdtwv[i]=x/(double)ii;
		   else
			   kdtwv[i]=-DBL_MAX;
		}
	}
	double max=DBL_MAX;
	int imax=0;
	for(unsigned int i=0; i<GramMat.nr(); i++){
		if(max<kdtwv[i] && clusterid[i]==nc){
			imax=i;
			max=kdtwv[i];
		}
	}
return imax;
}

/*
 * Returns the set of DTW medoids corresponding to each of the <nclust> clusters for each of the <nclass> categories.
 * std::vector<UniformTimeseries> dsTrain: the set of labelled time series
 * int ncat: number of categories (labels)
 * int nclust: number of cluster per category.
 */
std::vector<UniformTimeseries>  kmedoids_dtw(std::vector<UniformTimeseries> dsTrain, int ncat, int nclust){
std::vector<UniformTimeseries>  dsCat, ds;
  for(int nc=1; nc<=ncat; nc++){
	cout << "n class =" << nc << endl;
	dsCat=UniformTimeseries::extractTSfromLabelDS(dsTrain, nc);
	matrix<double> gramDTW=ComputeDTWGramMatrix(dsCat, nc);
	double **distDTW=vecmat2doubleptr(gramDTW);

	int npass=30, ifound, imedoid;
	int *clusterid = (int *)calloc(gramDTW.nr(), sizeof(int));
	double *errors = (double *)calloc(nclust, sizeof(double));
	kmedoids (nclust, gramDTW.nr(), distDTW, npass, clusterid, errors, &ifound);

	for (int k=0; k<nclust; k++){
		cout << "cluster_" << k <<endl;
		imedoid=getMedoid_min(gramDTW, clusterid, k);
		ds.push_back(dsCat[imedoid]);
	}
  }
return ds;
}

/*
 * Returns the set of KDTW medoids corresponding to the <nclust> clusters for each of the <nclass> categories.
 * std::vector<UniformTimeseries> dsTrain: the set of labelled time series
 * int ncat: number of categories (labels)
 * int nclust: number of cluster per category.
 * double sigma: parameter of the KDTW kernel
 */
std::vector<UniformTimeseries>  kmedoids_kdtw(std::vector<UniformTimeseries> dsTrain, int ncat, int nclust, double sigma){
std::vector<UniformTimeseries>  dsCat, ds;
  for(int nc=1; nc<=ncat; nc++){
	cout << "n class =" << nc << endl;
	dsCat=UniformTimeseries::extractTSfromLabelDS(dsTrain, nc);
	matrix<double> gramKDTW=ComputeKDTWGramMatrix(dsCat, sigma);

	for(unsigned int i=0; i<gramKDTW.nr(); i++)
		for(unsigned int j=0; j<gramKDTW.nc(); j++)
			gramKDTW(i,j) = - log(gramKDTW(i,j)+DBL_MIN);
	double **distKDTW=vecmat2doubleptr(gramKDTW);
	
	int npass=30, ifound, imedoid;
	int *clusterid = (int *)calloc(gramKDTW.nr(), sizeof(int));
	double *errors = (double *)calloc(nclust, sizeof(double));
	kmedoids(nclust, gramKDTW.nr(), distKDTW, npass, clusterid, errors, &ifound);

	for (int k=0; k<nclust; k++){
		cout << "cluster_" << k <<endl;
		imedoid=getMedoid_min(gramKDTW, clusterid, k);
		ds.push_back(dsCat[imedoid]);
	}
  }
return ds;
}
