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
#include "kcentroids.h"
using namespace std;

/*
 * Returns the set of centroids for pairs of time series taken two by two in dataset <ds>.
 * std::vector<UniformTimeseries> ds: a set of time series
 * double sigma: the parameter ok KDTW kernel.
 * int lab: the category associated to each of the time series in <ds>
 */
std::vector<UniformTimeseries> kadd2by2div2(std::vector<UniformTimeseries> ds, double sigma, int lab){
	std::vector<UniformTimeseries> ds2;
	UniformTimeseries med;
	matrix<double> ama;
	int ii=0;
	bool comp=false;
	int LL=ds.size();
	if(ds.size()==1){
		ds2.push_back(ds[0]);
	}
	else{
		int nbMax=LL/2;
		ii=0;
		int i=0;
		for (i=0; (i<LL&&ii<nbMax); i++){
			if(!comp){
				comp=true;
				med=ds[i];
			} else{
				ama=KDTW().kernel_AMA(med, ds[i], sigma);
				med = KDTW().add_kernelT(med, ds[i], ama, 1, 2);
				med.label=lab;
				ds2.push_back(med);
				cout <<  "norm= " << med.norml1() << endl;
				ii++;
				comp=false;
			}
		}

		while(i<LL&&ii<nbMax){
			ds2.push_back(ds[i]);
			i++;
		}
	}
	return ds2;
}

/*
 * Returns the set of centroids evaluated for all cluster in each category, one centroid per cluster (hence nclust.ncat centroids in total)
 * std::vector<UniformTimeseries> dsTrain: the set of labelled time series
 * int ncat: number of categories (labels)
 * int nclust: number of cluster per category.
 */
std::vector<UniformTimeseries>  kcentroid_kdtw_pwa(std::vector<UniformTimeseries> dsTrain, int ncat, int nclust, double sigma){
	std::vector<UniformTimeseries> centroid, dsCentroid;
	std::vector<UniformTimeseries> dsCat, ds, ds2;
	bool change=true;
	for(int nc=1; nc<=ncat; nc++){
		cout << "nc=" << nc << endl;
		dsCat=UniformTimeseries::extractTSfromLabelDS(dsTrain, nc);
		matrix<double> gramKDTW=ComputeKDTWGramMatrix(dsCat, sigma);
		for(unsigned int i=0; i<gramKDTW.nr(); i++)
			for(unsigned int j=0; j<gramKDTW.nc(); j++)
				gramKDTW(i,j) = -log(gramKDTW(i,j)+1e-300); // transform the gram matrix in a dissimilarity (distance like) matrix.
		double **distKDTW=vecmat2doubleptr(gramKDTW);

		int npass=30, ifound, N;
		int *clusterid = (int *)calloc(gramKDTW.nr(), sizeof(int));
		double *errors = (double *)calloc(nclust, sizeof(double));
		unsigned int n=nclust;
		if(n>dsCat.size())
			n=dsCat.size();
		if(nclust>1)
			kmedoids(n, gramKDTW.nr(), distKDTW, npass, clusterid, errors, &ifound);
		else
			for(unsigned int i=0; i<gramKDTW.nr(); i++)
				clusterid[i]=0;
		int ii=0;
		change=true;
		double inertia=0;
		while (ii< 10 && change){
			if(centroid.size()>0)
				centroid.clear();
			for (int k=0; k<nclust; k++){
				cout << "cluster_" << k << endl;
				ds.clear();
				for(unsigned int i=0; i<gramKDTW.nr(); i++){
					if(clusterid[i]==k){
						ds.push_back(dsCat[i]);
						cout <<"  " << i ;//<< " clusterid_" << clusterid[i];
					}
				}
				cout << endl;
				N=ds.size();
				ds2=ds;
				while(N>1){
					cout << "N=" << N << endl;
					ds2=kadd2by2div2(ds2, sigma, nc);
					N=ds2.size();
				}
				if(ds2.size()==0){
					cout << "kcentroid.cpp erreur ds2 line 307" << endl;
				}
				else{
					centroid.push_back(ds2[0]);
				}
			}
			change=false;
			cout << "Change=" << change << " with inertia=" << inertia <<  endl;
			ii++;
		}
		for(unsigned int i=0; i<centroid.size(); i++){
			centroid[i].label=nc;
			dsCentroid.push_back(centroid[i]);
		}
		dsCat.clear();
		free(clusterid);
	}
	return dsCentroid;
}

/*
 * Returns the DBA average (centroid estimate) of a set of times series <ds> containing a single label.
 * UniformTimeseries averageSequence: the current (initial) average time series.
 * std::vector<UniformTimeseries> ds: the set of time series to average.
 */
UniformTimeseries getDTWCentroid(UniformTimeseries averageSequence, std::vector<UniformTimeseries> ds){
	DBA *dba=new(DBA);
	averageSequence=dba->computeAverage(averageSequence, ds);
	averageSequence.label=ds[0].label;
	delete dba;
	return averageSequence;
}

/*
 * Returns the KDBA average (centroid estimate) of a set of times series <ds> containing a single label.
 * UniformTimeseries average: the current (initial) average.
 * std::vector<UniformTimeseries> ds: the set of time series to average.
 * double sigma: the parameter of the KDTW kernel.
 */
UniformTimeseries getKDTWCentroid(UniformTimeseries average, std::vector<UniformTimeseries> ds, double sigma){
	average=KDTW().add_kernelDBA(average, ds, sigma, ds[0].label);
	average.label=ds[0].label;
	return average;
}

/*
 * Returns the "inertia" (averaged sum of squared distances) of time series <ctrd> relatively to dataset <ds>.
 * UniformTimeseries ctrd: an average estimate.
 * std::vector<UniformTimeseries> ds: the set of time series to average.
 */
double getInertiaDTW(UniformTimeseries ctrd, std::vector<UniformTimeseries> ds){
	double inertia=0;
	for(unsigned int i=0; i< ds.size(); i++){
		inertia+=DTW().dist_squared(ctrd, ds[i]);
	}
	return(inertia/ds.size());
}


std::vector<UniformTimeseries>  kcentroid_idba(std::vector<UniformTimeseries> dsTrain, int nclass, int nclust){
	std::vector<UniformTimeseries> centroid, dsCentroid;
	std::vector<UniformTimeseries> dsCat, ds;
	UniformTimeseries ctrd, ctrdp;
	bool change=true;
	int imed;
	for(int nc=1; nc<=nclass; nc++){
		cout << "nc=" << nc << endl;
		dsCat=UniformTimeseries::extractTSfromLabelDS(dsTrain, nc);
		matrix<double> gramDTW=ComputeDTWGramMatrix(dsCat);

		double **distDTW=vecmat2doubleptr(gramDTW);

		int npass=30, ifound;
		int *clusterid = (int *)calloc(gramDTW.nr(), sizeof(int));
		double *errors = (double *)calloc(nclust, sizeof(double));
		kmedoids (nclust, gramDTW.nr(), distDTW, npass, clusterid, errors, &ifound);
		int ii=0;
		change=true;
		double inertia=0;
		double inertiap=DBL_MAX;
		while (ii< 10 && change){
			if(centroid.size()>0)
				centroid.clear();
			for (int k=0; k<nclust; k++){
				cout << "cluster_" << k <<endl;
				ds.clear();
				for(unsigned int i=0; i<gramDTW.nr(); i++){
					if(clusterid[i]==k){
						ds.push_back(dsCat[i]);
						cout <<"  " << i ;//<< " clusterid_" << clusterid[i];
					}
				}
				cout << endl;
				imed=getMedoid_min(gramDTW, clusterid, k);
				ctrd=getDTWCentroid(dsCat[imed],ds);
				inertia=getInertiaDTW(ctrd, ds);
				int ii=0;
				cout << "inertia=";
				while(ii<10 && inertia<inertiap){
					cout << inertia << ",";
					ctrdp=ctrd;
					inertiap=inertia;
					ctrd=getDTWCentroid(ctrdp,ds);
					inertia=getInertiaDTW(ctrd, ds);
					ii++;
				}
				cout << endl;
				ctrdp.label=nc;
				centroid.push_back(ctrdp);
			}
			change=false;
			cout << "Change=" << change << " with inertia=" << inertia <<  endl;
			ii++;
		}
		for(unsigned int i=0; i<centroid.size(); i++)
			dsCentroid.push_back(centroid[i]);
		dsCat.clear();
	}
	return dsCentroid;
}

/*
 * Returns the "inertia" (averaged sum of kernel values KDTW(ctrd,.)) of time series <ctrd> relatively to dataset <ds>.
 * UniformTimeseries ctrd: an average estimate.
 * std::vector<UniformTimeseries> ds: the set of time series to average.
 * double sigma: the parameter of the KDTW kernel.
 */
double getInertiaKDTW(UniformTimeseries ctrd, std::vector<UniformTimeseries> ds, double sigma){
	double inertia=0;
	double nrm=sqrt(KDTW().kernel(ctrd, ctrd, sigma));
	for(unsigned int i=0; i< ds.size(); i++){
		inertia+=KDTW().kernel(ctrd, ds[i], sigma)/nrm;
	}
	return(inertia/ds.size());
}

std::vector<UniformTimeseries>  kcentroid_ikdba(std::vector<UniformTimeseries> dsTrain, int nclass, int nclust, double sigma){
	std::vector<UniformTimeseries> centroid, dsCentroid;
	std::vector<UniformTimeseries> dsCat, ds;
	UniformTimeseries ctrd, ctrdp;
	bool change=true;
	int imed;
	for(int nc=1; nc<=nclass; nc++){
		cout << "nc=" << nc << endl;
		dsCat=UniformTimeseries::extractTSfromLabelDS(dsTrain, nc);
		matrix<double> gramKDTW=ComputeKDTWGramMatrix(dsCat, sigma);
		for(unsigned int i=0; i<gramKDTW.nr(); i++)
			for(unsigned int j=0; j<gramKDTW.nc(); j++)
				gramKDTW(i,j) = -log(gramKDTW(i,j)+1e-300);
		double **distKDTW=vecmat2doubleptr(gramKDTW);

		int npass=30, ifound;
		int *clusterid = (int *)calloc(gramKDTW.nr(), sizeof(int));
		double *errors = (double *)calloc(nclust, sizeof(double));
		kmedoids (nclust, gramKDTW.nr(), distKDTW, npass, clusterid, errors, &ifound);
		int ii=0;
		change=true;
		double inertia=0;
		double inertiap=-DBL_MAX;
		while (ii< 10 && change){
			if(centroid.size()>0)
				centroid.clear();
			for (int k=0; k<nclust; k++){
				cout << "cluster_" << k <<endl;
				ds.clear();
				for(unsigned int i=0; i<gramKDTW.nr(); i++){
					if(clusterid[i]==k){
						ds.push_back(dsCat[i]);
						cout <<"  " << i ;//<< " clusterid_" << clusterid[i];
					}
				}
				cout << endl;
				imed=getMedoid_min(gramKDTW, clusterid, k);
				ctrd=getKDTWCentroid(dsCat[imed],ds, sigma);
				inertia=getInertiaKDTW(ctrd, ds, sigma);
				int ii=0;
				cout << "inertia=";
				cout << inertia << ",";
				while(ii<10 && inertia>inertiap){
					ctrdp=ctrd;
					inertiap=inertia;
					ctrd=getKDTWCentroid(ctrdp,ds, sigma);
					inertia=getInertiaKDTW(ctrd, ds,sigma);
					cout << inertia << ",";
					ii++;
				}
				cout << endl;
				ctrdp.label=nc;
				centroid.push_back(ctrdp);
			}
			change=false;
			cout << "Change=" << change << " with inertia=" << inertiap <<  endl;
			ii++;
		}
		for(unsigned int i=0; i<centroid.size(); i++)
			dsCentroid.push_back(centroid[i]);
		dsCat.clear();
	}
	return dsCentroid;
}

