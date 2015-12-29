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
#include <stdexcept>
#include <limits>
#include <iostream>
#include <fstream>
#include "KDTW.h"
#include "DistanceL2.h"

using namespace std;

static const double INFTY = numeric_limits<double>::infinity();

// Singleton for use throughout the whole program.
std::shared_ptr<KDTW> KDTW_DIST = std::shared_ptr<KDTW>(new KDTW());

/*
 * kernel based distance (similar to the Euclidean distance) between two time series.
 */
double KDTW::kdistance(const UniformTimeseries &ts1,
                            const UniformTimeseries &ts2,
                            double sigma){
	return sqrt(kernel(ts1, ts1, sigma) + (kernel(ts2, ts2, sigma) -2.0*kernel(ts1, ts2, sigma)));
}

/*
 * local distance kernel matrix resulting from the aligment of two time series.
 */
matrix<double>  KDTW::kdistance_mat(const UniformTimeseries &ts1,
                            const UniformTimeseries &ts2,
                            double sigma){
	return kernel_mat(ts1, ts1, sigma) + kernel_mat(ts2, ts2, sigma) - kernel_mat(ts1, ts2, sigma)*2.0;
}

/*
 * kernel (KDTW) value for a pair of time series, as defined in
 * Marteau, P.F., Gibet, S., On Recursive Edit Distance Kernels with Application to Time Series Classification,
 * IEEE Trans Neural Netw Learn Syst 2015 Jun 25;26(6):1121-33.
 */

double KDTW::kernel(const UniformTimeseries &ts1,
                            const UniformTimeseries &ts2,
                            double sigma)
{
    const int l1 = ts1.size();
    const int l2 = ts2.size();
    const int dim=ts1.data[0].size();
    sigma=sigma*dim;

    if (l1 == 0 || l2 == 0)
        throw std::runtime_error("time series cannot be empty");

    int lmax=l1;
    int lmin=l2;
    if(lmax<l2){
	lmax=l2;
	lmin=l1;
	}
    size_t size = lmax;
    auto mat_Rx = vector<double>(size, 0.0);
    vector<vector<double> > mat_kxy(l1+1, std::vector<double>(l2+1));
    vector<vector<double> > mat_kxx(l1+1, std::vector<double>(l2+1));

    for (int i = 0; i < lmax; i++) {
	if(i<lmin){
		mat_kxy[i][0]=0; mat_kxy[0][i]=0;
		mat_kxx[i][0]=0; mat_kxx[0][i]=0;
	    const Point &p1 = ts1.data[i];
        const Point &p2 = ts2.data[i];
		mat_Rx[i]=exp(-DistanceL2().dist_squared(p1, p2)/sigma);
		}
	else if(i<l1){
		mat_kxy[i][0]=0;
		mat_kxx[i][0]=0;
	    const Point &p1 = ts1.data[i];
        const Point &p2 = ts2.data[l2-1];
		mat_Rx[i]=exp(-DistanceL2().dist_squared(p1, p2)/sigma);
		}
	else if (i<l2){
		mat_kxy[0][i]=0;
		mat_kxx[0][i]=0;
	    const Point &p1 = ts1.data[l1-1];
        const Point &p2 = ts2.data[i];
		mat_Rx[i]=exp(-DistanceL2().dist_squared(p1, p2)/sigma);
		}
    }
	
    double factor=1.0/3.0;
    mat_kxy[0][0] = 1.0;
    mat_kxx[0][0] = 1.0;

    //---------
    const Point &p2 = ts2.data[0];
    for (int i1 = 1; i1 <= l1; i1++) {
    	const Point &p1 = ts1.data[i1-1];
    	mat_kxy[i1][0]=mat_kxy[i1-1][0]*exp(-DistanceL2().dist_squared(p1, p2)/sigma);
    	mat_kxx[i1][0]=mat_kxx[i1-1][0]*mat_Rx[i1-1];
    }
    const Point &p1 = ts1.data[0];
    for (int i2 = 1; i2 <= l2; i2++) {
    	const Point &p2 = ts2.data[i2-1];
    	mat_kxy[0][i2]=mat_kxy[0][i2-1]*exp(-DistanceL2().dist_squared(p1, p2)/sigma);
       	mat_kxx[0][i2]=mat_kxx[0][i2-1]*mat_Rx[i2-1];
    }
    //---------

    for (int i1 = 1; i1 <= l1; i1++) {
        for (int i2 = 1; i2 <= l2; i2++) {
            const Point &p1 = ts1.data[i1-1];
            const Point &p2 = ts2.data[i2-1];
            const double cost = factor*exp(-DistanceL2().dist_squared(p1, p2)/sigma);
            const double choicexx1 = mat_kxx[i1-1][i2]*mat_Rx[i1-1];
            const double choicexx2 = mat_kxx[i1][i2-1]*mat_Rx[i2-1];
            mat_kxy[i1][i2] = cost*(mat_kxy[i1-1][i2] + mat_kxy[i1][i2-1] + mat_kxy[i1-1][i2-1]);
	    if(i1==i2){
	    	const double choicexx3 = cost*mat_kxx[i1-1][i2-1];
	    	mat_kxx[i1][i2] = factor*(choicexx1 + choicexx2) + choicexx3;
	        }
		else
			mat_kxx[i1][i2] = factor*(choicexx1 + choicexx2);
        }

    }

    double result = mat_kxy[l1][l2]+mat_kxx[l1][l2];
    return result;
}

/*
 * local kernel matrix resulting from the aligment of two time series.
 */
matrix<double> KDTW::kernel_mat(const UniformTimeseries &ts1,
                            const UniformTimeseries &ts2,
                            double sigma)
{
    const int l1 = ts1.size();
    const int l2 = ts2.size();
    const int dim=ts1.data[0].size();
    sigma=sigma*dim;

    if (l1 == 0 || l2 == 0)
        throw std::runtime_error("timeseries cannot be empty");

    int lmax=l1;
    int lmin=l2;
    if(lmax<l2){
	lmax=l2;
	lmin=l1;
	}
    size_t size = lmax;
    auto mat_Rx = vector<double>(size, 0.0);
    matrix<double> mat_kxy(l1+1, l2+1);
    matrix<double> mat_kxx(l1+1, l2+1);
    matrix<double> mat_out(l1+1, l2+1);

    for (int i = 0; i < lmax; i++) {
	if(i<lmin){
		mat_kxy(i,0)=0; mat_kxy(0,i)=0;
		mat_kxx(i,0)=0; mat_kxx(0,i)=0;
	    const Point &p1 = ts1.data[i];
        const Point &p2 = ts2.data[i];
		mat_Rx[i]=exp(-DistanceL2().dist_squared(p1, p2)/sigma);
		}
	else if(l1==lmax){
		mat_kxy(i,0)=0;
		mat_kxx(i,0)=0;
	    const Point &p1 = ts1.data[i];
        const Point &p2 = ts2.data[l2-1];
		mat_Rx[i]=exp(-DistanceL2().dist_squared(p1, p2)/sigma);
		}
	else {
		mat_kxy(0,i)=0;
		mat_kxx(0,i)=0;
	    const Point &p1 = ts1.data[l1-1];
        const Point &p2 = ts2.data[i];
		mat_Rx[i]=exp(-DistanceL2().dist_squared(p1, p2)/sigma);
		}
    }

    double factor=1.0/3.0;
    mat_kxy(0,0) = 1.0;
    mat_kxx(0,0) = 1.0;

    //--------
    const Point &p2 = ts2.data[0];
    for (int i1 = 1; i1 <= l1; i1++) {
    	const Point &p1 = ts1.data[i1-1];
    	mat_kxy(i1,0)=mat_kxy(i1-1,0)*exp(-DistanceL2().dist_squared(p1, p2)/sigma);
    	mat_kxx(i1,0)=mat_kxx(i1-1,0)*mat_Rx[i1-1];
    }
    const Point &p1 = ts1.data[0];
    for (int i2 = 1; i2 <= l2; i2++) {
    	const Point &p2 = ts2.data[i2-1];
    	mat_kxy(0,i2)=mat_kxy(0,i2-1)*exp(-DistanceL2().dist_squared(p1, p2)/sigma);
    	mat_kxx(0,i2)=mat_kxx(0,i2-1)*mat_Rx[i2-1];
    }
    //--------

    for (int i1 = 1; i1 <= l1; i1++) {
        for (int i2 = 1; i2 <= l2; i2++) {
            const Point &p1 = ts1.data[i1-1];
            const Point &p2 = ts2.data[i2-1];
            const double cost = factor*exp(-DistanceL2().dist_squared(p1, p2)/sigma);

            mat_kxy(i1,i2) = cost*(mat_kxy(i1-1,i2) + mat_kxy(i1,i2-1) + mat_kxy(i1-1,i2-1));

            mat_kxx(i1,i2) = factor*(mat_kxx(i1-1,i2)*mat_Rx[i1] + mat_kxx(i1,i2-1)*mat_Rx[i2]);
	        if(i1==i2){
		       mat_kxx(i1,i2) += cost*mat_kxx(i1-1,i2-1);
	        }
	    mat_out(i1,i2)=mat_kxy(i1,i2);
        }
    }

    return mat_out;

}

double _max(double a, double b){
	if(a<b)
		return b;
	else
		return a;
}

/*
 * AMA (Alignment kernel Matrix Average) between two time series.
 */
matrix<double> KDTW::kernel_AMA(UniformTimeseries ts1, UniformTimeseries ts2, double sigma)
{
	UniformTimeseries tsout;
	matrix<double> matk, matkr;
	matk = KDTW().kernel_mat(ts1, ts2, sigma);
	UniformTimeseries ts1r = ts1.reverse();
	UniformTimeseries ts2r = ts2.reverse();
	matkr = KDTW().kernel_mat(ts1r, ts2r, sigma);
	int l1=matk.nr();
	int l2=matk.nc();
	for(int i1=1; i1<l1; i1++){
		for(int j=1; j<l2; j++){
			matk(i1,j)= matk(i1,j) * matkr(l1-i1-1,l2-j-1);
		}
	}
    return matk;
}

/*
 * Kernelized time elastic centroid of a pair of time series.
 */
UniformTimeseries KDTW::averagePairT(UniformTimeseries &ts1, UniformTimeseries &ts2, double sigma)
{
	UniformTimeseries tsout;
	matrix<double> ama=KDTW().kernel_AMA(ts1, ts2, sigma);
	tsout = KDTW().add_kernelT(ts1, ts2, ama, 1, 2);

    return tsout;
}

/*
 * Barycenter of a vector of Points.
 */
Point barycenter(vector<Point> tab) {
	if (tab.size() < 1) {
		cout << "RuntimeException('empty double tab');" << endl;
	}

	unsigned int dim = tab[0].size();
	Point sum(dim);
	for (unsigned int i=0; i< dim; i++)
		sum[i]=0;
	for (unsigned int i=0; i< tab.size(); i++) {
		for(unsigned int j=0; j< dim; j++)
			sum[j] += (tab[i][j]);
	}
	for (unsigned int i=0; i< dim; i++)
		sum[i]=sum[i]/tab.size();
	return sum;
}

/*
 * Evaluate the KDBA centroid (kernelized version of the DBA algorithm.
 */
UniformTimeseries KDTW::add_kernelDBA(const UniformTimeseries C, std::vector<UniformTimeseries> sequences,
		double sigma, double lab)
{
double _dbl_min=1e-300;
UniformTimeseries ts2;
std::vector<std::vector<Point>> tupleAssociation;
tupleAssociation.resize(C.size());
const int l1 = C.size();
int lmax=l1;
const int dim=C.data[0].size();

for(unsigned int n=0; n<sequences.size(); n++){
	ts2=sequences[n];
	matrix<double> mat=kernel_AMA(C,ts2,sigma);
    const int l2 = ts2.size();

    sigma=sigma*dim;

    if (l1 == 0 || l2 == 0)
        throw std::runtime_error("timeseries cannot be empty");

    lmax=l1;
    if(lmax<l2){
	lmax=l2;
	}

    double z;
    UniformTimeseries ts22 = UniformTimeseries(lmax, dim);
    std::vector<double> normzj(lmax);
    for(int i=0; i<lmax; i++){
    	for(int k=0; k<dim; k++){
    		ts22.data[i][k]=0;
    	}
    	normzj[i]=0;
     }
     for(int i=0; i<l1; i++){
    	    for(int j=0; j<l2; j++){
    	    	z=mat(i+1,j+1)+_dbl_min;
     	    	for(int k=0; k<dim; k++){
      	    	  ts22.data[i][k]+=(ts2.data[j][k])*z;
    	    	}
  	    	    normzj[i]+=z;
    	    }
            for(int k=0; k<dim; k++){
            	ts22.data[i][k]=ts22.data[i][k]/normzj[i];
            }
        	tupleAssociation[i].push_back(ts22.data[i]);
    }
}
UniformTimeseries ts11=UniformTimeseries();
for(int i=0; i<lmax; i++){
    ts11.data.push_back(barycenter(tupleAssociation[i]));
}

ts11.label=lab;
return ts11;
}

/*
 * Kernelized centroid of a pair of time series. Version with floor and cell time approximations
// *** requires the AMA matrix in input and two scalar values used for the averaging
 */
UniformTimeseries KDTW::add_kernelT(const UniformTimeseries &ts1, const UniformTimeseries &ts2,
		matrix<double> ama, double scal, double gscal)
{
	double _dbl_min=1e-300;
    const int l1 = ts1.size();
    const int l2 = ts2.size();
    const int dim=ts1.data[0].size();

    if (l1 == 0 || l2 == 0)
        throw std::runtime_error("timeseries cannot be empty");

    double alpha;
    int lmax=l1;
    if(lmax<l2){
	lmax=l2;
	}

    double z;
    UniformTimeseries ts11 = UniformTimeseries(lmax, dim);
    UniformTimeseries ts22 = UniformTimeseries(lmax, dim);
    std::vector<double> normzi(lmax);
    std::vector<double> normzj(lmax);
    for(int i=0; i<lmax; i++){
    	for(int k=0; k<dim; k++){
    		ts11.data[i][k]=0;
    		ts22.data[i][k]=0;
    	}
    	normzi[i]=_dbl_min;
    	normzj[i]=_dbl_min;
     }

     for(int i=0; i<lmax; i++){

    	if(i<l1){
    	    const Point &p1 = ts1.data[i];
    	    for(int j=0; j<l2; j++){
                const Point &p2 = ts2.data[j];
    	    	z=ama(i+1,j+1)+_dbl_min;
    	    	alpha = ((float)(i+j))/2.0 - floor((i+j)/2);

    	    	for(int k=0; k<dim; k++){
    	    	  ts11.data[floor((float)(i+j)/2.0)][k]+=alpha*(scal*p1[k]+p2[k])*z;
    	    	  ts11.data[ceil((float)(i+j)/2.0)][k]+=(1-alpha)*(scal*p1[k]+p2[k])*z;
    	    	}
  	    	    normzi[floor((float)(i+j)/2.0)]+=alpha*z;
  	    	    normzi[ceil((float)(i+j)/2.0)]+=(1-alpha)*z;
    	    }
    	}
    	if(i<l2){
    	    const Point &p2 = ts2.data[i];
    	    for(int j=0; j<l1; j++){
                const Point &p1 = ts1.data[j];
    	    	z=ama(j+1,i+1)+_dbl_min;
    	    	alpha = ((float)(i+j))/2.0 - floor((i+j)/2);

    	    	for(int k=0; k<dim; k++){
      	    	  ts22.data[floor((float)(i+j)/2.0)][k]+=alpha*(scal*p1[k]+p2[k])*z;
      	    	  ts22.data[ceil((float)(i+j)/2.0)][k]+=(1-alpha)*(scal*p1[k]+p2[k])*z;
    	    	}
  	    	    normzj[floor((float)(i+j)/2.0)]+=alpha*z;
  	    	    normzj[ceil((float)(i+j)/2.0)]+=(1-alpha)*z;
    	    }
    	}
    }

    for(int i=0; i<lmax; i++){
     	for(int k=0; k<dim; k++)
			ts11.data[i][k]=(ts11.data[i][k]/normzi[i]+ts22.data[i][k]/normzj[i])/(gscal)/2;
    }

    return ts11;

}

/*
 * Kernelized centroid of a pair of time series. Version with round time approximation
// *** requires the AMA matrix in input and two scalar values used for the averaging
 */
UniformTimeseries KDTW::add_kernelTbis(const UniformTimeseries &ts1, const UniformTimeseries &ts2,
		matrix<double> ama, double scal, double gscal)
{
	double _dbl_min=1e-300;
    const int l1 = ts1.size();
    const int l2 = ts2.size();
    const int dim=ts1.data[0].size();

    if (l1 == 0 || l2 == 0)
        throw std::runtime_error("timeseries cannot be empty");

    int lmax=l1;
    if(lmax<l2){
	lmax=l2;
	}

    std::vector<double> normzi(lmax);
    std::vector<double> normzj(lmax);
    double z;
    UniformTimeseries ts11 = UniformTimeseries(lmax, dim);
    UniformTimeseries ts22 = UniformTimeseries(lmax, dim);

    for(int i=0; i<lmax; i++){
    	for(int k=0; k<dim; k++){
    		ts11.data[i][k]=0;
    		ts22.data[i][k]=0;
    	}
    	normzi[i]=0.0;
    	normzj[i]=0.0;
     }

    for(int i=0; i<lmax; i++){

    	if(i<l1){
    	    const Point &p1 = ts1.data[i];
    	    for(int j=0; j<l2; j++){
                const Point &p2 = ts2.data[j];
    	    	z=ama(i+1,j+1)+_dbl_min;
    	    	for(int k=0; k<dim; k++){
    	    	  ts11.data[round((float)(i+j)/2.0)][k]+=(scal*p1[k]+p2[k])*z;
    	    	}
  	    	    normzi[round((float)(i+j)/2.0)]+=z;
    	    }
    	}
    	if(i<l2){
    	    const Point &p2 = ts2.data[i];
    	    for(int j=0; j<l1; j++){
                const Point &p1 = ts1.data[j];
    	    	z=ama(j+1,i+1)+_dbl_min;
    	    	for(int k=0; k<dim; k++){
    	   	    	  ts22.data[round((float)(i+j)/2.0)][k]+=(scal*p1[k]+p2[k])*z;
    	    	    	}
    	  	     normzj[round((float)(i+j)/2.0)]+=z;
    	    }
    	}

    }
    for(int i=0; i<lmax; i++){
     	for(int k=0; k<dim; k++)
    			ts11.data[i][k]=(ts11.data[i][k]/normzi[i]+ts22.data[i][k]/normzj[i])/gscal/2;
    }
    return ts11;
}

