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
#include <memory>
#include <iostream>
#include <fstream>
#include <sstream>
#include "UniformTimeseries.h"

using namespace std;

/*
 * Constructors
 */
UniformTimeseries::UniformTimeseries(){
data={};
}

UniformTimeseries::UniformTimeseries(const std::vector<Point> &data)
: data(data)
{label=0;}
UniformTimeseries::UniformTimeseries(const std::vector<Point> &data, double lab)
: data(data)
{label=lab; subject=0;}
UniformTimeseries::UniformTimeseries(const std::vector<Point> &data, int _subject, double lab)
: data(data)
{label=lab; subject=_subject;}


UniformTimeseries::UniformTimeseries(int L, int D){
for(int i=0; i<L; i++){
	Point p={};
	for(int k=0; k<D; k++)
		p.push_back(0.0);
	data.push_back(p);
   }
}

/*
 * Time reversal of a time series
 */
UniformTimeseries UniformTimeseries::reverse(){
int L=this->data.size();
int D=this->data[0].size();
vector<Point> _data;
for(int i=0; i<L; i++){
	Point p;
	for(int k=0; k<D; k++)
		p.push_back(this->data[L-i-1][k]);
	_data.push_back(p);
   }
return UniformTimeseries(_data);
}

/*
 * Print the content of a time series
 */
void UniformTimeseries::print(){
   for(unsigned  int j=0; j < data.size(); j++){
        for(unsigned  int k=0; k < data[j].size(); k++){
        	std::cout << data[j][k] << " ";
	}
    std::cout << std::endl;
  }
}

/*
 * Sum the samples of a time series. The sum is a vector of values according to the multidimensionality of the time series
 */
std::vector<double> UniformTimeseries::sumSamples(){
   std::vector<double> sum (data[0].size(),0);
   for(unsigned  int j=0; j < data.size(); j++){
        for(unsigned  int k=0; k < data[j].size(); k++){
        	sum[k]+=data[j][k];
	}
  }
   return sum;
}

/*
 * Multiply all the sample of a time series by a scalar value
 */
void UniformTimeseries::scalMult(double scal){
   for(unsigned  int j=0; j < data.size(); j++){
        for(unsigned  int k=0; k < data[j].size(); k++){
        	data[j][k]=scal*data[j][k];
	}
  }
}

/*
 * Computes the min sample value along dimension k
 */
double UniformTimeseries::minOnDim(int k){
   double min=DBL_MAX, x;
   for(unsigned  int j=0; j < data.size(); j++){
        	x=data[j][k];
        	if(x<min)
        		min=x;
  }
  return min;
}

/*
 * Computes the max sample value along dimension k
 */
double UniformTimeseries::maxOnDim(int k){
   double max=-DBL_MAX, x;
   for(unsigned  int j=0; j < data.size(); j++){
        	x=data[j][k];
        	if(x>max)
        		max=x;
  }
  return max;
}


/*
 * Computes the Euclidean norm of a time series
 */
double UniformTimeseries::norm(){
   double norm=0;
   for(unsigned  int j=0; j < data.size(); j++){
        for(unsigned  int k=0; k < data[j].size(); k++){
        	norm+=data[j][k]*data[j][k];
	}
  }
   return sqrt(norm);
}

/*
 * Computes the L1 norm of a time series
 */
double UniformTimeseries::norml1(){
   double norm=0;
   for(unsigned  int j=0; j < data.size(); j++){
        for(unsigned  int k=0; k < data[j].size(); k++){
        	norm+=fabs(data[j][k]);
	}
  }
   return (norm);
}


/*
 * Load a time series from a file complying with the UCR file format (unidimensional time series)
 */
std::vector<UniformTimeseries> UniformTimeseries::parse_ucr_file(
        std::string &filename,
        int limit) // `limit` parses only a finite number of lines
{
    auto results = std::vector<UniformTimeseries>{};

    std::string line;
    std::ifstream infile (filename.c_str(), std::ifstream::in);
    if (!infile.good()) {
        std::cerr << "file not found: " << filename << std::endl;
    }
    int count = 0;
    while (std::getline(infile, line)) {
        double label;
        auto &&sline = std::stringstream{line};
        sline >> label;
        double d;
        auto ts = std::vector<std::vector<double>>{};
        while (sline >> d) {
            const auto point = std::vector<double>{d};
            ts.push_back(point);
        }
        ts.shrink_to_fit();

        results.push_back(UniformTimeseries(ts, label));

        count++;
        if (limit > 0 && count >= limit) break;
    }

    results.shrink_to_fit();
    infile.close();
    return results;
}

/*
 * Load a time series from a file complying with the MDTS file format (multidimensional time series)
 */
std::vector<UniformTimeseries> UniformTimeseries::parse_mdts_file(
        std::string &filename,
        int dim,
        int limit, // `limit` parses only a finite number of lines
        const bool withSubject)
{
    auto results = std::vector<UniformTimeseries>{};

    std::string line;
    std::ifstream infile (filename.c_str(), std::ifstream::in);
    if (!infile.good()) {
       std::cerr << "file not found: " <<  filename;
    }
    int count = 0;
    while (std::getline(infile, line)) {
        double label;
        int subject=0;
        auto &&sline = std::stringstream{line};
        if(withSubject)
        	sline >> subject;
        sline >> label;
        double idx, d;
        char c;
        auto ts = std::vector<std::vector<double>>{};
        while (sline >> idx && sline >> c && sline >> d) {
	    int i=0;
	    std::vector<double> point = std::vector<double>{};
	    while(i< dim && sline >> idx && sline >> c && sline >> d){
            	point.push_back(d);
 		i++;
	    }
	    ts.push_back(point);
        }
        ts.shrink_to_fit();

        results.push_back(UniformTimeseries(ts, subject, label));

        count++;
        if (limit > 0 && count >= limit) break;
    }
    infile.close();
    results.shrink_to_fit();
    return results;
}


/*
 * Extract from the <ds> data set  (vector of uniform time series) a subset containing all the time series tagged with label <lab>
 */
std::vector<UniformTimeseries> UniformTimeseries::extractTSfromLabelDS(std::vector<UniformTimeseries> ds, int lab){
    auto results = std::vector<UniformTimeseries>{};
    for(unsigned int i=0; i < ds.size(); i++){
    	if(ds[i].label==lab)
        	 results.push_back(ds[i]);
     }
    return results;
}


/*
 * save a time series on file <fname>
 */
void UniformTimeseries::saveFile(std::string fname){
ofstream fout;
fout.open(fname);
for(unsigned int j=0; j < data.size(); j++){
     for(unsigned int k=0; k < data[j].size(); k++){
    	 fout << data[j][k] << " ";
     }
     fout << endl;
 }
fout.close();
}

/*
 * save dataset <ds> of time series on file <fname>, according to "SVM" format
 */
void UniformTimeseries::saveUniformDatasetSVM(std::string fname, std::vector<UniformTimeseries> ds){
std::ofstream fout;
fout.open(fname.c_str());
int kk;
  for (unsigned  int i = 0 ; i < ds.size(); i++ ){
     fout << ds[i].subject << " " << ds[i].label << " ";
     kk=0;
     for(unsigned int j=0; j < ds[i].data.size(); j++){
	    fout << kk << ":" <<j << " "; kk++;
        for(unsigned int k=0; k < ds[i].data[j].size(); k++){
		    fout << kk << ":" << ds[i].data[j][k] << " "; kk++;
	}
     }
     fout << std::endl; fout.flush();
  }
fout.close();
}

/*
 * save dataset <ds> of time series on file <fname>, according to "UCR" format
 */
void UniformTimeseries::saveUniformDatasetUCR(std::string fname, std::vector<UniformTimeseries> ds){
std::ofstream fout;
fout.open(fname.c_str());

  for (unsigned  int i = 0 ; i < ds.size(); i++ ){
     fout << ds[i].label << " ";
     for(unsigned int j=0; j < ds[i].data.size(); j++){
        for(unsigned int k=0; k < ds[i].data[j].size(); k++){
		    fout << ds[i].data[j][k] << " ";
        }
     }
     fout << std::endl; fout.flush();
  }
fout.close();
}

/*
 * save dataset <ds> of time series on file <fname>, in comumn.
 */
void UniformTimeseries::saveUniformDatasetInColumns(std::string fname, std::vector<UniformTimeseries> ds){
std::ofstream fout;
fout.open(fname.c_str());
  for (unsigned int i = 0 ; i < ds.size(); i++ ){
     //fout << ds[i].label << " ";
     for(unsigned  int j=0; j < ds[i].data.size(); j++){
        for(unsigned  int k=0; k < ds[i].data[j].size(); k++){
		    fout << ds[i].data[j][k] << " ";
        }
        fout << endl;
     }
     fout << std::endl; fout.flush();
  }
fout.close();
}

/*
 * Relabel  dataset <ds> of time series from 1 to #categories
 */
std::vector<UniformTimeseries>  UniformTimeseries::relabelDS(std::vector<UniformTimeseries>& ds, std::map<int,int>& tab){
	unsigned int LL=ds.size();
	int l=1;//tab.size();
	for (unsigned int i=0; i<LL; i++){
		if(tab.find(ds[i].label) == tab.end()){
			tab[(int)ds[i].label]=l;
			ds[i].label=l;
			l++;
		}
		else
			ds[i].label=tab[(int)ds[i].label];
	}
	return ds;
}

