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
#include <iostream>
#include <cfloat>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <map>
#include "UniformTimeseries.h"
#include "DTW.h"
#include "GramMatrix.h"
#include "cluster.h"
#include "kmedoids.h"
#include "kcentroids.h"
#include "DBA.h"

#define ALL 0
#define DTW_MED 1
#define iDBA 2
#define KDTW_MED 3
#define iKDBA 4
#define pKDTW_PWA 5

int main(int argc, char **argv)
{
	// This example program produce 1 medoid/centroid for each class of
	// time series contained into the dataset:
	// - medoids are obtained either through DTW or KDTW mesures
	// - centroid estimates are obtained through iDBA, iKDBA or pKDTW_PWA methods described in
	//   https://hal.archives-ouvertes.fr/hal-01155134
	srand (time(NULL));
	std::string path="../DATASETS/", fileName;
    ofstream fout;
	string name_dataset = "CBF";
    double sigma=1;
    int nclass=1.0;
    int nclust =1; // here we consider a single cluser
	int type=ALL;
	std::map<int,int> tabLab;
	int dim=1;
	int FileType=_UCR;

	int i=1;
	while(i<argc){
		if(strcmp(argv[i],"-h")==0){
			cout << "./eKats -sigma <sigma> -path <PATH>  -DS <dataset>  -TYPE <DTW_MED|KDTW_MED|iDBA|pKDTW-PWA|iKDBA|ALL>" << endl;
			cout << "example: ./eKats -sigma 1 -PATH ../DATASETS/ -DS CBF -TYPE ALL" << endl;
			cout << "Options: " << endl;
			cout << "-sigma <sigma> : value for sigma (sigma = 1.0/nu in the papers) for KDTW_MED, pKDTW-PWA and iKDBA" << endl;
			cout << "-path <PATH> : path to locate the data set" << endl;
			cout << "-DS <dataset> : name of the data set. <dataset>_TRAIN and <dataset>_TEST should exist in /PATH/<dataset>/ directory in the UCR or MDTS format" << endl;
			cout << "-TYPE <DTW_MED|KDTW_MED|iDBA|pKTDW-PWA|iKDBA|ALL> is the tipe of centroid or methoid method" << endl;
			cout << "-MDTS : the dataset file is in the MDTS format." << endl;
			cout << "-UCR : the dataset file is in the UCR format." << endl;
			exit(0);
		}
		else if(strcmp(argv[i],"-sigma")==0){ // sigma value for KDTW medoid or iKDBA/p-KDTW-PWA centroid methods
			sigma=atof(argv[i+1]);
			i+=2;
		}
		else if(strcmp(argv[i],"-DS")==0){ // dataset name
			name_dataset=(argv[i+1]);
			i+=2;
		}
		else if(strcmp(argv[i],"-PATH")==0){ // path to the dataset
			path=(argv[i+1]);
			i+=2;
		}
		else if(strcmp(argv[i],"-DIM")==0){
			dim=atoi(argv[i+1]);
			i+=2;
		}
		else if(strcmp(argv[i],"-MTDS")==0){
			FileType=_MTDS;
			i+=1;
		}
		else if(strcmp(argv[i],"-UCR")==0){
			FileType=_UCR;
			i+=1;
		}
		else if(strcmp(argv[i],"-TYPE")==0){
			if(strcmp(argv[i+1],"pKDTW-PWA")==0) // probabilistic KDTW-PWA centroid
				type=pKDTW_PWA;
			else if(strcmp(argv[i+1],"KDTW_MED")==0) // KDTW Medoid
				type=KDTW_MED;
			else if(strcmp(argv[i+1],"DTW_MED")==0) // DTW MEdoid
				type=DTW_MED;
			else if(strcmp(argv[i+1],"iDBA")==0) // iterative DBA centroid
				type=iDBA;
			else if(strcmp(argv[i+1],"iKDBA")==0) // iterative KDBA centroid
				type=iKDBA;
			else if(strcmp(argv[i+1],"ALL")==0) // All methods
				type=ALL;
			else{
						cout << "unknown centroid/medoid type: " << argv[i+1] << endl;
						exit(0);
					}
			i+=2;
		}
		else{
		   cout << "unknown argument: " << argv[i] << endl;
		   i++;
		}
	}

	std::string file_train;
	file_train.assign(path);
	file_train.append(name_dataset);
	std::vector<UniformTimeseries> dsTRAIN, ds2, ds;

	std::vector<UniformTimeseries> centroid;
	std::vector<UniformTimeseries> dsCat;

	if(FileType == _UCR)
		dsTRAIN  = UniformTimeseries::parse_ucr_file(file_train);
	else
		 dsTRAIN  = UniformTimeseries::parse_mdts_file(file_train, dim, 0, false);

	if(dsTRAIN.size()==0){
		exit(0);
	}
	dsTRAIN=UniformTimeseries::relabelDS(dsTRAIN, tabLab);
	UniformTimeseries::saveUniformDatasetUCR(file_train.append(".rel"), dsTRAIN);

	std::cout << "Processing dataset " <<name_dataset << endl;
	std::cout <<"# Time series: " << dsTRAIN.size() << std::endl;
	std::cout << "Time series length: " << dsTRAIN[0].size() << std::endl;
	nclass=tabLab.size();
	std::cout << "nb CLASS: " << nclass << std::endl;

	// DTW Medoids
	if(type==ALL || type==DTW_MED){
	  cout << "Running DTW MEDOID .." << endl;
	  ds = kmedoids_dtw(dsTRAIN, nclass, nclust);
	  for(unsigned int n=0; n<ds.size(); n++){
	    fileName.assign("../octave/DTW_medoid_");
	    fileName.append(name_dataset);
	    fileName.append("_");
	    fileName.append(std::to_string(n+1));
	    fileName.append(".dat");
	    ds[n].saveFile(fileName);
	    }
	  ds.clear();
	}
	// iDBA iterative DBA
	if(type==ALL || type==iDBA){
		cout << "Running DTW iDBA CENTROID .." << endl;
		ds = kcentroid_idba(dsTRAIN, nclass, nclust);
		for(unsigned int n=0; n<ds.size(); n++){
			fileName.assign("../octave/iDBA_centroid_");
			fileName.append(name_dataset);
			fileName.append("_");
			fileName.append(std::to_string(n+1));
			fileName.append(".dat");
			ds[n].saveFile(fileName);
		}
		ds.clear();
	}
	// KDTW Medoids
	if(type==ALL || type==KDTW_MED){
		cout << "Running KDTW MEDOID .." << endl;
		ds = kmedoids_kdtw(dsTRAIN, nclass, nclust, sigma);
		for(unsigned int n=0; n<ds.size(); n++){
			fileName.assign("../octave/KDTW_medoid_");
			fileName.append(name_dataset);
			fileName.append("_");
			fileName.append(std::to_string(n+1));
			fileName.append(".dat");
			ds[n].saveFile(fileName);
		}
		ds.clear();
	}
	// pKDTW-PWA: probabilistic KDTW Pairwise agregate
	if(type==ALL || type==pKDTW_PWA){
		cout << "Running pKDTW-PWA CENTROID .." << endl;
		ds = kcentroid_kdtw_pwa(dsTRAIN, nclass, nclust, sigma);
		for(unsigned int n=0; n<ds.size(); n++){
			fileName.assign("../octave/pKDTW_PWA_centroid_");
			fileName.append(name_dataset);
			fileName.append("_");
			fileName.append(std::to_string(n+1));
			fileName.append(".dat");
			ds[n].saveFile(fileName);
		}
		ds.clear();
	}
	// iKDBA: iterative Kernel DBA
	if(type==ALL || type==iKDBA){
		cout << "Running iKDBA CENTROID .." << endl;
		ds = kcentroid_ikdba(dsTRAIN, nclass, nclust, sigma);
		for(unsigned int n=0; n<ds.size(); n++){
			fileName.assign("../octave/iKDBA_centroid_");
			fileName.append(name_dataset);
			fileName.append("_");
			fileName.append(std::to_string(n+1));
			fileName.append(".dat");
			ds[n].saveFile(fileName);
		}
		ds.clear();
	}
    return 0;
}
