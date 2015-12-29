/*************************************************************************
eKATS1.0 (source code generated 2015-12-15)
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
#include "utils.h"

std::vector<int> shuffleIndex(int N){
	std::map<int,int> tab;
    int alea;
    auto results = std::vector<int>{};
    int n=0;
    while (n<N){
      alea = (int) rand() % N;
	  if(tab.find(alea) == tab.end()){
		  results.push_back(alea);
		  tab[alea] = n;
		  n++;
	  }
    }
    return results;
}

void saveFileMatrix(std::string fName, std::vector<std::vector<double>> mat){
	std::ofstream fout;
	fout.open(fName.c_str());
	for(unsigned int i=0; i< mat.size(); i++){
		for(unsigned int j=0; j< mat[i].size(); j++){
			fout << mat[i][j] << " ";
		}
		fout << std::endl;
	}
	fout.close();
}

