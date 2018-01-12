#include <iostream>
#include <fstream>
#include <cmath>
#include <list>
#include <algorithm>
#include <vector>
#include <string>
#include "atom.h"
double convert(std::string input){
	double result,x1,x2;
  x1=std::stof(input.substr(0,input.find_first_of(" ")));
	x2=std::stof(input.substr(input.find_first_of(" ")+1));
	result=(x2-x1);
	return result;
};
double average(std::list<double> &input){
	double sum=0;
	for(std::list<double>::iterator a=input.begin();a!=input.end();a++){
		sum=sum+*a;
	}
	return sum/input.size();
};
atom& operator <<(atom& input,std::fstream& file){
			std::string line;
			getline(file,line);
			std::vector<double> p(3,0);
			std::string substemp,stringtoken;
			size_t i=0;
			while(stringtoken!=line){
			stringtoken=line.substr(0,line.find_first_of(" "));
			line=line.substr(line.find_first_of(" ")+1);
			if(i>2) break;
			p[i]=std::stof(stringtoken);
			i++;
			};
			input.position=p;
			return input;
		};
//compute the vector from A to B where A is the origin and B is the goal.
std::vector<double> dist(std::vector<double> A,std::vector<double> B,std::vector<double> P){
	//we need the periodically length of the crystal.
	std::vector<double> len(3,0);
	double tempdata;
	for(size_t i=0;i<3;i++){
		tempdata=B[i]-A[i];
		tempdata=(tempdata/P[i]-round(tempdata/P[i]))*P[i];
		len[i]=tempdata;
	}
	return len;
};
double distance(std::vector<double> A){
	double sum=0;
	for(size_t i=0;i<3;i++){
		sum=sum+A[i]*A[i];
	}
	return sqrt(sum);
};
atom::atom(double elect,std::string place){
			std::vector<double> p(3,0);
			std::string substemp,stringtoken;
			size_t i=0;
			while(stringtoken!=place){
			stringtoken=place.substr(0,place.find_first_of(" "));
			place=place.substr(place.find_first_of(" ")+1);
			if(i>2) break;
			p[i]=std::stof(stringtoken);
			i++;
			};
			atom(elect,p);
	  };
std::vector<int> changeindex(int index,int cell){
	std::vector<int> a(3,0);
	a[2]=floor(index/(cell*cell));
	index=index-a[2]*cell*cell;
	a[1]=floor(index/cell);
	a[0]=index-cell*a[1];
	return a;
}
std::vector<int> findneighbor_oxy(int index,int cell){
	std::vector<int> a(6,0);
	std::vector<int> temp;
	temp=changeindex(index,cell);
	a[0]=temp[0]+temp[1]*cell+temp[2]*cell*cell+0*cell*cell*cell;
	a[1]=temp[0]+temp[1]*cell+((temp[2]+1)%cell)*cell*cell+0*cell*cell*cell;
	a[2]=temp[0]+temp[1]*cell+temp[2]*cell*cell+1*cell*cell*cell;
	a[3]=temp[0]+((temp[1]+1)%cell)*cell+temp[2]*cell*cell+1*cell*cell*cell;
	a[4]=temp[0]+temp[1]*cell+temp[2]*cell*cell+2*cell*cell*cell;
	a[5]=(temp[0]+1)%cell+temp[1]*cell+temp[2]*cell*cell+2*cell*cell*cell;
	return a;
}
std::vector<double>& operator +=(std::vector<double>& A,std::vector<double>& B){
	for(size_t i=0;i<A.size();i++){
		A[i]=A[i]+B[i];	
	}
	return A;
}
