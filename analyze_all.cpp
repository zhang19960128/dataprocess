#include <iostream>
#include <fstream>
#include <list>
#include <algorithm>
#include <vector>
#include <string>
#include "atom.h"
int main(){
	std::fstream filein;
	std::string templine;
	filein.open("/global/cscratch1/sd/jiahaoz/datatest/T5/dump.xyz",std::fstream::in);
	std::string patternla_constant="ITEM: BOX BOUNDS pp pp pp";
	std::string pattern_atoms="ITEM: ATOMS x y z";
	std::vector<std::string> p_string(3);
	std::vector<double> p(3,0);
	std::list<double> pxall;
	std::list<double> pyall;
	std::list<double> pzall;
	std::vector<double> temp_displace(3,0.0);
	std::vector<double> displace(3,0.0);
	std::list<double> displace_x;
	std::list<double> displace_y;
	std::list<double> displace_z;
	int cell=10;
	std::vector<int> temp;
	std::vector<atom> Ba(cell*cell*cell);//on the corner of the cell.
	std::vector<atom> Ti(cell*cell*cell);//in the cell of the cell.
	std::vector<atom> Oxy(3*cell*cell*cell);//on the face of the cell.
	size_t count=0;
	//read the data from the txt files.
	for(std::string line;getline(filein,line);){
		if(line.find(patternla_constant)!=std::string::npos){
			//find the crystal constant in the file and then starting processing them.
			for(size_t i=0;i<3;i++){
				getline(filein,p_string[i]);
				p[i]=convert(p_string[i]);//this also store the periodically length of the crystal.
			}
			pxall.push_back(p[0]/cell);
			pyall.push_back(p[1]/cell);
			pzall.push_back(p[2]/cell);
		}
		if(line.find(pattern_atoms)!=std::string::npos){
			//find the items in the file and starting processing them.
			for(std::vector<atom>::iterator a=Ba.begin();a!=Ba.end();a++){
				(*a)<<filein;
				(*a).setcharge(1.34730);
			}
			for(std::vector<atom>::iterator a=Ti.begin();a!=Ti.end();a++){
				(*a)<<filein;
				(*a).setcharge(1.28905);
			}	
			for(std::vector<atom>::iterator a=Oxy.begin();a!=Oxy.end();a++){
				(*a)<<filein;
				(*a).setcharge(-0.87878);
			}
			//specify the index for Ti
			for(size_t i=0;i<cell*cell*cell;i++){
				//the index for the neighbor oxygen should be (0,n1,n2,n3),(0,n1+1,n2,n3),(1,n1,n2,n3),(1,n1,n2+1,n3),(2,n1,n2,n3),(2,n1,n2,n3+1)
				temp=findneighbor_oxy(i,cell);
				for(size_t k=0;k<3;k++){
					displace[k]=0.0;
				}
				for(std::vector<int>::iterator a=temp.begin();a!=temp.end();a++){
				temp_displace=dist(Oxy[*a].getposition(),Ti[i].getposition(),p);
				displace+=temp_displace;
				}
				displace_x.push_back(displace[0]/6);
				displace_y.push_back(displace[1]/6);
				displace_z.push_back(displace[2]/6);
			}
		}
	}
	//starting calculating the lattice constant.
	std::fstream fileout;
	fileout.open("result.txt",std::fstream::out);
	fileout<<"the average lattice constant is:"<<std::endl;
	fileout<<average(pxall)<<" "<<average(pyall)<<" "<<average(pzall)<<std::endl;
	fileout<<"the average displacement is:"<<std::endl;
	fileout<<average(displace_x)<<" "<<average(displace_y)<<" "<<average(displace_z)<<std::endl;
	//starting calculating the displacement of the material;
}
