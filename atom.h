#ifndef atom_h
#define atom_h
#include <iostream>
#include <fstream>
#include <list>
#include <array>
#include <cmath>
#include <algorithm>
#include <vector>
#include <string>
class atom{
	public:
		atom()=default;
		atom(double elect,std::vector<double> p):charge(elect),position(p){};
		atom(double elect,std::string place);
		atom(atom& B):position(B.position),charge(B.charge){
		};
		atom& operator=(atom& A){
			position=A.position;
			charge=A.charge;
			tick=A.tick;
			return *this;
		}
	  friend atom& operator <<(atom&,std::fstream&);
		friend void sorted_ba(std::vector<atom> &Ba,std::vector<double> p,int cell);
		void setcharge(double elect){
			charge=elect;
		}
		double getcharge(){
			return charge;
		}
		std::vector<double> getposition(){
			return position;
		}
		int gettick(){
			return tick;
		}
		~atom(){
			position.clear();
		};
	private:
		std::vector<double> position;
		double charge;
		int tick;
};
double convert(std::string input);
double average(std::list<double>& input);
std::vector<double> dist(std::vector<double>,std::vector<double>,std::vector<double>&);
std::vector<double> polar(atom& A,atom& B,std::vector<double>& p);
double distance(std::vector<double> A);
std::vector<int> changeindex(int,int);
std::vector<int> findneighbor_oxy(int,int);
std::vector<int> findneighbor_ba(int,int);
std::vector<double>& operator +=(std::vector<double>&,std::vector<double>&);
std::vector<double>& operator /=(std::vector<double>&,double);
std::vector<double>& operator *=(std::vector<double>&,double);
#endif
