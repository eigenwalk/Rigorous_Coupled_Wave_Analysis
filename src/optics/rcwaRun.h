#ifndef RCWA_RUN_H
#define RCWA_RUN_H

#include <armadillo> 
#include <iostream>
#include <math.h>
#include <string>
#include <map>
#include <filesystem>
//#include "manager/jobManager.h"
#include "structures/inputStructure.h"

namespace fs = std::filesystem;
using namespace std;
using namespace arma; 

class RCWA 
{
    private:

    public:
		RCWA();
		RCWA(Inputs inpts):
			m_input(inpts){};
		~RCWA();
        //auto init(const YAML::Node& in)->void;
        auto init()->bool;
        auto start()->void;
        auto saveResult()->void;
		//auto transferMatrix(double& wave, double& angle)->void;

		template <typename T>
		auto stonSplit(string& str)->vector<T>;
		auto stringSplit(string& str)->vector<string>;
		auto fillVoxel(vector<int>& v, mat& M, int iy)->void;
		template <typename T>
		auto checkvLength(vector<T>& v, int& s)->bool;
		auto loadNK(string& mat)->void;
		auto loadVoxel()->bool;

	public:
		Inputs m_input;
		arma::vec m_waves;
		arma::vec m_polar_angles;
		arma::vec m_azi_angles;
		vector<arma::mat> m_stack;
		map<string, arma::vec> m_n;
		map<string, arma::vec> m_k;
		map<string, arma::cx_vec> m_eps;
		map<string, fs::path> m_material;
		map<double, map<double, vector<double>>> m_Rs;
		map<double, map<double, vector<double>>> m_Rp;
		map<double, map<double, vector<double>>> m_R;
		
		float m_unit;
};

#endif
