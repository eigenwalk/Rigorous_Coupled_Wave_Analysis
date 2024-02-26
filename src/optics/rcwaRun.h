#ifndef RCWA_RUN_H
#define RCWA_RUN_H

#include <armadillo> 
#include <iostream>
#include <math.h>
#include <string>
#include <map>
#include <filesystem>
#include "structures/inputStructure.h"
#include "light.h" 
#include "scatterMatrix.h" 

namespace fs = std::filesystem;
using namespace std;
using namespace arma; 

class RCWA 
{
    private:

    public:
		RCWA();
		RCWA(Input input):
			m_input(input){};
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
		auto scatteringMatrix(int& i, int& j, int& k)->void;
		auto getToeplitzMatrix(cx_mat& Er_toep, int& idx)->void;
		auto setToeplitzMatrix(int& wid)->void;
		auto setNormalVectorField(cx_mat& M, int i)->map<string, mat>;
		auto getEigen(cx_mat& W, cx_mat& V, cx_mat& Eval, string mode, int idx)->void;
		auto getAB(cx_mat& A, cx_mat& B, cx_mat& Wi, cx_mat& Vi)->void;
		auto getABX(cx_mat& A, cx_mat& B, cx_mat& X, cx_mat& Wi, cx_mat& Vi, cx_mat& W_air, cx_mat& V_air, cx_mat& Eval, double& k0, double& thk)->void;
		auto getSmat_r(cx_mat& S11, cx_mat& S12, cx_mat& S21, cx_mat& S22, cx_mat& A, cx_mat& B)->void;
		auto getSmat_t(cx_mat& S11, cx_mat& S12, cx_mat& S21, cx_mat& S22, cx_mat& A, cx_mat& B)->void;
		auto getSmat_l(cx_mat& S11, cx_mat& S12, cx_mat& S21, cx_mat& S22, cx_mat& A, cx_mat& B, cx_mat& X)->void;
		auto setZero(cx_mat& M, double decimal)->void;
		
		auto saveArma(cx_mat& M, string file, int id)->void;
		auto saveArma(mat& M, string file, int id)->void;
		auto fftshift2d(const cx_mat& X)->cx_mat;

		template <typename T>
		auto censhift2d(const T& X)->T;
		auto findNVb(map<pair<int, int>, pair<double, double>>& nvb, mat& Nx, mat& Ny)->void;
		auto updateNV(map<pair<int, int>, pair<double, double>>& nvb, mat& Nx, mat& Ny)->void;
		auto normalizeNV(mat& Nx, mat& Ny)->void;

	public:
		Input m_input;
		Light* m_light;
		arma::vec m_waves;
		arma::vec m_polar_angles;
		arma::vec m_azi_angles;
		vector<arma::mat> m_stack;
		vector<arma::cx_mat> m_toep;
		vector<arma::cx_mat> m_toep_i;
		vector<arma::cx_mat> m_toep_delta;
		vector<arma::cx_mat> m_Nxx;
		vector<arma::cx_mat> m_Nxy;
		vector<arma::cx_mat> m_Nyy;
		double m_xLen;
		double m_yLen;
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
