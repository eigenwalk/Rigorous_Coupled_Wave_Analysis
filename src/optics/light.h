#ifndef LIGHT_H
#define LIGHT_H

#include <vector>
#include <armadillo> 
#include "structures/inputStructure.h"

using namespace std;
using namespace arma;

class Light 
{
  private:

  public:
	  Light();
	  Light(Input& input);
	  ~Light();
	  auto setKvector(const double& wav, const double& polar, const double& azi)->void;

  public:
	  Input* m_input;
	  int m_2nhx;
	  int m_2nhy;
	  int m_nHxy;
	  double m_Lx;
	  double m_Ly;
	  vec m_ivec;
	  vec m_jvec;
	  double m_k0;
	  vec m_ks;
	  vec m_nz;
	  vec m_TE;
	  vec m_TM;
	  vec m_P; 
	  vec m_Px_inc; 
	  vec m_Py_inc; 
	  vec m_Pxy_inc; 
	  vec m_delta_inc;
	  cx_vec m_kx;
	  cx_vec m_ky;
	  cx_vec m_kz;
	  cx_mat m_Kx;
	  cx_mat m_Ky;
	  cx_mat m_Kz;
	  cx_mat m_Kz_t;
	  cx_mat m_Kz_r;
	  
};

#endif
