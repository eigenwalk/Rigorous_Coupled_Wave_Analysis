#include "light.h" 
	
#define PI 3.1415926535897
#define RAD PI/180

Light::Light(Input& input){
  m_input = &input;
  m_ks = vec(3);
}

auto Light::setKvector(const double& wav, const double& polar, const double& azi)->void
{
  m_2nhx = 2 * m_input->hx + 1;
  m_2nhy = 2 * m_input->hy + 1;
  m_nHxy = m_2nhx * m_2nhy;
  m_ivec = linspace<vec>(-m_input->hx, m_input->hx, m_2nhx);
  m_jvec = linspace<vec>(-m_input->hy, m_input->hy, m_2nhy);
  m_Lx = m_input->dxy * m_input->nx;
  m_Ly = m_input->dxy * m_input->ny;


  // m_k0 to be used in (z .= k0 * z ). Because Kx/Ky is normalized with k0
  m_k0 = 2 * PI / wav;   
  m_ks(0) = sin(polar * RAD) * cos(azi * RAD);
  m_ks(1) = sin(polar * RAD) * sin(azi * RAD);
  m_ks(2) = cos(polar * RAD);

  vector<double> ny = {0, 1, 0};
  vector<double> nz = {0, 0, 1};
  m_nz = vec(nz);

  // If polar == 0 (normal incident), TE hss only Ey component 
  if (polar == 0) m_TE = vec(ny);
  else m_TE = cross(m_nz, m_ks);
  double normTE = norm(m_TE);
  m_TE = m_TE / normTE;
  m_TM = cross(m_ks, m_TE);
  double normTM = norm(m_TM);
  m_TM = m_TM / normTM;

  m_delta_inc = vec(m_nHxy, fill::zeros);
  m_delta_inc(int(m_nHxy/2)) = 1;
  m_P = m_input->pol * m_TE + (1 - m_input->pol) * m_TM;
  double normP = norm(m_P);
  m_Px_inc = m_P(0) * m_delta_inc;
  m_Py_inc = m_P(1) * m_delta_inc;
  m_Pxy_inc = join_cols(m_Px_inc, m_Py_inc);

  // m_P is The Incident Light Vector!
  m_P = m_P / normP;
  m_Kx = cx_mat(m_nHxy, m_nHxy, fill::zeros);
  m_Ky = cx_mat(m_nHxy, m_nHxy, fill::zeros); 
  m_Kz = cx_mat(m_nHxy, m_nHxy, fill::zeros); 
  m_Kz_t = cx_mat(m_nHxy, m_nHxy, fill::zeros); 
  m_Kz_r = cx_mat(m_nHxy, m_nHxy, fill::zeros); 

  for(int i = 0; i < m_2nhx; ++i){
    for(int j = 0; j < m_2nhy; ++j){
      int pos = i * m_2nhx + j;
	  m_Kx(pos, pos) = m_ks(0) - (2 * PI * m_ivec(i) / m_Lx);
	  m_Ky(pos, pos) = m_ks(1) - (2 * PI * m_jvec(j) / m_Ly);
	  m_Kz(pos, pos) = sqrt(1 - pow(m_ks(0),2) - pow(m_ks(1),2));
	  m_Kz_r(pos, pos) = sqrt(1 - pow(m_ks(0),2) - pow(m_ks(1),2));
	  m_Kz_t(pos, pos) = sqrt(1 - pow(m_ks(0),2) - pow(m_ks(1),2));
    }
  }
}

Light::Light(){
}

Light::~Light(){
}



