#include "rcwaRun.h" 
#include <omp.h> 
#include <vector> 
#include <cmath>
#include <complex>
#include <sstream>
#include <fstream>


	
auto RCWA::init()->bool
{
// Load light and make optic parameters
  double wave_s = m_input.wavelength[0];
  double wave_f = m_input.wavelength[1];
  double dwave = m_input.wavelength[2];
  double polar_angle_s = m_input.polar_angles[0];
  double polar_angle_f = m_input.polar_angles[1];
  double polar_dangle = m_input.polar_angles[2];
  double azi_angle_s = m_input.azi_angles[0];
  double azi_angle_f = m_input.azi_angles[1];
  double azi_dangle = m_input.azi_angles[2];
  int wave_n = (wave_f - wave_s) / dwave;
  int polar_angle_n = (polar_angle_f - polar_angle_s) / polar_dangle;
  int azi_angle_n = (azi_angle_f - azi_angle_s) / azi_dangle;
  m_waves = linspace<vec>(wave_s, wave_f, (wave_n+1));
  m_polar_angles = linspace<vec>(polar_angle_s, polar_angle_f, (polar_angle_n+1));
  m_azi_angles = linspace<vec>(azi_angle_s, azi_angle_f, (azi_angle_n+1));
  
  // Open NK lib folder and make dictionary
  fs::path cur_path = fs::current_path();
  std::string path = cur_path.generic_string() + "/" + m_input.nklib;
  for (const auto & files : fs::directory_iterator(path)){
	fs::path file = files.path(); 
	std::string filename = file.stem().generic_string();
	m_material[filename] = file;
  }

  // Check Structure First
  loadVoxel();
  
  // Then Check if NK is_exist
  set<std::string> mat_set;
  for (int i = 0; i < m_input.num_mater; ++i){
    string mater = m_input.mater[i];
    if(mat_set.find(mater) != mat_set.end()) continue;
    mat_set.insert(mater);
    if (m_material.find(mater) == m_material.end()){
      cout << "[Error] No Material is found in the nklib folder.. MAT: " << mater << endl;
      return false;
    }
    else{
      cout << "[INFO] Loadng NK values for material: " << mater << endl;
      loadNK(mater);
    }
  }
  return true;
}

auto RCWA::loadVoxel()->bool 
{
  cout << "[INFO] Loading structure file: " << m_input.voxel_fullpath << endl;
  ifstream vfile(m_input.voxel_fullpath);
  int idx = 0;
  int idx2 = 0;
  mat mat_lay;
  if (vfile.is_open() ) {
    while(vfile){
      std::string str;
  	  getline(vfile, str);
  	  if (idx == 0){
  	    // Pasing structure info
  	    auto mat_v = stringSplit(str); 
        m_input.nx = stoi(mat_v[1]);
        m_input.ny = stoi(mat_v[3]);
  	    m_input.dxy = stod(mat_v[5]);
  	    m_input.dz = stod(mat_v[7]);
        m_input.nz = stoi(mat_v[9]);
  	    m_stack.resize(m_input.nz);
		m_xLen = m_input.nx * m_input.dxy;
		m_yLen = m_input.ny * m_input.dxy;
  	    cout << "[INFO] Parsing voxel structure nx: " << m_input.nx << " ny: " << m_input.ny << " nz: " << m_input.nz << " dxy: " << m_input.dxy << " dz: " << m_input.dz << endl;
      }
      else if (idx == 1){
  	    // Pasing Materials index
  	    auto mat_v = stringSplit(str); 
  	    int loop = mat_v.size()/2;
  	    m_input.num_mater = loop;
  	    for(int i = 0; i < loop; ++i){
  	      m_input.mater[stoi(mat_v[2*i])] = mat_v[2*i + 1];
          cout << "[INFO] Parsing Material : " << stoi(mat_v[2*i]) << " " << m_input.mater[stoi(mat_v[2*i])] << endl;
        }
  	  }
      else if (idx == 2){
  	    if (m_input.dz < 0){
          // Adaptive dz for structure. vector<double> dz
          auto mat_v = stringSplit(str); 
          int loop = mat_v.size() - 1;
          m_input.dz_v.resize(loop);
          for(int i = 0; i < loop; ++i){
            m_input.dz_v[i] = stod(mat_v[i + 1]);
            cout << "[INFO] Parsing mesh_z info: " << i << " " << m_input.dz_v[i] << endl;
          }
          if (m_input.dz_v.size() != m_input.nz){
            cout << "[ERROR] number of layers != number of mesh_z. Simulation stopped.." << endl;
            return false;
          }
        }
        else{
          // Constant dz for all structure. double dz
          m_input.dz_v.resize(1);
  		  m_input.dz_v[0] = m_input.dz;
  	    }
  	  }
      else if (idx == 2 * idx2 + 3){
  			mat_lay = mat(m_input.nx, m_input.ny);
  			auto mat_v = stringSplit(str); 
  			cout << "[INFO] Voxel data is under construction for " << mat_v[0] << endl;
  			if (idx2 == m_input.nz) continue;
  			++idx2; 
      }
      else{
        // Voxel Data Parsing here
        auto vox_v = stonSplit<int>(str); 
        if(!checkvLength<int>(vox_v, m_input.nx)) return false;

        fillVoxel(vox_v, mat_lay, 0);
        for(int j = 1; j < m_input.ny; ++j){
          getline(vfile, str);
          vox_v = stonSplit<int>(str); 
          if(!checkvLength<int>(vox_v, m_input.nx)) return false;
          fillVoxel(vox_v, mat_lay, j);
        }
        cout << "debug: " << idx2 << endl;
        m_stack[idx2-1] = mat_lay;
      }
      ++idx;
  	}
  }
  return true;
}

template <typename T>
auto RCWA::checkvLength(vector<T>& v, int& s)->bool 
{
  if (v.size() != s){
    cout << "[ERROR] length of voxel.x != nx in the voxel file. Simulation stopped.." << endl;
    return false;
  }
  return true;
}

auto RCWA::fillVoxel(vector<int>& v, mat& M, int iy)->void 
{
  for (int ix = 0; ix < v.size(); ++ix) M(iy, ix) = v[ix];
}

auto RCWA::loadNK(std::string& mat)->void 
{
  vector<double> w0(5000);
  vector<double> n0(5000);
  vector<double> k0(5000);
  ifstream myfile(m_material[mat]);
  std::string str;
  int idx = 0;
  if (myfile.is_open() ) {
    while(myfile){
      getline(myfile, str);
      if (idx>0){
        auto mat_v = stonSplit<double>(str); 
        w0[idx-1] = mat_v[0];
        n0[idx-1] = mat_v[1];
        k0[idx-1] = mat_v[2];
      }
      ++idx;
    }
  }
  w0.resize(idx);
  n0.resize(idx);
  k0.resize(idx);
  
  // Interpolation
  vec vec_w0(w0);
  vec vec_n0(n0);
  vec vec_k0(k0);
  interp1(vec_w0, vec_n0, m_waves, m_n[mat], "linear");
  interp1(vec_w0, vec_k0, m_waves, m_k[mat], "linear");
  cx_vec N0(m_n[mat], -1*m_k[mat]);
  
  m_eps[mat] = pow(N0, 2);
  //m_eps[mat] = cx_vec(pow(m_n[mat], 2)-pow(m_k[mat],2), -2 * m_n[mat] * m_k[mat]);
}

auto RCWA::stringSplit(std::string& str)->vector<string> 
{
  stringstream ss(str);
  std::string val;
  vector<string> res(1000);
  int idx = 0;
  while(getline(ss, val, '\t')){
    res[idx] = val;
    ++idx;
  }
  res.resize(idx);
  return res;
}

template <typename T>
auto RCWA::stonSplit(std::string& str)->vector<T> 
{
  stringstream ss(str);
  std::string val;
  vector<T> res(5000);
  int idx = 0;
  while(getline(ss, val, '\t')){
    T num; 
	if (typeid(T) == typeid(int)) num = stoi(val);
    else if (typeid(T) == typeid(float)) num = stof(val);
    else if (typeid(T) == typeid(double)) num = stod(val);
    res[idx] = num;
    ++idx;
  }
  res.resize(idx);
  return res;
}

auto RCWA::start()->void
{
  cout << "[INFO] RCWA started over wavelength" << endl;
  for(int k = 0; k < m_azi_angles.size(); ++k){
    cout << " " << endl;
    double azi_ang = m_azi_angles[k];
    for(int j = 0; j < m_polar_angles.size(); ++j){
      cout << "Polar Angles: " << m_polar_angles[j] << ", Azimuthal Angles: " << m_azi_angles[k] << endl;
      // For R and T
      double polar_ang = m_polar_angles[j];
      m_Rs[azi_ang][polar_ang].resize(m_waves.size());
      m_Rp[azi_ang][polar_ang].resize(m_waves.size());
      m_R[azi_ang][polar_ang].resize(m_waves.size());

      //cout << "Wave[nm] Rs Rp R(total) Ts Tp T(Total) As Ap A(Total)" << endl;
      //#pragma omp parallel for
      for (int i = 0; i < m_waves.size(); ++i){
        double wav = m_waves[i];
        m_light = new Light(m_input);
		m_light->setKvector(wav, polar_ang, azi_ang);
		//m_light->m_Kx.print("Kx");
		//m_light->m_Ky.print("Ky");

		setToeplitzMatrix(i);
        scatteringMatrix(i, j, k);	
        //redhefferProduct(i, j, k);	
		//calcReflection(i, j, k);
      }
    }
  }
}


auto RCWA::scatteringMatrix(int& i, int& j, int& k)->void
{
  cx_mat Eigvec_W_air;
  cx_mat Eigvec_V_air;
  cx_mat Eigval_air;
  cx_mat A_air;
  cx_mat B_air;
  cx_mat S11_air;
  cx_mat S12_air;
  cx_mat S21_air;
  cx_mat S22_air;
  getEigen(Eigvec_W_air, Eigvec_V_air, Eigval_air, "air", 0);
  getAB(A_air, B_air, Eigvec_W_air, Eigvec_V_air);
  getSmat_r(S11_air, S12_air, S21_air, S22_air, A_air, B_air);

  // Layers Loop!
  //for ()
}

auto RCWA::getSmat_r(cx_mat& S11, cx_mat& S12, cx_mat& S21, cx_mat& S22, cx_mat& A, cx_mat& B)->void
{
  S11 = -inv(A) * B;
  S12 = 2 * inv(A);
  S21 = 0.5 * (A - B * inv(A) * B);
  S22 = B * inv(A);
}

auto RCWA::getSmat_t(cx_mat& S11, cx_mat& S12, cx_mat& S21, cx_mat& S22, cx_mat& A, cx_mat& B)->void
{
  S11 = inv(B) * A;
  S12 = 0.5 * (A - B * inv(A) * B);
  S22 = 2 * inv(A);
  S22 = -inv(A) * B;
}

auto RCWA::getSmat_l(cx_mat& S11, cx_mat& S12, cx_mat& S21, cx_mat& S22, cx_mat& A, cx_mat& B, cx_mat& X)->void
{
  cx_mat XBA_X = X * B * inv(A) * X;
  cx_mat XBA_XB = XBA_X * B;
  cx_mat XBA_XA = XBA_X * A;
  cx_mat C = A - XBA_XB;
  S11 = inv(C) * XBA_XA - B;
  S12 = (inv(C) * X) * (A - B * inv(A) * B);
  S21 = S12;
  S22 = S11;
}

auto RCWA::getAB(cx_mat& A, cx_mat& B, cx_mat& Wi, cx_mat& Vi)->void
{
  A = inv(Wi) * Wi + inv(Vi) * Vi;
  B = inv(Wi) * Wi - inv(Vi) * Vi;
  setZero(A, 1e-10);
  setZero(B, 1e-10);
}

auto RCWA::getABX(cx_mat& A, cx_mat& B, cx_mat& X, cx_mat& Wi, cx_mat& Vi, cx_mat& W_air, cx_mat& V_air, cx_mat& Eval, double& k0, double& thk)->void
{
  A = inv(Wi) * Wi + inv(Vi) * Vi;
  B = inv(Wi) * Wi - inv(Vi) * Vi;
  X = expmat(-Eval * m_light->m_k0 * thk);
  setZero(A, 1e-10);
  setZero(B, 1e-10);
  setZero(X, 1e-10);
}

auto RCWA::getEigen(cx_mat& W, cx_mat& V, cx_mat& Eval, string mode, int idx)->void
{
  int matxy = m_light->m_nHxy;
  cx_mat I = cx_mat(matxy, matxy, fill::eye);
  cx_mat P = cx_mat(2 * matxy, 2 * matxy, fill::zeros);
  cx_mat Q = cx_mat(2 * matxy, 2 * matxy, fill::zeros);
  cx_mat Er_toep   = cx_mat(matxy, matxy, fill::eye);
  cx_mat Er_toep_i = inv(Er_toep);
  if (mode != "air") getToeplitzMatrix(Er_toep, idx);

  P.submat(0, 0, matxy-1, matxy-1) = m_light->m_Kx * Er_toep_i * m_light->m_Ky;
  P.submat(0, matxy, matxy-1, 2*matxy-1) = I - m_light->m_Kx * Er_toep_i * m_light->m_Kx;
  P.submat(matxy, 0, 2*matxy-1, matxy-1) = m_light->m_Ky * Er_toep_i * m_light->m_Ky - 1 * I;
  P.submat(matxy, matxy, 2*matxy-1, 2*matxy-1) = -1 * m_light->m_Ky * Er_toep_i * m_light->m_Kx;

  Q.submat(0, 0, matxy-1, matxy-1) = m_light->m_Kx * m_light->m_Ky;
  Q.submat(0, matxy, matxy-1, 2*matxy-1) = Er_toep - m_light->m_Kx * m_light->m_Kx;
  Q.submat(matxy, 0, 2*matxy-1, matxy-1) = m_light->m_Ky * m_light->m_Ky - 1 * Er_toep;
  Q.submat(matxy, matxy, 2*matxy-1, 2*matxy-1) = -1 * m_light->m_Ky * m_light->m_Kx;

  cx_mat PQ = P * Q;
  setZero(PQ, 1e-10);
  cx_vec eval;
  eig_gen(eval, W, PQ);

  eval = sqrt(eval);
  Eval = diagmat(eval);
  cx_mat Eval_i = inv(Eval);
  V = Q * W * Eval_i;
}

auto RCWA::setZero(cx_mat& M, double decimal)->void
{
  mat M_abs = abs(M);
  uvec indices = find(M_abs < decimal);
  M.elem(indices).zeros();
}

auto RCWA::getToeplitzMatrix(cx_mat& Er_toep, int& idx)->void
{

}

auto RCWA::setToeplitzMatrix(int& wid)->void
{
  vector<cx_mat> voxel(m_input.nz);
  // For layers
  int matxy = m_light->m_nHxy;
  m_toep = vector<arma::cx_mat>(m_input.nz);
  m_toep_i = vector<arma::cx_mat>(m_input.nz);
  m_toep_delta = vector<arma::cx_mat>(m_input.nz);
  m_Nxx = vector<arma::cx_mat>(m_input.nz);
  m_Nxy = vector<arma::cx_mat>(m_input.nz);
  m_Nyy = vector<arma::cx_mat>(m_input.nz);
  for(int i = 0; i < m_input.nz;  ++i){
    mat Mat = m_stack[i];
	voxel[i] = cx_mat(m_input.ny, m_input.nx, fill::zeros);
	m_toep[i] = cx_mat(matxy, matxy, fill::zeros);
	m_toep_i[i] = cx_mat(matxy, matxy, fill::zeros);
	m_toep_delta[i] = cx_mat(matxy, matxy, fill::zeros);
	m_Nxx[i] = cx_mat(matxy, matxy, fill::zeros);
	m_Nxy[i] = cx_mat(matxy, matxy, fill::zeros);
	m_Nyy[i] = cx_mat(matxy, matxy, fill::zeros);
	for(int iy = 0; iy < m_input.ny; ++iy){
      for(int ix = 0; ix < m_input.nx; ++ix){
          int matidx = Mat(iy, ix);
		  string mat = m_input.mater[matidx];
		  complex<double> eps = m_eps[mat](wid); 
		  voxel[i](iy, ix) = eps; 
	  }
	}

	auto N = setNormalVectorField(voxel[i], i);
	//Nx = N["Nx"]
	//Ny = N["Ny"]

    cx_mat eps_fft = fft2(voxel[i])/(voxel[i].n_rows * voxel[i].n_cols);
    cx_mat eps_fft_i = fft2(1/voxel[i])/(voxel[i].n_rows * voxel[i].n_cols);
    for(int i1 = 0; i1 < m_light->m_2nhx; ++i1){
      for(int j1 = 0; j1 < m_light->m_2nhy; ++j1){
        int I1 = i1 * m_light->m_2nhy + j1;
        for(int i2 = 0; i2 < m_light->m_2nhx; ++i2){
          for(int j2 = 0; j2 < m_light->m_2nhy; ++j2){
          int J1 = i2 * m_light->m_2nhy + j2;
		  int I2 = (voxel[i].n_cols + i1 - i2)%voxel[i].n_cols;
		  int J2 = (voxel[i].n_rows + j1 - j2)%voxel[i].n_rows;
		  m_toep[i](J1, I1) = eps_fft(J2, I2);
		  m_toep_i[i](J1, I1) = eps_fft_i(J2, I2);
		  //m_Nxx[i](J1, I1) = nxx(J2, I2);
		  //m_Nxy[i](J1, I1) = nxy(J2, I2);
		  //m_Nyy[i](J1, I1) = nyy(J2, I2);
		  }
		}
	  }
	}
	// Inverse Rule
	m_toep_i[i] = inv(m_toep_i[i]);
    m_toep_delta[i] = m_toep[i] - m_toep_i[i];
	saveArma(voxel[i], "eps", i);
	saveArma(eps_fft, "eps_fft", i);
	saveArma(m_toep[i], "eps_toep", i);
	saveArma(m_toep_i[i], "eps_toep_i", i);
	saveArma(m_toep_delta[i], "eps_toep_delta", i);

	//eps_xx = toep - delta * Nxx
	//eps_xy = - delta * Nxy
	//eps_yx = - delta * Nxy
	//eps_yy = toep - delta * Nyy
	//eps_zz = toep
  }	
}

auto RCWA::setNormalVectorField(cx_mat& M, int i)->map<string, mat>
{
  // Normal vector method for the RCWA with automated vector field generation
  // by Peter Gotz et al. (2008) Optical Express

  map<string, mat> Norm;
  mat R = abs(M);
  int y_num = R.n_rows;
  int x_num = R.n_cols;

  if(R.max() == R.min()){
    cout << "R is constant" << endl;
	mat nx(y_num, x_num, fill::zeros);
	mat ny(y_num, x_num, fill::zeros);
	Norm["Nx"] = nx;
	Norm["Ny"] = ny;
    return Norm;
  }

  vec vx = linspace<vec>(-1, 1, (x_num));
  vec gaus1 = normpdf(vx, 0, 1.0/(x_num));
  double mvx = max(gaus1);
  gaus1 = gaus1/mvx;

  //vec gaus2 = vx % gaus1;
  vec gaus2(x_num-1, fill::zeros);
  for (int i = 0; i < x_num; ++i){
    gaus2 = gaus1 % vx / abs(vx);
  }
  gaus2.replace(datum::nan, 0);

  mat gcon  = circ_toeplitz(gaus1).t();
  mat gcon2 = circ_toeplitz(gaus2).t();

  gcon = censhift2d<mat>(gcon);
  gcon2 = censhift2d<mat>(gcon2);

  mat nx = gcon * R * gcon2.t();
  mat ny = gcon2 * R * gcon;

//  mat Anxy = sqrt(nx%nx + ny%ny);
//  nx /= Anxy; 
//  ny /= Anxy; 
//  nx.replace(datum::nan, 0);
//  ny.replace(datum::nan, 0);

  // save normal vector at boundaries
  map<pair<int, int>, pair<double, double>> nvb;
  findNVb(nvb, nx, ny);
  updateNV(nvb, nx, ny);
  normalizeNV(nx, ny);

  Norm["Nx"] = nx;
  Norm["Ny"] = ny;

  //saveArma(M, "gcon_M", i);
  //saveArma(gcon, "gcon", i);
  //saveArma(gcon2, "gcon2", i);
  saveArma(nx, "normal_x", i);
  saveArma(ny, "normal_y", i);

  return Norm;
}

auto RCWA::normalizeNV(mat& Nx, mat& Ny)->void
{
  
  int rows = Nx.n_rows;
  int cols = Nx.n_cols;
  for(int j = 0; j < rows; ++j){
    for(int i = 0; i < cols; ++i){
      double ds = sqrt(pow(Nx(j, i), 2) + pow(Ny(j, i), 2));
	  Nx(j, i) /= ds;
	  Ny(j, i) /= ds;
	}
  }
}

auto RCWA::updateNV(map<pair<int, int>, pair<double, double>>& nvb, mat& Nx, mat& Ny)->void
{
  int rows = Nx.n_rows;
  int cols = Nx.n_cols;
  for(auto& mp: nvb){
    auto key = mp.first;
    auto val = mp.second;
    int ix = key.first;
    int iy = key.second;
    double nx = val.first;
    double ny = val.second;
    for(int j = 0; j < rows; ++j){
      for(int i = 0; i < cols; ++i){
  	    if (i == ix and j == iy) continue;
          double rxy = 1.0 * pow(i - ix, 2) + pow(j - iy, 2);
          double nxx = (nx / rxy);
          double nyy = (ny / rxy);
		  //cout << "rxy, nx, ny, nxx, nyy: " << rxy << " " << nx << " " << ny << " " << nxx << " " << nyy << endl;
  	      Nx(j, i) += nxx;
      	  Ny(j, i) += nyy;
        }
	}
  }
}



auto RCWA::findNVb(map<pair<int, int>, pair<double, double>>& nvb, mat& Nx, mat& Ny)->void
{
  int rows = Nx.n_rows;
  int cols = Nx.n_cols;
  int num_nvc = 0;
  for(int j = 0; j < rows; ++j){
    for(int i = 0; i < cols; ++i){
      double nx = Nx(j, i);
      double ny = Ny(j, i);
	  if (abs(nx) > 0.005 or abs(ny) > 0.005){
        pair ixy(i, j);
		pair nxy(nx, ny);
        nvb[ixy] = nxy;
		num_nvc++;
	  }
	}
  } 
}

auto RCWA::fftshift2d(const cx_mat& X)->cx_mat
{
  int rows = X.n_rows;
  int cols = X.n_cols;
  int midRow = rows / 2;
  int midCol = cols / 2;

  cx_mat Res(rows, cols, fill::zeros);
  Res.submat(0, 0, midRow - 1, midCol - 1) = X.submat(midRow, midCol, rows - 1, cols - 1);
  Res.submat(midRow, 0, rows - 1, midCol - 1) = X.submat(0, midCol, midRow - 1, cols - 1);
  Res.submat(0, midCol, midRow - 1, cols - 1) = X.submat(midRow, 0, rows - 1, midCol - 1);
  Res.submat(midRow, midCol, rows - 1, cols - 1) = X.submat(0, 0, midRow - 1, midCol - 1);
  return Res;
}

template <typename T>
auto RCWA::censhift2d(const T& X)->T
{
  int rows = X.n_rows;
  int cols = X.n_cols;
  int midRow = rows / 2;
  int midCol = cols / 2;

  T Res(rows, cols, fill::zeros);
  Res.submat(0, 0, midRow - 1, midCol - 1) = X.submat(0, midCol, midRow - 1, cols - 1);
  Res.submat(midRow, 0, rows - 1, midCol - 1) = X.submat(midRow, midCol, rows - 1, cols - 1);
  Res.submat(0, midCol, midRow - 1, cols - 1) = X.submat(0, 0, midRow - 1, midCol - 1);
  Res.submat(midRow, midCol, rows - 1, cols - 1) = X.submat(midRow, 0, rows - 1, midCol - 1);
  
  return Res;
}

auto RCWA::saveArma(cx_mat& M, string file, int id)->void
{
  mat A = real(M);
  mat B = imag(M);
  saveArma(A, file+"_re", id);
  saveArma(B, file+"_im", id);
}

auto RCWA::saveArma(mat& M, string file, int id)->void
{
  int nr = M.n_rows;
  int nc = M.n_cols;
  string file1 = file + "_" + to_string(id) + ".txt";
  fs::path cur_path = fs::current_path();
  string path = cur_path.generic_string() + "/tmp/";
  if (!fs::is_directory(path)) fs::create_directory(path);
  string outfile1 = path + file1; 
  cout << "[INFO] Saving tmp results.. file: "  << endl;
  ofstream newFile1;
  newFile1.open(outfile1);
  for(int r = 0; r < nr; ++r){
    for(int c = 0; c < nc; ++c){
 	  newFile1 << real(M(r, c));
  	  if (c == nc-1) continue;
      else{
	    newFile1 << '\t';
	  }
    }
    if (r == nr-1) continue;
    else{
	  newFile1 << '\n';
    }
  }
  newFile1.close();
}


RCWA::RCWA(){
}

RCWA::~RCWA(){
	delete m_light;
}

