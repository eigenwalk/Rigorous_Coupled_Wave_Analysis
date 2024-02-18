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
    T num; if (typeid(T) == typeid(int)) num = stoi(val);
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
		m_light->m_Kx.print("Kx");
		m_light->m_Ky.print("Ky");
        //toeplitzMatrix(i, j, k);
        //scatteringMatrix(i, j, k);	
		//calcReflection(i, j, k);
      }
    }
  }
}



RCWA::RCWA(){
}

RCWA::~RCWA(){
	delete m_light;
}

