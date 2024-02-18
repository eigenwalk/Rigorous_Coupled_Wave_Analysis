#ifndef INPUTS_H
#define INPUTS_H

#include <string>
#include <vector>
#include <set>

using namespace std;

struct Input{
	// To be deleted
	vector<string> stack;
	vector<float> thick;
	string substrate;
	//

	map<int, string> mater;
	int nx;
	int ny;
	int nz;
	double dxy;
	double dz;
	vector<double> dz_v;
	

	string voxlib;
	string voxel;
	string voxel_fullpath;
	string nklib;
	vector<float> wavelength;
	vector<float> polar_angles;
	vector<float> azi_angles;
	int num_mater;
	float pol;
	float dop;
	int hx;
	int hy;
};

#endif
