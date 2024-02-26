#include <string>
#include <filesystem>
#include "jobManager.h"

auto JobManager::init()->void
{
	cout << "[INFO] JobManager initialization..." << endl;
    YAML::Node config = YAML::LoadFile(m_inputfile);
	Input input;
	if (inputmanager(input, config["Configuration"])){
		m_rcwa = new RCWA(input);
	}
	else{
		cout << "[Error] Input size is different.." << endl;
	}
}

auto JobManager::start()->void
{
    cout << "[INFO] JobManager started.." << endl;
	m_rcwa->init();
	m_rcwa->start();
	//m_rcwa->saveResult();
}

auto JobManager::makeReport()->void
{
    cout << "[INFO] Optics Transfer Matrix Method Simulation Done Successfully.." << endl;
}

auto JobManager::printModelInfo()->void
{
}

auto JobManager::inputmanager(Input& input, const YAML::Node& in)->bool
{
	if (structureUpdate(input, in))
	{
		opticsUpdate(input, in);
		return true;
	}
	return false;
}

auto JobManager::structureUpdate(Input& input, const YAML::Node& in)->bool
{
	fs::path cur_path = fs::current_path();
	input.voxlib = in["Structure"]["folder"].as<std::string>();
	input.voxel = in["Structure"]["voxel_file"].as<std::string>();
	input.nvmode = in["Structure"]["NV_method"].as<bool>();
	input.voxel_fullpath = cur_path.generic_string() + "/" + input.voxlib + "/" + input.voxel;

	//if (input.stack.size() != input.thick.size()){
	if (!fs::exists(input.voxel_fullpath)){
		cout << "[Error] No Voxel file found in the folder: " << input.voxel_fullpath << endl;
		return false;
	}
	return true;
}

auto JobManager::opticsUpdate(Input& input, const YAML::Node& in)->void
{
	input.nklib = in["Optics"]["nklib"].as<std::string>();
	input.wavelength = in["Optics"]["wavelength"].as<vector<float>>();
	input.polar_angles = in["Optics"]["polar_angle"].as<vector<float>>();
	input.azi_angles = in["Optics"]["azimuthal_angle"].as<vector<float>>();
	input.pol = in["Optics"]["polarization"].as<float>();
	input.hx = in["Optics"]["harmony_x"].as<int>();
	input.hy = in["Optics"]["harmony_y"].as<int>();
	//input.dop = in["Optics"]["dop"].as<float>();
}


JobManager::~JobManager()
{
	delete m_rcwa;
}
