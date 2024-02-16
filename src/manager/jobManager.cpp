#include <string>
#include <filesystem>
#include "jobManager.h"

auto JobManager::init()->void
{
	cout << "[INFO] JobManager initialization..." << endl;
    YAML::Node config = YAML::LoadFile(m_inputfile);
	Inputs inpts;
	if (inputmanager(inpts, config["Configuration"])){
		m_rcwa = new RCWA(inpts);
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

auto JobManager::inputmanager(Inputs& inpts, const YAML::Node& in)->bool
{
	if (structureUpdate(inpts, in))
	{
		opticsUpdate(inpts, in);
		return true;
	}
	return false;
}

auto JobManager::structureUpdate(Inputs& inpts, const YAML::Node& in)->bool
{
	fs::path cur_path = fs::current_path();
	inpts.voxlib = in["Structure"]["folder"].as<std::string>();
	inpts.voxel = in["Structure"]["voxel_file"].as<std::string>();
	inpts.voxel_fullpath = cur_path.generic_string() + "/" + inpts.voxlib + "/" + inpts.voxel;

	//if (inpts.stack.size() != inpts.thick.size()){
	if (!fs::exists(inpts.voxel_fullpath)){
		cout << "[Error] No Voxel file found in the folder: " << inpts.voxel_fullpath << endl;
		return false;
	}
	return true;
}

auto JobManager::opticsUpdate(Inputs& inpts, const YAML::Node& in)->void
{
	inpts.nklib = in["Optics"]["nklib"].as<std::string>();
	inpts.wavelength = in["Optics"]["wavelength"].as<vector<float>>();
	inpts.polar_angles = in["Optics"]["polar_angle"].as<vector<float>>();
	inpts.azi_angles = in["Optics"]["azimuthal_angle"].as<vector<float>>();
	inpts.pol = in["Optics"]["polarization"].as<float>();
	//inpts.dop = in["Optics"]["dop"].as<float>();
}


JobManager::~JobManager()
{
	delete m_rcwa;
}
