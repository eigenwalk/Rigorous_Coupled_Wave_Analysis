#ifndef JOB_MANAGER_H
#define JOB_MANAGER_H

#define YAML_CPP_API 

#include <iostream>
#include <math.h>
#include <string>
#include <yaml-cpp/yaml.h>
#include "optics/rcwaRun.h"
#include "structures/inputStructure.h"

namespace fs = std::filesystem;
using namespace std;

class JobManager
{
    private:

    public:
        //JobManager(int id, string model){
        JobManager(std::string inputfile){
            m_inputfile = inputfile;
        }
		~JobManager();
        auto init()->void;
        auto start()->void;
        auto makeReport()->void;
        auto printModelInfo()->void;
		auto inputmanager(Inputs& inpts, const YAML::Node& in)->bool;
		auto structureUpdate(Inputs& inpts, const YAML::Node& in)->bool;
		auto opticsUpdate(Inputs& inpts, const YAML::Node& in)->void;
	public:
		RCWA* m_rcwa;
		std::string m_inputfile;
		
};

#endif

