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
        JobManager(std::string inputfile){
            m_inputfile = inputfile;
        }
		~JobManager();
        auto init()->void;
        auto start()->void;
        auto makeReport()->void;
        auto printModelInfo()->void;
		auto inputmanager(Input& input, const YAML::Node& in)->bool;
		auto structureUpdate(Input& input, const YAML::Node& in)->bool;
		auto opticsUpdate(Input& input, const YAML::Node& in)->void;
	public:
		unique_ptr<RCWA> m_rcwa;
		std::string m_inputfile;
		
};

#endif

