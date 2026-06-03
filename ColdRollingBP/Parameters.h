#pragma once

#include <vector>
#include <string>
#include <map>
#include <algorithm>

struct Inputs {
	int coilInPool;
	int target_512;
	int target_008;
	int planUp;
	int planLower;
	double alpha_1;
	double alpha_2;
};

class Parameters
{
public:
	Parameters() = default;

	Parameters(Inputs inputs) {
		//set 
		inputSet = inputs;
		setValues();
	}

	~Parameters() = default;

	void setMaxPlans(int _K) {
		K = _K;
	}

	void setValues() {
		N = 4 * inputSet.coilInPool;
		S = inputSet.coilInPool;

		flow_mass = {
			std::pair<std::string, double>("C512", inputSet.target_512),
			std::pair<std::string, double>("C008", inputSet.target_008),
		};
		planUp = inputSet.planUp;
		planLower = inputSet.planLower;
		MassFlow[0] = inputSet.target_512;
		MassFlow[1] = inputSet.target_008;

		alpha1 = inputSet.alpha_1;
		alpha2 = inputSet.alpha_2;
	}

	Inputs inputSet;

	int C = 4;
	int L = 2;      
	int K = 0;       
	int N = 4 * inputSet.coilInPool;
	int S = inputSet.coilInPool;

	// the connection relationship of steel mark
	int steelArr[4][4] =
	{
			{ 1, 1, 1, 0 },
			{ 1, 1, 1, 0 },
			{ 1, 1, 1, 1 },
			{ 0, 0, 1, 1 }
	};
	// for obtaining the connection pool
	std::map<std::string, int> chsType = {
			{ "C502016", 0 },
			{ "C502017", 1 },
			{ "C502018", 2 },
			{ "C502019", 3 },
	};

	/**
	 * @brief for math model and master model, the first for C512
	*/
	int MassFlow[2] = { inputSet.target_512, inputSet.target_008 };

	// for sub-problem to get weight of a target flow
	std::vector<std::string> target_flow = { "C512", "C008" };

	/**
	 * @brief for local search
	*/
	std::map<std::string, double> flow_mass = {
			std::pair<std::string, double>("C512", inputSet.target_512),
			std::pair<std::string, double>("C008", inputSet.target_008),
	};

	std::vector<std::string> chs_type = { "C502016", "C502017", "C502018", "C502019" };

	std::vector<std::string> chs_connect = { "C502016_C502017", "C502017_C502018", "C502018_C502016",
											 "C502018_C502019" };

	/**
	 * @brief for local search and math model
	*/
	double planUp = inputSet.planUp;
	double planLower = inputSet.planLower;


	double setUpCost = 500;
	// the coefficient of objective
	double alpha1 = 1;
	double alpha2 = 0.5;

	//parameters of solver
	bool isShowLog = false;
	int _timeLimit = 3600;
	std::string _logFile;

	double epsilon = 1e-4;
};

namespace Config {
	constexpr int MAX_COILS = 100;
	constexpr int BUCKET_SIZE = 10;
}
