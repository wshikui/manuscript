#pragma once

#include "Algorithm.h"
#include "Column.h"

class InitColumns
{
public:
	InitColumns() = default;

	~InitColumns() = default;

	InitColumns(const std::vector<Material*>& _materials, Problem* _problem, Parameters* _parameters);

	// actual production schedule
	void getSolutions();

	std::vector<Column*> getInitPool();

	bool isFeasible();

	// print initial column information and the cost
	void getInitCost();

private:
	std::vector<Material*> materials;
	std::vector<Column*> initPool;

	// materials in each pool : not use
	std::map<std::string, std::vector<Material*>> chsMaterials;

	// materials in each pool and each flow unit, this local search need
	std::map<std::string, std::map<std::string, std::vector<Material*>>> chsMaterial;

	// connection pool name -> its materials
	std::map<std::string, std::vector<Material*>> getConnectPool(const std::string& _poolName);

	// this problem : describe the all materials
	Problem* mpProblem{};

	Parameters* parameters{};

	std::map<std::string, double> currentFlowTarget;

	//-----------------private method----------------------//
	static double totalWeight(std::vector<Material*>& materialSet);

	bool isFinished(std::map<std::string, double>& _currentFlowTarget);

	bool isFinishedAll(std::map<std::string, double>& _currentFlowTarget);

	void flowPool(std::map<std::string, double>& _currentFlowTarget);

	std::string finder(std::map<std::string, double>& _currentFlowTarget);

	void join(std::string& pool, std::vector<Material*>& selected);

	std::vector<std::string> getConnectNames(std::string& _poolName);

	static void
		mergePool(std::vector<Material*>& _merge, std::vector<Material*>& _original, std::vector<Material*>& _substitute);

	Solution* solver(std::string& _poolName, std::vector<Material*>& _poolMaterials);
};


