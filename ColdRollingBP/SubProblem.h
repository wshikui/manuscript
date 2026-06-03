#pragma once
#include "BaseModel.h"
#include "MasterProblem.h"
#include "Labelling.h"
#include <unordered_set>
#include <bitset>
#include <algorithm>

class SubProblem : public BaseModel
{
public:
	SubProblem() = default;
	~SubProblem() override = default;
	
	SubProblem(MasterProblem* _masterProblem, std::vector<int>& vertexWithZero,
		std::map<std::string, std::vector<int>>& _poolIndex,
		std::map<std::string, std::vector<std::bitset<Config::MAX_COILS>>>& _pool_coil_ngsets,
		Problem* _problem, Parameters* _parameters);

	SubProblem(MasterProblem* _masterProblem, std::vector<int>& vertexWithZero, 
		std::map<int, std::vector<int>> _forbidCoils,
		const std::vector<std::vector<int>>& _subset,
		std::map<std::string, std::vector<int>>& _poolIndex,
		std::map<std::string, std::vector<std::bitset<Config::MAX_COILS>>>& _pool_coil_ngsets,
		Problem* _problem, Parameters* _parameters);


	std::vector<Column>& extractColumns();
	std::vector<Column>& extractColumnsLabel();

	std::vector<Column>& extractColumnsSR();
	std::vector<Column>& extractColumnsLabelSR();

	void updateSubset(std::vector<int> s);

	static int countContainedElements(const std::vector<int>& mainVec, const std::vector<int>& searchVec);

	int solveNum = 0;
	void reSetSubNum();

private:
	// get dual variables from Master problem and get coil attribute
	MasterProblem* rmp{};
	Problem* problem{};
	Parameters* parameters{};

	//store the new columns form the sub-problem
	std::vector<Column> columns;

	std::map<std::string, std::vector<int>> poolIndex;
	std::map<std::string, std::vector<std::bitset<Config::MAX_COILS>>> pool_coil_ngsets;

	std::map<int, std::vector<int>> forbidCoils;
	std::vector<int> vertexZero{};
	std::vector<std::vector<int>> subset{};

	std::vector<std::string> conPoolNames = { "C502016_17", "C502017_18", "C502016_18","C502018_19" };

	/**
	 * get the weight of target flow in this plan (no longer needed)
	 * @param index coil index vector
	 * @param l the l-th target flow
	 * @return
	 */
	double getPlanWeight(std::vector<int>& index, int l);

	/**
	 * split all the coils into single or connect pool
	 * @param _poolIndex map -> key: pool name; value: coil index in this pool
	 */
	void splitPoolIndex(std::map<std::string, std::vector<int>>& _poolIndex);

	/**
	 * add vertex into column according to the sequence of cplex result
	 * @param _cplex
	 * @param x
	 * @param _column
	 * @param _index
	 */
	void getRoute(IloCplex& _cplex, IloArray<IloNumVarArray>& x, Column& _column, const std::vector<int>& _index, IloNumVarArray& y);


	/**
	 * @brief add vertex into column and calculate the reduced cost and objective value
	 * @param _cplex
	 * @param x
	 * @param _column
	 * @param _index
	 * @param y
	 * @param _subset
	 */
	void getRouteSR(IloCplex& _cplex, IloArray<IloNumVarArray>& x, Column& _column, const std::vector<int>& _index, IloNumVarArray& y,
		std::vector<std::vector<int>>& _subset);

	/**
	 * solve_lp the sub-problem according to this pool
	 * @param _name the pool name
	 * @param index the coil index in this pool
	 */
	bool solverPool(const std::string& _name, std::vector<int>& index);
	bool solverPoolLabel(const std::string& _name, std::vector<int>& index);

	bool solverPoolLabelSR(const std::string& _name, std::vector<int>& index);
	bool solverPoolSR(const std::string& _name, std::vector<std::vector<int>>& _subset, std::vector<int>& index);

	bool solverPoolLabelSR_ng_routes(const std::string& _name, std::vector<int>& index);
};

