#pragma once

#include "SubProblem.h"
#include "StopWatch.h"
#include "BranchRule.h"
#include "SPmodel.h"

class ColumnGeneration
{
public:
	ColumnGeneration() = default;

	~ColumnGeneration();

	//for ng routes(lack branching for successive coils), without SR inequality
	ColumnGeneration(std::vector<Column>& columnPool, std::vector<int>& vertexEquality, std::vector<int>& vertexZero,
		std::map<std::string, std::vector<int>>& poolIndex, 
		std::map<std::string, std::vector<std::bitset<Config::MAX_COILS>>> pool_coil_ngsets,
		Parameters* _parameters, Problem* _problem);

	//all branching rule with SR
	ColumnGeneration(std::vector<Column>& columnPool, std::vector<int>& vertexEquality, std::vector<int>& vertexZero,
		const std::vector<std::vector<int>>& subset, 
		std::map<int, std::vector<int>> _forbidCoils,
		std::map<std::string, std::vector<int>>& poolIndex,
		std::map<std::string, std::vector<std::bitset<Config::MAX_COILS>>> pool_coil_ngsets,
		Parameters* _parameters, Problem* _problem);

	bool execute();
	bool executeAndCutR_depth(unsigned int node_number, int depth);


	void getBasicColumns(std::vector<std::pair<Column, double>>& pool);
	void setBasicColumns();

	/**
	 * @brief get all no duplicate columns in RMP
	 * @return
	*/
	void updateColumnPool(std::vector<Column>& pool);

	void getTime(double& valMP, double& valSp);
	double getSolValue();

	MasterProblem* getRMPmodel();

	//std::map<int, std::vector<std::string>> getPoolInfo() const;

	int solveSub = 0;

private:
	MasterProblem* rmp{};
	SubProblem* subProblem{};

	Problem* problem{};
	Parameters* parameters{};

	// all columns in the solved restricted master problem
	std::vector<Column> columnPoolCG;
	std::vector<std::pair<Column, double>> basicColumns;

	//TODO recode the pool name in the process
	//std::map<int, std::vector<std::string>> poolInfo;

	double timeMP{};
	double timeSP{};

	void getResult();

	// generate new column based RMP's dual variables
	bool generateColumns();
	bool generateColumnsAndSR();

	bool isSameColumn(const Column& column1, const Column& column2);
	void remove_duplicate_column();
	void printPool();

};

