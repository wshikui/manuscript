#pragma once
#include "BaseModel.h"
#include "Column.h"
#include <cmath>

class MasterProblem : public BaseModel
{
public:
	MasterProblem() = default;

	~MasterProblem() override;

	/**
	 * @brief create RMP based on initial columns, decision variable is linear
	 * @param _columns
	 * @param _problem
	 * @param _parameters
	*/
	MasterProblem(std::vector<Column*>& _columns, Problem* _problem, Parameters* _parameters);

	/**
	 * @brief create RMP based on initial columns, branch rule, decision variable is linear
	 * @param localPool
	 * @param vertexWithEquality
	 * @param _problem
	 * @param _parameters
	*/
	MasterProblem(std::vector<Column>& localPool, std::vector<int>& vertexWithEquality, Problem* _problem, Parameters* _parameters);

	// MasterProblem(std::vector<Column>& localPool, std::vector<int>& vertexWithEquality, Problem* _problem, Parameters* _parameters, int cuts);

	/**
	 * @brief create RMP based on initial columns, branch rule and variable type
	 * @param localPool
	 * @param vertexWithEquality
	 * @param _problem
	 * @param _parameters
	 * @param isInteger : is decision variable integer or linear
	*/
	MasterProblem(std::vector<Column>& localPool, std::vector<int>& vertexWithEquality, Problem* _problem, Parameters* _parameters, bool isInteger);

	/**
	 * @brief create RMP based on branch rule and SR inequality
	 * @param localPool
	 * @param vertexWithEquality
	 * @param subSet
	 * @param _problem
	 * @param _parameters
	*/
	MasterProblem(std::vector<Column>& localPool, std::vector<int>& vertexWithEquality, const std::vector<std::vector<int>>& _subSet,
		Problem* _problem, Parameters* _parameters);

	// get dual variables, CPLEX has the same function
	void getDual(IloNumArray& variables);

	// add columns
	void addColumn(const Column& _column);

	//statistical information 

	//get the number of variables in the i-th constraint
	int getVarNum(int i);

	//get the coefficient and value of decision variable including the i-th coil
	std::tuple<std::vector<double>, std::vector<double>> getCoefAndValue(int i);

	//get the weight in inequality of plan number
	double getWeightRH(int i);

	//get the number of SR including the coil
	double getNumSRofCoil(int i);

	//number of basic columns in coil i's constraint
	double getBCconstraint(int i);

	//the constraint is binding constraint
	std::tuple<double, double, double, double, double> isBinding(int i);


	/**
	 * @brief add column from sub-problem to RMP consdier SR inequality
	 * @param _column
	*/
	void addColumnSR(const Column& _column);


	double getDualVariable(int i);


	void printResult(int iter);

	IloNumVarArray& getY();
	IloRangeArray& getCons();

	//get the basic columns in RMP
	std::vector<std::pair<Column, double>>& getBasicColumns();

	std::vector<std::vector<int>>& getSubSet();

	void updateSubset(std::vector<int> s);


private:
	IloNumVarArray y;
	IloRangeArray cons;
	IloObjective obj;
	IloNumArray upLimits;
	IloNumArray downLimits;

	Problem* problem{};

	Parameters* parameters{};

	//S (maybe increase in iteration process), basic columns for solving separation problem
	std::vector<std::vector<int>> subSet;
	std::vector<Column> columnPool;
	std::vector<std::pair<Column, double>> basicColumns;

	int** h{};
	void initMatrixH();
	double getTargetMass();
	int getMaxPlanNum();

	/**
	 * get the target flow mass of the initial columns from heuristic algorithm
	 * @param _columns : initial columns vector
	 * @param r : the index of column in vector
	 * @param l : the index of target flow
	 * @return : mass
	 */
	double getInitColumnMass(std::vector<Column*>& _columns, int r, int l);
	double getInitColumnMass(std::vector<Column>& _columns, int r, int l);

	// rmp at 
	void initModel(std::vector<Column*>& _columns, Problem* _problem);
	void initMasterModel(std::vector<Column>& localPool, std::vector<int>& vertexWithEquality);
	void initMasterModel(std::vector<Column>& localPool, std::vector<int>& vertexWithEquality, bool integer);

	//initialize the RMP(variable, objective£¬ and constraint), add the rule of branching and SR inequality
	void initMasterModel(std::vector<Column>& localPool, std::vector<int>& vertexWithEquality,
		const std::vector<std::vector<int>>& subSet);

};
