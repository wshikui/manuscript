#pragma once
#include "Material.h"
#include "Solution.h"
#include "Parameters.h"
#include "Tool.h"
//#include <numeric>

class Problem
{
public:
	Problem() = default;

	explicit Problem(const std::vector<Material*>& materials);

	~Problem();

	int getVertices() const;

	void initialization();

	// for integer programing, real coils
	double getCost(int prev, int post) const;

	void printCostMatrix() const;

	void printCostMatrixPool(const std::string& name);

	void printWeight() const;

	// for column : print information: include the virtual vertex


	double getWeight(int i) const;

	double getOutWidth(int i) const;

	std::string getFlow(int i) const;

	//create an initial solution by out width
	Solution* createInitialSolution() const;

	bool feasible(int prev, int post);

	bool feasible_info(int prev, int post);

	//for local search algorithm
	void evaluate(Solution* solution, const std::map<std::string, double>& _currentFlowTarget);

	void printSolution(Solution* solution);

	//for label setting 
	std::map<int, std::vector<int>> forbidVertex{};

private:

	std::vector<Material*> _materials;
	int vertices = 0;

	// 类似无穷大的数
	double M = 10000;

	//钢种跳跃规则
	std::map<std::string, std::string> steel_rules;

	// flow unit target
	Parameters parameter;

	//厚度、宽度跳跃值 do not consider in or out
	double up_width_limit{};
	double down_width_limit{};
	double up_thick_limit{};
	double down_thick_limit{};

	//一维数组
	double* in_weight{};
	double* in_width{};
	double* out_width{};
	double* in_thick{};
	double* out_thick{};

	//二维数组 (两个coil之间的跳跃值)
	double** in_width_matrix{};
	double** out_width_matrix{};
	double** in_thick_matrix{};
	double** out_thick_matrix{};
	bool** steel_type_matrix{};
	double** penalty_matrix{};

	bool** feasible_matrix{};

	// material's flow
	std::vector<std::string> materialsFlow;

	// private method
	double in_width_jump(int prev, int post);

	double out_width_jump(int prev, int post);

	double in_thick_jump(int prev, int post);

	double out_thick_jump(int prev, int post);

	bool steel_type_jump(const std::string& s_prev, const std::string& s_post);

	double get_penalty(int prev, int post);

	void getForbidVertex();


};
