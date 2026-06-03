#pragma once

#include "BaseModel.h"
#include "Problem.h"
#include "Column.h"

class MathModel : public BaseModel
{
public:
	MathModel() = default;

	~MathModel() override;

	MathModel(const Parameters& parameters, std::vector<Material*>& materials);

	void printResult(Problem* _problem);

	void printResultR(IloCplex& cplex, IloArray<IloArray<IloNumVarArray>>& x);

	const IloArray<IloArray<IloNumVarArray>>& getValues();

	[[nodiscard]] double getObjective(const Parameters& parameters) const;

private:

	// 模型变量 此时的变量尚未和env关联起来 在initial model时 关联env 并指定数组大小 类型
	IloArray<IloArray<IloNumVarArray>> X;
	IloArray<IloNumVarArray> Z;
	IloArray<IloNumVarArray> Y;
	IloArray<IloNumVarArray> U;

	int N{};
	int C{};
	int L{};
	int K{};
	int S{};

	// 含虚拟板卷在内的全部重量
	std::vector<double> massA;
	// 含虚拟板卷在内的流向矩阵
	int** flowArr{};
	// 含虚拟板卷在内的各个板卷的cost矩阵
	double** pena{};
	// 各流向的需求量
	int MassFlow[2] = { 200, 150 };
	//开机成本
	int set_up = 500;

	[[maybe_unused]] int INFIN = 10000;

	static void getMass(const std::vector<Material*>& coils, std::vector<double>& _massAll);

	int** getFlowArr(const std::vector<Material*>& coils) const;

	double** getPena(const std::vector<Material*>& coils);

	// 打印信息
	void printFlowArr(int** _flowArr) const;

	void printPena(double** _pena) const;

	void printMassAll();

	void initModel(const Parameters& parameters, std::vector<Material*>& materials);

	void initModel_S(const Parameters& parameters, std::vector<Material*>& materials);

	static std::map<std::string, std::vector<Material*>> selectChs(const std::vector<Material*>& coils, int m, const Parameters& parameters);
};

