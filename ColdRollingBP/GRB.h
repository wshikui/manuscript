#pragma once
#include "gurobi_c++.h"
#include "Problem.h"
#include "Column.h"
#include <vector>

class SolverModel
{
public:
	SolverModel(const Parameters& parameters, std::vector<Material*>& materials,
		GRBEnv& env, GRBModel& model);

	void optimize(GRBModel& model);
	void print_result(GRBModel& model);
	void print_variables(GRBModel& model, Problem* _problem, Parameters parameters);

private:

	//GRBEnv env;
	//GRBModel model;

	std::vector<std::vector<std::vector<GRBVar>>> X;
	std::vector<std::vector<GRBVar>> Y;
	std::vector<std::vector<GRBVar>> Z;
	std::vector<std::vector<GRBVar>> U;

	int N{};
	int C{};
	int L{};
	int K{};
	int S{};

	std::vector<double> massA;
	int** flowArr{};
	double** pena{};
	int MassFlow[2] = { 200, 150 };
	int set_up = 500;
	[[maybe_unused]] int INFIN = 10000;

	void get_all_mass(const std::vector<Material*>& coils);
	int** getFlowArr(const std::vector<Material*>& coils) const;
	double** getPena(const std::vector<Material*>& coils);

	void printFlowArr(int** _flowArr) const;
	void printPena(double** _pena) const;
	void printMassAll();


	static std::map<std::string, std::vector<Material*>> selectChs(const std::vector<Material*>& coils, int m, const Parameters& parameters);

	//
	void initializeVariables(GRBModel& model);
	void initiallizeCons(GRBModel& model, std::vector<Material*>& materials, const Parameters& parameters);
	void init_model(const Parameters& parameters, std::vector<Material*>& materials, GRBModel& model);
};