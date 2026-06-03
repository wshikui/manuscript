#include "GRB.h"


SolverModel::SolverModel(const Parameters& parameters, std::vector<Material*>& materials, GRBEnv& env, GRBModel& model)
{
	C = parameters.C;
	L = parameters.L;
	K = parameters.K;
	N = parameters.N;
	S = parameters.S;

	init_model(parameters, materials, model);
}

void SolverModel::get_all_mass(const std::vector<Material*>& coils)
{
	massA.emplace_back(0);
	for (const auto& item : coils)
	{
		massA.emplace_back(item->getValue(Attribute::IN_MAT_WT));
	}
}

int** SolverModel::getFlowArr(const std::vector<Material*>& coils) const
{
	int** _flowArr = new int* [N + 1];
	for (int i = 0; i < N + 1; i++)
	{
		_flowArr[i] = new int[L];
	}
	for (int i = 0; i < L; i++)
	{
		_flowArr[0][i] = 0;
	}
	for (int i = 1; i < N + 1; i++)
	{
		if (coils[i - 1]->getStr(Attribute::FLOW) == "C512")
		{
			_flowArr[i][0] = 1;
			_flowArr[i][1] = 0;
		}
		else
		{
			_flowArr[i][0] = 0;
			_flowArr[i][1] = 1;
		}
	}
	return _flowArr;
}

double** SolverModel::getPena(const std::vector<Material*>& coils)
{
	auto* problem = new Problem(coils);

	pena = new double* [N + 1];
	for (int i = 0; i < N + 1; i++)
	{
		pena[i] = new double[N + 1];
		if (i == 0)
		{
			for (int j = 0; j < N + 1; j++)
			{
				if (i == j)
				{
					pena[i][j] = INFIN;
				}
				else
				{
					pena[i][j] = 0;
				}
			}
		}
		else
		{
			pena[i][0] = 0;
			for (int j = 0; j < N; j++)
			{
				pena[i][j + 1] = problem->getCost(i - 1, j);
			}
		}
	}

	for (int i = 0; i < N + 1; i++)
	{
		for (int j = 0; j < N + 1; j++)
		{
			if (pena[i][j] > INFIN)
			{
				pena[i][j] = INFIN;
			}
		}
	}
	delete problem;

	return pena;
}

void SolverModel::printFlowArr(int** _flowArr) const
{
	int width = 8;
	std::cout << "_____coils flow unit_____" << std::endl;
	std::cout << std::setw(width) << std::left << "index" << "|"
		<< std::setw(width) << std::left << "C512" << "|"
		<< std::setw(width) << std::left << "C008" << "|" << std::endl;
	for (int i = 0; i < N + 1; i++)
	{
		std::cout << std::setw(width) << std::left << i << "|";
		for (int j = 0; j < L; j++)
		{
			std::cout << std::setw(width) << std::left << _flowArr[i][j] << "|";
		}
		std::cout << std::endl;
	}
}

void SolverModel::printPena(double** _pena) const
{
	int width = 8;
	std::cout << "______Math model penalty matrix include unreal coil______" << std::endl;
	for (int i = 0; i < N + 1; i++)
	{
		for (int j = 0; j < N + 1; j++)
		{
			std::cout << std::setw(width) << std::left << _pena[i][j] << "|";
		}
		std::cout << std::endl;
	}
}

void SolverModel::printMassAll()
{
	std::cout << "_____coils weight_____" << std::endl;
	for (const auto& item : massA)
	{
		std::cout << item << std::endl;
	}
}

std::map<std::string, std::vector<Material*>> SolverModel::selectChs(const std::vector<Material*>& coils, int m, const Parameters& parameters)
{
	std::map<std::string, std::vector<Material*>> chs_connect;

	// Parameters parameters;

	int numCoilCHS = parameters.S;
	int numChs = parameters.C;
	//std::cout << "numCoilCHS:" << numCoilCHS << " numChs" << numChs << std::endl;

	std::vector<Material*> coils_chs; // m池的材料 m池索引加上池内材料数量
	if (m * numCoilCHS + numCoilCHS <= static_cast<int>(coils.size()))
	{
		coils_chs.assign(coils.begin() + m * numCoilCHS, coils.begin() + m * numCoilCHS + numCoilCHS);
	}
	else
	{
		std::cout << "select coils error" << std::endl;
	}

	//获取连接关系
	std::vector<Material*> coils_no_join;
	int* arr = new int[parameters.C];
	// int arr[] = {};
	for (int i = 0; i < parameters.C; i++)
	{
		arr[i] = parameters.steelArr[m][i];
		//std::cout << "select chs: " << m << " with " << i << ", " << arr[i] << std::endl;
	}
	for (int i = 0; i < parameters.C; i++)
	{
		if (arr[i] == 0)
		{
			//i 对应不可连接的池的索引
			std::vector<Material*> coils_temp;
			coils_temp.assign(coils.begin() + i * numCoilCHS, coils.begin() + i * numCoilCHS + numCoilCHS);
			for (const auto& item : coils_temp)
			{
				coils_no_join.emplace_back(item);
			}
		}
	}

	// itself
	chs_connect["it"] = coils_chs;
	// can not connect
	chs_connect["no"] = coils_no_join;
	delete[] arr;
	return chs_connect;
}


void SolverModel::initializeVariables(GRBModel& model)
{
	X.resize(N + 1, std::vector<std::vector<GRBVar>>(N + 1, std::vector<GRBVar>(K)));
	Y.resize(N + 1, std::vector<GRBVar>(K));
	Z.resize(C, std::vector<GRBVar>(K));
	U.resize(N + 1, std::vector<GRBVar>(K));

	for (int i = 0; i < N + 1; i++) {
		for (int j = 0; j < N + 1; j++) {
			for (int k = 0; k < K; k++) {
				X[i][j][k] = model.addVar(0, 1, 0.0, GRB_BINARY,
					"x_" + std::to_string(i) + "_" + std::to_string(j) + "_" + std::to_string(k));
			}
		}
	}
	for (int i = 0; i < N + 1; i++) {
		for (int k = 0; k < K; k++) {
			Y[i][k] = model.addVar(0, 1, 0.0, GRB_BINARY,
				"y_" + std::to_string(i) + "_" + std::to_string(k));
		}
	}

	for (int c = 0; c < C; c++) {
		for (int k = 0; k < K; k++) {
			Z[c][k] = model.addVar(0, 1, 0.0, GRB_BINARY,
				"z_" + std::to_string(c) + "_" + std::to_string(k));
		}
	}

	for (int i = 0; i < N + 1; i++) {
		for (int k = 0; k < K; k++) {
			U[i][k] = model.addVar(1, N, 0.0, GRB_INTEGER,
				"u_" + std::to_string(i) + "_" + std::to_string(k));
		}
	}
	model.update();
}

void SolverModel::initiallizeCons(GRBModel& model, std::vector<Material*>& materials, const Parameters& parameters)
{
	GRBLinExpr Cons0;
	for (int i = 0; i < N + 1; i++)
	{
		for (int k = 0; k < K; k++)
		{
			Cons0 += X[i][i][k];
		}
	}
	model.addConstr(Cons0 == 0);

	//约束1 根据连接表进行池间热卷连接约束
	for (int k = 0; k < K; k++)
	{
		for (int m = 0; m < C; m++)
		{
			std::map<std::string, std::vector<Material*>> chs_m = selectChs(materials, m, parameters);
			if (!chs_m["no"].empty())
			{
				GRBLinExpr Cons1;
				for (int i = 0; i < chs_m["it"].size(); i++)
				{
					for (int j = 0; j < chs_m["no"].size(); j++)
					{
						Cons1 += X[chs_m["it"][i]->index][chs_m["no"][j]->index][k];
					}
				}
				model.addConstr(Cons1 == 0);
			}
		}
	}
	for (int k = 0; k < K; k++)
	{
		for (int m = 0; m < C; m++)
		{
			std::map<std::string, std::vector<Material*>> chs_m = selectChs(materials, m, parameters);
			if (!chs_m["no"].empty())
			{
				GRBLinExpr Cons1_2;
				for (int i = 0; i < chs_m["no"].size(); i++)
				{
					for (int j = 0; j < chs_m["it"].size(); j++)
					{
						Cons1_2 += X[chs_m["no"][i]->index][chs_m["it"][j]->index][k];
					}
				}
				model.addConstr(Cons1_2 == 0);
			}
		}
	}

	//约束2 计划内池数量
	for (int k = 0; k < K; k++)
	{
		GRBLinExpr Cons2;
		for (int c = 0; c < C; c++)
		{
			Cons2 += Z[c][k];
		}
		model.addConstr(Cons2 <= 2);
	}
	//约束3 池不选中 其中材料不选中 不针对全部y_ik约束
	for (int k = 0; k < K; k++)
	{
		for (int c = 0; c < C; c++)
		{
			std::map<std::string, std::vector<Material*>> chs_c = selectChs(materials, c, parameters);
			GRBLinExpr Cons3;
			for (int i = 0; i < chs_c["it"].size(); i++)
			{
				Cons3 += Y[chs_c["it"][i]->index][k];
			}
			model.addConstr(Cons3 <= N * Z[c][k]);
			model.addConstr(Cons3 >= Z[c][k]);
		}
	}
	//约束4 池不选中 其中材料没有连接
	for (int k = 0; k < K; k++)
	{
		for (int c = 0; c < C; c++)
		{
			GRBLinExpr Cons4;
			std::map<std::string, std::vector<Material*>> chs_c = selectChs(materials, c, parameters);
			for (int i = 0; i < chs_c["it"].size(); i++)
			{
				for (int j = 0; j < N + 1; j++)
				{
					if (j != chs_c["it"][i]->index)
					{
						Cons4 += X[chs_c["it"][i]->index][j][k];
					}
				}
			}
			model.addConstr(Cons4 <= N * Z[c][k]);
			model.addConstr(Cons4 >= Z[c][k]);
		}
	}
	for (int k = 0; k < K; k++)
	{
		for (int c = 0; c < C; c++)
		{
			GRBLinExpr Cons4_2;
			std::map<std::string, std::vector<Material*>> chs_c = selectChs(materials, c, parameters);
			for (int j = 0; j < chs_c["it"].size(); j++)
			{
				for (int i = 0; i < N + 1; i++)
				{
					if (i != chs_c["it"][j]->index)
					{
						Cons4_2 += X[i][chs_c["it"][j]->index][k];
					}
				}
			}
			model.addConstr(Cons4_2 <= N * Z[c][k]);
			model.addConstr(Cons4_2 >= Z[c][k]);
		}
	}
	//约束5 板卷至多属于一个计划 虚拟板卷除外
	for (int i = 1; i < N + 1; i++)
	{
		GRBLinExpr Cons5;
		for (int k = 0; k < K; k++)
		{
			Cons5 += Y[i][k];
		}
		model.addConstr(Cons5 <= 1);
	}
	//约束6 流平衡和子环消除
	for (int i = 0; i < N + 1; i++)
	{
		for (int k = 0; k < K; k++)
		{
			GRBLinExpr Cons6;
			for (int j = 0; j < N + 1; j++)
			{
				if (j != i)
				{
					Cons6 += X[i][j][k];
				}
			}
			model.addConstr(Cons6 == Y[i][k]);

			GRBLinExpr Cons6_2;
			for (int j = 0; j < N + 1; j++)
			{
				if (j != i)
				{
					Cons6_2 += X[j][i][k];
				}
			}
			model.addConstr(Cons6 == Cons6_2);
			model.addConstr(Cons6_2 == Y[i][k]);
		}
	}

	for (int k = 0; k < K; k++)
	{
		for (int i = 1; i < N + 1; i++)
		{
			for (int j = 1; j < N + 1; j++)
			{
				if (j != i)
				{
					GRBLinExpr Cons6_5;
					Cons6_5 = U[i - 1][k] - U[j - 1][k] + N * X[i][j][k];
					model.addConstr(Cons6_5 <= N - 1);
				}
			}
		}
	}
	//约束7 计划容量限制
	for (int k = 0; k < K; k++)
	{
		GRBLinExpr Cons7;
		for (int i = 1; i < N + 1; i++)
		{
			Cons7 += Y[i][k] * massA[i];
		}
		model.addConstr(Cons7 <= parameters.planUp * Y[0][k]);
		model.addConstr(Cons7 >= parameters.planLower * Y[0][k]);
	}

	//约束8 目标流向
	for (int l = 0; l < L; l++)
	{
		GRBLinExpr Cons8;
		for (int i = 1; i < N + 1; i++)
		{
			for (int k = 0; k < K; k++)
			{
				Cons8 += Y[i][k] * massA[i] * flowArr[i][l];
			}
		}
		model.addConstr(Cons8 >= MassFlow[l]);
	}
	model.update();

	GRBLinExpr obj;
	// first part
	for (int k = 0; k < K; k++)
	{
		for (int i = 0; i < N + 1; i++)
		{
			for (int j = 0; j < N + 1; j++)
			{
				if (j != i)
				{
					obj += pena[i][j] * X[i][j][k];
				}
			}
		}
	}
	obj = obj * parameters.alpha1;
	for (int k = 0; k < K; k++)
	{
		for (int j = 1; j < N + 1; j++)
		{
			obj += parameters.alpha2 * X[0][j][k] * set_up;
		}
	}
	model.setObjective(obj, GRB_MINIMIZE);
	model.update();
}



void SolverModel::init_model(const Parameters& parameters, std::vector<Material*>& materials, GRBModel& model)
{
	//model.set(GRB_StringAttr_ModelName, "GRB_model_schedule");

	get_all_mass(materials);
	flowArr = getFlowArr(materials);
	//printFlowArr(flowArr);

	for (int i = 0; i < L; i++)
	{
		MassFlow[i] = parameters.MassFlow[i];
	}
	pena = getPena(materials);
	//printPena(pena);

	initializeVariables(model);

	initiallizeCons(model, materials, parameters);
}

void SolverModel::optimize(GRBModel& model)
{
	model.set(GRB_DoubleParam_TimeLimit, 3600.0);
	model.update();
	model.optimize();
}

void SolverModel::print_result(GRBModel& model)
{
	int status = model.get(GRB_IntAttr_Status);
	//std::cout << "obj_value:" << model.get(GRB_DoubleAttr_ObjVal) << std::endl;
	std::cout << "lower:" << model.get(GRB_DoubleAttr_ObjBound) << std::endl;
	std::cout << "upper:" << model.get(GRB_DoubleAttr_ObjVal) << std::endl;
	std::cout << "Gap:" << model.get(GRB_DoubleAttr_MIPGap) * 100 << "%" << std::endl;
	std::cout << "time:" << model.get(GRB_DoubleAttr_Runtime) << std::endl;
}


void SolverModel::print_variables(GRBModel& model, Problem* _problem, Parameters parameters)
{
	std::vector<Column> columns;

	//Parameters parameters;
	for (int k = 0; k < K; k++) {
		for (int c = 0; c < C; c++) {
			int zck = 0;
			if (Z[c][k].get(GRB_DoubleAttr_X) > 1e-5) {
				zck = 1;
			}
			//std::cout << "Z[" << c << "][" << k << "]" << zck << " ";
		}
		//std::cout << std::endl;
	}
	for (int k = 0; k < K; k++) {
		for (int i = 0; i < 1 + N; i++) {
			//std::cout << "Y[" << i << "][" << k << "]" << Y[i][k].get(GRB_DoubleAttr_X) << " ";
		}
	}
	//std::cout << std::endl;

	std::cout << "_____________MathModel result_____________" << std::endl;
	for (int k = 0; k < K; k++)
	{
		// get chs names
		std::string chsName;
		for (int c = 0; c < C; c++)
		{
			if (Z[c][k].get(GRB_DoubleAttr_X) > 0.9)
			{
				chsName += parameters.chs_type[c] + "_";
			}
		}
		std::vector<int> temp;
		std::cout << "plan_" << k + 1 << " -> " << chsName << std::endl;

		// get vertices in this k route
		int begin = -1;
		for (int j = 0; j < N + 1; j++)
		{
			if (X[0][j][k].get(GRB_DoubleAttr_X) > 0.99)
			{
				begin = j;
				break;
			}
		}
		if (begin == -1)
		{
			std::cout << "this " << k + 1 << "plan has no coils" << std::endl;
			continue;
		}
		std::vector<int> route = { begin };
		for (int i = 0; i < N + 1; i++)
		{
			for (int j = 0; j < N + 1; j++)
			{
				if (X[begin][j][k].get(GRB_DoubleAttr_X) > 0.99)
				{
					begin = j;
					if (begin == 0)
					{
						break;
					}
					route.emplace_back(j);
					break;
				}
			}
			if (begin == 0)
			{
				break;
			}
		}

		for (int i = 0; i < N + 1; i++)
		{
			for (int j = 0; j < N + 1; j++)
			{
				if (X[i][j][k].get(GRB_DoubleAttr_X) > 0.99)
				{
					std::cout << "X[" << i << "][" << j << "][" << k << "]" << "\t";
					temp.emplace_back(i);
				}
			}
		}
		std::cout << std::endl;

		Column column(_problem);
		column.setName(chsName);
		for (const auto& item : route)
		{
			column.addVertex(item);
		}
		column.calculateCost();
		column.printRoute();
		columns.emplace_back(column);

		std::cout << "--------------------------------------------------" << std::endl;
	}
	double cost = 0.0;
	for (const auto& item : columns)
	{
		cost += item.getCost();
	}
	cost = cost * parameters.alpha1;
	//	for (const auto& item: columns)
	//	{
	//		cost += (parameters.planUp - item.getDemand()) * 0.7;
	//	}
	cost += int(columns.size()) * set_up * parameters.alpha2;
	std::cout << "gurobi column cost -> " << cost << std::endl;
}
