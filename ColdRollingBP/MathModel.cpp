//
// Created by wshikui on 2023/9/18.
//
#include "MathModel.h"

MathModel::~MathModel()
{
	for (int i = 0; i < N + 1; i++)
	{
		delete[] flowArr[i];
		delete[] pena[i];
	}
	delete[] flowArr;
	delete[] pena;

	mEnv.end();
}

std::map<std::string, std::vector<Material*>> MathModel::selectChs(const vector<Material*>& coils, int m, const Parameters& parameters)
{
	std::map<std::string, std::vector<Material*>> chs_connect;
	// Parameters parameters;

	int numCoilCHS = parameters.S;
	int numChs = parameters.C;

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

void MathModel::getMass(const std::vector<Material*>& coils, std::vector<double>& _massAll)
{
	_massAll.emplace_back(0);
	for (const auto& item : coils)
	{
		_massAll.emplace_back(item->getValue(Attribute::IN_MAT_WT));
	}
}

// 生成包括虚拟板卷在内的各板卷流向矩阵 第一位为C512 第二位为C008
int** MathModel::getFlowArr(const std::vector<Material*>& coils) const
{
	// N + 1 * L
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

//获取包括虚拟板卷在内的全部板卷的cost矩阵
double** MathModel::getPena(const std::vector<Material*>& coils)
{
	// 新建问题类 获取真实板卷之间的cost
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

	// 是否需要统一无穷大数值
	// 统一
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

void MathModel::printFlowArr(int** _flowArr) const
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

void MathModel::printMassAll()
{
	std::cout << "_____coils weight_____" << std::endl;
	for (const auto& item : massA)
	{
		std::cout << item << std::endl;
	}
}

void MathModel::initModel(const Parameters& parameters, vector<Material*>& materials)
{
	// 对model重命名
	mModel.setName("Mathematical model");

	getMass(materials, massA);
	//printMassAll();

	flowArr = getFlowArr(materials);
	//printFlowArr(flowArr);

	for (int i = 0; i < L; i++)
	{
		MassFlow[i] = parameters.MassFlow[i];
	}
	pena = getPena(materials);
	//printPena(pena);

	try
	{
		// X
		X = IloArray<IloArray<IloNumVarArray>>(mEnv, N + 1);
		for (int i = 0; i < N + 1; i++)
		{
			X[i] = IloArray<IloNumVarArray>(mEnv, N + 1);
			for (int j = 0; j < N + 1; j++)
			{
				X[i][j] = IloNumVarArray(mEnv, K);
				for (int k = 0; k < K; k++)
				{
					X[i][j][k] = IloNumVar(mEnv, 0, 1, ILOINT);
				}
			}
		}

		// 定义Z变量
		Z = IloArray<IloNumVarArray>(mEnv, C);
		for (int c = 0; c < C; c++)
		{
			Z[c] = IloNumVarArray(mEnv, K);
			for (int k = 0; k < K; k++)
			{
				Z[c][k] = IloNumVar(mEnv, 0, 1, ILOINT);
			}
		}

		// 定义Y
		Y = IloArray<IloNumVarArray>(mEnv, N + 1);
		for (int i = 0; i < N + 1; i++)
		{
			Y[i] = IloNumVarArray(mEnv, K);
			for (int k = 0; k < K; k++)
			{
				Y[i][k] = IloNumVar(mEnv, 0, 1, ILOINT);
			}
		}

		// 定义U
		U = IloArray<IloNumVarArray>(mEnv, N + 1);
		for (int i = 0; i < N; i++)
		{
			U[i] = IloNumVarArray(mEnv, K);
			for (int k = 0; k < K; k++)
			{
				U[i][k] = IloNumVar(mEnv, 1, N, ILOINT);
			}
		}

		// 约束1
		IloExpr Cons0(mEnv);
		for (int i = 0; i < N + 1; i++)
		{
			for (int k = 0; k < K; k++)
			{
				Cons0 += X[i][i][k];
			}
		}
		mModel.add(Cons0 == 0);
		Cons0.end();

		//约束1 根据连接表进行池间热卷连接约束
		for (int k = 0; k < K; k++)
		{
			for (int m = 0; m < C; m++)
			{
				std::map<std::string, std::vector<Material*>> chs_m = selectChs(materials, m, parameters);
				if (!chs_m["no"].empty())
				{
					IloExpr Cons1(mEnv);
					for (int i = 0; i < chs_m["it"].size(); i++)
					{
						for (int j = 0; j < chs_m["no"].size(); j++)
						{
							Cons1 += X[chs_m["it"][i]->index][chs_m["no"][j]->index][k];
						}
					}
					mModel.add(Cons1 == 0);
					Cons1.end();
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
					IloExpr Cons1_2(mEnv);
					for (int i = 0; i < chs_m["no"].size(); i++)
					{
						for (int j = 0; j < chs_m["it"].size(); j++)
						{
							Cons1_2 += X[chs_m["no"][i]->index][chs_m["it"][j]->index][k];
						}
					}
					mModel.add(Cons1_2 == 0);
					Cons1_2.end();
				}
			}
		}

		//约束2 计划内池数量
		for (int k = 0; k < K; k++)
		{
			IloExpr Cons2(mEnv);
			for (int c = 0; c < C; c++)
			{
				Cons2 += Z[c][k];
			}
			mModel.add(Cons2 <= 2);
			Cons2.end();
		}

		//约束3 池不选中 其中材料不选中 不针对全部y_ik约束
		for (int k = 0; k < K; k++)
		{
			for (int c = 0; c < C; c++)
			{
				std::map<std::string, std::vector<Material*>> chs_c = selectChs(materials, c, parameters);
				IloExpr Cons3(mEnv);
				for (int i = 0; i < chs_c["it"].size(); i++)
				{
					Cons3 += Y[chs_c["it"][i]->index][k];
				}
				mModel.add(Cons3 <= N * Z[c][k]);
				mModel.add(Cons3 >= Z[c][k]);
				Cons3.end();
			}
		}

		//约束4 池不选中 其中材料没有连接
		for (int k = 0; k < K; k++)
		{
			for (int c = 0; c < C; c++)
			{
				IloExpr Cons4(mEnv);
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
				mModel.add(Cons4 <= N * Z[c][k]);
				mModel.add(Cons4 >= Z[c][k]);
				Cons4.end();
			}
		}
		for (int k = 0; k < K; k++)
		{
			for (int c = 0; c < C; c++)
			{
				IloExpr Cons4_2(mEnv);
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
				mModel.add(Cons4_2 <= N * Z[c][k]);
				mModel.add(Cons4_2 >= Z[c][k]);
				Cons4_2.end();
			}
		}

		//约束5 板卷至多属于一个计划 虚拟板卷除外
		for (int i = 1; i < N + 1; i++)
		{
			IloExpr Cons5(mEnv);
			for (int k = 0; k < K; k++)
			{
				Cons5 += Y[i][k];
			}
			mModel.add(Cons5 <= 1);
			Cons5.end();
		}

		//约束6 流平衡和子环消除
		for (int i = 0; i < N + 1; i++)
		{
			for (int k = 0; k < K; k++)
			{
				IloExpr Cons6(mEnv);
				for (int j = 0; j < N + 1; j++)
				{
					if (j != i)
					{
						Cons6 += X[i][j][k];
					}
				}
				mModel.add(Cons6 == Y[i][k]);

				IloExpr Cons6_2(mEnv);
				for (int j = 0; j < N + 1; j++)
				{
					if (j != i)
					{
						Cons6_2 += X[j][i][k];
					}
				}
				mModel.add(Cons6 == Cons6_2);
				mModel.add(Cons6_2 == Y[i][k]);
				Cons6.end();
				Cons6_2.end();
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
						IloExpr Cons6_5(mEnv);
						Cons6_5 = U[i - 1][k] - U[j - 1][k] + N * X[i][j][k];
						mModel.add(Cons6_5 <= N - 1);
						Cons6_5.end();
					}
				}
			}
		}

		//约束7 计划容量限制
		for (int k = 0; k < K; k++)
		{
			IloExpr Cons7(mEnv);
			for (int i = 1; i < N + 1; i++)
			{
				Cons7 += Y[i][k] * massA[i];
			}
			mModel.add(Cons7 <= parameters.planUp * Y[0][k]);
			mModel.add(Cons7 >= parameters.planLower * Y[0][k]);
			Cons7.end();
		}

		//约束8 目标流向
		for (int l = 0; l < L; l++)
		{
			IloExpr Cons8(mEnv);
			for (int i = 1; i < N + 1; i++)
			{
				for (int k = 0; k < K; k++)
				{
					Cons8 += Y[i][k] * massA[i] * flowArr[i][l];
				}
			}
			mModel.add(Cons8 >= MassFlow[l]);
			Cons8.end();
		}

		//目标函数
		IloExpr obj(mEnv);

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

		// second part
//		for (int k = 0; k < K; k++)
//		{
//			obj += 0.7 * parameters.planUp;
//			for (int i = 1; i < N + 1; i++)
//			{
//				obj -= 0.7 * Y[i][k] * massA[i];
//			}
//		}

		// third part
		for (int k = 0; k < K; k++)
		{
			for (int j = 1; j < N + 1; j++)
			{
				obj += parameters.alpha2 * X[0][j][k] * set_up;
			}
		}
		mModel.add(IloMinimize(mEnv, obj));
		obj.end();
	}
	catch (IloException& ex)
	{
		cerr << "InitModel Error -> " << ex << endl;
	}
	catch (...)
	{
		cerr << "Error when extractColumn" << endl;
	}

}

void MathModel::initModel_S(const Parameters& parameters, std::vector<Material*>& materials)
{
	// 对model重命名
	mModel.setName("Mathematical model");

	getMass(materials, massA);
	printMassAll();

	flowArr = getFlowArr(materials);
	printFlowArr(flowArr);

	for (int i = 0; i < L; i++)
	{
		MassFlow[i] = parameters.MassFlow[i];
	}
	pena = getPena(materials);
	printPena(pena);

	try
	{
		// X
		X = IloArray<IloArray<IloNumVarArray>>(mEnv, N + 1);
		for (int i = 0; i < N + 1; i++)
		{
			X[i] = IloArray<IloNumVarArray>(mEnv, N + 1);
			for (int j = 0; j < N + 1; j++)
			{
				X[i][j] = IloNumVarArray(mEnv, K);
				for (int k = 0; k < K; k++)
				{
					X[i][j][k] = IloNumVar(mEnv, 0, 1, ILOINT);
				}
			}
		}

		// 定义Z变量
		Z = IloArray<IloNumVarArray>(mEnv, C);
		for (int c = 0; c < C; c++)
		{
			Z[c] = IloNumVarArray(mEnv, K);
			for (int k = 0; k < K; k++)
			{
				Z[c][k] = IloNumVar(mEnv, 0, 1, ILOINT);
			}
		}

		// 定义Y
		Y = IloArray<IloNumVarArray>(mEnv, N + 1);
		for (int i = 0; i < N + 1; i++)
		{
			Y[i] = IloNumVarArray(mEnv, K);
			for (int k = 0; k < K; k++)
			{
				Y[i][k] = IloNumVar(mEnv, 0, 1, ILOINT);
			}
		}

		// 定义U
		U = IloArray<IloNumVarArray>(mEnv, N + 1);
		for (int i = 0; i < N; i++)
		{
			U[i] = IloNumVarArray(mEnv, K);
			for (int k = 0; k < K; k++)
			{
				U[i][k] = IloNumVar(mEnv, 1, N, ILOINT);
			}
		}

		// 约束1
		IloExpr Cons0(mEnv);
		for (int i = 0; i < N + 1; i++)
		{
			for (int k = 0; k < K; k++)
			{
				Cons0 += X[i][i][k];
			}
		}
		mModel.add(Cons0 == 0);
		Cons0.end();

		//约束1 根据连接表进行池间热卷连接约束
		for (int k = 0; k < K; k++)
		{
			for (int m = 0; m < C; m++)
			{
				std::map<std::string, std::vector<Material*>> chs_m = selectChs(materials, m, parameters);
				if (!chs_m["no"].empty())
				{
					IloExpr Cons1(mEnv);
					for (int i = 0; i < chs_m["it"].size(); i++)
					{
						for (int j = 0; j < chs_m["no"].size(); j++)
						{
							Cons1 += X[chs_m["it"][i]->index][chs_m["no"][j]->index][k];
						}
					}
					mModel.add(Cons1 == 0);
					Cons1.end();
				}
			}
		}

		//this constraint is redundant

		//for (int k = 0; k < K; k++)
		//{
		//	for (int m = 0; m < C; m++)
		//	{
		//		std::map<std::string, std::vector<Material*>> chs_m = selectChs(materials, m);
		//		if (!chs_m["no"].empty())
		//		{
		//			IloExpr Cons1_2(mEnv);
		//			for (int i = 0; i < chs_m["no"].size(); i++)
		//			{
		//				for (int j = 0; j < chs_m["it"].size(); j++)
		//				{
		//					Cons1_2 += X[chs_m["no"][i]->index][chs_m["it"][j]->index][k];
		//				}
		//			}
		//			mModel.add(Cons1_2 == 0);
		//			Cons1_2.end();
		//		}
		//	}
		//}

		//约束2 计划内池数量
		for (int k = 0; k < K; k++)
		{
			IloExpr Cons2(mEnv);
			for (int c = 0; c < C; c++)
			{
				Cons2 += Z[c][k];
			}
			mModel.add(Cons2 <= 2);
			Cons2.end();
		}

		//约束3 池不选中 其中材料不选中 不针对全部y_ik约束
		for (int k = 0; k < K; k++)
		{
			for (int c = 0; c < C; c++)
			{
				std::map<std::string, std::vector<Material*>> chs_c = selectChs(materials, c, parameters);
				IloExpr Cons3(mEnv);
				for (int i = 0; i < chs_c["it"].size(); i++)
				{
					Cons3 += Y[chs_c["it"][i]->index][k];
				}
				mModel.add(Cons3 <= N * Z[c][k]);
				//mModel.add(Cons3 >= Z[c][k]);
				Cons3.end();
			}
		}

		//约束4 池不选中 其中材料没有连接
		//for (int k = 0; k < K; k++)
		//{
		//	for (int c = 0; c < C; c++)
		//	{
		//		IloExpr Cons4(mEnv);
		//		std::map<std::string, std::vector<Material*>> chs_c = selectChs(materials, c);
		//		for (int i = 0; i < chs_c["it"].size(); i++)
		//		{
		//			for (int j = 0; j < N + 1; j++)
		//			{
		//				if (j != chs_c["it"][i]->index)
		//				{
		//					Cons4 += X[chs_c["it"][i]->index][j][k];
		//				}
		//			}
		//		}
		//		mModel.add(Cons4 <= N * Z[c][k]);
		//		mModel.add(Cons4 >= Z[c][k]);
		//		Cons4.end();
		//	}
		//}
		//for (int k = 0; k < K; k++)
		//{
		//	for (int c = 0; c < C; c++)
		//	{
		//		IloExpr Cons4_2(mEnv);
		//		std::map<std::string, std::vector<Material*>> chs_c = selectChs(materials, c);
		//		for (int j = 0; j < chs_c["it"].size(); j++)
		//		{
		//			for (int i = 0; i < N + 1; i++)
		//			{
		//				if (i != chs_c["it"][j]->index)
		//				{
		//					Cons4_2 += X[i][chs_c["it"][j]->index][k];
		//				}
		//			}
		//		}
		//		mModel.add(Cons4_2 <= N * Z[c][k]);
		//		mModel.add(Cons4_2 >= Z[c][k]);
		//		Cons4_2.end();
		//	}
		//}

		//约束5 板卷至多属于一个计划 虚拟板卷除外
		for (int i = 1; i < N + 1; i++)
		{
			IloExpr Cons5(mEnv);
			for (int k = 0; k < K; k++)
			{
				Cons5 += Y[i][k];
			}
			mModel.add(Cons5 <= 1);
			Cons5.end();
		}

		//约束6 流平衡和子环消除
		for (int i = 0; i < N + 1; i++)
		{
			for (int k = 0; k < K; k++)
			{
				IloExpr Cons6(mEnv);
				for (int j = 0; j < N + 1; j++)
				{
					if (j != i)
					{
						Cons6 += X[i][j][k];
					}
				}
				mModel.add(Cons6 == Y[i][k]);

				IloExpr Cons6_2(mEnv);
				for (int j = 0; j < N + 1; j++)
				{
					if (j != i)
					{
						Cons6_2 += X[j][i][k];
					}
				}
				mModel.add(Cons6 == Cons6_2);
				mModel.add(Cons6_2 == Y[i][k]);
				Cons6.end();
				Cons6_2.end();
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
						IloExpr Cons6_5(mEnv);
						Cons6_5 = U[i - 1][k] - U[j - 1][k] + N * X[i][j][k];
						mModel.add(Cons6_5 <= N - 1);
						Cons6_5.end();
					}
				}
			}
		}

		//约束7 计划容量限制
		for (int k = 0; k < K; k++)
		{
			IloExpr Cons7(mEnv);
			for (int i = 1; i < N + 1; i++)
			{
				Cons7 += Y[i][k] * massA[i];
			}
			mModel.add(Cons7 <= parameters.planUp * Y[0][k]);
			mModel.add(Cons7 >= parameters.planLower * Y[0][k]);
			Cons7.end();
		}

		//约束8 目标流向
		for (int l = 0; l < L; l++)
		{
			IloExpr Cons8(mEnv);
			for (int i = 1; i < N + 1; i++)
			{
				for (int k = 0; k < K; k++)
				{
					Cons8 += Y[i][k] * massA[i] * flowArr[i][l];
				}
			}
			mModel.add(Cons8 >= MassFlow[l]);
			Cons8.end();
		}

		//目标函数
		IloExpr obj(mEnv);

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

		// second part
//		for (int k = 0; k < K; k++)
//		{
//			obj += 0.7 * parameters.planUp;
//			for (int i = 1; i < N + 1; i++)
//			{
//				obj -= 0.7 * Y[i][k] * massA[i];
//			}
//		}

		// third part
		for (int k = 0; k < K; k++)
		{
			for (int j = 1; j < N + 1; j++)
			{
				obj += parameters.alpha2 * X[0][j][k] * set_up;
			}
		}
		mModel.add(IloMinimize(mEnv, obj));
		obj.end();
	}
	catch (IloException& ex)
	{
		cerr << "InitModel Error -> " << ex << endl;
	}
	catch (...)
	{
		cerr << "Error when extractColumn" << endl;
	}
}

MathModel::MathModel(const Parameters& parameters, vector<Material*>& materials) : BaseModel()
{
	// 创建默认环境
	createModel();

	//初始化计划数 热卷数等常数
	C = parameters.C;
	L = parameters.L;
	K = parameters.K;
	N = parameters.N;
	S = parameters.S;

	// 初始化模型 添加变量和约束
	initModel(parameters, materials);
	//initModel_S(parameters, materials);
}

void MathModel::printPena(double** _pena) const
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

void MathModel::printResult(Problem* _problem)
{
	std::vector<Column> columns;

	Parameters parameters;

	std::cout << std::endl;
	std::cout << "_____________MathModel result_____________" << std::endl;
	for (int k = 0; k < K; k++)
	{
		// get chs names
		std::string chsName;
		for (int c = 0; c < C; c++)
		{
			if (cplex.getValue(Z[c][k]) > 0.99)
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
			if (cplex.getValue(X[0][j][k]) > 0.99)
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
				if (cplex.getValue(X[begin][j][k]) > 0.99)
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
				if (cplex.getValue(X[i][j][k]) > 0.99)
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
	std::cout << "cplex column cost -> " << cost << std::endl;
	// printResultR(cplex, X);

	std::cout << "-----------------------------------------" << std::endl;
	std::cout << "cplex_status : " << cplex.getStatus() << std::endl;
	std::cout << "cplex_gap : " << cplex.getMIPRelativeGap() << std::endl;
	std::cout << "cplex_upper : " << cplex.getObjValue() << std::endl;
	std::cout << "cplex_lower : " << cplex.getBestObjValue() << std::endl;
	std::cout << "cplex_cal_gap : " << (cplex.getObjValue() - cplex.getBestObjValue()) / cplex.getBestObjValue() << std::endl;
}

const IloArray<IloArray<IloNumVarArray>>& MathModel::getValues()
{
	return X;
}

double MathModel::getObjective(const Parameters& parameters) const
{
	double obj = 0.0;
	for (int k = 0; k < K; k++)
	{
		for (int i = 0; i < N + 1; i++)
		{
			for (int j = 0; j < N + 1; j++)
			{
				if (j != i)
				{
					obj += pena[i][j] * cplex.getValue(X[i][j][k]);
				}
			}
		}
	}
	obj = obj * 0.3;

	for (int k = 0; k < K; k++)
	{
		obj += 0.7 * parameters.planUp;
		for (int i = 1; i < N + 1; i++)
		{
			obj -= 0.7 * cplex.getValue(Y[i][k]) * massA[i];
		}
	}

	for (int k = 0; k < K; k++)
	{
		for (int j = 1; j < N + 1; j++)
		{
			obj += 0.7 * cplex.getValue(X[0][j][k]) * set_up;
		}
	}
	return obj;
}

void MathModel::printResultR(IloCplex& cplex, IloArray<IloArray<IloNumVarArray>>& x)
{
	double obj = 0;
	for (int k = 0; k < K; k++)
	{
		double a = 0;
		for (int i = 0; i < N + 1; i++)
		{
			for (int j = 0; j < N + 1; j++)
			{
				if (j != i && cplex.getValue(x[i][j][k]) > 0)
				{
					a += pena[i][j];
				}
			}
		}
		std::cout << a << std::endl;
	}

	int a = 7;
	for (int i = 0; i < N + 1; i++)
	{
		for (int j = 0; j < N + 1; j++)
		{
			if (cplex.getValue(x[i][j][6]) > 0)
			{
				std::cout << i << " " << j << std::endl;
			}
		}
	}

	for (int k = 0; k < K; k++)
	{
		for (int j = 1; j < N + 1; j++)
		{
			if (cplex.getValue(x[0][j][k]) > 0)
			{
				obj += 0.7 * 1000 * cplex.getValue(x[0][j][k]);
			}
		}
	}
	std::cout << obj << std::endl;
}

