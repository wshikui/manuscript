//
// Created by wshikui on 2023/9/16.
//
#include "BaseModel.h"

//BaseModel::BaseModel(IloEnv& env, IloModel& model)
//{
//	mEnv = env;
//	mModel = model;
//	mStatus = -1;
//	cplex = IloCplex(mModel);
//}

int BaseModel::getStatusI() const
{
	return mStatus;
}

void BaseModel::update()
{
	//在创建的基础模型上添加新的变量 更新模型
}

void BaseModel::createModel()
{
	// 默认创建线性规划模型时调用 模型命名为BaseModel
	try
	{
		mEnv = IloEnv();
		mModel = IloModel(mEnv);
		mModel.setName("BaseModel");
	}
	catch (IloException& exception)
	{
		std::cerr << "Error: " << exception << std::endl;
	}
	catch (...)
	{
		std::cerr << "Unknown Error when create model" << std::endl;
	}
}

IloModel& BaseModel::getModel()
{
	return mModel;
}

IloEnv& BaseModel::getEnv()
{
	return mEnv;
}

IloCplex& BaseModel::getSolver()
{
	return cplex;
}

bool BaseModel::optimize()
{
	try
	{
		cplex = IloCplex(mModel);
		std::string pathLP = "Log/";
		pathLP += mModel.getName();
		pathLP += ".lp";

		//cplex.exportModel(pathLP.c_str());

		//not print any log from cplex
		cplex.setParam(IloCplex::Param::MIP::Display, 0);
		std::ofstream nullStream;
		nullStream.open("NUL");
		cplex.setOut(nullStream);

		Parameters _solver;
		_solver._logFile = "Log/";
		_solver._timeLimit = 3600;
		_solver.isShowLog = false;
		cplex.setParam(IloCplex::Param::TimeLimit, _solver._timeLimit);
		// cplex.setParam(IloCplex::EpGap, 0.01);
		if (!cplex.solve())
		{
			//mEnv.error() << "Failed!" << std::endl;
			//mEnv.out() << "Solution status: " << cplex.getStatus() << std::endl;

			mEnv.end();
			return false;
		}
		//mEnv.out() << "Solution status: " << cplex.getStatus() << std::endl;
		//mEnv.out() << "Solution value:  " << cplex.getObjValue() << std::endl;

		mStatus = cplex.getStatus();
		return true;
	}
	catch (IloException& exception)
	{
		cerr << "Error: " << exception << std::endl;
	}
	catch (...)
	{
		std::cerr << "Unknown Error when optimize Base model" << std::endl;
	}
	return false;
}

bool BaseModel::optimize(const Parameters& _solver)
{

	try
	{
		cplex = IloCplex(mModel);
		// IloCplex cplex(mModel);
		cplex.setParam(IloCplex::Param::TimeLimit, _solver._timeLimit);

		if (_solver.isShowLog)
		{
			std::string logFilePath = _solver._logFile + mModel.getName() + ".txt";
			std::ofstream logfile(logFilePath.c_str());
			cplex.setOut(logfile);
			logfile.close();
		}
		if (!cplex.solve())
		{
			mEnv.error() << "Failed!" << std::endl;
			mEnv.out() << "Solution status: " << cplex.getStatus() << std::endl;
			mEnv.end();
			return false;
		}
		mEnv.out() << "Solution status: " << cplex.getStatus() << std::endl;
		mEnv.out() << "Solution value:  " << cplex.getObjValue() << std::endl;
		mStatus = cplex.getStatus();
		mEnv.end();
		return true;
	}
	catch (IloException& exception)
	{
		cerr << "Error: " << exception << std::endl;
	}
	catch (...)
	{
		std::cerr << "Unknown Error when optimize Base model" << std::endl;
	}

	return false;
}

double BaseModel::getObj() const
{
	double obj = cplex.getObjValue();
	return obj;
}



