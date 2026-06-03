#pragma once

#include <iostream>
#include "Parameters.h"

#ifndef IL_STD
#define IL_STD
#endif

#include <ilcplex/ilocplex.h>

ILOSTLBEGIN


class BaseModel
{
public:
	BaseModel() = default;
//	BaseModel(IloEnv& env, IloModel& model);

	virtual ~BaseModel() = default;

	int getStatusI() const;

	IloModel& getModel();

	IloEnv& getEnv();

	IloCplex& getSolver();

	void createModel();

	void update();

	bool optimize();

	bool optimize(const Parameters& _solver);

	double getObj() const;

protected:
	int mStatus{ -1 };
	IloEnv mEnv;
	IloModel mModel;
	IloCplex cplex;
};

