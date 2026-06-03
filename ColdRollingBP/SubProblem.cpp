//
// Created by wshikui on 2023/10/18.
//
#include "SubProblem.h"

SubProblem::SubProblem(MasterProblem* _masterProblem, std::vector<int>& vertexWithZero,
	std::map<std::string, std::vector<int>>& _poolIndex, 
	std::map<std::string, std::vector<std::bitset<Config::MAX_COILS>>>& _pool_coil_ngsets, Problem* _problem, Parameters* _parameters)
{
	rmp = _masterProblem;
	problem = _problem;
	parameters = _parameters;
	vertexZero = vertexWithZero;

	poolIndex = _poolIndex;
	pool_coil_ngsets = _pool_coil_ngsets;
}

SubProblem::SubProblem(MasterProblem* _masterProblem, std::vector<int>& vertexWithZero, 
	std::map<int, std::vector<int>> _forbidCoils,
	const std::vector<std::vector<int>>& _subset, 
	std::map<std::string, std::vector<int>>& _poolIndex, 
	std::map<std::string, std::vector<std::bitset<Config::MAX_COILS>>>& _pool_coil_ngsets, 
	Problem* _problem, Parameters* _parameters)
{
	rmp = _masterProblem;
	problem = _problem;
	parameters = _parameters;
	vertexZero = vertexWithZero;
	subset = _subset;
	poolIndex = _poolIndex;
	pool_coil_ngsets = _pool_coil_ngsets;
	forbidCoils = _forbidCoils;
}

double SubProblem::getPlanWeight(vector<int>& index, int l)
{
	double weight = 0.0;
	for (const auto& item : index)
	{
		if (problem->getFlow(item) == parameters->target_flow[l])
		{
			weight += problem->getWeight(item);
		}
	}
	return weight;
}

void SubProblem::splitPoolIndex(std::map<std::string, std::vector<int>>& _poolIndex)
{
	// index of single pool, include the virtual vertex for get coil's index
	for (int c = 0; c < parameters->C; c++)
	{
		std::string name = parameters->chs_type[c];
		std::vector<int> index = { 0 };
		for (int i = 0; i < parameters->S; i++)
		{
			index.emplace_back(c * parameters->S + 1 + i);
		}
		_poolIndex[name] = index;
	}

	// index of connect pool
	std::string conName1 = "C502016_17";
	std::vector<int> index1 = { 0 };
	index1.insert(index1.end(), _poolIndex["C502016"].begin() + 1, _poolIndex["C502016"].end());
	index1.insert(index1.end(), _poolIndex["C502017"].begin() + 1, _poolIndex["C502017"].end());
	_poolIndex[conName1] = index1;

	std::string conName2 = "C502017_18";
	std::vector<int> index2 = { 0 };
	index2.insert(index2.end(), _poolIndex["C502017"].begin() + 1, _poolIndex["C502017"].end());
	index2.insert(index2.end(), _poolIndex["C502018"].begin() + 1, _poolIndex["C502018"].end());
	_poolIndex[conName2] = index2;

	std::string conName3 = "C502016_18";
	std::vector<int> index3 = { 0 };
	index3.insert(index3.end(), _poolIndex["C502016"].begin() + 1, _poolIndex["C502016"].end());
	index3.insert(index3.end(), _poolIndex["C502018"].begin() + 1, _poolIndex["C502018"].end());
	_poolIndex[conName3] = index3;

	std::string conName4 = "C502018_19";
	std::vector<int> index4 = { 0 };
	index4.insert(index4.end(), _poolIndex["C502018"].begin() + 1, _poolIndex["C502018"].end());
	index4.insert(index4.end(), _poolIndex["C502019"].begin() + 1, _poolIndex["C502019"].end());
	_poolIndex[conName4] = index4;
}

void
SubProblem::getRoute(IloCplex& _cplex, IloArray<IloNumVarArray>& x, Column& _column, const std::vector<int>& _index,
	IloNumVarArray& y)
{
	// get the vertex connect with virtual vertex
	int begin = -1;
	int n = static_cast<int>(x.getSize());

	for (int i = 1; i < n; i++)
	{
		if (_cplex.getValue(x[0][i]) > 0.99)
		{
			begin = static_cast<int>(i);
			break;
		}
	}
	if (begin == -1)
	{
		std::cout << "route begin == -1" << std::endl;
	}
	// add vertex into column
	std::vector<int> route = { _index[begin] };
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			if (_cplex.getValue(x[begin][j]) > 0.99)
			{
				begin = j;
				if (begin == 0)
				{
					break;
				}
				route.emplace_back(_index[j]);
				break;
			}
		}
		if (begin == 0)
		{
			break;
		}
	}

	for (auto& item : route)
	{
		_column.addVertex(item);
	}
	_column.calculateCost();
	_column.printRoute();

	// calculate cost of column according to \sum_{x_{ij}c_{ij}}
	double reduceCost;
	double cost = 0.0;
	for (IloInt i = 1; i < x.getSize(); i++)
	{
		for (IloInt j = 1; j < x.getSize(); j++)
		{
			if (i != j && _cplex.getValue(x[i][j]) > 0.99)
			{
				cost += problem->getCost(_index[i] - 1, _index[j] - 1);
			}
		}
	}
	reduceCost = cost * parameters->alpha1 + parameters->alpha2 * parameters->setUpCost;
	// calculate the reduce cost of this column
	for (IloInt i = 1; i < y.getSize(); i++)
	{
		if (_cplex.getValue(y[i]) > 0.99)
		{
			reduceCost -= rmp->getDualVariable(_index[i] - 1);
		}
	}
	for (int l = 0; l < parameters->L; l++)
	{
		double mass = 0.0;
		for (IloInt i = 1; i < y.getSize(); i++)
		{
			if (_cplex.getValue(y[i]) > 0.99)
			{
				int h_il = problem->getFlow(_index[i]) == parameters->target_flow[l] ? 1 : 0;
				mass += h_il * problem->getWeight(_index[i]);
			}
		}
		reduceCost -= mass * rmp->getDualVariable(problem->getVertices() + l);
	}
	reduceCost -= rmp->getDualVariable(problem->getVertices() + parameters->L);
	reduceCost -= rmp->getDualVariable(problem->getVertices() + parameters->L + 1);
	std::cout << "cost -> " << cost << std::endl;
	std::cout << "reduced cost -> " << reduceCost << std::endl;
	std::cout << "objective -> " << _cplex.getObjValue() << std::endl;
}

void SubProblem::getRouteSR(IloCplex& _cplex, IloArray<IloNumVarArray>& x, Column& _column, const std::vector<int>& _index, IloNumVarArray& y,
	std::vector<std::vector<int>>& _subset)
{
	// get the vertex connect with virtual vertex
	int begin = -1;
	int n = static_cast<int>(x.getSize());

	for (int i = 1; i < n; i++)
	{
		if (_cplex.getValue(x[0][i]) > 0.99)
		{
			begin = static_cast<int>(i);
			break;
		}
	}
	if (begin == -1)
	{
		std::cout << "route begin == -1" << std::endl;
	}
	// add vertex into column
	std::vector<int> route = { _index[begin] };
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			if (_cplex.getValue(x[begin][j]) > 0.99)
			{
				begin = j;
				if (begin == 0)
				{
					break;
				}
				route.emplace_back(_index[j]);
				break;
			}
		}
		if (begin == 0)
		{
			break;
		}
	}

	for (auto& item : route)
	{
		_column.addVertex(item);
	}
	_column.calculateCost();
	//_column.printRoute();

	// calculate cost of column according to \sum_{x_{ij}c_{ij}}
	double reduceCost;
	double cost = 0.0;
	for (IloInt i = 1; i < x.getSize(); i++)
	{
		for (IloInt j = 1; j < x.getSize(); j++)
		{
			if (i != j && _cplex.getValue(x[i][j]) > 0.99)
			{
				cost += problem->getCost(_index[i] - 1, _index[j] - 1);
			}
		}
	}
	reduceCost = cost * parameters->alpha1 + parameters->alpha2 * parameters->setUpCost;
	// calculate the reduce cost of this column
	for (IloInt i = 1; i < y.getSize(); i++)
	{
		if (_cplex.getValue(y[i]) > 0.99)
		{
			reduceCost -= rmp->getDualVariable(_index[i] - 1);
		}
	}
	for (int l = 0; l < parameters->L; l++)
	{
		double mass = 0.0;
		for (IloInt i = 1; i < y.getSize(); i++)
		{
			if (_cplex.getValue(y[i]) > 0.99)
			{
				int h_il = problem->getFlow(_index[i]) == parameters->target_flow[l] ? 1 : 0;
				mass += h_il * problem->getWeight(_index[i]);
			}
		}
		reduceCost -= mass * rmp->getDualVariable(problem->getVertices() + l);
	}
	reduceCost -= rmp->getDualVariable(problem->getVertices() + parameters->L);
	reduceCost -= rmp->getDualVariable(problem->getVertices() + parameters->L + 1);

	//add SR dual values
	for (int s = 0; s < subset.size(); s++) {
		std::vector<int> S = subset[s];
		if (countContainedElements(_column.getRoute(), S) >= 2) {
			reduceCost -= rmp->getDualVariable(problem->getVertices() + parameters->L + 2 + s);
		}
	}

	//std::cout << "cost -> " << cost << std::endl;
	//std::cout << "reduced cost -> " << reduceCost << std::endl;
	//std::cout << "objective -> " << _cplex.getObjValue() << std::endl;
}

bool SubProblem::solverPool(const std::string& _name, std::vector<int>& index)
{
	createModel();

	std::string name = "SubModel_";
	name += _name;
	mModel.setName(name.c_str());

	// notice : the parameter "index" include virtual coil
	int numCoils = static_cast<int>(index.size()) - 1;

	try
	{
		IloArray<IloNumVarArray> X(mEnv, numCoils + 1);
		IloNumVarArray Y(mEnv, numCoils + 1);
		IloNumVarArray U(mEnv, numCoils);

		for (int i = 0; i < numCoils + 1; i++)
		{
			X[i] = IloNumVarArray(mEnv, numCoils + 1);
			for (int j = 0; j < numCoils + 1; j++)
			{
				X[i][j] = IloNumVar(mEnv, 0, 1, ILOINT);
			}
		}
		for (int i = 0; i < numCoils + 1; i++)
		{
			Y[i] = IloNumVar(mEnv, 0, 1, ILOINT);
		}
		for (int i = 0; i < numCoils; i++)
		{
			U[i] = IloNumVar(mEnv, 1, numCoils, ILOINT);
		}

		IloExpr cons0(mEnv);
		for (int i = 0; i < numCoils + 1; i++)
		{
			cons0 += X[i][i];
		}
		mModel.add(cons0 == 0);
		cons0.end();

		for (int i = 0; i < numCoils + 1; i++)
		{
			IloExpr cons1_1(mEnv);
			IloExpr cons1_2(mEnv);
			for (int j = 0; j < numCoils + 1; j++)
			{
				if (j != i)
				{
					cons1_1 += X[j][i];
				}
			}
			for (int j = 0; j < numCoils + 1; j++)
			{
				if (j != i)
				{
					cons1_2 += X[i][j];
				}
			}
			mModel.add(cons1_1 == cons1_2);
			mModel.add(cons1_2 == Y[i]);
			cons1_1.end();
			cons1_2.end();
		}
		// flow constraint
		for (int i = 1; i < numCoils + 1; i++)
		{
			for (int j = 1; j < numCoils + 1; j++)
			{
				if (j != i)
				{
					IloExpr cons3(mEnv);
					cons3 = U[i - 1] - U[j - 1] + numCoils * X[i][j];
					mModel.add(cons3 <= numCoils - 1);
					cons3.end();
				}
			}
		}
		// route weight constraint: exclude virtual coil
		IloExpr cons4(mEnv);
		for (int i = 1; i < numCoils + 1; i++)
		{
			cons4 += Y[i] * problem->getWeight(index[i]);
		}
		mModel.add(cons4 <= parameters->planUp);
		mModel.add(cons4 >= parameters->planLower);
		cons4.end();

		// no selecting some coils based branch rule
		if (!vertexZero.empty()) {
			for (const auto& coilIndex : vertexZero) {
				if (std::find(index.begin(), index.end(), coilIndex) != index.end()) {
					auto it = std::find(index.begin(), index.end(), coilIndex);
					// get the no select coil index in this "pool index"
					int noSelectIndex = static_cast<int>(std::distance(index.begin(), it));
					IloExpr consNoSelect(mEnv);
					consNoSelect = Y[noSelectIndex];
					mModel.add(consNoSelect == 0);
					consNoSelect.end();
				}
			}
		}

		// objective, need the dual variables from RMP (address)
		IloExpr obj(mEnv);
		for (int i = 1; i < numCoils + 1; i++)
		{
			for (int j = 1; j < numCoils + 1; j++)
			{
				if (j != i)
				{
					obj += X[i][j] * problem->getCost(index[i] - 1, index[j] - 1) * parameters->alpha1;
				}
			}
		}
		obj += parameters->alpha2 * parameters->setUpCost;
		for (int i = 1; i < numCoils + 1; i++)
		{
			obj -= Y[i] * rmp->getDualVariable(index[i] - 1);
		}
		for (int l = 0; l < parameters->L; l++)
		{
			for (int i = 1; i < numCoils + 1; i++)
			{
				int h_il = problem->getFlow(index[i]) == parameters->target_flow[l] ? 1 : 0;
				obj -= Y[i] * h_il * problem->getWeight(index[i]) * rmp->getDualVariable(parameters->N + l);
			}
		}
		mModel.add(IloMinimize(mEnv, obj));
		obj.end();

		cplex = IloCplex(mModel);
		cplex.setParam(IloCplex::Param::MIP::Display, 0);
		std::ofstream nullStream;
		nullStream.open("NUL");
		cplex.setOut(nullStream);
		if (!cplex.solve())
		{
			mEnv.error() << "Failed!" << std::endl;
			mEnv.out() << "Solution status: " << cplex.getStatus() << std::endl;
			mEnv.end();
			return false;
		}

		// create candidate column
		Column column(problem);
		getRoute(cplex, X, column, index, Y);
		double rc = cplex.getObjValue() - rmp->getDualVariable(problem->getVertices() + parameters->L)
			- rmp->getDualVariable(problem->getVertices() + parameters->L + 1);
		column.setRC(rc);
		columns.emplace_back(column);
		mEnv.end();
		return true;
	}
	catch (IloException& ex)
	{
		cerr << "single pool Problem solve_lp Error -> " << ex << endl;
	}
	catch (...)
	{
		cerr << "Error when single pool extractColumn" << endl;
	}
	return false;
}

bool SubProblem::solverPoolLabel(const std::string& _name, std::vector<int>& index)
{

	//initialize a begin label, notice the index set has the begin index
	std::vector<int> allIndex;
	allIndex.insert(allIndex.begin(), index.begin() + 1, index.end());
	allIndex.emplace_back(parameters->N + 1);

	//remove coils that are not selected due to branch rule
	for (const auto& item : vertexZero) {
		// item in all coil index of pool
		if (std::find(allIndex.begin(), allIndex.end(), item) != allIndex.end()) {
			allIndex.erase(std::remove(allIndex.begin(), allIndex.end(), item), allIndex.end());
		}
	}

	ElementaryLabel start_label(allIndex);

	//extend function
	LabelExtender extender(vertexZero, problem, rmp, parameters);

	LabellingAlgorithm<ElementaryLabel, LabelExtender> alg;
	std::vector<Column> newColumns = alg.solveR(0, allIndex.back(), start_label, extender, problem);


	//add new column to pool : the column with the highest mass or reduced cost
	if (!newColumns.empty()) {
		double mass = 0;
		double cost = 0;
		int index = 0;
		int index_cost = 0;
		for (int i = 0; i < newColumns.size(); i++) {
			if (newColumns[i].getRC() < cost) {
				cost = newColumns[i].getRC();
				index_cost = i;
			}
		}
		if (newColumns[index_cost].getRC() >= 0) {
			return false;
		}
		columns.emplace_back(newColumns[index_cost]);
		newColumns[index_cost].printRoute();

		return true;
	}
	return false;
}

bool SubProblem::solverPoolLabelSR(const std::string& _name, std::vector<int>& index)
{
	std::vector<int> allIndex;
	allIndex.insert(allIndex.begin(), index.begin() + 1, index.end());
	allIndex.emplace_back(parameters->N + 1);

	//remove coils that are not selected due to branch rule
	for (const auto& item : vertexZero) {
		// item in all coil index of pool
		if (std::find(allIndex.begin(), allIndex.end(), item) != allIndex.end()) {
			allIndex.erase(std::remove(allIndex.begin(), allIndex.end(), item), allIndex.end());
		}
	}

	//the start label with SR inequality
	int constraintsNum = parameters->N + parameters->L + 2;
	ElementaryLabelSR start_label(allIndex, subset, rmp, constraintsNum);

	//extend function
	LabelExtender extender(vertexZero, problem, rmp, parameters);

	//test dominated

	//auto next_1_label = extender(start_label, 1, parameters->N + 1);
	//auto next_2_label = extender(start_label, 2, parameters->N + 1);
	//auto next_2_1_label = extender(next_2_label.value(), 1, parameters->N + 1);
	//if (next_3_label < next_2_3_label) {
	//	std::cout << "debug" << std::endl;
	//}

	//label setting algorithm
	LabellingAlgorithm<ElementaryLabelSR, LabelExtender> alg;

	std::vector<Column> newColumns = alg.solveR(0, allIndex.back(), start_label, extender, problem);

	//add new column to pool : the column with the highest mass or reduced cost
	if (!newColumns.empty()) {
		double mass = 0;
		double cost = 0;
		int index = 0;
		int index_cost = 0;
		for (int i = 0; i < newColumns.size(); i++) {
			if (newColumns[i].getRC() < cost) {
				cost = newColumns[i].getRC();
				index_cost = i;
			}
		}
		if (newColumns[index_cost].getRC() >= 0) {
			return false;
		}
		newColumns[index_cost].setName(_name);
		columns.emplace_back(newColumns[index_cost]);
		return true;
	}
	return false;
}

bool SubProblem::solverPoolSR(const std::string& _name, std::vector<std::vector<int>>& _subset, std::vector<int>& index)
{
	createModel();

	std::string name = "SubModel_sr_";
	name += _name;
	mModel.setName(name.c_str());

	int numCoils = static_cast<int>(index.size()) - 1;
	int numS = static_cast<int>(_subset.size());

	try {
		IloArray<IloNumVarArray> X(mEnv, numCoils + 1);
		IloNumVarArray Y(mEnv, numCoils + 1);
		IloNumVarArray U(mEnv, numCoils);
		IloNumVarArray Z(mEnv, numS);

		for (int i = 0; i < numCoils + 1; i++)
		{
			X[i] = IloNumVarArray(mEnv, numCoils + 1);
			for (int j = 0; j < numCoils + 1; j++)
			{
				X[i][j] = IloNumVar(mEnv, 0, 1, ILOINT);
			}
		}
		for (int i = 0; i < numCoils + 1; i++)
		{
			Y[i] = IloNumVar(mEnv, 0, 1, ILOINT);
		}
		for (int i = 0; i < numCoils; i++)
		{
			U[i] = IloNumVar(mEnv, 1, numCoils, ILOINT);
		}
		if (!_subset.empty())
		{
			for (int i = 0; i < numS; i++)
			{
				Z[i] = IloNumVar(mEnv, 0, 1, ILOINT);
			}
		}


		IloExpr cons0(mEnv);
		for (int i = 0; i < numCoils + 1; i++)
		{
			cons0 += X[i][i];
		}
		mModel.add(cons0 == 0);
		cons0.end();

		for (int i = 0; i < numCoils + 1; i++)
		{
			IloExpr cons1_1(mEnv);
			IloExpr cons1_2(mEnv);
			for (int j = 0; j < numCoils + 1; j++)
			{
				if (j != i)
				{
					cons1_1 += X[j][i];
				}
			}
			for (int j = 0; j < numCoils + 1; j++)
			{
				if (j != i)
				{
					cons1_2 += X[i][j];
				}
			}
			mModel.add(cons1_1 == cons1_2);
			mModel.add(cons1_2 == Y[i]);
			cons1_1.end();
			cons1_2.end();
		}
		// flow constraint
		for (int i = 1; i < numCoils + 1; i++)
		{
			for (int j = 1; j < numCoils + 1; j++)
			{
				if (j != i)
				{
					IloExpr cons3(mEnv);
					cons3 = U[i - 1] - U[j - 1] + numCoils * X[i][j];
					mModel.add(cons3 <= numCoils - 1);
					cons3.end();
				}
			}
		}
		// route weight constraint: exclude virtual coil
		IloExpr cons4(mEnv);
		for (int i = 1; i < numCoils + 1; i++)
		{
			cons4 += Y[i] * problem->getWeight(index[i]);
		}
		mModel.add(cons4 <= parameters->planUp);
		mModel.add(cons4 >= parameters->planLower);
		cons4.end();

		// no selecting some coils based branch rule
		if (!vertexZero.empty()) {
			for (const auto& coilIndex : vertexZero) {
				if (std::find(index.begin(), index.end(), coilIndex) != index.end()) {
					auto it = std::find(index.begin(), index.end(), coilIndex);
					// get the no select coil index in this "pool index"
					int noSelectIndex = static_cast<int>(std::distance(index.begin(), it));
					IloExpr consNoSelect(mEnv);
					consNoSelect = Y[noSelectIndex];
					mModel.add(consNoSelect == 0);
					consNoSelect.end();
				}
			}
		}
		//SR inequality, the sequence of index in S should be ascending
		for (int s = 0; s < numS; s++) {
			std::vector<int> S = _subset[s];
			if (countContainedElements(index, S) >= 2) {
				for (int i = 0; i < 3; i++) {
					for (int j = i + 1; j < 3; j++) {
						// i and j exist in this index set
						if (std::count(index.begin(), index.end(), S[i]) > 0 &&
							std::count(index.begin(), index.end(), S[j]) > 0) {
							IloExpr conSR(mEnv);
							int coil_i = S[i];
							int coil_j = S[j];
							auto iti = std::find(index.begin(), index.end(), coil_i);
							auto itj = std::find(index.begin(), index.end(), coil_j);
							int index_i = static_cast<int>(std::distance(index.begin(), iti));
							int index_j = static_cast<int>(std::distance(index.begin(), itj));
							conSR = Y[index_i] + Y[index_j] - 1;
							mModel.add(conSR <= Z[s]);
							conSR.end();
						}
					}
				}
			}
		}

		// objective, need the dual variables from RMP (address)
		// the first part : x_{ij} don't include virtual coil
		IloExpr obj(mEnv);
		for (int i = 1; i < numCoils + 1; i++)
		{
			for (int j = 1; j < numCoils + 1; j++)
			{
				if (j != i)
				{
					obj += X[i][j] * problem->getCost(index[i] - 1, index[j] - 1) * parameters->alpha1;
				}
			}
		}
		obj += parameters->alpha2 * parameters->setUpCost;
		// notice: get the dual variable according to the coil index, no the index "i"
		for (int i = 1; i < numCoils + 1; i++)
		{
			obj -= Y[i] * rmp->getDualVariable(index[i] - 1);
		}
		for (int l = 0; l < parameters->L; l++)
		{
			for (int i = 1; i < numCoils + 1; i++)
			{
				int h_il = problem->getFlow(index[i]) == parameters->target_flow[l] ? 1 : 0;
				obj -= Y[i] * h_il * problem->getWeight(index[i]) * rmp->getDualVariable(parameters->N + l);
			}
		}
		//may exist error
		for (int s = 0; s < numS; s++) {
			obj -= Z[s] * rmp->getDualVariable(problem->getVertices() + parameters->L + 2 + s);
		}

		mModel.add(IloMinimize(mEnv, obj));
		obj.end();

		cplex = IloCplex(mModel);
		cplex.setParam(IloCplex::Param::MIP::Display, 0);
		std::ofstream nullStream;
		nullStream.open("NUL");
		cplex.setOut(nullStream);

		if (!cplex.solve())
		{
			mEnv.error() << "Failed!" << std::endl;
			mEnv.out() << "Solution status: " << cplex.getStatus() << std::endl;
			mEnv.end();
			return false;
		}

		// create candidate column
		Column column(problem);

		//add vertex into the column and calculate the RC
		getRouteSR(cplex, X, column, index, Y, subset);

		double rc = cplex.getObjValue() - rmp->getDualVariable(problem->getVertices() + parameters->L)
			- rmp->getDualVariable(problem->getVertices() + parameters->L + 1);
		column.setRC(rc);
		column.setName(_name);
		columns.emplace_back(column);
		mEnv.end();
		return true;

	}
	catch (IloException& ex)
	{
		cerr << "single pool Problem solve_lp Error -> " << ex << endl;
	}
	catch (...)
	{
		cerr << "Error when single pool extractColumn" << endl;
	}

	return false;
}

bool SubProblem::solverPoolLabelSR_ng_routes(const std::string& _name, std::vector<int>& index)
{
	std::vector<int> allIndex;
	allIndex.insert(allIndex.begin(), index.begin() + 1, index.end());
	allIndex.emplace_back(parameters->N + 1);

	for (const auto& item : vertexZero) {
		// item in all coil index of pool
		if (std::find(allIndex.begin(), allIndex.end(), item) != allIndex.end()) {
			allIndex.erase(std::remove(allIndex.begin(), allIndex.end(), item), allIndex.end());
		}
	}

	//the start label with SR inequality
	int constraintsNum = parameters->N + parameters->L + 2;
	LabelSR start_label(parameters, allIndex, subset, rmp, constraintsNum);
	start_label.compatibleVertex = allIndex;
	start_label.ng_memory.reset();
	start_label.ng_memory.set(0);

	//extend function
	LabelExtender extender(allIndex, vertexZero, forbidCoils , poolIndex[_name], pool_coil_ngsets[_name], problem, rmp, parameters);
	LabellingAlgorithm<LabelSR, LabelExtender> alg;
	std::vector<Column> newColumns = alg.solveR(0, allIndex.back(), start_label, extender, problem);

	//add new column to pool
	if (!newColumns.empty()) {
		double mass = 0;
		double cost = 0;
		int index = 0;
		int index_cost = 0;
		for (int i = 0; i < newColumns.size(); i++) {
			if (newColumns[i].getRC() < cost) {
				cost = newColumns[i].getRC();
				index_cost = i;
			}
		}
		if (newColumns[index_cost].getRC() >= 0) {
			return false;
		}
		newColumns[index_cost].setName(_name);
		columns.emplace_back(newColumns[index_cost]);
		return true;
	}
	return false;
}

void SubProblem::reSetSubNum()
{
	solveNum = 0;
}

int SubProblem::countContainedElements(const std::vector<int>& mainVec, const std::vector<int>& searchVec)
{
	std::unordered_set<int> searchSet(searchVec.begin(), searchVec.end());

	return std::count_if(mainVec.begin(), mainVec.end(), [&searchSet](int elem) {
		return searchSet.count(elem);
		});
}


std::vector<Column>& SubProblem::extractColumns()
{
	columns.clear();

	//the flag of solving merged pool
	std::map<std::string, bool> flags = {
		std::pair<std::string, bool>("C502016_17", false),
		std::pair<std::string, bool>("C502017_18", false),
		std::pair<std::string, bool>("C502016_18", false),
		std::pair<std::string, bool>("C502018_19", false)
	};

	//wether need to solve
	std::map<std::string, bool> flagsSingel = {
		std::pair<std::string, bool>("C502016", true),
		std::pair<std::string, bool>("C502017", true),
		std::pair<std::string, bool>("C502018", true),
		std::pair<std::string, bool>("C502019", true)
	};

	//for each merged pool firstly
	for (const auto& name : conPoolNames) {
		bool flag = solverPool(name, poolIndex[name]);
		flags[name] = flag;
	}

	//solve the single pool based on merged pool's flag
	if (!flags["C502016_17"]) {
		flagsSingel["C502016"] = false;
		flagsSingel["C502017"] = false;
	}
	if (!flags["C502017_18"]) {
		flagsSingel["C502017"] = false;
		flagsSingel["C502018"] = false;
	}
	if (!flags["C502016_18"]) {
		flagsSingel["C502018"] = false;
		flagsSingel["C502016"] = false;
	}
	if (!flags["C502018_19"]) {
		flagsSingel["C502019"] = false;
		flagsSingel["C502018"] = false;
	}

	for (const auto& item : columns) {
		std::string mergedPool = item.getName();
		std::cout << mergedPool << std::endl;
		item.printRoute();
		std::string pool1 = mergedPool.substr(0, 7);
		std::string pool2 = mergedPool.substr(0, 5) + mergedPool.substr(8, 10);

		std::vector<int> indexes = item.getRoute();
		std::vector<int> indexPool1 = poolIndex[pool1];
		std::vector<int> indexPool2 = poolIndex[pool2];
		bool isFromPool1 = true;
		bool isFromPool2 = true;

		for (const auto& index : indexes) {
			if (std::find(indexPool1.begin(), indexPool1.end(), index) == indexPool1.end()) {
				isFromPool1 = false;
			}
			if (std::find(indexPool2.begin(), indexPool2.end(), index) == indexPool2.end()) {
				isFromPool2 = false;
			}
		}
		if (isFromPool1) {
			flagsSingel[pool1] = false;
		}
		if (isFromPool2) {
			flagsSingel[pool2] = false;
		}
	}

	for (const auto& item : flagsSingel) {
		if (item.second) {
			bool flag = solverPool(item.first, poolIndex[item.first]);
		}
	}

	//remove the duplicate columns
	std::vector<Column> columnsNew;
	for (const auto& item : columns) {
		if (columnsNew.empty()) {
			columnsNew.emplace_back(item);
			continue;
		}
		auto it = std::find_if(columnsNew.begin(), columnsNew.end(),
			[&](const Column& c_in_new_pool) {
				return item.getRoute() == c_in_new_pool.getRoute();
			});
		if (it != columnsNew.end()) {
			continue;
		}
		columnsNew.emplace_back(item);
	}
	columns = columnsNew;

	return columns;
}

std::vector<Column>& SubProblem::extractColumnsLabel()
{
	columns.clear();

	//the flag of solving merged pool
	std::map<std::string, bool> flags = {
		std::pair<std::string, bool>("C502016_17", false),
		std::pair<std::string, bool>("C502017_18", false),
		std::pair<std::string, bool>("C502016_18", false),
		std::pair<std::string, bool>("C502018_19", false)
	};

	//wether need to solve
	std::map<std::string, bool> flagsSingel = {
		std::pair<std::string, bool>("C502016", true),
		std::pair<std::string, bool>("C502017", true),
		std::pair<std::string, bool>("C502018", true),
		std::pair<std::string, bool>("C502019", true)
	};

	//for each merged pool firstly
	for (const auto& name : conPoolNames) {
		auto start = std::chrono::high_resolution_clock::now();
		bool flag = solverPoolLabel(name, poolIndex[name]);
		auto end = std::chrono::high_resolution_clock::now();
		std::cout << "solver pool time = " << std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count() << std::endl;
		flags[name] = flag;
	}

	//solve the single pool based on merged pool's flag
	if (!flags["C502016_17"]) {
		flagsSingel["C502016"] = false;
		flagsSingel["C502017"] = false;
	}
	if (!flags["C502017_18"]) {
		flagsSingel["C502017"] = false;
		flagsSingel["C502018"] = false;
	}
	if (!flags["C502016_18"]) {
		flagsSingel["C502018"] = false;
		flagsSingel["C502016"] = false;
	}
	if (!flags["C502018_19"]) {
		flagsSingel["C502019"] = false;
		flagsSingel["C502018"] = false;
	}

	for (const auto& item : columns) {
		std::string mergedPool = item.getName();
		std::cout << mergedPool << std::endl;
		item.printRoute();
		std::string pool1 = mergedPool.substr(0, 7);
		std::string pool2 = mergedPool.substr(0, 5) + mergedPool.substr(8, 10);

		std::vector<int> indexes = item.getRoute();
		std::vector<int> indexPool1 = poolIndex[pool1];
		std::vector<int> indexPool2 = poolIndex[pool2];
		bool isFromPool1 = true;
		bool isFromPool2 = true;

		for (const auto& index : indexes) {
			if (std::find(indexPool1.begin(), indexPool1.end(), index) == indexPool1.end()) {
				isFromPool1 = false;
			}
			if (std::find(indexPool2.begin(), indexPool2.end(), index) == indexPool2.end()) {
				isFromPool2 = false;
			}
		}
		if (isFromPool1) {
			flagsSingel[pool1] = false;
		}
		if (isFromPool2) {
			flagsSingel[pool2] = false;
		}
	}

	for (const auto& item : flagsSingel) {
		if (item.second) {
			bool flag = solverPoolLabel(item.first, poolIndex[item.first]);
		}
	}

	std::vector<Column> columnsNew;
	for (const auto& item : columns) {
		if (columnsNew.empty()) {
			columnsNew.emplace_back(item);
			continue;
		}
		auto it = std::find_if(columnsNew.begin(), columnsNew.end(),
			[&](const Column& c_in_new_pool) {
				return item.getRoute() == c_in_new_pool.getRoute();
			});
		if (it != columnsNew.end()) {
			continue;
		}
		columnsNew.emplace_back(item);
	}
	columns = columnsNew;

	return columns;
}

std::vector<Column>& SubProblem::extractColumnsSR()
{
	columns.clear();

	//the flag of solving merged pool
	std::map<std::string, bool> flags = {
		std::pair<std::string, bool>("C502016_17", false),
		std::pair<std::string, bool>("C502017_18", false),
		std::pair<std::string, bool>("C502016_18", false),
		std::pair<std::string, bool>("C502018_19", false)
	};

	//wether need to solve
	std::map<std::string, bool> flagsSingel = {
		std::pair<std::string, bool>("C502016", true),
		std::pair<std::string, bool>("C502017", true),
		std::pair<std::string, bool>("C502018", true),
		std::pair<std::string, bool>("C502019", true)
	};

	//for each merged pool firstly
	for (const auto& name : conPoolNames) {
		auto start = std::chrono::high_resolution_clock::now();
		bool flag = solverPoolSR(name, subset, poolIndex[name]);
		auto end = std::chrono::high_resolution_clock::now();
		//std::cout << "solver pool time = " << std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count() << std::endl;
		flags[name] = flag;
		solveNum++;
	}

	////solve the single pool based on merged pool's flag
	if (!flags["C502016_17"]) {
		flagsSingel["C502016"] = false;
		flagsSingel["C502017"] = false;
	}
	if (!flags["C502017_18"]) {
		flagsSingel["C502017"] = false;
		flagsSingel["C502018"] = false;
	}
	if (!flags["C502016_18"]) {
		flagsSingel["C502018"] = false;
		flagsSingel["C502016"] = false;
	}
	if (!flags["C502018_19"]) {
		flagsSingel["C502019"] = false;
		flagsSingel["C502018"] = false;
	}

	//the columns form merged pool are same as signal pool
	for (const auto& item : columns) {
		std::string mergedPool = item.getName();
		//std::cout << mergedPool << std::endl;
		//item.printRoute();
		std::string pool1 = mergedPool.substr(0, 7);
		std::string pool2 = mergedPool.substr(0, 5) + mergedPool.substr(8, 10);

		std::vector<int> indexes = item.getRoute();
		std::vector<int> indexPool1 = poolIndex[pool1];
		std::vector<int> indexPool2 = poolIndex[pool2];
		bool isFromPool1 = true;
		bool isFromPool2 = true;

		for (const auto& index : indexes) {
			if (std::find(indexPool1.begin(), indexPool1.end(), index) == indexPool1.end()) {
				isFromPool1 = false;
			}
			if (std::find(indexPool2.begin(), indexPool2.end(), index) == indexPool2.end()) {
				isFromPool2 = false;
			}
		}
		if (isFromPool1) {
			flagsSingel[pool1] = false;
		}
		if (isFromPool2) {
			flagsSingel[pool2] = false;
		}
	}

	for (const auto& item : flagsSingel) {
		if (item.second) {
			bool flag = solverPoolSR(item.first, subset, poolIndex[item.first]);
			solveNum++;
		}
	}

	//remove the duplicate columns
	std::vector<Column> columnsNew;
	for (const auto& item : columns) {
		if (columnsNew.empty()) {
			columnsNew.emplace_back(item);
			continue;
		}
		auto it = std::find_if(columnsNew.begin(), columnsNew.end(),
			[&](const Column& c_in_new_pool) {
				return item.getRoute() == c_in_new_pool.getRoute();
			});
		if (it != columnsNew.end()) {
			continue;
		}
		columnsNew.emplace_back(item);
	}
	columns = columnsNew;

	return columns;
}

std::vector<Column>& SubProblem::extractColumnsLabelSR()
{
	columns.clear();

	//the flag of solving merged pool
	std::map<std::string, bool> flags = {
		std::pair<std::string, bool>("C502016_17", false),
		std::pair<std::string, bool>("C502017_18", false),
		std::pair<std::string, bool>("C502016_18", false),
		std::pair<std::string, bool>("C502018_19", false)
	};

	//wether need to solve
	std::map<std::string, bool> flagsSingel = {
		std::pair<std::string, bool>("C502016", true),
		std::pair<std::string, bool>("C502017", true),
		std::pair<std::string, bool>("C502018", true),
		std::pair<std::string, bool>("C502019", true)
	};

	//for each merged pool firstly
	for (const auto& name : conPoolNames) {
		auto start = std::chrono::high_resolution_clock::now();
		//bool flag = solverPoolLabelSR(name, poolIndex[name]);
		bool flag = solverPoolLabelSR_ng_routes(name, poolIndex[name]);
		auto end = std::chrono::high_resolution_clock::now();
		//std::cout << "solver pool time = " << std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count() << std::endl;
		flags[name] = flag;
		solveNum++;
	}

	//solve the single pool based on merged pool's flag
	if (!flags["C502016_17"]) {
		flagsSingel["C502016"] = false;
		flagsSingel["C502017"] = false;
	}
	if (!flags["C502017_18"]) {
		flagsSingel["C502017"] = false;
		flagsSingel["C502018"] = false;
	}
	if (!flags["C502016_18"]) {
		flagsSingel["C502018"] = false;
		flagsSingel["C502016"] = false;
	}
	if (!flags["C502018_19"]) {
		flagsSingel["C502019"] = false;
		flagsSingel["C502018"] = false;
	}

	//the columns form merged pool are same as signal pool
	for (const auto& item : columns) {
		std::string mergedPool = item.getName();
		//std::cout << mergedPool << std::endl;
		// item.printRoute();
		std::string pool1 = mergedPool.substr(0, 7);
		std::string pool2 = mergedPool.substr(0, 5) + mergedPool.substr(8, 10);

		std::vector<int> indexes = item.getRoute();
		std::vector<int> indexPool1 = poolIndex[pool1];
		std::vector<int> indexPool2 = poolIndex[pool2];
		bool isFromPool1 = true;
		bool isFromPool2 = true;

		for (const auto& index : indexes) {
			if (std::find(indexPool1.begin(), indexPool1.end(), index) == indexPool1.end()) {
				isFromPool1 = false;
			}
			if (std::find(indexPool2.begin(), indexPool2.end(), index) == indexPool2.end()) {
				isFromPool2 = false;
			}
		}
		if (isFromPool1) {
			flagsSingel[pool1] = false;
		}
		if (isFromPool2) {
			flagsSingel[pool2] = false;
		}
	}

	for (const auto& item : flagsSingel) {
		if (item.second) {
			//bool flag = solverPoolLabelSR(item.first, poolIndex[item.first]);
			bool flag = solverPoolLabelSR_ng_routes(item.first, poolIndex[item.first]);
			solveNum++;
		}
	}

	//remove the duplicate columns
	std::vector<Column> columnsNew;
	for (const auto& item : columns) {
		if (columnsNew.empty()) {
			columnsNew.emplace_back(item);
			continue;
		}
		auto it = std::find_if(columnsNew.begin(), columnsNew.end(),
			[&](const Column& c_in_new_pool) {
				return item.getRoute() == c_in_new_pool.getRoute();
			});
		if (it != columnsNew.end()) {
			continue;
		}
		columnsNew.emplace_back(item);
	}
	columns = columnsNew;

	return columns;
}

void SubProblem::updateSubset(std::vector<int> s)
{
	for (const auto& vec : subset) {
		if (std::equal(vec.begin(), vec.end(), s.begin(), s.end())) {
			return;
		}
	}
	subset.emplace_back(s);
}




