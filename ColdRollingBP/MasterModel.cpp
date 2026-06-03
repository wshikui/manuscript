//
// Created by wshikui on 2023/10/4.
//
#include "MasterProblem.h"

MasterProblem::MasterProblem(std::vector<Column*>& _columns, Problem* _problem, Parameters* _parameters) : BaseModel()
{
	problem = _problem;
	parameters = _parameters;

	// initial the target flow matrix (h_{il}: (N + 1) * L)
	initMatrixH();

	createModel();
	initModel(_columns, _problem);
}

MasterProblem::~MasterProblem()
{
	for (int i = 0; i < problem->getVertices() + 1; i++)
	{
		delete[] h[i];
	}
	delete[] h;
	mEnv.end();
}

void MasterProblem::initMatrixH()
{
	h = new int* [problem->getVertices() + 1];
	for (int i = 0; i < problem->getVertices() + 1; i++)
	{
		h[i] = new int[parameters->L];
		if (i == 0)
		{
			for (int l = 0; l < parameters->L; l++)
			{
				h[i][l] = 1;
			}
			continue;
		}
		if (problem->getFlow(i) == "C512")
		{
			h[i][0] = 1;
			h[i][1] = 0;
		}
		else
		{
			h[i][0] = 0;
			h[i][1] = 1;
		}
	}
}

double MasterProblem::getTargetMass()
{
	double mass = 0.0;
	for (int i = 0; i < parameters->L; i++)
	{
		mass += parameters->MassFlow[i];
	}
	return mass;
}

double MasterProblem::getInitColumnMass(vector<Column*>& _columns, int r, int l)
{
	double mass = 0.0;
	for (int i = 1; i < problem->getVertices() + 1; i++)
	{
		if (_columns[r]->contains(i))
		{
			mass += h[i][l] * problem->getWeight(i);
		}
	}
	return mass;
}

double MasterProblem::getInitColumnMass(vector<Column>& _columns, int r, int l)
{
	double mass = 0.0;
	for (int i = 1; i < problem->getVertices() + 1; i++)
	{
		if (_columns[r].contains(i))
		{
			mass += h[i][l] * problem->getWeight(i);
		}
	}
	return mass;
}

void MasterProblem::initModel(vector<Column*>& _columns, Problem* _problem)
{
	mModel.setName("MasterModel");

	int numVariables = static_cast<int>(_columns.size());
	int numConstraint = _problem->getVertices() + parameters->L + 1 + 1;

	try
	{
		y = IloNumVarArray(mEnv);

		upLimits = IloNumArray(mEnv, numConstraint);
		downLimits = IloNumArray(mEnv, numConstraint);

		// for all vertices exclude virtual vertex 0
		// the lower limit of first |V| constraints should be determined
		for (int i = 0; i < _problem->getVertices(); i++)
		{
			upLimits[i] = 1;
			downLimits[i] = -IloInfinity;
		}
		//for target flow
		for (int i = 0; i < parameters->L; i++)
		{
			upLimits[_problem->getVertices() + i] = IloInfinity;
			downLimits[_problem->getVertices() + i] = parameters->MassFlow[i];
		}
		// for two constraints about plan numbers
		double massAll = getTargetMass();
		double numPlans = std::ceil(massAll / parameters->planUp);
		upLimits[_problem->getVertices() + parameters->L] = IloInfinity;
		downLimits[_problem->getVertices() + parameters->L] = numPlans;

		int maxPlans = getMaxPlanNum();
		upLimits[_problem->getVertices() + parameters->L + 1] = maxPlans;
		downLimits[_problem->getVertices() + parameters->L + 1] = -IloInfinity;

		cons = IloAdd(mModel, IloRangeArray(mEnv, downLimits, upLimits));
		obj = IloAdd(mModel, IloMinimize(mEnv));

		// add variables
		for (int r = 0; r < numVariables; r++)
		{
			IloNumColumn constraintColumn;
			constraintColumn = obj(
				_columns[r]->getCost() * parameters->alpha1 + parameters->alpha2 * parameters->setUpCost);
			for (int i = 0; i < _problem->getVertices(); i++)
			{
				constraintColumn += cons[i](_columns[r]->contains(i + 1));
			}
			for (int i = 0; i < parameters->L; i++)
			{
				double columnMass = getInitColumnMass(_columns, r, i);
				constraintColumn += cons[_problem->getVertices() + i](columnMass);
			}
			constraintColumn += cons[_problem->getVertices() + parameters->L](1);
			constraintColumn += cons[_problem->getVertices() + parameters->L + 1](1);
			y.add(IloNumVar(constraintColumn, 0, 1, ILOFLOAT));
		}
	}
	catch (IloException& ex)
	{
		cerr << "Master Problem InitModel Error -> " << ex << endl;
	}
	catch (...)
	{
		cerr << "Error when Master Problem extractColumn" << endl;
	}
}

/*
 * for all vertices exclude virtual vertex 0
 * the lower limit of first |V| constraints should be determined
 * the index represents the order of constraints, not the index of real coil
 */
void MasterProblem::initMasterModel(vector<Column>& localPool, vector<int>& vertexWithEquality)
{
	mModel.setName("MasterModel");

	int numVariables = static_cast<int>(localPool.size());
	int numConstraint = problem->getVertices() + parameters->L + 1 + 1;

	try
	{
		y = IloNumVarArray(mEnv);

		upLimits = IloNumArray(mEnv, numConstraint);
		downLimits = IloNumArray(mEnv, numConstraint);

		// the first N constraints for selecting vertex
		for (int i = 0; i < problem->getVertices(); i++)
		{
			if (std::find(vertexWithEquality.begin(), vertexWithEquality.end(), i + 1) != vertexWithEquality.end())
			{
				upLimits[i] = 1;
				downLimits[i] = 1;
				continue;
			}
			upLimits[i] = 1;
			downLimits[i] = -IloInfinity;
		}
		//for target flow
		for (int i = 0; i < parameters->L; i++)
		{
			upLimits[problem->getVertices() + i] = IloInfinity;
			downLimits[problem->getVertices() + i] = parameters->MassFlow[i];
		}
		// for two constraints about plan numbers
		double massAll = getTargetMass();
		double numPlans = std::ceil(massAll / parameters->planUp);
		upLimits[problem->getVertices() + parameters->L] = IloInfinity;
		downLimits[problem->getVertices() + parameters->L] = numPlans;
		//downLimits[problem->getVertices() + parameters->L] = -IloInfinity;

		int maxPlans = getMaxPlanNum();
		upLimits[problem->getVertices() + parameters->L + 1] = maxPlans;
		downLimits[problem->getVertices() + parameters->L + 1] = -IloInfinity;

		cons = IloAdd(mModel, IloRangeArray(mEnv, downLimits, upLimits));
		obj = IloAdd(mModel, IloMinimize(mEnv));

		// add variables
		for (int r = 0; r < numVariables; r++)
		{
			IloNumColumn constraintColumn;
			constraintColumn = obj(
				localPool[r].getCost() * parameters->alpha1 + parameters->alpha2 * parameters->setUpCost);
			for (int i = 0; i < problem->getVertices(); i++)
			{
				constraintColumn += cons[i](localPool[r].contains(i + 1));
			}
			for (int i = 0; i < parameters->L; i++)
			{
				double columnMass = getInitColumnMass(localPool, r, i);
				constraintColumn += cons[problem->getVertices() + i](columnMass);
			}
			constraintColumn += cons[problem->getVertices() + parameters->L](1);
			constraintColumn += cons[problem->getVertices() + parameters->L + 1](1);
			y.add(IloNumVar(constraintColumn, 0, 1, ILOFLOAT));
		}
	}
	catch (IloException& ex)
	{
		cerr << "Master Problem InitModel Error -> " << ex << endl;
	}
	catch (...)
	{
		cerr << "Error when Master Problem extractColumn" << endl;
	}
}


void MasterProblem::initMasterModel(std::vector<Column>& localPool, std::vector<int>& vertexWithEquality, bool integer)
{
	mModel.setName("MasterModel_MIP_UP_Bound");

	int numVariables = static_cast<int>(localPool.size());
	int numConstraint = problem->getVertices() + parameters->L + 1 + 1;

	try
	{
		y = IloNumVarArray(mEnv);

		upLimits = IloNumArray(mEnv, numConstraint);
		downLimits = IloNumArray(mEnv, numConstraint);

		// the first N constraints for selecting vertex
		for (int i = 0; i < problem->getVertices(); i++)
		{
			if (std::find(vertexWithEquality.begin(), vertexWithEquality.end(), i + 1) != vertexWithEquality.end())
			{
				upLimits[i] = 1;
				downLimits[i] = 1;
				continue;
			}
			upLimits[i] = 1;
			downLimits[i] = -IloInfinity;
		}
		//for target flow
		for (int i = 0; i < parameters->L; i++)
		{
			upLimits[problem->getVertices() + i] = IloInfinity;
			downLimits[problem->getVertices() + i] = parameters->MassFlow[i];
		}
		// for two constraints about plan numbers
		double massAll = getTargetMass();
		double numPlans = std::ceil(massAll / parameters->planUp);
		upLimits[problem->getVertices() + parameters->L] = IloInfinity;
		downLimits[problem->getVertices() + parameters->L] = numPlans;

		int maxPlans = getMaxPlanNum();
		upLimits[problem->getVertices() + parameters->L + 1] = maxPlans;
		downLimits[problem->getVertices() + parameters->L + 1] = -IloInfinity;

		cons = IloAdd(mModel, IloRangeArray(mEnv, downLimits, upLimits));
		obj = IloAdd(mModel, IloMinimize(mEnv));

		// add variables
		for (int r = 0; r < numVariables; r++)
		{
			IloNumColumn constraintColumn;
			constraintColumn = obj(
				localPool[r].getCost() * parameters->alpha1 + parameters->alpha2 * parameters->setUpCost);
			for (int i = 0; i < problem->getVertices(); i++)
			{
				constraintColumn += cons[i](localPool[r].contains(i + 1));
			}
			for (int i = 0; i < parameters->L; i++)
			{
				double columnMass = getInitColumnMass(localPool, r, i);
				constraintColumn += cons[problem->getVertices() + i](columnMass);
			}
			constraintColumn += cons[problem->getVertices() + parameters->L](1);
			constraintColumn += cons[problem->getVertices() + parameters->L + 1](1);
			y.add(IloNumVar(constraintColumn, 0, 1, (integer ? ILOINT : ILOFLOAT)));
		}
	}
	catch (IloException& ex)
	{
		cerr << "Master Problem InitModel Error -> " << ex << endl;
	}
	catch (...)
	{
		cerr << "Error when Master Problem extractColumn" << endl;
	}
}

void MasterProblem::initMasterModel(std::vector<Column>& localPool, std::vector<int>& vertexWithEquality, const std::vector<std::vector<int>>& subSet)
{
	mModel.setName("MasterModelWithSR");

	int numVariables = static_cast<int>(localPool.size());
	int numSR = static_cast<int>(subSet.size());
	int numConstraint = problem->getVertices() + parameters->L + 1 + 1;
	numConstraint += numSR; // the number of SR inequality

	try {
		y = IloNumVarArray(mEnv);

		upLimits = IloNumArray(mEnv, numConstraint);
		downLimits = IloNumArray(mEnv, numConstraint);

		// the first N constraints for selecting vertex
		for (int i = 0; i < problem->getVertices(); i++)
		{
			if (std::find(vertexWithEquality.begin(), vertexWithEquality.end(), i + 1) != vertexWithEquality.end())
			{
				upLimits[i] = 1;
				downLimits[i] = 1;
				continue;
			}
			upLimits[i] = 1;
			downLimits[i] = -IloInfinity;
		}
		//for target flow
		for (int i = 0; i < parameters->L; i++)
		{
			upLimits[problem->getVertices() + i] = IloInfinity;
			downLimits[problem->getVertices() + i] = parameters->MassFlow[i];
		}
		// for two constraints about plan numbers
		double massAll = getTargetMass();
		double numPlans = std::ceil(massAll / parameters->planUp);
		upLimits[problem->getVertices() + parameters->L] = IloInfinity;
		downLimits[problem->getVertices() + parameters->L] = numPlans;
		//downLimits[problem->getVertices() + parameters->L] = -IloInfinity;

		int maxPlans = getMaxPlanNum();
		upLimits[problem->getVertices() + parameters->L + 1] = maxPlans;
		downLimits[problem->getVertices() + parameters->L + 1] = -IloInfinity;

		//for SR, if subset is empty, skip modeling the SR
		for (int i = 0; i < numSR; i++) {
			upLimits[problem->getVertices() + parameters->L + 2 + i] = 1;
			downLimits[problem->getVertices() + parameters->L + 2 + i] = -IloInfinity;
		}

		cons = IloAdd(mModel, IloRangeArray(mEnv, downLimits, upLimits));
		obj = IloAdd(mModel, IloMinimize(mEnv));

		// add variables
		for (int r = 0; r < numVariables; r++)
		{
			IloNumColumn constraintColumn;
			constraintColumn = obj(
				localPool[r].getCost() * parameters->alpha1 + parameters->alpha2 * parameters->setUpCost);
			for (int i = 0; i < problem->getVertices(); i++)
			{
				constraintColumn += cons[i](localPool[r].contains(i + 1));
			}
			for (int i = 0; i < parameters->L; i++)
			{
				double columnMass = getInitColumnMass(localPool, r, i);
				constraintColumn += cons[problem->getVertices() + i](columnMass);
			}
			constraintColumn += cons[problem->getVertices() + parameters->L](1);
			constraintColumn += cons[problem->getVertices() + parameters->L + 1](1);
			//SR inequality
			for (int i = 0; i < numSR; i++) {
				int num = 0;
				for (const auto& item : subSet[i]) {
					if (localPool[r].contains(item)) {
						num++;
					}
				}
				if (num >= 2) {
					constraintColumn += cons[problem->getVertices() + parameters->L + 2 + i](1); //2503
				}
			}
			y.add(IloNumVar(constraintColumn, 0, 1, ILOFLOAT));
		}
	}
	catch (IloException& ex)
	{
		cerr << "Master Problem InitModel Error -> " << ex << endl;
	}
	catch (...)
	{
		cerr << "Error when Master Problem extractColumn" << endl;
	}


}

/*
 * need to determine the cost by x{ij} or vertices in this column
 * need to get original objective for c{r}y{r} to c{ij}x{ij} and e{i}y{i}
 */
void MasterProblem::addColumn(const Column& _column)
{
	IloNumColumn cplexColumn;
	cplexColumn = obj(_column.getCost() * parameters->alpha1 + parameters->alpha2 * parameters->setUpCost);
	// for real coils
	for (int i = 0; i < problem->getVertices(); i++)
	{
		cplexColumn += cons[i](_column.contains(i + 1));
	}
	// for target flow
	for (int l = 0; l < parameters->L; l++)
	{
		double mass = 0.0;
		// all real coils
		for (int i = 1; i < problem->getVertices() + 1; i++)
		{
			if (_column.contains(i) && problem->getFlow(i) == parameters->target_flow[l])
			{
				mass += problem->getWeight(i);
			}
		}
		cplexColumn += cons[problem->getVertices() + l](mass);
	}
	cplexColumn += cons[problem->getVertices() + parameters->L](1);

	//TODO 2503
	cplexColumn += cons[problem->getVertices() + parameters->L + 1](1);
	y.add(IloNumVar(cplexColumn, 0, 1, ILOFLOAT));
}

void MasterProblem::getDual(IloNumArray& variables)
{
	int numConstraint = problem->getVertices() + parameters->L + 1 + 1;
	for (int i = 0; i < numConstraint; i++)
	{
		variables.add(cplex.getDual(cons[i]));
	}
}

int MasterProblem::getVarNum(int i)
{
	IloExpr constraint = cons[i].getExpr();
	int numVariables = 0;
	for (IloExpr::LinearIterator it = constraint.getLinearIterator(); it.ok(); ++it) {
		numVariables++;
	}
	return numVariables;
}

std::tuple<std::vector<double>, std::vector<double>> MasterProblem::getCoefAndValue(int i)
{
	std::vector<double> coefs;
	std::vector<double> yr;
	for (int j = 0; j < columnPool.size(); j++) {
		//get the j-th coefficient and value in the objective function
		//columnPool[j].printRoute();
		if (columnPool[j].contains(i)) {
			for (IloExpr::LinearIterator it = obj.getLinearIterator(); it.ok(); ++it) {
				if (it.getVar().getId() == y[j].getId()) {
					coefs.emplace_back(it.getCoef());
					yr.emplace_back(cplex.getValue(y[j]));
					break;
				}
			}
		}
	}

	//for (IloExpr::LinearIterator it = obj.getLinearIterator(); it.ok(); ++it) {
	//	if (it.getVar().getId() == y[i].getId()) {
	//		coefs.emplace_back(it.getCoef());
	//		yr.emplace_back(cplex.getValue(y[i]));
	//	}
	//}
	return std::make_tuple(coefs, yr);
}

double MasterProblem::getWeightRH(int i)
{
	double result = 0;
	for (int j = 0; j < columnPool.size(); j++) {
		if (columnPool[j].contains(i)) {
			result += cplex.getValue(y[j]) * problem->getWeight(i);
		}
	}
	if (problem->getFlow(i) == "C512") {
		result = result / downLimits[problem->getVertices()];
	}
	else {
		result = result / downLimits[problem->getVertices() + 1];
	}
	return result;
}

double MasterProblem::getNumSRofCoil(int i)
{
	int num = 0;
	for (const auto& S : subSet) {
		for (const auto& index : S) {
			if (index == i) {
				num++;
				break;
			}
		}
	}
	if (subSet.size() == 0) {
		return 0;
	}
	double result = static_cast<double>(num / subSet.size());
	return result;
}

double MasterProblem::getBCconstraint(int i)
{
	int num = 0;
	IloExpr constraint_coil = cons[i].getExpr();
	for (IloExpr::LinearIterator it = constraint_coil.getLinearIterator(); it.ok(); ++it) {
		if (cplex.getValue(it.getVar()) != 0) {
			num++;
		}
	}

	return static_cast<double>(num / parameters->K);
}

std::tuple<double, double, double, double, double> MasterProblem::isBinding(int i)
{
	// i - the index of coil selection constraint

	std::vector<int> nums;
	//the number of columns in the i-th constraint
	//int numVar = getVarNum(i);
	IloExpr cons_coil = cons[i].getExpr();

	for (int j = 0; j < cons.getSize(); j++) {
		// a binding constraint
		if (cplex.getDual(cons[j]) != 0) {
			int num = 0;
			IloExpr binding = cons[j].getExpr();
			//the columns in the constraint appears in binding constraint
			for (IloExpr::LinearIterator it = cons_coil.getLinearIterator(); it.ok(); ++it) {
				for (IloExpr::LinearIterator itb = binding.getLinearIterator(); itb.ok(); ++itb) {
					if (itb.getVar().getId() == it.getVar().getId()) {
						num++;
						break;
					}
				}
			}
			nums.emplace_back(num);
		}
	}
	double count = nums.size();
	if (count == 0) {
		return std::make_tuple(static_cast<double>(count / parameters->N * 2), 0, 0, 0, 0);
	}
	else {
		//normalization
		std::vector<double> _nums;
		double _count = static_cast<double>(count / parameters->N * 2);
		for (const auto& item : nums) {
			_nums.emplace_back(static_cast<double>(item) / (parameters->N * 2));
		}
		double max = *(std::max_element(_nums.begin(), _nums.end()));
		double min = *(std::min_element(_nums.begin(), _nums.end()));
		double mean = std::accumulate(_nums.begin(), _nums.end(), 0.0) / static_cast<double>(_nums.size());
		double variance = Tool::calculateVariance(_nums);
		return std::make_tuple(_count, max, min, mean, variance);
	}
}

void MasterProblem::addColumnSR(const Column& _column)
{
	IloNumColumn cplexColumn;
	cplexColumn = obj(_column.getCost() * parameters->alpha1 + parameters->alpha2 * parameters->setUpCost);
	// for real coils
	for (int i = 0; i < problem->getVertices(); i++)
	{
		//cplexColumn += cons[i](_column.contains(i + 1));
		cplexColumn += cons[i](_column.contain_num(i + 1));
	}
	// for target flow
	for (int l = 0; l < parameters->L; l++)
	{
		double mass = 0.0;
		// all real coils
		for (int i = 1; i < problem->getVertices() + 1; i++)
		{
			if (_column.contains(i) && problem->getFlow(i) == parameters->target_flow[l])
			{
				//mass += problem->getWeight(i);
				mass += (problem->getWeight(i) * _column.contain_num(i));
			}
		}
		cplexColumn += cons[problem->getVertices() + l](mass);
	}
	//for vehicle number
	cplexColumn += cons[problem->getVertices() + parameters->L](1);
	cplexColumn += cons[problem->getVertices() + parameters->L + 1](1);
	//for SR 
	//for (int i = 0; i < subSet.size(); i++) {
	//	int num = 0;
	//	std::set<int> record;
	//	for (const auto& item : subSet[i]) {
	//		if (_column.contains(item)) {
	//			record.insert(item);
	//		}
	//	}
	//	num = record.size();
	//	if (num >= 2) {
	//		cplexColumn += cons[problem->getVertices() + parameters->L + 2 + i](1);
	//	}
	//}


	for (int i = 0; i < subSet.size(); i++) {
		int num = 0;
		std::set<int> record;
		for (const auto& item : subSet[i]) {
			if (_column.contains(item)) {
				record.insert(item);
			}
		}
		num = record.size();
		if (num >= 2) {
			cplexColumn += cons[problem->getVertices() + parameters->L + 2 + i](1);
		}
	}

	//for (int i = 0; i < subSet.size(); i++) {
	//	int num = 0;
	//	for (const auto& item : subSet[i]) {
	//		if (_column.contains(item)) {
	//			num++;
	//		}
	//	}
	//	if (num >= 2) {
	//		cplexColumn += cons[problem->getVertices() + parameters->L + 2 + i](1);
	//	}
	//}

	y.add(IloNumVar(cplexColumn, 0, 1, ILOFLOAT));
	columnPool.emplace_back(_column);
}

double MasterProblem::getDualVariable(int i)
{
	return cplex.getDual(cons[i]);
}

int MasterProblem::getMaxPlanNum()
{
	double mass = 0.0;
	for (int i = 0; i < problem->getVertices(); i++)
	{
		mass += problem->getWeight(i + 1);
	}
	int num = static_cast<int>(std::floor(mass / parameters->planLower));
	return num;
}

void MasterProblem::printResult(int iter)
{
	bool isMore = true;
	if (isMore) {
		std::cout << "____________" << iter << " variable Y____________" << std::endl;
		for (IloInt i = 0; i < y.getSize(); i++)
		{
			std::cout << "y[" << i << "] : " << cplex.getValue(y[i]) << std::endl;
		}
		std::cout << "____________dual variable____________" << std::endl;
		int index = 0;

		for (int i = 0; i < problem->getVertices(); i++)
		{
			if (i % parameters->S == 0 && i != 0)
			{
				std::cout << "---------------" << "CHS -> " << i / parameters->S << "--------------" << std::endl;
			}
			std::cout << "vertex_" << i + 1 << " : " << cplex.getDual(cons[i]) << std::endl;
		}
		std::cout << "C512 constraint : " << cplex.getDual(cons[problem->getVertices() + 0]) << std::endl;
		std::cout << "C008 constraint : " << cplex.getDual(cons[problem->getVertices() + 1]) << std::endl;
		std::cout << "plan lower : " << cplex.getDual(cons[problem->getVertices() + parameters->L]) << std::endl;
		std::cout << "plan upper : " << cplex.getDual(cons[problem->getVertices() + parameters->L + 1]) << std::endl;
		if (!subSet.empty()) {
			for (int i = 0; i < subSet.size(); i++) {
				std::cout << cplex.getDual(cons[problem->getVertices() + parameters->L + 2 + i]) << std::endl;
			}
		}
		std::cout << "____________" << iter << "  iteration -> objective = " << cplex.getObjValue() << std::endl;

		int numPlans = 0;
		double massC512 = 0.0;
		double massC008 = 0.0;

		std::cout << "print the plan info" << std::endl;
		for (IloInt i = 0; i < y.getSize(); i++) {
			if (cplex.getValue(y[i]) > 0) {
				numPlans++;
				//std::cout << "plan " << i << " : " << columnPool[i].getDemand() << std::endl;
				for (int j = 1; j < problem->getVertices() + 1; j++)
				{
					if (columnPool[i].contains(j) == 1)
					{
						if (problem->getFlow(j) == "C512")
						{
							massC512 += problem->getWeight(j) * cplex.getValue(y[i]);
						}
						else
						{
							massC008 += problem->getWeight(j) * cplex.getValue(y[i]);
						}
					}
				}
			}
		}
		std::cout << "num of plan : " << numPlans << std::endl;
		std::cout << "mass : " << massC512 << " " << massC008 << std::endl;
	}
	else {
		std::cout << "____________" << iter << "  iteration -> objective = " << cplex.getObjValue() << std::endl;

	}

}

IloNumVarArray& MasterProblem::getY()
{
	return y;
}

IloRangeArray& MasterProblem::getCons()
{
	return cons;
}

std::vector<std::pair<Column, double>>& MasterProblem::getBasicColumns()
{
	basicColumns.clear();

	for (int i = 0; i < y.getSize(); i++) {
		if (cplex.getValue(y[i]) > 0) {
			basicColumns.emplace_back(std::make_pair(columnPool[i], cplex.getValue(y[i])));
		}
	}
	return basicColumns;
}

std::vector<std::vector<int>>& MasterProblem::getSubSet()
{
	return subSet;
}

void MasterProblem::updateSubset(std::vector<int> s)
{
	for (const auto& vec : subSet) {
		if (std::equal(vec.begin(), vec.end(), s.begin())) {
			return;
		}
	}
	//if this S is new, save it
	subSet.emplace_back(s);
}


//IloObjective& MasterProblem::getObj()
//{
// contradiction with CPLEX API
//	return obj;
//}


MasterProblem::MasterProblem(vector<Column>& localPool, vector<int>& vertexWithEquality, Problem* _problem,
	Parameters* _parameters)
{
	problem = _problem;
	parameters = _parameters;

	// initialize h_{il}
	initMatrixH();

	createModel();

	initMasterModel(localPool, vertexWithEquality);
	//initMasterModel_no_cut(localPool, vertexWithEquality);
}

//MasterProblem::MasterProblem(std::vector<Column>& localPool, std::vector<int>& vertexWithEquality, Problem* _problem, Parameters* _parameters, int cuts)
//{
//	problem = _problem;
//	parameters = _parameters;
//
//	// initialize h_{il}
//	initMatrixH();
//
//	createModel();
//
//	if (cuts == 0) {
//		initMasterModel_no_cut(localPool, vertexWithEquality);
//	}
//	else {
//		initMasterModel(localPool, vertexWithEquality);
//	}
//}

MasterProblem::MasterProblem(std::vector<Column>& localPool, std::vector<int>& vertexWithEquality, Problem* _problem, Parameters* _parameters, bool isInteger)
{
	problem = _problem;
	parameters = _parameters;

	// initialize h_{il}
	initMatrixH();

	createModel();

	initMasterModel(localPool, vertexWithEquality, true);
}

MasterProblem::MasterProblem(std::vector<Column>& localPool, std::vector<int>& vertexWithEquality,
	const std::vector<std::vector<int>>& _subSet,
	Problem* _problem, Parameters* _parameters)
{
	problem = _problem;
	parameters = _parameters;

	subSet = _subSet;
	columnPool = localPool;

	// initialize h_{il}
	initMatrixH();

	createModel();

	initMasterModel(localPool, vertexWithEquality, _subSet);
}

