#include "ColumnGeneration.h"

ColumnGeneration::ColumnGeneration(std::vector<Column>& columnPool, std::vector<int>& vertexEquality, 
	std::vector<int>& vertexZero, std::map<std::string, std::vector<int>>& poolIndex, 
	std::map<std::string, std::vector<std::bitset<Config::MAX_COILS>>> pool_coil_ngsets, Parameters* _parameters, Problem* _problem)
{
	parameters = _parameters;
	problem = _problem;

	rmp = new MasterProblem(columnPool, vertexEquality, problem, parameters);
	subProblem = new SubProblem(rmp, vertexZero, poolIndex, pool_coil_ngsets, problem, parameters);

	for (const auto& item : columnPool)
	{
		columnPoolCG.emplace_back(item);
	}
}


ColumnGeneration:: ColumnGeneration(std::vector<Column>& columnPool, std::vector<int>& vertexEquality, std::vector<int>& vertexZero, 
	const std::vector<std::vector<int>>& subset, 
	std::map<int, std::vector<int>> _forbidCoils,
	std::map<std::string, std::vector<int>>& poolIndex, 
	std::map<std::string, std::vector<std::bitset<Config::MAX_COILS>>> pool_coil_ngsets, 
	Parameters* _parameters, Problem* _problem)
{
	parameters = _parameters;
	problem = _problem;

	//modify the RMP based on SR
	rmp = new MasterProblem(columnPool, vertexEquality, subset, problem, parameters);
	subProblem = new SubProblem(rmp, vertexZero, _forbidCoils, subset, poolIndex, pool_coil_ngsets, problem, parameters);

	for (const auto& item : columnPool)
	{
		columnPoolCG.emplace_back(item);
	}
}

ColumnGeneration::~ColumnGeneration()
{
	delete rmp;
	delete subProblem;
}

bool ColumnGeneration::generateColumns()
{
	bool flag = false;
	// this column's reduced cost maybe positive
	auto startSP = std::chrono::high_resolution_clock::now();

	//get new columns from sub-problem
	//std::vector<Column> columns = subProblem->extractColumns();
	std::vector<Column> columns = subProblem->extractColumnsLabel();

	auto endSP = std::chrono::high_resolution_clock::now();
	timeSP += std::chrono::duration_cast<std::chrono::duration<double>>(endSP - startSP).count();

	std::cout << "___________________generateColumn___________________ " << std::endl;
	for (const auto& item : columns) {
		if (item.getRC() < -parameters->epsilon) {
			item.printRoute();
			rmp->addColumn(item);
			columnPoolCG.emplace_back(item);
			flag = true;
		}
	}

	std::cout << "_________generate columns end and sub-problem was solver in " <<
		std::chrono::duration_cast<std::chrono::duration<double>>(endSP - startSP).count() << std::endl;

	return flag;
}

bool ColumnGeneration::generateColumnsAndSR()
{
	bool flag = false;
	auto startSP = std::chrono::high_resolution_clock::now();

	//get new columns from sub-problem with SR constraint
	//std::vector<Column> columns = subProblem->extractColumnsSR();

	std::vector<Column> columns = subProblem->extractColumnsLabelSR();

	//the number of solving the sub-problems
	solveSub += subProblem->solveNum;
	subProblem->reSetSubNum();

	auto endSP = std::chrono::high_resolution_clock::now();
	timeSP += std::chrono::duration_cast<std::chrono::duration<double>>(endSP - startSP).count();

	//std::cout << "___________________generateColumn in sub-problems___________________ " << std::endl;

	int newColumns = 0;
	for (const auto& item : columns) {
		if (item.getRC() < -parameters->epsilon) {
			//std::cout << "reduced cost -> " << item.getRC() << std::endl;
			//item.printRoute();
			
			//change this function when consider SR
			newColumns++;
			rmp->addColumnSR(item);
			columnPoolCG.emplace_back(item);
			flag = true;
		}
	}
	
	//std::cout << "create " << newColumns << " columns in subproblem" << std::endl;

	return flag;
}

bool ColumnGeneration::isSameColumn(const Column& column1, const Column& column2)
{
	std::vector<int> route1 = column1.getRoute();
	std::vector<int> route2 = column2.getRoute();
	if (route1 == route2) {
		return true;
	}
	else {
		return false;
	}
}

void ColumnGeneration::remove_duplicate_column()
{
	std::vector<Column> newPool;
	for (const auto& column : columnPoolCG) {
		if (newPool.empty()) {
			newPool.emplace_back(column);
			continue;
		}
		auto nc = std::find_if(newPool.begin(), newPool.end(),
			[&](const Column& c_in_new_pool) {
				return column.getRoute() == c_in_new_pool.getRoute();
			});

		if (nc != newPool.end()) {
			continue;
		}
		newPool.emplace_back(column);
	}
	//modify the column pool in CG
	columnPoolCG = newPool;
}

void ColumnGeneration::printPool()
{
	for (const auto& column : columnPoolCG) {
		column.printRoute();
	}
}

bool ColumnGeneration::execute()
{
	// the termination condition considers running time and reduce cost
	Stopwatch stopwatch("ColumnGeneration_execute");
	bool solveRMP;
	int num = 0;
	do
	{
		auto startMP = std::chrono::high_resolution_clock::now();
		solveRMP = rmp->optimize();
		auto endMP = std::chrono::high_resolution_clock::now();

		timeMP += std::chrono::duration_cast<std::chrono::duration<double>>(endMP - startMP).count();

		// std::string modelName = "../../../../Log/MasterModel_" + std::to_string(num) + ".lp";
		// rmp->getSolver().exportModel(modelName.c_str());

		if (solveRMP) {
			//rmp->printResult(num); 
		}
		num++;
		if (num > 600) {
			std::cout << "may error in column generation" << std::endl;
		}
	} while (solveRMP &&
		stopwatch.elapsed() < parameters->_timeLimit &&
		generateColumns());

	std::cout << "column generation time : " << stopwatch.elapsed() << std::endl;
	std::cout << "column generation end" << std::endl;

	if (!solveRMP) {
		//std::cout << "column generation -> RMP is failed" << std::endl;
		return false;
	}

	getResult();
	setBasicColumns();
	remove_duplicate_column();

	return true;
}

bool ColumnGeneration::executeAndCutR_depth(unsigned int node_number, int depth)
{

	Stopwatch stopwatch("ColumnGeneration_execute");
	bool solveRMP;
	int num = 0;
	do
	{
		auto startMP = std::chrono::high_resolution_clock::now();
		solveRMP = rmp->optimize();
		auto endMP = std::chrono::high_resolution_clock::now();

		timeMP += std::chrono::duration_cast<std::chrono::duration<double>>(endMP - startMP).count();

		// std::string modelName = "../../../../Log/MasterModel_" + std::to_string(num) + ".lp";
		// rmp->getSolver().exportModel(modelName.c_str());

		//solve the separation problem based on current columnsCG___
		if (solveRMP && depth < 5) {
			//separation
			std::vector<std::pair<Column, double>> basicColumns = rmp->getBasicColumns();
			SPmodel spModel(basicColumns);
			//bool flagSR = spModel.getS(basicColumns);
			bool flagSR = spModel.getSenum(basicColumns);
			if (flagSR) {
				//change RMP(address) and resolve it
				IloExpr conSP(rmp->getEnv());

				IloRangeArray cons = rmp->getCons();
				IloNumVarArray y = rmp->getY();
				for (int i = 0; i < rmp->getY().getSize(); i++) {
					int num = 0;
					for (const auto& item : spModel.S) {
						IloRange vertexCons = cons[item - 1];
						IloExpr expr(rmp->getEnv());
						expr = vertexCons.getExpr();
						for (IloExpr::LinearIterator it = expr.getLinearIterator(); it.ok(); ++it) {
							if (it.getVar().getId() == y[i].getId()) {
								num++;
								expr.end();
								break;
							}
						}
						expr.end();
					}
					if (num >= 2) {
						conSP += y[i];
					}
				}
				//update the cons(ILoRangeArray)
				IloRange rowConstraint(rmp->getEnv(), -IloInfinity, conSP, 1);
				//the cons in rmp have changed
				cons.add(rowConstraint);
				rmp->getModel().add(rowConstraint);
				conSP.end();

				solveRMP = rmp->optimize();
				//std::cout << rmp->getModel().getName() << std::endl;
				//std::cout << rmp->getCons().getSize() << std::endl;

				//update all S for RMP and sub-problem
				rmp->updateSubset(spModel.S);
				subProblem->updateSubset(spModel.S);
			}
		}
		if (solveRMP) {
			//rmp->printResult(num);
		}
		num++;
		if (num > 60) {
			std::cout << "maybe error in cg" << std::endl;
			solveRMP = true;
			break;
		}
	} while (solveRMP &&
		stopwatch.elapsed() < parameters->_timeLimit &&
		generateColumnsAndSR());


	//std::cout << "column generation time:" << stopwatch.elapsed() << std::endl;
	//std::cout << "column generation end" << std::endl;
	//std::cout << "the number of solving sub-problems:" << solveSub << std::endl;

	if (!solveRMP) {
		//std::cout << "column generation -> RMP is failed" << std::endl;
		return false;
	}

	getResult();
	setBasicColumns();
	remove_duplicate_column();

	return true;
}


void ColumnGeneration::getResult()
{
	IloNumVarArray y = rmp->getY();
	IloCplex cplex = rmp->getSolver();

	// mass of target unit in columns
	double massC512 = 0.0;
	double massC008 = 0.0;
	// for all columns
	for (int i = 0; i < columnPoolCG.size(); i++)
	{
		if (cplex.getValue(y[i]) > 0)
		{
			// for all real vertices in a column
			for (int j = 1; j < problem->getVertices() + 1; j++)
			{
				if (columnPoolCG[i].contains(j) == 1)
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
	//std::cout << "column generation mass result: " << massC512 << " " << massC008 << std::endl;
	//std::cout << "______________the plans information of column generation______________" << std::endl;
	int numPlans = 0;
	for (int i = 0; i < columnPoolCG.size(); i++)
	{
		if (cplex.getValue(y[i]) > 0)
		{
			double rc = 0;
			double weightC512 = 0;
			double weightC008 = 0;
			std::vector<int> route = columnPoolCG[i].getRoute();

			//std::cout << "y[" << i << "] = " << cplex.getValue(y[i]) << std::endl;

			//columnPoolCG[i].printRoute();

			rc += (parameters->alpha1 * columnPoolCG[i].getCost() + parameters->alpha2 * parameters->setUpCost);
			for (const auto& item : route) {
				rc -= rmp->getDualVariable(item - 1);
				if (problem->getFlow(item) == "C512") {
					weightC512 += problem->getWeight(item);
				}
				else {
					weightC008 += problem->getWeight(item);
				}
			}
			rc -= rmp->getDualVariable(problem->getVertices()) * weightC512;
			rc -= rmp->getDualVariable(problem->getVertices() + 1) * weightC008;
			rc -= rmp->getDualVariable(problem->getVertices() + parameters->L);
			rc -= rmp->getDualVariable(problem->getVertices() + parameters->L + 1);
			//std::cout << "actural reduced cost = " << rc << std::endl;
			numPlans++;
		}
	}
	//std::cout << "number of plans -> " << numPlans << std::endl;

	//std::cout << "______________vertex constraint______________" << std::endl;
	// for all real vertices
	bool isVertex = true;
	for (int i = 1; i < problem->getVertices() + 1; i++)
	{
		double result = 0;
		for (int r = 0; r < columnPoolCG.size(); r++)
		{
			if (cplex.getValue(y[r]) > 0)
			{
				result += columnPoolCG[r].contains(i) * cplex.getValue(y[r]);
			}
		}
		//std::cout << "vertex -> " << i << " " << result << std::endl;
		if (result > 1 + 1e-5)
		{
			isVertex = false;
		}
	}
	// check the constraints in master problem
	//std::cout << "___________check the feasible of Column generation___________" << std::endl;
	//if (massC008 >= parameters->planLower)
	//{
	//	std::cout << "target weight constraint of C008 -> " << 1 << std::endl;
	//}
	//else
	//{
	//	std::cout << "target weight constraint of C008 -> " << 0 << std::endl;
	//}
	//if (massC512 >= parameters->planLower)
	//{
	//	std::cout << "target weight constraint of C512 -> " << 1 << std::endl;
	//}
	//else
	//{
	//	std::cout << "target weight constraint of C512 -> " << 0 << std::endl;
	//}
	//std::cout << "real vertex constraint -> " << isVertex << std::endl;

	// calculate the objective of column generation: cost and penalty of plans number
	// the getCost function of class column returns coil connection cost

	//Add check to reduced cost
	double cost = 0.0;
	for (int i = 0; i < columnPoolCG.size(); i++)
	{
		if (cplex.getValue(y[i]) > 0)
		{
			cost += (parameters->alpha1 * columnPoolCG[i].getCost() + parameters->alpha2 * parameters->setUpCost) *
				cplex.getValue(y[i]);
			
		}
	}

	//std::cout << "objective of column generation: " << cost << std::endl;

}


void ColumnGeneration::getBasicColumns(std::vector<std::pair<Column, double>>& pool)
{
	pool.clear();
	pool.insert(pool.end(), basicColumns.begin(), basicColumns.end());
}

void ColumnGeneration::setBasicColumns()
{
	IloNumVarArray y = rmp->getY();
	IloCplex cplex = rmp->getSolver();
	for (int i = 0; i < y.getSize(); i++) {
		double value = cplex.getValue(y[i]);
		if (value > 0) {
			basicColumns.emplace_back(std::make_pair(columnPoolCG[i], value));
		}
	}

}


void ColumnGeneration::updateColumnPool(std::vector<Column>& pool)
{
	pool.clear();
	pool.insert(pool.end(), columnPoolCG.begin(), columnPoolCG.end());
}

void ColumnGeneration::getTime(double& valMP, double& valSp)
{
	valMP = timeMP;
	valSp = timeSP;
}

double ColumnGeneration::getSolValue()
{
	double value = rmp->getObj();
	return value;
}

MasterProblem* ColumnGeneration::getRMPmodel()
{
	return rmp;
}

//std::map<int, std::vector<std::string>> ColumnGeneration::getPoolInfo() const
//{
//	return poolInfo;
//}



