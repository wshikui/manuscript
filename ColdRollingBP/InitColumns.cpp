//
// Created by wshikui on 2023/10/3.
//

#include "InitColumns.h"

InitColumns::InitColumns(const std::vector<Material*>& _materials, Problem* _problem, Parameters* _parameters)
{
	// all materials
	materials = _materials;
	mpProblem = _problem;
	parameters = _parameters;

	// single pool materials
	for (const auto& item : _materials)
	{
		chsMaterials[item->getStr(Attribute::CHS_TYPE)].emplace_back(item);
	}
}

std::map<std::string, std::vector<Material*>> InitColumns::getConnectPool(const std::string& _poolName)
{
	int index = -1;
	for (const auto& item : parameters->chsType)
	{
		if (item.first == _poolName)
		{
			index = item.second;
			break;
		}
	}
	if (index == -1)
	{
		std::cout << "InitColumns : Error in get Connect pool materials" << std::endl;
	}

	std::map<std::string, std::vector<Material*>> connectPoolsMaterials;
	for (int i = 0; i < parameters->C; i++)
	{
		int flag = parameters->steelArr[index][i];
		if (flag == 1 && i != index)
		{
			int _conIndex = i;
			std::string conName = parameters->chs_type[_conIndex];
			std::string _newName = parameters->chs_type[index] + "_" + parameters->chs_type[_conIndex];
			connectPoolsMaterials[_newName].insert(connectPoolsMaterials[_newName].end(),
				chsMaterials[_poolName].begin(), chsMaterials[_poolName].end());
			connectPoolsMaterials[_newName].insert(connectPoolsMaterials[_newName].end(),
				chsMaterials[conName].begin(), chsMaterials[conName].end());
		}

	}
	return connectPoolsMaterials;
}


std::vector<Column*> InitColumns::getInitPool()
{
	return initPool;
}

double InitColumns::totalWeight(std::vector<Material*>& materialSet)
{
	if (materialSet.empty())
	{
		std::cout << "WARNING -> InitColumn -> totalWeight -> input empty " << std::endl;
		return 0.0;
	}
	double weights = 0.0;
	for (const auto& item : materialSet)
	{
		weights += item->getValue(Attribute::IN_MAT_WT);
	}
	return weights;
}

bool InitColumns::isFinished(std::map<std::string, double>& _currentFlowTarget)
{
	if (_currentFlowTarget.empty())
	{
		return true;
	}

	double current = 0.0;
	for (const auto& item : initPool)
	{
		current += item->getDemand();
	}

	double total = 0.0;
	for (int i = 0; i < parameters->L; i++)
	{
		total += parameters->MassFlow[i];
	}

	return current > total;
}

bool InitColumns::isFinishedAll(std::map<std::string, double>& _currentFlowTarget)
{
	if (_currentFlowTarget.empty())
	{
		return true;
	}
	if (initPool.size() < 8) {
		return false;
	}

	std::map<std::string, double> current_mass;
	for (const auto& item : parameters->flow_mass) {
		std::string flow = item.first;
		double value = 0;

		//get the weight of the flow in all columns
		for (const auto& col : initPool) {
			for (int i = 1; i < mpProblem->getVertices() + 1; i++) {
				if (col->contains(i)) {
					if (mpProblem->getFlow(i) == flow) {
						value += mpProblem->getWeight(i);
					}
				}
			}
		}
		current_mass[flow] = value;
	}
	for (const auto& item : parameters->flow_mass) {
		std::string flow = item.first;
		if (item.second > current_mass[flow]) {
			return false;
		}
	}
	return true;
}

void InitColumns::flowPool(std::map<std::string, double>& _currentFlowTarget)
{
	// update each running, get coils according to the current target flow
	chsMaterial.clear();

	for (int i = 0; i < mpProblem->getVertices(); i++)
	{
		// coil is not selected
		if (!materials[i]->getStatus())
		{
			std::string chs = materials[i]->getStr(Attribute::CHS_TYPE);
			std::string flow = materials[i]->getStr(Attribute::FLOW);
			if (_currentFlowTarget.count(flow) == 0)
			{
				continue;
			}
			chsMaterial[chs][flow].emplace_back(materials[i]);
		}
	}
}

std::string InitColumns::finder(std::map<std::string, double>& _currentFlowTarget)
{
	// before running  "finder", need update chsMaterial
	std::string result;
	double minValue = -1;

	for (auto& chs : chsMaterial)
	{
		double values = 0;
		for (auto& flowUnit : chs.second)
		{
			values += totalWeight(flowUnit.second) / (_currentFlowTarget[flowUnit.first] + 0.00001);
		}
		values /= static_cast<int>(_currentFlowTarget.size());
		if (values > minValue)
		{
			minValue = values;
			result = chs.first;
		}
	}
	return result;
}

void InitColumns::join(std::string& pool, std::vector<Material*>& selected)
{
	// get materials (no selected) from all flow unit in the pool
	// the sequence of coils in ¡°selected¡± is different from the materials
	// because the coils is split by flow and chs type
	for (const auto& item : chsMaterial[pool])
	{
		for (const auto& coil : item.second)
		{
			if (!coil->getStatus())
			{
				selected.emplace_back(coil);
			}
		}
	}
}

Solution* InitColumns::solver(std::string& _poolName, std::vector<Material*>& _poolMaterials)
{
	auto comparator = new Comparator(0, 1);
	auto criterion = new TerminationCriterion("MaxGeneration", 0, 200);
	auto problem = new Problem(_poolMaterials);
	int depth = std::min(static_cast<int>(mpProblem->getVertices() * 0.1), 5);

	auto algorithm = new Algorithm(problem, criterion, depth, comparator);
	algorithm->execute(currentFlowTarget);

	//std::cout << "InitColumn_CHS_name -> " << _poolName << std::endl;
	
	Solution* solution = algorithm->getSolutions(parameters->planLower);

	delete algorithm;
	delete problem;
	delete criterion;
	delete comparator;

	return solution;
}

std::vector<std::string> InitColumns::getConnectNames(std::string& _poolName)
{
	std::vector<std::string> result;

	int index = -1;
	for (const auto& item : parameters->chsType)
	{
		if (item.first == _poolName)
		{
			index = item.second;
			break;
		}
	}
	if (index == -1)
	{
		std::cout << "InitColumns : Error in get Connect pool materials" << std::endl;
	}

	for (int i = 0; i < parameters->C; i++)
	{
		int flag = parameters->steelArr[index][i];
		if (flag == 1 && i != index)
		{
			int _conIndex = i;
			std::string conName = parameters->chs_type[_conIndex];
			result.emplace_back(conName);
		}
	}
	return result;
}

void InitColumns::mergePool(std::vector<Material*>& _merge, std::vector<Material*>& _original,
	std::vector<Material*>& _substitute)
{
	_merge.insert(_merge.end(), _original.begin(), _original.end());
	_merge.insert(_merge.end(), _substitute.begin(), _substitute.end());
}

void InitColumns::getSolutions()
{
	currentFlowTarget = parameters->flow_mass;

	// this finished condition maybe not meet each target unit
	while (!isFinishedAll(currentFlowTarget))
	{
		// fill the pool based on current flow target and chs type and coil's status
		flowPool(currentFlowTarget);

		// choose pool which has most materials for the target flow unit
		std::string poolName = finder(currentFlowTarget);

		std::vector<Material*> selected;
		join(poolName, selected);

		double singleCapacity = 0.0;
		Solution* preferred = nullptr;
		if (totalWeight(selected) < parameters->planLower)
		{
			std::cout << "InitColumns -> " << poolName << " -> total weight < plan lower limit " << std::endl;
		}
		else
		{
			preferred = solver(poolName, selected);
			singleCapacity = preferred->get_objective(2);
		}
		if (singleCapacity < parameters->planLower)
		{
			std::vector<std::string> connectPoolNames = getConnectNames(poolName);

			std::vector<Solution*> solutionsCon;
			std::vector<std::string> poolCon;
			std::vector<std::vector<Material*>> coilSetsCon;

			// for all pools which can be connected
			for (auto& name : connectPoolNames)
			{
				std::vector<Material*> substitute;
				join(name, substitute);
				if (substitute.empty())
				{
					std::cout << "InitColumns -> get columns -> "
						<< poolName << " connect " << "name" << " -> still empty" << std::endl;
					continue;
				}
				std::vector<Material*> merge;
				mergePool(merge, selected, substitute);
				std::string newName = poolName;
				newName += "_";
				newName += name;
				Solution* preferredOther = solver(newName, merge);
				if (preferredOther != nullptr)
				{
					solutionsCon.emplace_back(preferredOther);
					poolCon.emplace_back(name);
					coilSetsCon.emplace_back(merge);
				}
			}
			// get a better solution from  solutionsCon
			double maxWeight = -1;
			int maxIndex = -1;
			int index = 0;
			for (const auto& item : solutionsCon)
			{
				double weight = item->get_objective(2);
				if (weight > maxWeight)
				{
					maxWeight = weight;
					maxIndex = index;
				}
				index++;
			}
			if (maxIndex != -1)
			{
				if (preferred != nullptr)
				{
					delete preferred;
					preferred = nullptr;
				}
				poolName += "_" + poolCon[maxIndex];
				selected = coilSetsCon[maxIndex];
				preferred = new Solution(*(solutionsCon[maxIndex]));
			}
			for (auto& item : solutionsCon)
			{
				delete item;
			}
		}

		// if single and connect pool can not generate solution break;
		if (preferred == nullptr)
		{
			break;
		}
		int num = 0;
		double _total = 0.0;
		for (int i = 0; i < preferred->get_number_of_variables() - preferred->get_number_of_violated_constraints(); i++)
		{
			_total += selected[preferred->get_variable(i)]->getValue(Attribute::IN_MAT_WT);
			// _total += materials[preferred->get_variable(i)]->getValue(Attribute::IN_MAT_WT);
			if (_total > parameters->planUp)
			{
				break;
			}
			num++;
		}
		if (_total < parameters->planLower)
		{
			for (int i = 0; i < num; i++)
			{
				int loc = preferred->get_variable(i);
				selected[loc]->setStatus();
			}
			continue;
		}
		auto column = new Column(mpProblem);
		column->setName(poolName);
		for (int i = 0; i < num; i++)
		{
			int loc = preferred->get_variable(i);
			int index = selected[loc]->index;	// notice: the index of coil
			selected[loc]->setStatus();
			column->addVertex(index);
			currentFlowTarget[selected[loc]->getStr(Attribute::FLOW)] -= selected[loc]->getValue(Attribute::IN_MAT_WT);
		}
		column->calculateCost();
		initPool.emplace_back(column);

		delete preferred;
		auto iter = currentFlowTarget.begin();
		while (iter != currentFlowTarget.end())
		{
			if (iter->second <= 0)
			{
				currentFlowTarget.erase(iter++);
			}
			else
			{
				iter++;
			}
		}
	}
}

bool InitColumns::isFeasible()
{
	double massC008 = 0.0;
	double massC512 = 0.0;
	bool isFeasible = false;
	for (const auto& item : initPool)
	{
		for (int i = 1; i < mpProblem->getVertices() + 1; i++)
		{
			if (item->contains(i))
			{
				if (materials[i - 1]->getStr(Attribute::FLOW) == "C008")
				{
					massC008 += materials[i - 1]->getValue(Attribute::IN_MAT_WT);
				}
				else if (materials[i - 1]->getStr(Attribute::FLOW) == "C512")
				{
					massC512 += materials[i - 1]->getValue(Attribute::IN_MAT_WT);
				}
			}
		}
	}
	for (const auto& item : initPool)
	{
		double mass = 0.0;
		for (int i = 1; i < mpProblem->getVertices() + 1; i++)
		{
			if (item->contains(i))
			{
				mass += materials[i - 1]->getValue(Attribute::IN_MAT_WT);
			}
		}
		if (mass < parameters->planLower)
		{
			std::cout << item->getName() << " : mass < plan's lower limit" << std::endl;
		}
	}
	std::cout << "mass_C512: " << massC512 << " mass_C008: " << massC008 << std::endl;

	// maybe neither unit meets the demand
	if (massC008 < parameters->flow_mass["C008"])
	{
		std::cout << "C008 mass is no meet" << std::endl;
		isFeasible = false;
		return isFeasible;
	}
	if (massC512 < parameters->flow_mass["C512"])
	{
		std::cout << "C512 mass is no meet" << std::endl;
		isFeasible = false;
		return isFeasible;
	}
	isFeasible = true;
	return isFeasible;
}

void InitColumns::getInitCost()
{
	std::cout << "print the initial solutions" << std::endl;
	for (const auto& item : initPool)
	{
		item->printRoute();
		std::cout << "_________________________" << std::endl;
	}

	bool flag = isFeasible();
	std::cout << "this initial column pool is feasible : " << flag << std::endl;

	// get objective of initial columns in math model
	double initCost = 0.0;
	for (const auto& item : initPool)
	{
		initCost += item->getCost();
	}
	initCost = initCost * parameters->alpha1;
	initCost += int(initPool.size()) * parameters->setUpCost * parameters->alpha2;
	std::cout << "init cost -> " << initCost << std::endl;
}

