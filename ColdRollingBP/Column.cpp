//
// Created by wshikui on 2023/9/20.
//

#include <utility>

#include "Column.h"

Column::Column(Problem* _problem)
{
	mCost = 0;
	mDemand = 0;
	problem = _problem;
	mContained = std::vector<int>(_problem->getVertices() + 1, 0);

	// include the virtual vertex and coils vertex
	// actually, the vertices in column don't include virtual vertex
	mRoute.reserve(_problem->getVertices() + 1);

	numRealCoils = _problem->getVertices();
}

Column::Column(Column* _column)
{
	if (_column != nullptr)
	{
		this->name = _column->name;
		this->mRoute = _column->mRoute;
		this->mCost = _column->mCost;
		this->mContained = _column->mContained;
		this->mDemand = _column->mDemand;
		this->problem = _column->problem;
		numRealCoils = _column->numRealCoils;
	}
}

Column::Column(Problem* _problem, std::vector<int>& _route)
{
	problem = _problem;
	mContained = std::vector<int>(_problem->getVertices() + 1, 0);

	for (const auto& item : _route) {
		mRoute.push_back(item);
		mContained[item] = 1;
	}
	double cost = 0;
	for (int i = 0; i < _route.size() - 1; i++) {
		cost += problem->getCost(mRoute[i] - 1, mRoute[i + 1] - 1);
	}
	mCost = cost;
	numRealCoils = _problem->getVertices();
}

std::vector<int>::const_iterator Column::begin() const
{
	return std::cbegin(mRoute);
}

std::vector<int>::const_iterator Column::end() const
{
	return std::cend(mRoute);
}

bool Column::contains(const int i) const
{
	return mContained[i];
}

int Column::contain_num(int i) const
{
	int num = std::count(mRoute.begin(), mRoute.end(), i);

	return num;
}

int Column::getNbVertices() const
{
	return static_cast<int>(mRoute.size());
}

double Column::getCost() const
{
	return mCost;
}

double Column::getDemand() const
{
	return mDemand;
}

bool Column::addVertex(const int i)
{
	if (contains(i) == 1)
	{
		//std::cout << "add duplicate coils to the column" << std::endl;
		has_cycle = true;
		//return false;
	}

	// after a vertex is added in this route, calculate the current demand
	// after all vertices are added in this route, calculate the cost of this route
	if (i == 0)
	{
		mDemand += 0;
	}
	else
	{
		mDemand += problem->getWeight(i);
	}
	mRoute.push_back(i);
	return mContained[i] = 1;
}

void Column::setCost(const double cost)
{
	mCost = cost;
}

void Column::printRoute() const
{
	std::cout << "__________route info : " << name << "__________" << std::endl;

	int width = 8;

	std::cout << std::setw(width) << std::left << "index" << "|"
		<< std::setw(width) << std::left << "weight" << "|"
		<< std::setw(width) << std::left << "outWidth" << "|"
		<< std::setw(width) << std::left << "flow unit" << "|"
		<< std::endl;

	for (const auto& item : mRoute)
	{
		std::cout << std::setw(width) << std::left << item << "|"
			<< std::setw(width) << std::left << problem->getWeight(item) << "|"
			<< std::setw(width) << std::left << problem->getOutWidth(item) << "|"
			<< std::setw(width) << std::left << problem->getFlow(item) << "|"
			<< std::endl;
	}
	std::cout << "cost -> " << mCost << "  total weight -> " << mDemand << " reduced cost -> " << reduceCost 
		<< " cycle -> " << has_cycle << std::endl;
	//std::cout << "cost -> " << mCost << "  total weight -> " << mDemand << std::endl;
}

void Column::setName(std::string _name)
{
	this->name = std::move(_name);
}

std::string Column::getName() const
{
	return name;
}

void Column::calculateCost()
{
	double _cost = 0.0;
	for (int i = 0; i < mRoute.size() - 1; i++)
	{
		_cost += problem->getCost(mRoute[i] - 1, mRoute[i + 1] - 1);
	}
	mCost = _cost;
}

void Column::setRC(double _reduceCost)
{
	reduceCost = _reduceCost;
}

double Column::getRC() const
{
	return reduceCost;
}

std::vector<int> Column::getRoute() const
{
	return mRoute;
}

bool Column::column_has_cycle() const
{
	return has_cycle;
}

std::map<int, std::vector<int>> Column::visited_coils_with_predecessors() const
{
	std::map<int, std::vector<int>> visited;

	//the predecessor of the first real coil is 0.
	visited[mRoute[0]] = { 0 };

	for (size_t i = 1; i < mRoute.size(); i++) {
		//key is the item in mRoute.
		int current_coil = mRoute[i];
		int predecessor_coil = mRoute[i - 1];
		visited[current_coil].push_back(predecessor_coil);
	}
	return visited;
}

std::optional<std::pair<int, int>> Column::common_coil_visited_from_two_different_predecessors(const Column& other_col) const
{
	std::map<int, std::vector<int>> this_coils_col = this->visited_coils_with_predecessors();
	std::map<int, std::vector<int>> other_coils_col = other_col.visited_coils_with_predecessors();

	for (const auto& kv : this_coils_col) {
		if (other_coils_col.find(kv.first) == other_coils_col.end()) { continue; }

		const int& successor = kv.first;
		const std::vector<int>& predecessors = kv.second;

		const std::vector<int>& other_predecessors = other_coils_col[successor];

		for (const auto& item : predecessors) {
			if (std::find(other_predecessors.begin(), other_predecessors.end(), item) == other_predecessors.end()) {
				return std::make_pair(item, successor);
			}
		}
	}

	return std::optional<std::pair<int, int>>();
}





