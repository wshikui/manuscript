#pragma once
#include "Problem.h"
#include <optional>

class Column
{
public:
	Column() = default;

	~Column() = default;

	explicit Column(Problem* _problem);
	explicit Column(Column* _column);

	Column(Problem* _problem, std::vector<int>& _route);

	[[nodiscard]] std::vector<int>::const_iterator begin() const;
	[[nodiscard]] std::vector<int>::const_iterator end() const;

	bool contains(int i) const;

	int contain_num(int i) const;

	int getNbVertices() const;

	double getCost() const;
	double getDemand() const;

	bool addVertex(int i);

	void setCost(double cost);

	void printRoute() const;

	void setName(std::string _name);
	std::string getName() const;

	void calculateCost();
	void setRC(double _reduceCost);
	[[nodiscard]] double getRC() const;
	std::vector<int> getRoute() const;

	//maybe not use
	int numRealCoils = 0;

	bool column_has_cycle() const;
	std::map<int, std::vector<int>> visited_coils_with_predecessors() const;
	std::optional<std::pair<int, int>> common_coil_visited_from_two_different_predecessors(const Column& other_col) const;

private:
	double mCost{};
	double mDemand{};
	double reduceCost{};
	bool has_cycle = false;

	std::vector<int> mContained;
	std::vector<int> mRoute;

	Problem* problem{};
	std::string name;
};

