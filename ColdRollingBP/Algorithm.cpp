//
// Created by wshikui on 2023/9/21.
//
#include "Algorithm.h"

Comparator::Comparator(int _index1, int _index2)
{
	index1 = _index1;
	index2 = _index2;
}

// one > two : 1
int Comparator::compare(void* one, void* two) const
{
	if (one == nullptr)
	{
		std::cout << "when compare two solutions -> one solution is null" << std::endl;
		return -1;
	}
	else if (two == nullptr)
	{
		std::cout << "when compare two solutions -> two solution is null" << std::endl;
		return 1;
	}
	std::vector<int> indexObjective = { index1, index2 };

	// multi compare
	for (const auto& item : indexObjective)
	{
		double objective1 = ((Solution*)one)->get_objective(item);
		double objective2 = ((Solution*)two)->get_objective(item);
		if (objective1 == objective2)
		{
			continue;
		}
		else
		{
			if (objective1 < objective2)
			{
				return -1;
			}
			else
			{
				return 1;
			}
		}
	}
	// std::cout << "when compare two solution : solution_1 == solution_2" << std::endl;
	return 0;
}

Archive::Archive(int size, Comparator* _comparator)
{
	capacity = size;
	comparator = _comparator;
	this->solution_list = new std::vector<Solution*>();
}

Archive::~Archive()
{
	//clear items in vector

	for (auto& item : *solution_list)
	{
		delete item;
	}
	(*solution_list).clear();
	delete solution_list;
}

std::vector<Solution*>* Archive::get_solution_list() const
{
	return solution_list;
}

int Archive::index_worst()
{
	// -1 : no find worst solution
	if (solution_list->empty())
	{
		return -1;
	}

	int index = 0;
	int flag;
	Solution* worstSolution = (*solution_list)[0];
	Solution* candidateSolution;

	for (size_t i = 1; i < solution_list->size(); i++)
	{
		candidateSolution = (*solution_list)[i];
		flag = comparator->compare(worstSolution, candidateSolution);
		// candidateSolution is more worse with more smaller objection
		if (flag == 1)
		{
			index = (int)i;
			worstSolution = candidateSolution;
		}
	}
	return index;
}

bool Archive::add(Solution* solution)
{
	if (solution == nullptr) {
		return false;
	}
	// it was in solution list
	for (auto& item : *solution_list)
	{
		if (*item == *solution)
		{
			return false;
		}
	}

	if ((int)(*solution_list).size() < capacity)
	{
		solution_list->emplace_back(solution);
		return true;
	}
	else
	{
		int worst = index_worst();

		if (worst == 0) {
			//std::cout << "get worst solution: index 0" << std::endl;
		}

		int cmp = comparator->compare(solution, (*solution_list)[worst]);

		if (cmp == 1)
		{
			//maybe == 1, objective of one is more large.
			delete (*solution_list)[worst];
			(*solution_list).erase((*solution_list).begin() + worst);
			(*solution_list).push_back(solution);
			return true;
		}
	}
	// add failed
	return false;
}


int Archive::size() const
{
	return (int)solution_list->size();
}

Solution* Archive::get(int index)
{
	return (*solution_list)[index];
}

void Archive::sort()
{
	if (!solution_list->empty())
	{
		for (size_t i = 0; i < solution_list->size() - 1; i++)
		{
			for (size_t j = i + 1; j < solution_list->size(); j++)
			{
				if (comparator->compare(solution_list->at(i), solution_list->at(j)) == -1)
				{
					// solution i is batter, exchange i and j
					Tool::swap(solution_list->at(i), solution_list->at(j));
				}
			}
		}
	}
	else
	{
		std::cout << "Now the solution list is empty" << std::endl;
	}
}


Neighbor::Neighbor(int _depth)
{
	depth = _depth;
}

void* Neighbor::execute(void* object)
{
	auto* solution = (Solution*)object;

	int nVariable = solution->get_number_of_variables();

	int loc = nVariable - solution->get_number_of_violated_constraints();

	auto* solutionSet = new std::vector<Solution*>();

	if (loc == nVariable)
	{
		return solutionSet;
	}

	// Step2 插入
	for (int i = 0; i < depth; ++i)
	{
		// Step2.1 向前
		int pre = loc - (i + 1);
		if (pre > -1)
		{
			// 删除
			// the new solution will be deleted in ~Archive
			auto* ss1 = new Solution(*solution);
			for (int j = pre; j < nVariable - (i + 1); ++j)
			{
				ss1->set_variable(j, solution->get_variable(j + (i + 1)));
			}

			for (int j = 0; j < (i + 1); ++j)
			{
				ss1->set_variable(nVariable - (i + 1) + j, solution->get_variable(loc - (i + 1) + j));
			}

			solutionSet->push_back(ss1);

			// 插入
			auto* ss2 = new Solution(*solution);
			for (int j = pre; j < loc; ++j)
			{
				ss2->set_variable(j + 1, solution->get_variable(j));
			}

			ss2->set_variable(pre, solution->get_variable(loc));

			solutionSet->push_back(ss2);

		}

		int pos = loc + (i + 1);
		if (pos < nVariable)
		{
			auto* ss3 = new Solution(*solution);
			for (int j = pos; j < nVariable; ++j)
			{
				ss3->set_variable(j - (i + 1), solution->get_variable(j));
			}

			for (int j = 0; j < (i + 1); ++j)
			{
				ss3->set_variable(nVariable - (i + 1) + j, solution->get_variable(loc + j));
			}

			solutionSet->push_back(ss3);

			auto* ss4 = new Solution(*solution);
			for (int j = loc; j < pos; ++j)
			{
				ss4->set_variable(j, solution->get_variable(j + 1));
			}
			ss4->set_variable(pos, solution->get_variable(loc));
			solutionSet->push_back(ss4);
		}
	}

	return solutionSet;
}


TerminationCriterion::TerminationCriterion(const std::string& _condition, long _state, long _upper_limit)
{
	condition = _condition;
	state = _state;
	upper_limit = _upper_limit;
}

void TerminationCriterion::update()
{
	state += 1;
}

bool TerminationCriterion::fulfilled() const
{
	return state >= upper_limit;
}

const std::string& TerminationCriterion::get_condition() const
{
	return condition;
}

long TerminationCriterion::get_state() const
{
	return state;
}

long TerminationCriterion::get_upper_limit() const
{
	return upper_limit;
}


Algorithm::Algorithm(Problem* _problem, TerminationCriterion* _criterion, int depth, Comparator* _comparator)
{
	problem = _problem;
	criterion = _criterion;
	archive = new Archive(50, _comparator);
	neighbor = new Neighbor(depth);
}

Algorithm::~Algorithm()
{
	// use ~archive to delete the solutions
	delete archive;
	delete neighbor;
}

Solution* Algorithm::find_x0()
{
	// if x0 == nullptr, all solution in archive have been searched
	Solution* x0 = nullptr;
	int* perm = Permutation::perm(archive->size());
	for (int i = 0; i < archive->size(); ++i)
	{
		int index = perm[i];
		if (archive->get(index)->get_label() == "no")
		{
			//debug
			std::cout << "archive size = " << archive->size() <<
				", find index = " << index << std::endl;

			x0 = archive->get(index);
			break;
		}
	}

	delete[] perm;
	return x0;
}

Solution* Algorithm::find_x0_fix()
{
	Solution* x0 = nullptr;
	for (int i = 0; i < archive->size(); ++i) {
		if (archive->get(i)->get_label() == "no") {
			x0 = archive->get(i);
			break;
		}
	}

	return x0;
}

std::vector<Solution*>* Algorithm::execute(const std::map<std::string, double>& _currentFlowTarget)
{
	Solution* solution = problem->createInitialSolution();
	solution->set_label("no");
	problem->evaluate(solution, _currentFlowTarget);
	criterion->update();
	archive->add(solution);

	while (!criterion->fulfilled() && solution->get_number_of_violated_constraints() != 0)
	{
		Solution* x0 = find_x0_fix();
		if (x0 == nullptr)
		{
			break;
		}
		auto* a = new Solution(*x0);
		x0->set_label("yes");
		auto* N = static_cast<std::vector<Solution*>*> (neighbor->execute(a));

		for (size_t i = 0; i < N->size(); ++i)
		{
			problem->evaluate(N->at(i), _currentFlowTarget);
			criterion->update();
			bool flag = archive->add(N->at(i));
			if (!flag)
			{
				delete N->at(i);
			}
		}

		N->clear();
		delete N;
		delete a;

	}
	archive->sort();
	return archive->get_solution_list();
}

Solution* Algorithm::getBetterSolution(double lower, double upper)
{
	// must return a pointer
	Solution* preferred;
	double minObj = FLT_MAX;
	double maxObj = -FLT_MAX;
	int minIndex = 0;
	int maxIndex = 0;

	std::vector<Solution*>* solutions = archive->get_solution_list();

	std::cout << "solutions size -> " << solutions->size() << std::endl;
	for (int i = 0; i < solutions->size(); i++)
	{
		double obj = solutions->at(i)->get_objective(2);
		if ((obj > maxObj) && (obj > lower) && (obj < upper))
		{
			maxObj = obj;
			maxIndex = i;
		}
		if (obj < minObj)
		{
			minObj = obj;
			minIndex = i;
		}
	}
	if (maxObj < 0)
	{
		preferred = new Solution(*(*solutions)[minIndex]);
	}
	else
	{
		preferred = new Solution(*(*solutions)[maxIndex]);
	}
	return preferred;
}

Solution* Algorithm::getSolutions(double lower)
{
	//get a solution with heaviest weight
	std::vector<Solution*>* solutions = archive->get_solution_list();
	int index = 0;
	int weight = solutions->at(0)->get_objective(2);
	Solution* preferred{};

	for (int i = 1; i < solutions->size(); i++) {
		double obj = solutions->at(i)->get_objective(2);
		if ((obj >= lower) && obj > weight) {
			index = i;
		}
	}
	preferred = new Solution(*(*solutions)[index]);
	return preferred;
}



