#pragma once
#include "Problem.h"
//#include <algorithm>
#include <random>

// generate random array
class Permutation
{
public:
	static int* perm(int length)
	{
		int* result = new int[length];

		int* aux = new int[length];

		// 使用随机数生成器引擎
		std::random_device rd;
		std::mt19937 gen(rd());

		// First, create an array from 0 to length - 1.
		// Also is needed to create a random array of size length
		for (int i = 0; i < length; i++)
		{
			result[i] = i;
			std::uniform_int_distribution<> dis(0, length - 1);
			aux[i] = dis(gen);
		}

		// Sort the random array with effect in result, and then we obtain a
		// permutation array between 0 and length - 1
		std::sort(result, result + length, [&](int a, int b)
			{
				return aux[a] < aux[b];
			});

		delete[] aux;
		return result;
	}
};

// for compare two solution by maximizing objective, include __two__ objective function values
class Comparator
{
public:
	Comparator() = default;

	~Comparator() = default;

	Comparator(int _index1, int _index2);

	int compare(void* one, void* two) const;

private:
	// index of objective function
	int index1;
	int index2;
};


class Archive
{
public:
	Archive(int size, Comparator* _comparator);

	// delete solutions in archive, the solutions from neighbor
	virtual ~Archive();

	[[nodiscard]] std::vector<Solution*>* get_solution_list() const;

	bool add(Solution* solution);

	int size() const;

	Solution* get(int index);

	void sort();

private:

	// comparator for add better solution to solution list
	Comparator* comparator;

	// vector pointer
	std::vector<Solution*>* solution_list;

	int capacity;

	int index_worst();
};

class Neighbor
{
public:
	explicit Neighbor(int _depth);

	virtual void* execute(void* object);

private:
	int depth;
};

class TerminationCriterion
{
public:
	TerminationCriterion(const std::string& _condition, long _state, long _upper_limit);

	TerminationCriterion() = default;

	void update();

	bool fulfilled() const;

	const std::string& get_condition() const;

	long get_state() const;

	long get_upper_limit() const;

protected:
	std::string condition;

	long state{};

	long upper_limit{};
};

class Algorithm
{
public:
	Algorithm(Problem* _problem, TerminationCriterion* _criterion, int depth, Comparator* _comparator);

	virtual ~Algorithm();

	virtual std::vector<Solution*>* execute(const std::map<std::string, double>& _currentFlowTarget);

	// get a solution from archive
	Solution* getBetterSolution(double lower, double upper);

	//get solutions from archive
	Solution* getSolutions(double lower);

private:
	Problem* problem;

	TerminationCriterion* criterion;

	Archive* archive;

	Neighbor* neighbor;

	Solution* find_x0();

	Solution* find_x0_fix();
};


