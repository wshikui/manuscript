//
// Created by wshikui on 2023/9/21.
//
#include "Solution.h"

Solution::Solution(int number_of_variables, int number_of_objectives)
{
	this->number_of_variables = number_of_variables;
	this->number_of_objectives = number_of_objectives;
	this->number_of_violated_constraints = 0;
	this->overall_constraint_violation = 0.;

	variables = new int[number_of_variables];
	objectives = new double[number_of_objectives];
	label = "";
}

Solution::Solution(const Solution& solution)
{
	number_of_variables = solution.get_number_of_variables();
	number_of_objectives = solution.get_number_of_objectives();
	number_of_violated_constraints = solution.get_number_of_violated_constraints();
	overall_constraint_violation = solution.get_overall_constraint_violation();

	variables = new int[number_of_variables];
	objectives = new double[number_of_objectives];

	for (int i = 0; i < number_of_variables; ++i)
	{
		variables[i] = solution.get_variable(i);
	}

	for (int i = 0; i < number_of_objectives; ++i)
	{
		objectives[i] = solution.get_objective(i);
	}

	label = solution.get_label();
}

Solution::~Solution()
{
	//std::cout << variables.
	delete[] variables;
	delete[] objectives;
}

Solution* Solution::copy()
{
	return new Solution(*this);
}

int Solution::get_number_of_variables() const
{
	return number_of_variables;
}

int Solution::get_number_of_objectives() const
{
	return number_of_objectives;
}

int* Solution::get_variables() const
{
	return variables;
}

int Solution::get_variable(int index) const
{
	return variables[index];
}

double* Solution::get_objectives() const
{
	return objectives;
}

double Solution::get_objective(int index) const
{
	return objectives[index];
}

void Solution::set_variable(int index, int value)
{
	variables[index] = value;
}

void Solution::set_objective(int index, double value)
{
	objectives[index] = value;
}

int Solution::get_number_of_violated_constraints() const
{
	return number_of_violated_constraints;
}

void Solution::set_number_of_violated_constraints(int _number_of_violated_constraints)
{
	Solution::number_of_violated_constraints = _number_of_violated_constraints;
}

double Solution::get_overall_constraint_violation() const
{
	return overall_constraint_violation;
}

void Solution::set_overall_constraint_violation(double _overall_constraint_violation)
{
	Solution::overall_constraint_violation = _overall_constraint_violation;
}

const std::string& Solution::get_label() const
{
	return label;
}

void Solution::set_label(const std::string& _label)
{
	Solution::label = _label;
}


