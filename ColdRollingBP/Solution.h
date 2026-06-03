#pragma once
#include <ostream>
#include <string>

class Solution
{
public:
	Solution(int number_of_variables, int number_of_objectives);

	Solution(const Solution& solution);

	virtual ~Solution();

	virtual Solution* copy();

	int get_number_of_variables() const;

	int get_number_of_objectives() const;

	int get_number_of_violated_constraints() const;

	void set_number_of_violated_constraints(int _number_of_violated_constraints);

	double get_overall_constraint_violation() const;

	void set_overall_constraint_violation(double _overall_constraint_violation);

	int* get_variables() const;

	int get_variable(int index) const;

	double* get_objectives() const;

	double get_objective(int index) const;

	void set_variable(int index, int value);

	void set_objective(int index, double value);

	const std::string& get_label() const;

	void set_label(const std::string& _label);


	bool operator==(const Solution& solution)
	{
		bool flag = true;

		for (int i = 0; i < this->get_number_of_variables(); ++i)
		{
			if (this->variables[i] != solution.get_variable(i))
			{
				flag = false;
				break;
			}
		}

		return flag;
	}


private:
	int number_of_variables;

	int number_of_objectives;

	int number_of_violated_constraints;

	double overall_constraint_violation;

	//整数决策变量数组
	int* variables;

	//目标函数数组
	double* objectives;

	std::string label;
};
