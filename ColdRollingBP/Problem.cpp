//
// Created by wshikui on 2023/9/15.
//
#include "Problem.h"

Problem::Problem(const std::vector<Material*>& materials)
{
	_materials = materials;

	// the number of coils, exclude the virtual vertex
	vertices = static_cast<int>(materials.size());

	materialsFlow.resize(vertices);
	initialization();

	// create the one dimensional matrix
	in_weight = new double[vertices];
	in_width = new double[vertices];
	out_width = new double[vertices];
	in_thick = new double[vertices];
	out_thick = new double[vertices];

	// create the two-dimensional matrix
	in_width_matrix = new double* [vertices];
	out_width_matrix = new double* [vertices];
	in_thick_matrix = new double* [vertices];
	out_thick_matrix = new double* [vertices];
	steel_type_matrix = new bool* [vertices];
	penalty_matrix = new double* [vertices];
	feasible_matrix = new bool* [vertices];

	for (int i = 0; i < vertices; i++)
	{
		in_weight[i] = materials[i]->getValue(Attribute::IN_MAT_WT);
		in_width[i] = materials[i]->getValue(Attribute::IN_MAT_WIDTH);
		out_width[i] = materials[i]->getValue(Attribute::OUT_MAT_WIDTH);
		in_thick[i] = materials[i]->getValue(Attribute::IN_MAT_THICK);
		out_thick[i] = materials[i]->getValue(Attribute::OUT_MAT_THICK);

		materialsFlow[i] = materials[i]->getStr(Attribute::FLOW);

		// create the two-dimensional matrix
		in_width_matrix[i] = new double[vertices];
		out_width_matrix[i] = new double[vertices];
		in_thick_matrix[i] = new double[vertices];
		out_thick_matrix[i] = new double[vertices];
		steel_type_matrix[i] = new bool[vertices];
		penalty_matrix[i] = new double[vertices];
		feasible_matrix[i] = new bool[vertices];

		// get the jump values between two coils from all materials
		for (int j = 0; j < vertices; j++)
		{
			in_width_matrix[i][j] = in_width_jump(i, j);
			out_width_matrix[i][j] = out_width_jump(i, j);
			in_thick_matrix[i][j] = in_thick_jump(i, j);
			out_thick_matrix[i][j] = out_thick_jump(i, j);
			std::string steel_prev = materials[i]->getStr(Attribute::RULE_CODE);
			std::string steel_post = materials[j]->getStr(Attribute::RULE_CODE);
			steel_type_matrix[i][j] = steel_type_jump(steel_prev, steel_post);
			penalty_matrix[i][j] = get_penalty(i, j);
			feasible_matrix[i][j] = feasible(i, j);
		}
	}

	//after initialize the cost(penalty) matrix, we get the forbid vertex of every vertex
	getForbidVertex();
}

int Problem::getVertices() const
{
	return vertices;
}

double Problem::getCost(int prev, int post) const
{
	return penalty_matrix[prev][post];
}

Problem::~Problem()
{
	delete[] in_weight;

	delete[] in_width;
	delete[] out_width;
	delete[] in_thick;
	delete[] out_thick;

	for (int i = 0; i < vertices; i++)
	{
		delete[] in_width_matrix[i];
		delete[] out_width_matrix[i];
		delete[] in_thick_matrix[i];
		delete[] out_thick_matrix[i];
		delete[] steel_type_matrix[i];
		delete[] penalty_matrix[i];
	}

	delete[] in_width_matrix;
	delete[] out_width_matrix;
	delete[] in_thick_matrix;
	delete[] out_thick_matrix;
	delete[] steel_type_matrix;
	delete[] penalty_matrix;
}

double Problem::in_width_jump(int prev, int post)
{
	double temp = prev == post ? M : _materials[prev]->getValue(Attribute::IN_MAT_WIDTH) -
		_materials[post]->getValue(Attribute::IN_MAT_WIDTH);
	return temp;
}

double Problem::out_width_jump(int prev, int post)
{
	double temp = prev == post ? M : _materials[prev]->getValue(Attribute::OUT_MAT_WIDTH) -
		_materials[post]->getValue(Attribute::OUT_MAT_WIDTH);
	return temp;
}

double Problem::in_thick_jump(int prev, int post)
{
	double temp = prev == post ? M : _materials[prev]->getValue(Attribute::IN_MAT_THICK) -
		_materials[post]->getValue(Attribute::IN_MAT_THICK);
	return temp;
}

double Problem::out_thick_jump(int prev, int post)
{
	double temp = prev == post ? M : _materials[prev]->getValue(Attribute::OUT_MAT_THICK) -
		_materials[post]->getValue(Attribute::OUT_MAT_THICK);
	return temp;
}

bool Problem::steel_type_jump(const std::string& s_prev, const std::string& s_post)
{
	std::vector<std::string> prev_rules = Tool::super_split(s_prev, ",");
	std::vector<std::string> post_rules = Tool::super_split(s_post, ",");

	for (size_t i = 0; i < prev_rules.size(); ++i)
	{
		for (size_t j = 0; j < post_rules.size(); ++j)
		{
			std::vector<std::string> prev_i = Tool::super_split(prev_rules.at(i), "-");
			std::vector<std::string> post_j = Tool::super_split(post_rules.at(j), "-");

			if (prev_i[0] == post_j[0])
			{
				return true;
			}
			if (prev_i.size() > 1 && post_j.size() > 1)
			{
				if (steel_rules[prev_rules.at(i)] == post_rules.at(j) ||
					prev_rules.at(i) == steel_rules[post_rules.at(j)])
				{
					return true;
				}
			}
		}
	}
	return false;
}

void Problem::initialization()
{
	steel_rules["C502GZ1-2"] = "C502GZ2-1";
	steel_rules["C502GZ2-1"] = "C502GZ1-2";
	steel_rules["C502GZ2-3"] = "C502GZ3-2";
	steel_rules["C502GZ3-2"] = "C502GZ2-3";
	steel_rules["C502GZ3-4"] = "C502GZ4-3";
	steel_rules["C502GZ4-3"] = "C502GZ3-4";

	up_width_limit = 150;
	down_width_limit = 150;
	up_thick_limit = 99;
	down_thick_limit = 99;

	parameter = Parameters();
}

double Problem::get_penalty(int prev, int post)
{
	// TODO set the >= M to M, narrow the field
	double c_ij = 0;
	// penalty of IN WIDTH
	if (in_width_matrix[prev][post] >= 0 && in_width_matrix[prev][post] < down_width_limit)
	{
		c_ij += 0.2 * in_width_matrix[prev][post];
	}
	else if (in_width_matrix[prev][post] > -up_width_limit && in_width_matrix[prev][post] <= 0)
	{
		c_ij += 0.8 * abs(in_width_matrix[prev][post]);
	}
	else
	{
		c_ij += M;
	}

	// penalty of OUT WIDTH
	if (out_width_matrix[prev][post] >= 0 && out_width_matrix[prev][post] < down_width_limit)
	{
		c_ij += 0.3 * out_width_matrix[prev][post];
	}
	else if (out_width_matrix[prev][post] >= -up_width_limit && out_width_matrix[prev][post] < 0)
	{
		c_ij += 0.3 * abs(out_width_matrix[prev][post]);
	}
	else
	{
		c_ij += M;
	}

	// penalty of IN THICK
	if (in_thick_matrix[prev][post] >= 0 && in_thick_matrix[prev][post] < down_thick_limit)
	{
		c_ij += 0.3 * in_thick_matrix[prev][post];
	}
	else if (in_thick_matrix[prev][post] >= -up_thick_limit && in_thick_matrix[prev][post] < 0)
	{
		c_ij += 0.3 * abs(in_thick_matrix[prev][post]);
	}
	else
	{
		c_ij += M;
	}
	// penalty of OUT THICK
	if (out_thick_matrix[prev][post] >= 0 && out_thick_matrix[prev][post] < down_thick_limit)
	{
		c_ij += 0.3 * out_thick_matrix[prev][post];
	}
	else if (out_thick_matrix[prev][post] > -up_thick_limit && out_thick_matrix[prev][post] <= 0)
	{
		c_ij += 0.3 * abs(out_thick_matrix[prev][post]);
	}
	else
	{
		c_ij += M;
	}

	// penalty of STEEL GRADE
	if (steel_type_matrix[prev][post])
	{
		c_ij += 0;
	}
	else
	{
		c_ij += M;
	}

	//	if(c_ij > M)
	//	{
	//		c_ij = M;
	//	}

	return c_ij;
}

void Problem::getForbidVertex()
{
	//for every real vertex
	for (int i = 0; i < vertices; i++) {
		forbidVertex[i + 1];
		//get the connected coils with infinity penalty
		for (int j = 0; j < vertices; j++) {
			if (i != j) {
				if (penalty_matrix[i][j] >= M) {
					forbidVertex[i + 1].push_back(j + 1);
				}
			}
		}
	}

}

void Problem::printCostMatrix() const
{
	std::cout << "__________real coils cost matrix__________" << std::endl;

	int width = 8;
	for (int i = 0; i < vertices; i++)
	{
		for (int j = 0; j < vertices; j++)
		{
			std::cout << std::setw(width) << std::left << penalty_matrix[i][j] << "|";
		}
		std::cout << std::endl;
	}
}

void Problem::printCostMatrixPool(const std::string& name)
{
	// bool istest = steel_type_jump("C502GZ3,C502GZ3-2,C502GZ3-4", "C502GZ4,C502GZ4-3");

	std::cout << "______________print cost matrix of " << name << "______________" << std::endl;
	std::vector<int> indexPool_16(parameter.S), indexPool_17(parameter.S);
	std::vector<int> indexPool_18(parameter.S), indexPool_19(parameter.S);
	std::iota(indexPool_16.begin(), indexPool_16.end(), parameter.S * 0);
	std::iota(indexPool_17.begin(), indexPool_17.end(), parameter.S * 1);
	std::iota(indexPool_18.begin(), indexPool_18.end(), parameter.S * 2);
	std::iota(indexPool_19.begin(), indexPool_19.end(), parameter.S * 3);
	std::vector<int> index;
	if (name == "C502016") {
		index.insert(index.end(), indexPool_16.begin(), indexPool_16.end());
	}
	if (name == "C502017") {
		index.insert(index.end(), indexPool_17.begin(), indexPool_17.end());
	}
	if (name == "C502018") {
		index.insert(index.end(), indexPool_18.begin(), indexPool_18.end());
	}
	if (name == "C502019") {
		index.insert(index.end(), indexPool_19.begin(), indexPool_19.end());
	}
	if (name == "C502016_17") {
		index.insert(index.end(), indexPool_16.begin(), indexPool_16.end());
		index.insert(index.end(), indexPool_17.begin(), indexPool_17.end());
	}
	if (name == "C502017_18") {
		index.insert(index.end(), indexPool_17.begin(), indexPool_17.end());
		index.insert(index.end(), indexPool_18.begin(), indexPool_18.end());
	}
	if (name == "C502016_18") {
		index.insert(index.end(), indexPool_16.begin(), indexPool_16.end());
		index.insert(index.end(), indexPool_18.begin(), indexPool_18.end());
	}
	if (name == "C502018_19") {
		index.insert(index.end(), indexPool_18.begin(), indexPool_18.end());
		index.insert(index.end(), indexPool_19.begin(), indexPool_19.end());
	}
	int width = 8;
	int num = 0;
	for (int i = 0; i < index.size(); i++) {
		for (int j = 0; j < index.size(); j++) {
			std::cout << std::setw(width) << std::left << penalty_matrix[index[i]][index[j]] << "|";
			if (penalty_matrix[index[i]][index[j]] >= M) {
				num++;
				feasible_info(index[i], index[j]);
			}
		}
		std::cout << std::endl;
	}
	std::cout << "penelty coefficient infinity number ： " << num << std::endl;
}

void Problem::printWeight() const
{
	std::cout << "__________real coils weight__________" << std::endl;
	int width = 8;
	for (int i = 0; i < vertices; i++)
	{
		std::cout << std::setw(width) << std::left << i + 1 << "|" 
			<< std::setw(width) << std::left << in_weight[i] << "|"
			<< std::setw(width) << std::left << materialsFlow[i] << "|"
			<< std::endl;
	}
}

double Problem::getWeight(int i) const
{
	// 注意虚拟节点 0, 1, ... , vertices
	if (i > vertices)
	{
		//std::cout << "Problem : get vertex weight : index exceed " << std::endl;
		return 0;
	}
	if (i == 0)
	{
		return 0;
	}
	else
	{
		return in_weight[i - 1];
	}
}

double Problem::getOutWidth(int i) const
{
	if (i > vertices)
	{
		std::cout << "Problem : get vertex weight : index exceed " << std::endl;
	}
	if (i == 0)
	{
		return 0;
	}
	else
	{
		return out_width[i - 1];
	}
}

std::string Problem::getFlow(int i) const
{
	std::string unit = "virtual coil";
	if (i > vertices)
	{
		std::cout << "Problem : get vertex weight : index exceed " << std::endl;
	}
	if (i == 0)
	{
		return unit;
	}
	else
	{
		return materialsFlow[i - 1];
	}
}

Solution* Problem::createInitialSolution() const
{
	// delete this solution in algorithm (archive)
	auto* solution = new Solution(vertices, 3);

	// the solution include vertex 0, (0, 1, 2, ...)
	std::vector<int> indexCoils;
	indexCoils.resize(vertices);
	for (int i = 0; i < vertices; i++)
	{
		indexCoils[i] = i;
	}

	// sort in descending out width order
	for (int i = 0; i < vertices - 1; i++)
	{
		for (int j = 0; j < vertices - 1 - i; j++)
		{
			if (out_width[indexCoils[j]] < out_width[indexCoils[j + 1]])
			{
				std::swap(indexCoils[j], indexCoils[j + 1]);
			}
		}
	}
	for (int i = 0; i < vertices; i++)
	{
		solution->set_variable(i, indexCoils[i]);
	}
	return solution;
}

bool Problem::feasible(int prev, int post)
{
	bool flag1, flag2, flag3, flag4, flag5;

	//判断入口宽度是否可行
	flag1 = (in_width_matrix[prev][post] > 0 ?
		in_width_matrix[prev][post] < down_width_limit : abs(in_width_matrix[prev][post]) < up_width_limit);
	if (!flag1)
	{
		return false;
	}

	//出口宽度
	flag2 = (out_width_matrix[prev][post] > 0 ?
		out_width_matrix[prev][post] < down_width_limit : abs(out_width_matrix[prev][post]) < up_width_limit);
	if (!flag2)
	{
		return false;
	}

	//入口厚度
	flag3 = (in_thick_matrix[prev][post] > 0 ?
		in_thick_matrix[prev][post] < down_thick_limit : abs(in_thick_matrix[prev][post]) < up_thick_limit);
	if (!flag3)
	{
		return false;
	}

	flag4 = (out_thick_matrix[prev][post] > 0 ?
		out_thick_matrix[prev][post] < down_thick_limit : abs(out_thick_matrix[prev][post]) < up_thick_limit);
	if (!flag4)
	{
		return false;
	}

	flag5 = steel_type_matrix[prev][post];
	return flag5;
}

bool Problem::feasible_info(int prev, int post)
{
	bool flag1, flag2, flag3, flag4, flag5;

	//判断入口宽度是否可行
	flag1 = (in_width_matrix[prev][post] > 0 ?
		in_width_matrix[prev][post] < down_width_limit : abs(in_width_matrix[prev][post]) < up_width_limit);
	if (!flag1)
	{
		//std::cout << "in_width : infeasible" << std::endl;
		return false;
	}

	//出口宽度
	flag2 = (out_width_matrix[prev][post] > 0 ?
		out_width_matrix[prev][post] < down_width_limit : abs(out_width_matrix[prev][post]) < up_width_limit);
	if (!flag2)
	{
		//std::cout << "out_width : infeasible" << std::endl;
		return false;
	}

	//入口厚度
	flag3 = (in_thick_matrix[prev][post] > 0 ?
		in_thick_matrix[prev][post] < down_thick_limit : abs(in_thick_matrix[prev][post]) < up_thick_limit);
	if (!flag3)
	{
		//std::cout << "in_thick : infeasible" << std::endl;
		return false;
	}

	flag4 = (out_thick_matrix[prev][post] > 0 ?
		out_thick_matrix[prev][post] < down_thick_limit : abs(out_thick_matrix[prev][post]) < up_thick_limit);
	if (!flag4)
	{
		//std::cout << "out_thick : infeasible" << std::endl;
		return false;
	}

	flag5 = steel_type_matrix[prev][post];
	if (!flag5) {
		//std::cout << "steel grade : infeasible" << std::endl;
	}
	return flag5;
}

void Problem::evaluate(Solution* solution, const std::map<std::string, double>& _currentFlowTarget)
{
	int* x = solution->get_variables();

	// assign weight
	std::map<std::string, double> weights;
	double totalTarget = 0.0;
	for (const auto& item : _currentFlowTarget)
	{
		totalTarget += item.second;
	}
	for (const auto& item : _currentFlowTarget)
	{
		weights[item.first] = item.second / totalTarget;
	}

	// the weight of coils in each flow unit in this solution
	std::map<std::string, double> flows;
	for (const auto& item : parameter.flow_mass)
	{
		flows[item.first] = 0.0;
	}

	// number of feasible coils
	int count = 1;

	// number of feasible coils (out width)
	int countWidth = 0;

	// total weight in this solution
	double totalWeight = in_weight[x[0]];

	for (int i = 0; i < vertices - 1; i++)
	{
		if (feasible_matrix[x[i]][x[i + 1]])
		{
			totalWeight += in_weight[x[i + 1]];
			count++;
			flows[materialsFlow[x[i + 1]]] += in_weight[i + 1];
			if (out_width[x[i]] > out_width[x[i + 1]])
			{
				countWidth++;
			}
		}
		else
		{
			break;
		}
	}
	double flowRatio = 0.00001;
	for (const auto& item : flows)
	{
		flowRatio += item.second * weights[item.first];
	}
	solution->set_objective(0, countWidth);
	solution->set_objective(1, flowRatio);
	solution->set_objective(2, totalWeight);

	// index of coil : begin to violate out width constraint
	solution->set_number_of_violated_constraints(vertices - count);
}

void Problem::printSolution(Solution* solution)
{
	int width = 15;
	std::cout << "Algorithm___execute___print___solution" << std::endl;
	std::cout << std::setw(width) << std::left << "index" << "|"
		<< std::setw(width) << std::left << "In_width" << "|"
		<< std::setw(width) << std::left << "Out_width" << "|"
		<< std::setw(width) << std::left << "In_thick" << "|"
		<< std::setw(width) << std::left << "Out_thick" << "|"
		<< std::setw(width) << std::left << "Grade" << "|"
		<< std::setw(width) << std::left << "Flow" << "|"
		<< std::setw(width) << std::left << "NO" << "|"
		<< std::endl;
	for (int i = 0; i < vertices; i++)
	{
		std::cout << std::setw(width) << std::left << _materials[solution->get_variable(i)]->index << "|"
			<< std::setw(width) << std::left << in_width[solution->get_variable(i)] << "|"
			<< std::setw(width) << std::left << out_width[solution->get_variable(i)] << "|"
			<< std::setw(width) << std::left << in_thick[solution->get_variable(i)] << "|"
			<< std::setw(width) << std::left << out_thick[solution->get_variable(i)] << "|"
			<< std::setw(width) << std::left << _materials[solution->get_variable(i)]->getStr(Attribute::RULE_CODE) << "|"
			<< std::setw(width) << std::left << materialsFlow[solution->get_variable(i)] << "|"
			<< std::setw(width) << std::left << _materials[solution->get_variable(i)]->getStr(Attribute::MAT_NO) << "|"
			<< std::endl;
	}
}


