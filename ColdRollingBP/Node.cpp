//
// Created by wsk on 2023/11/29.
//

#include <utility>
#include "Node.h"

Node::Node(Problem* _pro, Parameters* _para, std::vector<int>& _vertices, std::vector<Column>* _pool, 
	std::vector<Column>& _localPool, std::optional<double> _fatherLb, std::map<std::string, std::vector<int>> _poolIndex, 
	std::map<std::string, std::vector<std::bitset<Config::MAX_COILS>>> _pool_coil_ngsets)
{
	pro = _pro;
	para = _para;

	verticesWithEquality = _vertices;
	forbid_coils = pro->forbidVertex;

	globalPool = _pool;
	localPool = _localPool;
	fatherLb = _fatherLb;
	depth = 0;
	name = "root";

	totalTimeSp = 0;
	totalTimeMP = 0;
	totalTime = 0;
	solValue = std::numeric_limits<double>::max();
	mipSolValue = std::numeric_limits<double>::max();

	this->cgSolve = true;

	poolIndex = _poolIndex;
	pool_coil_ngsets = _pool_coil_ngsets;
}


Node::Node(const Node& father_node, std::shared_ptr<BranchingRule> branchingRule, std::string _name)
{
	pro = father_node.pro;
	para = father_node.para;

	name = _name;
	fatherLb = father_node.solValue;
	depth = father_node.depth + 1;

	solValue = std::numeric_limits<double>::max();
	mipSolValue = std::numeric_limits<double>::max();
	totalTimeSp = 0;
	totalTimeMP = 0;
	totalTime = 0;

	/*
	 * equality constraints of father node
	 */
	verticesWithEquality = father_node.verticesWithEquality;
	vertexWithZero = father_node.vertexWithZero;
	branchRule = std::move(branchingRule);
	S = father_node.S;	

	globalPool = father_node.globalPool;
	localPool = father_node.localPool;

	poolIndex = father_node.poolIndex;
	pool_coil_ngsets = father_node.pool_coil_ngsets;

	//according to branching rule, remove the incompatible columns in RMP and coils in sub-pro
	remove_incompatible_columns();
	determine_equality_constraints();
	no_select_vertex();

	//change the forbid coils with branching on successive coils
	add_forbid_coils();

}

Node::Node(const Node& father_node)
{
	this->pro = father_node.pro;
	this->para = father_node.para;

	this->name = father_node.name;
	// the lower bound of this node from father node
	this->fatherLb = father_node.fatherLb;
	this->depth = father_node.depth;

	this->solValue = father_node.solValue;
	this->cgSolve = true;
	this->mipSolValue = father_node.mipSolValue;

	this->totalTime = father_node.totalTime;
	this->totalTimeMP = father_node.totalTimeMP;
	this->totalTimeSp = father_node.totalTimeSp;

	// the branch rule
	this->verticesWithEquality = father_node.verticesWithEquality;
	this->vertexWithZero = father_node.vertexWithZero;

	this->localPool = father_node.localPool;
	this->basicColumns = father_node.basicColumns;
	this->mipBaseColumns = father_node.mipBaseColumns;
}

Node::~Node()
{
	delete cg;
}

void Node::solver(unsigned int node_number)
{
	basicColumns = std::vector<std::pair<Column, double>>();

	auto nodeStart = std::chrono::high_resolution_clock::now();
	//cg = new ColumnGeneration(localPool, verticesWithEquality, vertexWithZero, para, pro);
	cg = new ColumnGeneration(localPool, verticesWithEquality, vertexWithZero, poolIndex, pool_coil_ngsets, para, pro);
	cgSolve = cg->execute();
	auto nodeEnd = std::chrono::high_resolution_clock::now();
	totalTime = std::chrono::duration_cast<std::chrono::duration<double>>(nodeEnd - nodeStart).count();

	if (!cgSolve) {
		return;
	}

	cg->updateColumnPool(localPool);
	cg->getBasicColumns(basicColumns);
	cg->getTime(totalTimeMP, totalTimeSp);
	solValue = cg->getSolValue();
}

void Node::solverAddSR(unsigned int node_number)
{
	basicColumns = std::vector<std::pair<Column, double>>();

	auto nodeStart = std::chrono::high_resolution_clock::now();

	cg = new ColumnGeneration(localPool, verticesWithEquality, vertexWithZero, S, forbid_coils, poolIndex, pool_coil_ngsets, para, pro);
	cgSolve = cg->executeAndCutR_depth(node_number, depth);

	//TODO recode the pool name during iteration process
	//std::map<int, std::vector<std::string>> temp = cg->getPoolInfo();
	//Tool::writeToFile(temp, "Log/pool_info_less_4_12.txt");


	auto nodeEnd = std::chrono::high_resolution_clock::now();
	totalTime = std::chrono::duration_cast<std::chrono::duration<double>>(nodeEnd - nodeStart).count();
	
	//std::cout << "the runtime at root node: " << totalTime << std::endl;

	if (!cgSolve) {
		return;
	}

	cg->updateColumnPool(localPool);
	cg->getBasicColumns(basicColumns);
	cg->getTime(totalTimeMP, totalTimeSp);
	solValue = cg->getSolValue();

	//update the S of this node
	MasterProblem* rmpModel = cg->getRMPmodel();
	std::vector<std::vector<int>> sInRMP = rmpModel->getSubSet();
	for (const auto& _s : sInRMP) {
		if (std::find(this->S.begin(), this->S.end(), _s) == this->S.end()) {
			this->S.emplace_back(_s);
		}
	}

}

//bool Node::solve_lp(std::vector<Column>& pool, std::vector<int>& vertexEquality)
//{
//	MasterProblem* linear = new MasterProblem(pool, vertexEquality, pro, para);
//	bool flag = linear->optimize();
//	delete linear;
//
//	return flag;
//}

bool Node::solve_mip(std::vector<Column>& pool)
{
	MasterProblem* integerRMP = new MasterProblem(pool, verticesWithEquality, pro, para, true);
	integerRMP->getModel().setName("solver_MIP_up");
	bool flag = integerRMP->optimize();

	if (flag) {
		mipSolValue = integerRMP->getObj();
		// update the mip basic columns
		IloNumVarArray y = integerRMP->getY();
		IloCplex cplex = integerRMP->getSolver();
		for (int i = 0; i < y.getSize(); i++) {
			double value = cplex.getValue(y[i]);
			if (value > 0) {
				mipBaseColumns.emplace_back(std::make_pair(pool[i], value));
			}
		}
	}

	delete integerRMP;
	return flag;
}

std::tuple<bool, double> Node::get_lower_sr(unsigned int node_number)
{
	double obj = 0;
	ColumnGeneration* cg_temp = new ColumnGeneration(localPool, verticesWithEquality, vertexWithZero, S, forbid_coils, poolIndex, pool_coil_ngsets, para, pro);
	// ColumnGeneration* cg_temp = new ColumnGeneration(localPool, verticesWithEquality, vertexWithZero, S, para, pro);
	bool flag = cg_temp->executeAndCutR_depth(node_number, depth);

	if (!flag) {
		obj = this->fatherLb.value() * 2;
	}
	else {
		obj = cg_temp->getSolValue();
	}
	delete cg_temp;

	return std::make_tuple(flag, obj);
}

bool Node::isFeasible() const
{
	if (!cgSolve) {
		return false;
	}
	bool isVertex = true;
	for (int i = 1; i < pro->getVertices() + 1; i++)
	{
		double result = 0;
		for (int r = 0; r < basicColumns.size(); r++)
		{
			result += basicColumns[r].first.contains(i) * basicColumns[r].second;
		}
		if (result >= 1 + 1e-5)
		{
			isVertex = false;
		}
	}

	bool flag1;
	bool flag2;
	double mass_512 = 0;
	double mass_008 = 0;
	for (const auto& item : basicColumns)
	{
		for (int i = 1; i < pro->getVertices() + 1; i++)
		{
			if (item.first.contains(i))
			{
				pro->getFlow(i) == "C512" ? mass_512 += pro->getWeight(i) : mass_008 += pro->getWeight(i);
			}
		}
	}
	flag1 = (mass_512 >= para->flow_mass["C512"]);
	flag2 = (mass_008 >= para->flow_mass["008"]);

	return isVertex && flag1 && flag2;
}

bool Node::has_fractional_solution() const
{
	return std::any_of(basicColumns.begin(), basicColumns.end(),
		[](const auto& cc) {
			return cc.second > 0.0 && cc.second < 1.0;
		});
}

bool Node::has_basic_column_with_cycles() const
{
	return std::any_of(basicColumns.begin(), basicColumns.end(),
		[](const auto& cc) {
			return cc.first.column_has_cycle();
		}
	);
}

bool Node::has_integer_coils_with_frac_sol() const
{
	bool is_all_coils_integer = true;
	for (int i = 1; i < pro->getVertices() + 1; i++) {
		double val = 0;
		for (const auto& item : basicColumns) {
			val = item.first.contains(i) * item.second;
		}
		if (val > 0 && val < 1) {
			is_all_coils_integer = false;
		}
	}
	if (is_all_coils_integer) {
		if (has_fractional_solution()) {
			return true;
		}
	}
	return false;
}

void Node::determine_equality_constraints()
{
	if (!branchRule)
	{
		return;
	}
	// for all real vertices (N)
	for (int i = 1; i < pro->getVertices() + 1; i++)
	{
		if (branchRule->should_row_be_equality(i))
		{
			verticesWithEquality.emplace_back(i);
		}
	}
}

// exclude some incompatible columns
void Node::remove_incompatible_columns()
{
	if (!branchRule)
	{
		return;
	}
	std::vector<Column> newLocalPool;
	for (const auto& item : localPool)
	{
		if (branchRule->is_column_compatible(item))
		{
			newLocalPool.emplace_back(item);
		}
	}
	localPool = newLocalPool;
}

// update vector: index of unselected coil
void Node::no_select_vertex()
{
	if (!branchRule)
	{
		return;
	}
	for (int i = 1; i < pro->getVertices() + 1; i++)
	{
		if (branchRule->remove_vertex(i))
		{
			vertexWithZero.emplace_back(i);
		}
	}
}

void Node::add_forbid_coils()
{
	if (!branchRule)
	{
		return;
	}
	for (int i = 1; i < pro->getVertices() + 1; i++) {
		for (int j = 1; j < pro->getVertices() + 1; j++) {
			//add forbid coils to each coil
			if (branchRule->remove_arc(i, j)) {
				forbid_coils[i].push_back(j);
			}
		}
	}
}


bool Node::hasDuplicates(const std::vector<std::vector<int>>& A)
{
	for (const auto& vec : A) {
		std::vector<int> sortedVec = vec;
		std::sort(sortedVec.begin(), sortedVec.end());
		for (size_t i = 1; i < sortedVec.size(); ++i) {
			if (sortedVec[i] == sortedVec[i - 1]) {
				return true;
			}
		}
	}
	return false;
}


