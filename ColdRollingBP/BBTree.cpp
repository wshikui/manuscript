//
// Created by wshikui on 2023/11/29.
//
#include "BBTree.h"

BBTree::BBTree(const std::string& _name, Problem* _problem, Parameters* _para, std::vector<Column*>& _initColumns) {
	name = _name;
	problem = _problem;
	para = _para;

	ub = std::numeric_limits<double>::max();
	lb = std::numeric_limits<double>::min();

	for (const auto& item : _initColumns) {
		Column column(item);
		globalPool->emplace_back(column);
	}

	std::vector<int> verticesEquality{};
	std::optional<double> father_lb;

	splitsplitPoolIndex();
	initialize_ng_sets();

	//at the root node, the verticesEquality is empty and
	//the forbidden coils of each coil are initialized in PROBLEM
	//Node* rootNode = new Node(problem, para, verticesEquality, globalPool, *globalPool, father_lb);
	Node* rootNode = new Node(problem, para, verticesEquality, globalPool, *globalPool, father_lb, poolIndex, pool_coil_ngsets);

	unexploredNodes.push(rootNode);

	// nodeAttainUb = rootNode;
	nodeBoundType = BoundType::FROM_LP;
	maxDepth = 0;
	elapsedTime = 0;
	total_time_on_master = 0;
	total_time_on_pricing = 0;
	//the root node
	nodesGenerated = 1;
}

BBTree::~BBTree() {
	delete globalPool;

	// delete nodeAttainUb;
	delete nodeAttainUb;
}

void BBTree::explore_tree() {
	// load the trained xgboost
	//BoosterHandle xg_model;
	//XGBoosterCreate(NULL, 0, &xg_model);
	//XGBoosterLoadModel(xg_model, "Model/xgboost_model_72.json");

	printHeader();
	auto nodeNumber = 0u;
	auto startTime = std::chrono::high_resolution_clock::now();

	//TODO maybe gap<0.1% can stop
	while (!unexploredNodes.empty())
	{
		//std::cerr << "Nodes in tree: " << unexploredNodes.size() << std::endl;
		auto current_node = unexploredNodes.top();
		unexploredNodes.pop(); // the address of current node still exists

		//current_node->solver(nodeNumber++);
		current_node->solverAddSR(nodeNumber++);

		if (current_node->depth > maxDepth) {
			maxDepth = current_node->depth;
		}

		if (!current_node->isFeasible()) {
			if (nodeNumber == 1u) { std::cout << "Root node infeasible" << std::endl; }
			update_lb(current_node, nodeNumber);
			delete current_node;
			gap = std::abs((ub - lb) / ub) * 100;
			continue;
		}

		if (current_node->solValue >= ub) {
			update_lb(current_node, nodeNumber);
			delete current_node;
			gap = std::abs((ub - lb) / ub) * 100;
			continue;
		}

		try_to_obtain_ub(current_node);

		if (current_node->has_fractional_solution() || current_node->has_basic_column_with_cycles() || 
			current_node->has_integer_coils_with_frac_sol()) {
			branch(current_node);
			//branch_strong(current_node, nodeNumber);
			//branch_strong_learn_MF(xg_model, current_node, nodeNumber);
		}

		update_lb(current_node, nodeNumber);
		double gap_node = std::abs((ub - current_node->solValue) / ub) * 100;
		gap = std::abs((ub - lb) / ub) * 100;

		if (nodeNumber == 1u) {
			gapRoot = gap_node;
		}
		total_time_on_master += current_node->totalTimeMP;
		total_time_on_pricing += current_node->totalTimeSp;

		auto curr_time = std::chrono::high_resolution_clock::now();
		auto el_time = std::chrono::duration_cast<std::chrono::duration<double>>(curr_time - startTime).count();
		printRow(current_node, gap_node, el_time);

		if (el_time > 3600) {
			std::cout << "over time limit" << std::endl;
			break;
		}
		if (gap < 0.1) {
			break;
		}

		// add time limit of this while
		delete current_node;
	}
	auto end_time = std::chrono::high_resolution_clock::now();
	elapsedTime = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - startTime).count();
	std::cout << "all run time of branch and price : " << elapsedTime << std::endl;
	std::cout << "all node solved: " << nodeNumber << std::endl;
	printResult();

	//save features of all nodes
	//saveFatures(all_node_features);
}

void BBTree::printHeader() {
	//column name

	std::cout << std::left
		<< std::setw(12) << "Unexplored"
		<< std::setw(10) << "Total"
		<< std::setw(12) << "LB at node"
		<< std::setw(10) << "LB best"
		<< std::setw(10) << "UB best"
		<< std::setw(12) << "Gap at node"
		<< std::setw(10) << "Gap best"
		<< std::setw(18) << "Time of cur MP"
		<< std::setw(16) << "Time of cur SP"
		<< std::setw(12) << "Time Total"
		<< std::setw(8) << "Depth"
		<< std::endl;
}

void BBTree::printRow(const Node* currentNode, double gapNode, double current_time)
{
	auto printUb = (ub < std::numeric_limits<double>::max() - 100);
	std::string printU = "inf";
	if (printUb) {
		printU = std::to_string(ub);
	}
	std::cout << std::setw(12) << std::left << unexploredNodes.size()
		<< std::setw(10) << std::left << nodesGenerated
		<< std::setw(12) << std::left << currentNode->solValue
		<< std::setw(10) << std::left << lb
		<< std::setw(10) << std::left << ub
		<< std::setw(12) << std::left << gapNode
		<< std::setw(10) << std::left << gap
		<< std::setw(18) << std::left << currentNode->totalTimeMP
		<< std::setw(16) << std::left << currentNode->totalTimeSp
		<< std::setw(12) << std::left << current_time
		<< std::setw(8) << std::left << currentNode->depth
		<< std::endl;
}

void BBTree::printResult()
{
	std::cout << "______________RESULT OF BRANCH AND PRICE______________" << std::endl;
	std::cout << "final GAP : " << gap << std::endl;
	std::cout << "lower bound : " << lb << std::endl;
	std::cout << "upper bound : " << ub << std::endl;

	if (nodeAttainUb != nullptr) {
		std::cout << "objective of MIP: " << nodeAttainUb->solValue << std::endl;
		std::cout << "objective of MP: " << nodeAttainUb->mipSolValue << std::endl;
		std::cout << "_____ SOLUTION _____" << std::endl;
		if (nodeBoundType == BoundType::FROM_MIP) {
			std::cout << "from mip" << std::endl;
			for (const auto& item : nodeAttainUb->mipBaseColumns) {
				item.first.printRoute();
			}
		}
		else {
			std::cout << "from lp" << std::endl;
			for (const auto& item : nodeAttainUb->basicColumns) {
				item.first.printRoute();
			}
		}
	}
}

void BBTree::try_to_obtain_ub(Node* current_node) {

	// solve RMP as an integer model

	std::vector<Column> feasibleColumns;
	for (const auto& item : current_node->localPool) {
		feasibleColumns.emplace_back(item);
	}

	if (current_node->has_fractional_solution()) {
		if (current_node->solve_mip(feasibleColumns)) {
			if (ub - current_node->mipSolValue > Node::cplexEpsilon) {
				// update the global upper bound
				ub = current_node->mipSolValue;
				// new one node for nodeAttainUb based on current node
				if (nodeAttainUb != nullptr) {
					delete nodeAttainUb;
					nodeAttainUb = nullptr;
				}
				nodeAttainUb = new Node(*current_node);
				nodeBoundType = BoundType::FROM_MIP;
			}
		}
		else {
			//std::cout << "MIP is feasible" << std::endl;
		}
	}
	else {
		if (ub - current_node->solValue > Node::cplexEpsilon) {
			// update the global upper bound
			ub = current_node->solValue;
			// new one node for nodeAttainUb based on current node
			if (nodeAttainUb != nullptr) {
				delete nodeAttainUb;
				nodeAttainUb = nullptr;
			}
			nodeAttainUb = new Node(*current_node);
			nodeBoundType = BoundType::FROM_LP;

		}
	}

}

void BBTree::update_lb(Node* current_node, unsigned int nodeNumber) {
	// Root node. Global LB is node's LB.
	if (nodeNumber <= 1u) {
		lb = current_node->solValue;
		return;
	}
	// Tree explored completely => Solution found! => LB = UB
	if (unexploredNodes.empty()) {
		lb = ub;
		return;
	}

	lb = *unexploredNodes.top()->fatherLb;
}

void BBTree::branch(Node* current_node) {
	branch_on_vertex_select(current_node) ||
	branch_on_successive_coils(current_node);
}

bool BBTree::branch_on_vertex_select(Node* current_node) {
	// the vertex's value is fractional and most close to 0.5
	double most_fractional_val = .5f;
	int index = -1;
	for (int i = 1; i < problem->getVertices() + 1; i++) {
		double val = 0;
		for (const auto& item : current_node->basicColumns) {
			auto& column = item.first;
			auto& coefficient = item.second;
			val += coefficient * column.contains(i);
		}
		if (val > 0 && val < 1) {
			if (std::fabs(val - .5f) < most_fractional_val) {
				most_fractional_val = std::fabs(val - .5f);
				index = i;
			}
		}
	}
	if (index == -1) {
		std::cout << "no find vertex with fractional flow" << std::endl;
		return false;
	}
	std::shared_ptr<BranchingRule> includeCoil = std::make_shared<IncludeCoil>(index);
	std::shared_ptr<BranchingRule> excludeCoil = std::make_shared<ExcludeCoil>(index);

	Node* includeNode = new Node(*current_node, includeCoil,
		current_node->name + "->select_" + std::to_string(index));
	Node* excludeNode = new Node(*current_node, excludeCoil,
		current_node->name + "->no_select_" + std::to_string(index));

	unexploredNodes.push(includeNode);
	unexploredNodes.push(excludeNode);
	nodesGenerated += 2;

	return true;
}

bool BBTree::branch_on_successive_coils(Node* current_node)
{
	double most_fractional_val = .5f;
	std::optional<std::pair<int, int>> most_fractioanl_successive_coils;

	//get the most fractional column with cycle or successive coils within other column
	for (auto it = current_node->basicColumns.begin(); it != current_node->basicColumns.end(); ++it) {
		const auto& col = it->first;
		const auto& coefficient = it->second;

		//in first iteration, not run this continue
		if (std::fabs(coefficient - 0.5) > most_fractional_val) { continue; }

		bool has_cycle = col.column_has_cycle();
		if (has_cycle) {
			auto visited = col.visited_coils_with_predecessors();
			for (const auto& coil_pred : visited) {
				if (coil_pred.second.size() > 1) {
					most_fractioanl_successive_coils = std::make_pair(coil_pred.second.back(), coil_pred.first);
					most_fractional_val = coefficient - 0.5;
					break;
				}
			}
		}
	}

	if (!most_fractioanl_successive_coils.has_value()) {
		most_fractional_val = .5f;
		for (auto it = current_node->basicColumns.begin(); it != current_node->basicColumns.end(); ++it) {
			const auto& col = it->first;
			const auto& coefficient = it->second;

			for (auto other_it = it + 1; other_it != current_node->basicColumns.end(); ++other_it) {
				const auto& other_sol = other_it->first;
				auto coil_succ = col.common_coil_visited_from_two_different_predecessors(other_sol);
				if (coil_succ.has_value()) {
					most_fractioanl_successive_coils = coil_succ;
					most_fractional_val = coefficient - 0.5;
					break;
				}
			}
		}
	}

	if (most_fractioanl_successive_coils.has_value()) {
		std::pair<int, int> branch_successive_coils = most_fractioanl_successive_coils.value();
		std::shared_ptr<BranchingRule> force_successive = std::make_shared<ForceConsecutive>(branch_successive_coils);
		std::shared_ptr<BranchingRule> forbid_successive = std::make_shared<ForbidConsecutive>(branch_successive_coils);
		Node* forceNode = new Node(*current_node, force_successive,
			current_node->name + "->force_" + std::to_string(branch_successive_coils.first) + "-" + std::to_string(branch_successive_coils.second));
		Node* forbidNode = new Node(*current_node, forbid_successive,
			current_node->name + "->forbid_" + std::to_string(branch_successive_coils.first) + "-" + std::to_string(branch_successive_coils.second));
		return true;
	}

	return  false;
}


void BBTree::branch_on_random(Node* current_node)
{
	std::vector<int> indexes;
	std::vector<double> coeff;
	for (int i = 1; i < problem->getVertices() + 1; i++) {
		double val = 0;
		for (const auto& item : current_node->basicColumns) {
			auto& column = item.first;
			auto& coefficient = item.second;
			val += coefficient * column.contains(i);
		}
		if (val > 0 && val < 1) {
			indexes.emplace_back(i);
			coeff.emplace_back(val);
		}
	}

	//int temp = min_element(coeff.begin(), coeff.end()) - coeff.begin();
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_int_distribution<> dis(0, indexes.size() - 1);
	int randomIndex = dis(gen);

	int index = indexes[randomIndex];

	//std::string node_name = current_node->name;
	//std::cout << node_name << std::endl;

	std::shared_ptr<BranchingRule> includeCoil = std::make_shared<IncludeCoil>(index);
	std::shared_ptr<BranchingRule> excludeCoil = std::make_shared<ExcludeCoil>(index);

	Node* includeNode = new Node(*current_node, includeCoil,
		current_node->name + "->select_" + std::to_string(index));
	Node* excludeNode = new Node(*current_node, excludeCoil,
		current_node->name + "->no_select_" + std::to_string(index));

	unexploredNodes.push(includeNode);
	unexploredNodes.push(excludeNode);
	nodesGenerated += 2;
}

void BBTree::branch_strong(Node* current_node, unsigned int nodeNumber)
{
	//get the feature of candidate variable
	std::map<int, std::map<std::string, double>> node_variables_features;
	//initilize_features(node_variables_features);

	MasterProblem* rmp = current_node->cg->getRMPmodel();

	//std::string modelName = "Model/MasterModel_SR_" + std::to_string(nodeNumber) + ".lp";
	//rmp->getSolver().exportModel(modelName.c_str());

	for (int i = 0; i < problem->getVertices(); i++) {
		//slack and ceil distances
		double val = 0;
		for (const auto& item : current_node->basicColumns) {
			auto& column = item.first;
			auto& coefficient = item.second;
			val += coefficient * column.contains(i + 1);
		}
		if (val <= 0 || val >= 1) {
			continue;
		}

		node_variables_features[i]["is_candidate"] = 1;
		node_variables_features[i]["cand_solfracs"] = val;
		node_variables_features[i]["cand_slack"] = 1 - val;
		//dual and degree
		node_variables_features[i]["dual_coil"] = rmp->getDualVariable(i) / (para->setUpCost * para->alpha2 + 200);
		node_variables_features[i]["coil_cons_degree"] = static_cast<double>(rmp->getVarNum(i)) / problem->getVertices();

		//objective coefficient
		auto coef_values = rmp->getCoefAndValue(i + 1);
		std::vector<double> coefs_ = std::get<0>(coef_values);
		std::vector<double> coefs;
		for (int j = 0; j < coefs_.size(); j++) {
			coefs.emplace_back(coefs_[j] / (para->setUpCost * para->alpha2 + 200));
		}
		std::vector<double> values = std::get<1>(coef_values);
		std::vector<double> product(coefs.size());
		std::transform(coefs.begin(), coefs.end(), values.begin(), product.begin(), std::multiplies<double>());

		node_variables_features[i]["coil_obj_count"] = static_cast<double>(coefs.size()) / (para->N * 2);
		node_variables_features[i]["coil_obj_max"] = *(std::max_element(coefs.begin(), coefs.end()));
		node_variables_features[i]["coil_obj_min"] = *(std::max_element(coefs.begin(), coefs.end()));
		node_variables_features[i]["coil_obj_mean"] = std::accumulate(coefs.begin(), coefs.end(), 0.0) / static_cast<double>(coefs.size());
		node_variables_features[i]["coil_obj_variance"] = Tool::calculateVariance(coefs);

		node_variables_features[i]["coil_y_count"] = static_cast<double>(values.size()) / (para->N * 2);
		node_variables_features[i]["coil_y_max"] = *(std::max_element(values.begin(), values.end()));
		node_variables_features[i]["coil_y_min"] = *(std::max_element(values.begin(), values.end()));
		node_variables_features[i]["coil_y_mean"] = std::accumulate(values.begin(), values.end(), 0.0) / static_cast<double>(values.size());
		node_variables_features[i]["coil_y_variance"] = Tool::calculateVariance(values);

		node_variables_features[i]["coil_product_mean"] = std::accumulate(product.begin(), product.end(), 0.0) / static_cast<double>(product.size());
		node_variables_features[i]["coil_product_max"] = *(std::max_element(product.begin(), product.end()));
		node_variables_features[i]["coil_product_min"] = *(std::min_element(product.begin(), product.end()));
		node_variables_features[i]["coil_product_variance"] = Tool::calculateVariance(product);

		//weight
		node_variables_features[i]["coil_weight"] = problem->getWeight(i + 1) / 50;
		if (problem->getFlow(i + 1) == "C512") {
			node_variables_features[i]["coil_weight_C512"] = problem->getWeight(i + 1) / para->MassFlow[0];
		}
		else {
			node_variables_features[i]["coil_weight_C008"] = problem->getWeight(i + 1) / para->MassFlow[1];
		}
		node_variables_features[i]["coil_weight_RH"] = rmp->getWeightRH(i + 1);

		//number in SR set
		node_variables_features[i]["coil_num_SR"] = rmp->getNumSRofCoil(i + 1);

		//number of basic columns in coil i's constraint
		node_variables_features[i]["coil_constraint_bc"] = rmp->getBCconstraint(i);

		//coil in the integer solution
		int _inIntegerSol = 0;
		bool isFound = false;
		for (const auto& item : nodeAttainUb->mipBaseColumns) {
			for (const auto& _coil : item.first.getRoute()) {
				if (_coil == i + 1) {
					_inIntegerSol = 1;
					isFound = true;
					break;
				}
			}
			if (isFound) {
				break;
			}
		}
		node_variables_features[i]["coil_in_up"] = _inIntegerSol;

		auto _binding = rmp->isBinding(i);
		node_variables_features[i]["coil_cons_bind_count"] = std::get<0>(_binding);
		node_variables_features[i]["coil_cons_bind_max"] = std::get<1>(_binding);
		node_variables_features[i]["coil_cons_bind_min"] = std::get<2>(_binding);
		node_variables_features[i]["coil_cons_bind_mean"] = std::get<3>(_binding);
		node_variables_features[i]["coil_cons_bind_variance"] = std::get<4>(_binding);
	}

	//get the varible with largest score
	double obj_up = 0;
	double obj_down = 0;
	double pseudo_cost = 0;

	int select_coil = -1;

	//label
	std::vector<int> cand_indexes;
	std::vector<double> strong_branch_scores;

	for (int i = 1; i < problem->getVertices() + 1; i++) {
		double val = 0;
		for (const auto& item : current_node->basicColumns) {
			auto& column = item.first;
			auto& coefficient = item.second;
			val += coefficient * column.contains(i);
		}

		if (val > 0 && val < 1) {
			std::shared_ptr<BranchingRule> includeCoil = std::make_shared<IncludeCoil>(i);
			std::shared_ptr<BranchingRule> excludeCoil = std::make_shared<ExcludeCoil>(i);
			Node* includeNode = new Node(*current_node, includeCoil,
				current_node->name + "->select_" + std::to_string(i));
			Node* excludeNode = new Node(*current_node, excludeCoil,
				current_node->name + "->no_select_" + std::to_string(i));
			auto result_up = includeNode->get_lower_sr(nodeNumber);
			auto result_down = excludeNode->get_lower_sr(nodeNumber);
			obj_up = std::get<1>(result_up) - current_node->solValue;
			obj_down = std::get<1>(result_down) - current_node->solValue;
			if (std::max(obj_up, 1e-6) * std::max(obj_down, 1e-6) > pseudo_cost) {
				pseudo_cost = std::max(obj_up, 1e-6) * std::max(obj_down, 1e-6);
				select_coil = i;
			}
			delete includeNode;
			delete excludeNode;

			cand_indexes.emplace_back(i - 1);
			strong_branch_scores.emplace_back(std::max(obj_up, 1e-6) * std::max(obj_down, 1e-6));
		}
	}
	for (int i = 0; i < cand_indexes.size(); i++) {
		if (cand_indexes[i] == select_coil - 1) {
			node_variables_features[cand_indexes[i]]["coil_label"] = 1;
			continue;
		}
		if (strong_branch_scores[i] > 0.8 * pseudo_cost) {
			node_variables_features[cand_indexes[i]]["coil_label"] = 1;
		}
	}
	//save the features of node
	all_node_features.emplace_back(node_variables_features);

	//get the final branch node
	if (select_coil == -1) {
		return;
	}

	std::shared_ptr<BranchingRule> includeCoil = std::make_shared<IncludeCoil>(select_coil);
	std::shared_ptr<BranchingRule> excludeCoil = std::make_shared<ExcludeCoil>(select_coil);

	Node* includeNode = new Node(*current_node, includeCoil,
		current_node->name + "->select_" + std::to_string(select_coil));
	Node* excludeNode = new Node(*current_node, excludeCoil,
		current_node->name + "->no_select_" + std::to_string(select_coil));

	unexploredNodes.push(includeNode);
	unexploredNodes.push(excludeNode);
	nodesGenerated += 2;
}

//void BBTree::branch_strong_learn(BoosterHandle pre_model, Node* current_node, unsigned int nodeNumber)
//{
//    const double epsilon = 1e-6;
//
//    std::map<int, std::map<std::string, double>> node_variables_features;
//    initilize_features(node_variables_features);
//
//    std::vector<int> can_indexes;
//    std::vector<float> features;
//
//    MasterProblem* rmp = current_node->cg->getRMPmodel();
//
//    for (int i = 0; i < problem->getVertices(); i++) {
//        //slack and ceil distances
//        double val = 0;
//        for (const auto& item : current_node->basicColumns) {
//            auto& column = item.first;
//            auto& coefficient = item.second;
//            val += coefficient * column.contains(i + 1);
//        }
//        if (val <= 0 || val >= 1) {
//            continue;
//        }
//        //the idnex of candidate coil
//        can_indexes.emplace_back(i + 1);
//
//        node_variables_features[i]["is_candidate"] = 1;
//        node_variables_features[i]["cand_solfracs"] = val;
//        node_variables_features[i]["cand_slack"] = 1 - val;
//        //dual and degree
//        node_variables_features[i]["dual_coil"] = rmp->getDualVariable(i) / (para->setUpCost * para->alpha2 + 200);
//        node_variables_features[i]["coil_cons_degree"] = static_cast<double>(rmp->getVarNum(i)) / problem->getVertices();
//
//        //objective coefficient
//        auto coef_values = rmp->getCoefAndValue(i + 1);
//        std::vector<double> coefs_ = std::get<0>(coef_values);
//        std::vector<double> coefs;
//        for (int j = 0; j < coefs_.size(); j++) {
//            coefs.emplace_back(coefs_[j] / (para->setUpCost * para->alpha2 + 200));
//        }
//        std::vector<double> values = std::get<1>(coef_values);
//        std::vector<double> product(coefs.size());
//        std::transform(coefs.begin(), coefs.end(), values.begin(), product.begin(), std::multiplies<double>());
//
//        node_variables_features[i]["coil_obj_count"] = static_cast<double>(coefs.size()) / (para->N * 2);
//        node_variables_features[i]["coil_obj_max"] = *(std::max_element(coefs.begin(), coefs.end()));
//        node_variables_features[i]["coil_obj_min"] = *(std::max_element(coefs.begin(), coefs.end()));
//        node_variables_features[i]["coil_obj_mean"] = std::accumulate(coefs.begin(), coefs.end(), 0.0) / static_cast<double>(coefs.size());
//        node_variables_features[i]["coil_obj_variance"] = Tool::calculateVariance(coefs);
//
//        node_variables_features[i]["coil_y_count"] = static_cast<double>(values.size()) / (para->N * 2);
//        node_variables_features[i]["coil_y_max"] = *(std::max_element(values.begin(), values.end()));
//        node_variables_features[i]["coil_y_min"] = *(std::max_element(values.begin(), values.end()));
//        node_variables_features[i]["coil_y_mean"] = std::accumulate(values.begin(), values.end(), 0.0) / static_cast<double>(values.size());
//        node_variables_features[i]["coil_y_variance"] = Tool::calculateVariance(values);
//
//        node_variables_features[i]["coil_product_mean"] = std::accumulate(product.begin(), product.end(), 0.0) / static_cast<double>(product.size());
//        node_variables_features[i]["coil_product_max"] = *(std::max_element(product.begin(), product.end()));
//        node_variables_features[i]["coil_product_min"] = *(std::min_element(product.begin(), product.end()));
//        node_variables_features[i]["coil_product_variance"] = Tool::calculateVariance(product);
//
//        //weight
//        node_variables_features[i]["coil_weight"] = problem->getWeight(i + 1) / 50;
//        if (problem->getFlow(i + 1) == "C512") {
//            node_variables_features[i]["coil_weight_C512"] = problem->getWeight(i + 1) / para->MassFlow[0];
//        }
//        else {
//            node_variables_features[i]["coil_weight_C008"] = problem->getWeight(i + 1) / para->MassFlow[1];
//        }
//        node_variables_features[i]["coil_weight_RH"] = rmp->getWeightRH(i + 1);
//
//        //number in SR set
//        node_variables_features[i]["coil_num_SR"] = rmp->getNumSRofCoil(i + 1);
//
//        //number of basic columns in coil i's constraint
//        node_variables_features[i]["coil_constraint_bc"] = rmp->getBCconstraint(i);
//
//        //coil in the integer solution
//        int _inIntegerSol = 0;
//        bool isFound = false;
//        for (const auto& item : nodeAttainUb->mipBaseColumns) {
//            for (const auto& _coil : item.first.getRoute()) {
//                if (_coil == i + 1) {
//                    _inIntegerSol = 1;
//                    isFound = true;
//                    break;
//                }
//            }
//            if (isFound) {
//                break;
//            }
//        }
//        node_variables_features[i]["coil_in_up"] = _inIntegerSol;
//
//        auto _binding = rmp->isBinding(i);
//        node_variables_features[i]["coil_cons_bind_count"] = std::get<0>(_binding);
//        node_variables_features[i]["coil_cons_bind_max"] = std::get<1>(_binding);
//        node_variables_features[i]["coil_cons_bind_min"] = std::get<2>(_binding);
//        node_variables_features[i]["coil_cons_bind_mean"] = std::get<3>(_binding);
//        node_variables_features[i]["coil_cons_bind_variance"] = std::get<4>(_binding);
//
//        //add variable's feature to vector
//        features.emplace_back(node_variables_features[i]["cand_solfracs"]);
//        features.emplace_back(node_variables_features[i]["cand_slack"]);
//        features.emplace_back(node_variables_features[i]["dual_coil"]);
//        features.emplace_back(node_variables_features[i]["coil_cons_degree"]);
//        features.emplace_back(node_variables_features[i]["coil_obj_count"]);
//        features.emplace_back(node_variables_features[i]["coil_obj_max"]);
//        features.emplace_back(node_variables_features[i]["coil_obj_min"]);
//        features.emplace_back(node_variables_features[i]["coil_obj_mean"]);
//        features.emplace_back(node_variables_features[i]["coil_obj_variance"]);
//        features.emplace_back(node_variables_features[i]["coil_y_count"]);
//        features.emplace_back(node_variables_features[i]["coil_y_max"]);
//        features.emplace_back(node_variables_features[i]["coil_y_min"]);
//        features.emplace_back(node_variables_features[i]["coil_y_mean"]);
//        features.emplace_back(node_variables_features[i]["coil_y_variance"]);
//        features.emplace_back(node_variables_features[i]["coil_product_mean"]);
//        features.emplace_back(node_variables_features[i]["coil_product_max"]);
//        features.emplace_back(node_variables_features[i]["coil_product_min"]);
//        features.emplace_back(node_variables_features[i]["coil_product_variance"]);
//        features.emplace_back(node_variables_features[i]["coil_weight"]);
//        features.emplace_back(node_variables_features[i]["coil_weight_C512"]);
//        features.emplace_back(node_variables_features[i]["coil_weight_C008"]);
//        features.emplace_back(node_variables_features[i]["coil_weight_RH"]);
//        features.emplace_back(node_variables_features[i]["coil_num_SR"]);
//        features.emplace_back(node_variables_features[i]["coil_constraint_bc"]);
//        features.emplace_back(node_variables_features[i]["coil_in_up"]);
//        features.emplace_back(node_variables_features[i]["coil_cons_bind_count"]);
//        features.emplace_back(node_variables_features[i]["coil_cons_bind_max"]);
//        features.emplace_back(node_variables_features[i]["coil_cons_bind_min"]);
//        features.emplace_back(node_variables_features[i]["coil_cons_bind_mean"]);
//        features.emplace_back(node_variables_features[i]["coil_cons_bind_variance"]);
//    }
//    std::vector<float> predictions;
//    bst_ulong out_len = can_indexes.size();
//    bst_ulong num_features = node_variables_features[0].size() - 2;
//
//    const float* out_result;
//    DMatrixHandle dtest;
//    XGDMatrixCreateFromMat(features.data(), out_len, num_features, -999.0f, &dtest);
//    XGBoosterPredict(pre_model, dtest, 0, 0, 0, &out_len, &out_result);
//
//    predictions.assign(out_result, out_result + out_len);
//
//    //Obtain the index with the highest probability value
//    int pre_index = std::distance(predictions.begin(), std::max_element(predictions.begin(), predictions.end()));
//
//    int select_coil = can_indexes[pre_index];
//    std::shared_ptr<BranchingRule> includeCoil = std::make_shared<IncludeCoil>(select_coil);
//    std::shared_ptr<BranchingRule> excludeCoil = std::make_shared<ExcludeCoil>(select_coil);
//
//    Node* includeNode = new Node(*current_node, includeCoil,
//        current_node->name + "->select_" + std::to_string(select_coil));
//    Node* excludeNode = new Node(*current_node, excludeCoil,
//        current_node->name + "->no_select_" + std::to_string(select_coil));
//
//    unexploredNodes.push(includeNode);
//    unexploredNodes.push(excludeNode);
//    nodesGenerated += 2;
//
//}

//void BBTree::branch_strong_learn_MF(BoosterHandle pre_model, Node* current_node, unsigned int nodeNumber)
//{
//    if (current_node->depth >= 3) {
//        double most_fractional_val = .5f;
//        int index = -1;
//        for (int i = 1; i < problem->getVertices() + 1; i++) {
//            double val = 0;
//            for (const auto& item : current_node->basicColumns) {
//                auto& column = item.first;
//                auto& coefficient = item.second;
//                val += coefficient * column.contains(i);
//            }
//            if (val > 0 && val < 1) {
//                if (std::fabs(val - .5f) < most_fractional_val) {
//                    most_fractional_val = std::fabs(val - .5f);
//                    index = i;
//                }
//            }
//        }
//        if (index == -1) {
//            std::cout << "no find vertex with fractional flow" << std::endl;
//            return;
//        }
//        std::shared_ptr<BranchingRule> includeCoil = std::make_shared<IncludeCoil>(index);
//        std::shared_ptr<BranchingRule> excludeCoil = std::make_shared<ExcludeCoil>(index);
//
//        Node* includeNode = new Node(*current_node, includeCoil,
//            current_node->name + "->select_" + std::to_string(index));
//        Node* excludeNode = new Node(*current_node, excludeCoil,
//            current_node->name + "->no_select_" + std::to_string(index));
//
//        unexploredNodes.push(includeNode);
//        unexploredNodes.push(excludeNode);
//        nodesGenerated += 2;
//
//        return;
//    }
//
//    const double epsilon = 1e-6;
//    std::map<int, std::map<std::string, double>> node_variables_features;
//    initilize_features(node_variables_features);
//
//    std::vector<int> can_indexes;
//    std::vector<float> features;
//
//    MasterProblem* rmp = current_node->cg->getRMPmodel();
//
//    for (int i = 0; i < problem->getVertices(); i++) {
//        //slack and ceil distances
//        double val = 0;
//        for (const auto& item : current_node->basicColumns) {
//            auto& column = item.first;
//            auto& coefficient = item.second;
//            val += coefficient * column.contains(i + 1);
//        }
//        if (val <= 0 || val >= 1) {
//            continue;
//        }
//        //the idnex of candidate coil
//        can_indexes.emplace_back(i + 1);
//
//        node_variables_features[i]["is_candidate"] = 1;
//        node_variables_features[i]["cand_solfracs"] = val;
//        node_variables_features[i]["cand_slack"] = 1 - val;
//        //dual and degree
//        node_variables_features[i]["dual_coil"] = rmp->getDualVariable(i) / (para->setUpCost * para->alpha2 + 200);
//        node_variables_features[i]["coil_cons_degree"] = static_cast<double>(rmp->getVarNum(i)) / problem->getVertices();
//
//        //objective coefficient
//        auto coef_values = rmp->getCoefAndValue(i + 1);
//        std::vector<double> coefs_ = std::get<0>(coef_values);
//        std::vector<double> coefs;
//        for (int j = 0; j < coefs_.size(); j++) {
//            coefs.emplace_back(coefs_[j] / (para->setUpCost * para->alpha2 + 200));
//        }
//        std::vector<double> values = std::get<1>(coef_values);
//        std::vector<double> product(coefs.size());
//        std::transform(coefs.begin(), coefs.end(), values.begin(), product.begin(), std::multiplies<double>());
//
//        node_variables_features[i]["coil_obj_count"] = static_cast<double>(coefs.size()) / (para->N * 2);
//        node_variables_features[i]["coil_obj_max"] = *(std::max_element(coefs.begin(), coefs.end()));
//        node_variables_features[i]["coil_obj_min"] = *(std::max_element(coefs.begin(), coefs.end()));
//        node_variables_features[i]["coil_obj_mean"] = std::accumulate(coefs.begin(), coefs.end(), 0.0) / static_cast<double>(coefs.size());
//        node_variables_features[i]["coil_obj_variance"] = Tool::calculateVariance(coefs);
//
//        node_variables_features[i]["coil_y_count"] = static_cast<double>(values.size()) / (para->N * 2);
//        node_variables_features[i]["coil_y_max"] = *(std::max_element(values.begin(), values.end()));
//        node_variables_features[i]["coil_y_min"] = *(std::max_element(values.begin(), values.end()));
//        node_variables_features[i]["coil_y_mean"] = std::accumulate(values.begin(), values.end(), 0.0) / static_cast<double>(values.size());
//        node_variables_features[i]["coil_y_variance"] = Tool::calculateVariance(values);
//
//        node_variables_features[i]["coil_product_mean"] = std::accumulate(product.begin(), product.end(), 0.0) / static_cast<double>(product.size());
//        node_variables_features[i]["coil_product_max"] = *(std::max_element(product.begin(), product.end()));
//        node_variables_features[i]["coil_product_min"] = *(std::min_element(product.begin(), product.end()));
//        node_variables_features[i]["coil_product_variance"] = Tool::calculateVariance(product);
//
//        //weight
//        node_variables_features[i]["coil_weight"] = problem->getWeight(i + 1) / 50;
//        if (problem->getFlow(i + 1) == "C512") {
//            node_variables_features[i]["coil_weight_C512"] = problem->getWeight(i + 1) / para->MassFlow[0];
//        }
//        else {
//            node_variables_features[i]["coil_weight_C008"] = problem->getWeight(i + 1) / para->MassFlow[1];
//        }
//        node_variables_features[i]["coil_weight_RH"] = rmp->getWeightRH(i + 1);
//
//        //number in SR set
//        node_variables_features[i]["coil_num_SR"] = rmp->getNumSRofCoil(i + 1);
//
//        //number of basic columns in coil i's constraint
//        node_variables_features[i]["coil_constraint_bc"] = rmp->getBCconstraint(i);
//
//        //coil in the integer solution
//        int _inIntegerSol = 0;
//        bool isFound = false;
//        for (const auto& item : nodeAttainUb->mipBaseColumns) {
//            for (const auto& _coil : item.first.getRoute()) {
//                if (_coil == i + 1) {
//                    _inIntegerSol = 1;
//                    isFound = true;
//                    break;
//                }
//            }
//            if (isFound) {
//                break;
//            }
//        }
//        node_variables_features[i]["coil_in_up"] = _inIntegerSol;
//
//        auto _binding = rmp->isBinding(i);
//        node_variables_features[i]["coil_cons_bind_count"] = std::get<0>(_binding);
//        node_variables_features[i]["coil_cons_bind_max"] = std::get<1>(_binding);
//        node_variables_features[i]["coil_cons_bind_min"] = std::get<2>(_binding);
//        node_variables_features[i]["coil_cons_bind_mean"] = std::get<3>(_binding);
//        node_variables_features[i]["coil_cons_bind_variance"] = std::get<4>(_binding);
//
//        //add variable's feature to vector
//        features.emplace_back(node_variables_features[i]["cand_solfracs"]);
//        features.emplace_back(node_variables_features[i]["cand_slack"]);
//        features.emplace_back(node_variables_features[i]["dual_coil"]);
//        features.emplace_back(node_variables_features[i]["coil_cons_degree"]);
//        features.emplace_back(node_variables_features[i]["coil_obj_count"]);
//        features.emplace_back(node_variables_features[i]["coil_obj_max"]);
//        features.emplace_back(node_variables_features[i]["coil_obj_min"]);
//        features.emplace_back(node_variables_features[i]["coil_obj_mean"]);
//        features.emplace_back(node_variables_features[i]["coil_obj_variance"]);
//        features.emplace_back(node_variables_features[i]["coil_y_count"]);
//        features.emplace_back(node_variables_features[i]["coil_y_max"]);
//        features.emplace_back(node_variables_features[i]["coil_y_min"]);
//        features.emplace_back(node_variables_features[i]["coil_y_mean"]);
//        features.emplace_back(node_variables_features[i]["coil_y_variance"]);
//        features.emplace_back(node_variables_features[i]["coil_product_mean"]);
//        features.emplace_back(node_variables_features[i]["coil_product_max"]);
//        features.emplace_back(node_variables_features[i]["coil_product_min"]);
//        features.emplace_back(node_variables_features[i]["coil_product_variance"]);
//        features.emplace_back(node_variables_features[i]["coil_weight"]);
//        features.emplace_back(node_variables_features[i]["coil_weight_C512"]);
//        features.emplace_back(node_variables_features[i]["coil_weight_C008"]);
//        features.emplace_back(node_variables_features[i]["coil_weight_RH"]);
//        features.emplace_back(node_variables_features[i]["coil_num_SR"]);
//        features.emplace_back(node_variables_features[i]["coil_constraint_bc"]);
//        features.emplace_back(node_variables_features[i]["coil_in_up"]);
//        features.emplace_back(node_variables_features[i]["coil_cons_bind_count"]);
//        features.emplace_back(node_variables_features[i]["coil_cons_bind_max"]);
//        features.emplace_back(node_variables_features[i]["coil_cons_bind_min"]);
//        features.emplace_back(node_variables_features[i]["coil_cons_bind_mean"]);
//        features.emplace_back(node_variables_features[i]["coil_cons_bind_variance"]);
//    }
//    std::vector<float> predictions;
//    bst_ulong out_len = can_indexes.size();
//    bst_ulong num_features = node_variables_features[0].size() - 2;
//
//    const float* out_result;
//    DMatrixHandle dtest;
//    XGDMatrixCreateFromMat(features.data(), out_len, num_features, -999.0f, &dtest);
//    XGBoosterPredict(pre_model, dtest, 0, 0, 0, &out_len, &out_result);
//
//    predictions.assign(out_result, out_result + out_len);
//
//    //Obtain the index with the highest probability value
//    int pre_index = std::distance(predictions.begin(), std::max_element(predictions.begin(), predictions.end()));
//
//    int select_coil = can_indexes[pre_index];
//    std::shared_ptr<BranchingRule> includeCoil = std::make_shared<IncludeCoil>(select_coil);
//    std::shared_ptr<BranchingRule> excludeCoil = std::make_shared<ExcludeCoil>(select_coil);
//
//    Node* includeNode = new Node(*current_node, includeCoil,
//        current_node->name + "->select_" + std::to_string(select_coil));
//    Node* excludeNode = new Node(*current_node, excludeCoil,
//        current_node->name + "->no_select_" + std::to_string(select_coil));
//
//    unexploredNodes.push(includeNode);
//    unexploredNodes.push(excludeNode);
//    nodesGenerated += 2;
//}
//
//void BBTree::initilize_features(std::map<int, std::map<std::string, double>>& _node_features)
//{
//    for (int i = 0; i < problem->getVertices(); i++) {
//        _node_features[i];
//        //the statics features maybe not necessary
//        //coefs(3) 
//
//        //_node_features[i]["coefs"] = 0; 
//        //_node_features[i]["coefs_pos"] = 0; 
//        //_node_features[i]["coefs_neg"] = 0; 
//
//        ////Numfeatures
//        //_node_features[i]["nnzrs"] = 0;
//        ////Stafeatures degrees(4)
//        //_node_features[i]["root_cdeg_mean"] = 0;
//        //_node_features[i]["root_cdeg_var"] = 0;
//        //_node_features[i]["root_cdeg_min"] = 0;
//        //_node_features[i]["root_cdeg_max"] = 0;
//
//        ////Stafeatures coeffs.(10)
//        //_node_features[i]["root_pcoefs_count"] = 0;
//        //_node_features[i]["root_pcoefs_var"] = 0;
//        //_node_features[i]["root_pcoefs_mean"] = 0;
//        //_node_features[i]["root_pcoefs_min"] = 0;
//        //_node_features[i]["root_pcoefs_max"] = 0;
//        //_node_features[i]["root_ncoefs_count"] = 0;
//        //_node_features[i]["root_ncoefs_var"] = 0;
//        //_node_features[i]["root_ncoefs_mean"] = 0;
//        //_node_features[i]["root_ncoefs_min"] = 0;
//        //_node_features[i]["root_ncoefs_max"] = 0;
//
//        //---featuresamic------------------
//        //Slafeaturesces(2)
//        _node_features[i]["is_candidate"] = 0;
//        //_node_features[i]["is_branch"] = 0;
//        _node_features[i]["cand_solfracs"] = 0;
//        _node_features[i]["cand_slack"] = 0;
//
//        //duafeaturesiables in this constraint
//        _node_features[i]["dual_coil"] = 0;
//        _node_features[i]["coil_cons_degree"] = 0;
//
//        //objfeatures
//        _node_features[i]["coil_obj_count"] = 0;
//        _node_features[i]["coil_obj_mean"] = 0;
//        _node_features[i]["coil_obj_max"] = 0;
//        _node_features[i]["coil_obj_min"] = 0;
//        _node_features[i]["coil_obj_variance"] = 0;
//
//        _node_features[i]["coil_y_count"] = 0;
//        _node_features[i]["coil_y_mean"] = 0;
//        _node_features[i]["coil_y_max"] = 0;
//        _node_features[i]["coil_y_min"] = 0;
//        _node_features[i]["coil_y_variance"] = 0;
//
//        _node_features[i]["coil_product_mean"] = 0;
//        _node_features[i]["coil_product_max"] = 0;
//        _node_features[i]["coil_product_min"] = 0;
//        _node_features[i]["coil_product_variance"] = 0;
//
//        //coifeatures
//        _node_features[i]["coil_weight"] = 0;
//        _node_features[i]["coil_weight_C512"] = 0;
//        _node_features[i]["coil_weight_C008"] = 0;
//        _node_features[i]["coil_weight_RH"] = 0;
//
//        //Psefeaturesranching on this variable, compute
//        //_node_features[i]["cand_obj_up"] = 0;
//        //_node_features[i]["cand_obj_down"] = 0;
//
//        //_node_features[i]["cand_ps_up"] = 0;
//        //_node_features[i]["cand_ps_down"] = 0;
//        //_node_features[i]["cand_ps_sum"] = 0;
//        //_node_features[i]["cand_ps_ratio"] = 0;
//        //_node_features[i]["cand_ps_product"] = 0;
//
//        //Inffeaturesics (4)
//        //_node_features[i]["cand_frac_up_infeas"] = 0;
//        //_node_features[i]["cand_frac_down_infeas"] = 0;
//        //_node_features[i]["infeasible_num"] = 0;
//
//        //Stafeatures degrees (7)
//        //_node_features[i]["cand_cdeg_mean"] = 0;
//        //_node_features[i]["cand_cdeg_var"] = 0;
//        //_node_features[i]["cand_cdeg_min"] = 0;
//        //_node_features[i]["cand_cdeg_max"] = 0;
//        //_node_features[i]["cand_cdeg_mean_ratio"] = 0;
//        //_node_features[i]["cand_cdeg_min_ratio"] = 0;
//        //_node_features[i]["cand_cdeg_max_ratio"] = 0;
//
//        //Minfeatures constraint coeffs. to RHS (4)
//        //_node_features[i]["prhs_ratio_max"] = -1;
//        //_node_features[i]["prhs_ratio_min"] = 1;
//        //_node_features[i]["nrhs_ratio_max"] = -1;
//        //_node_features[i]["nrhs_ratio_min"] = 1;
//
//        //Minfeaturesl coefficient ratios (8)
//        //_node_features[i]["ota_pp_max"] = 0;
//        //_node_features[i]["ota_pp_min"] = 1;
//        //_node_features[i]["ota_pn_max"] = 0;
//        //_node_features[i]["ota_pn_min"] = 1;
//        //_node_features[i]["ota_np_max"] = 0;
//        //_node_features[i]["ota_np_min"] = 1;
//        //_node_features[i]["ota_nn_max"] = 0;
//        //_node_features[i]["ota_nn_min"] = 1;
//
//        _node_features[i]["coil_num_SR"] = 0;
//        _node_features[i]["coil_constraint_bc"] = 0;
//        _node_features[i]["coil_in_up"] = 0;
//        _node_features[i]["coil_cons_bind_count"] = 0;
//        _node_features[i]["coil_cons_bind_max"] = 0;
//        _node_features[i]["coil_cons_bind_min"] = 0;
//        _node_features[i]["coil_cons_bind_mean"] = 0;
//        _node_features[i]["coil_cons_bind_variance"] = 0;
//
//        _node_features[i]["coil_label"] = 0;
//    }
//}

void BBTree::saveFatures(std::vector<std::map<int, std::map<std::string, double>>>& all_features)
{
	std::string filename = "Feature_R/coils_32_2_R.txt";

	std::ofstream outFile(filename);
	if (!outFile.is_open()) {
		std::cerr << "Failed to open file!" << std::endl;
		return;
	}
	//write columns
	outFile << "QID" << "\t" << "cand_solfracs" << "\t" << "cand_slack" << "\t"
		<< "dual_coil" << "\t" << "coil_cons_degree" << "\t" << "coil_obj_count" << "\t"
		<< "coil_obj_mean" << "\t" << "coil_obj_max" << "\t" << "coil_obj_min" << "\t"
		<< "coil_obj_variance" << "\t" << "coil_y_count" << "\t" << "coil_y_mean" << "\t"
		<< "coil_y_max" << "\t" << "coil_y_min" << "\t" << "coil_y_variance" << "\t"
		<< "coil_product_mean" << "\t" << "coil_product_max" << "\t" << "coil_product_min" << "\t"
		<< "coil_product_variance" << "\t" << "coil_weight" << "\t" << "coil_weight_C512" << "\t"
		<< "coil_weight_C008" << "\t" << "coil_weight_RH" << "\t" << "coil_num_SR" << "\t"
		<< "coil_constraint_bc" << "\t" << "coil_in_up" << "\t" << "coil_cons_bind_count" << "\t"
		<< "coil_cons_bind_max" << "\t" << "coil_cons_bind_min" << "\t" << "coil_cons_bind_mean" << "\t"
		<< "coil_cons_bind_variance" << "\t" << "coil_label" << "\t" << std::endl;

	for (int i = 0; i < all_features.size(); i++) {

		// the candidate variables' features in the i-th node
		for (int j = 0; j < all_features[i].size(); j++) {
			if (all_features[i][j]["is_candidate"] == 1) {
				outFile << i << "\t"
					<< all_features[i][j]["cand_solfracs"] << "\t"
					<< all_features[i][j]["cand_slack"] << "\t"
					<< all_features[i][j]["dual_coil"] << "\t"
					<< all_features[i][j]["coil_cons_degree"] << "\t"
					<< all_features[i][j]["coil_obj_count"] << "\t"
					<< all_features[i][j]["coil_obj_mean"] << "\t"
					<< all_features[i][j]["coil_obj_max"] << "\t"
					<< all_features[i][j]["coil_obj_min"] << "\t"
					<< all_features[i][j]["coil_obj_variance"] << "\t"
					<< all_features[i][j]["coil_y_count"] << "\t"
					<< all_features[i][j]["coil_y_mean"] << "\t"
					<< all_features[i][j]["coil_y_max"] << "\t"
					<< all_features[i][j]["coil_y_min"] << "\t"
					<< all_features[i][j]["coil_y_variance"] << "\t"
					<< all_features[i][j]["coil_product_mean"] << "\t"
					<< all_features[i][j]["coil_product_max"] << "\t"
					<< all_features[i][j]["coil_product_min"] << "\t"
					<< all_features[i][j]["coil_product_variance"] << "\t"
					<< all_features[i][j]["coil_weight"] << "\t"
					<< all_features[i][j]["coil_weight_C512"] << "\t"
					<< all_features[i][j]["coil_weight_C008"] << "\t"
					<< all_features[i][j]["coil_weight_RH"] << "\t"
					<< all_features[i][j]["coil_num_SR"] << "\t"
					<< all_features[i][j]["coil_constraint_bc"] << "\t"
					<< all_features[i][j]["coil_in_up"] << "\t"
					<< all_features[i][j]["coil_cons_bind_count"] << "\t"
					<< all_features[i][j]["coil_cons_bind_max"] << "\t"
					<< all_features[i][j]["coil_cons_bind_min"] << "\t"
					<< all_features[i][j]["coil_cons_bind_mean"] << "\t"
					<< all_features[i][j]["coil_cons_bind_variance"] << "\t"
					<< all_features[i][j]["coil_label"] << "\t"
					<< std::endl;
			}
		}
	}
}

void BBTree::splitsplitPoolIndex()
{
	for (int c = 0; c < para->C; c++)
	{
		std::string name = para->chs_type[c];
		std::vector<int> index = { 0 };
		for (int i = 0; i < para->S; i++)
		{
			index.emplace_back(c * para->S + 1 + i);
		}
		poolIndex[name] = index;
	}
	// index of connect pool
	std::string conName1 = "C502016_17";
	std::vector<int> index1 = { 0 };
	index1.insert(index1.end(), poolIndex["C502016"].begin() + 1, poolIndex["C502016"].end());
	index1.insert(index1.end(), poolIndex["C502017"].begin() + 1, poolIndex["C502017"].end());
	poolIndex[conName1] = index1;

	std::string conName2 = "C502017_18";
	std::vector<int> index2 = { 0 };
	index2.insert(index2.end(), poolIndex["C502017"].begin() + 1, poolIndex["C502017"].end());
	index2.insert(index2.end(), poolIndex["C502018"].begin() + 1, poolIndex["C502018"].end());
	poolIndex[conName2] = index2;

	std::string conName3 = "C502016_18";
	std::vector<int> index3 = { 0 };
	index3.insert(index3.end(), poolIndex["C502016"].begin() + 1, poolIndex["C502016"].end());
	index3.insert(index3.end(), poolIndex["C502018"].begin() + 1, poolIndex["C502018"].end());
	poolIndex[conName3] = index3;

	std::string conName4 = "C502018_19";
	std::vector<int> index4 = { 0 };
	index4.insert(index4.end(), poolIndex["C502018"].begin() + 1, poolIndex["C502018"].end());
	index4.insert(index4.end(), poolIndex["C502019"].begin() + 1, poolIndex["C502019"].end());
	poolIndex[conName4] = index4;
}

void BBTree::initialize_ng_sets()
{
	for (const auto& item : poolIndex) {
		int num_coils_pool = item.second.size();
		pool_coil_ngsets[item.first].resize(item.second.size());
		pool_coil_ngsets[item.first][0].reset();
		pool_coil_ngsets[item.first][0].set(0);

		for (int i = 1; i < num_coils_pool; i++) {
			pool_coil_ngsets[item.first][i].reset();
			pool_coil_ngsets[item.first][i].set(poolIndex[item.first][i]);
			std::vector<size_t> neighbors = get_k_nearest_neighbors(poolIndex[item.first][i], poolIndex[item.first], 6);
			for (const auto& index : neighbors) {
				pool_coil_ngsets[item.first][i].set(poolIndex[item.first][index]);
			}
			pool_coil_ngsets[item.first][i].reset(0);
		}
	}
	
	//std::string name = "C502016_18";
	//for (int i = 0; i < poolIndex[name].size(); i++) {
	//	std::cout << i << ": ";
	//	for (int j = 0; j < Config::MAX_COILS; j++) {
	//		if (pool_coil_ngsets[name][i][j] == 1) {
	//			std::cout << j << " ";
	//		}
	//	}
	//	std::cout << std::endl;
	//}
}

std::vector<size_t> BBTree::get_k_nearest_neighbors(int coil_index, std::vector<int>& all_coils, int k_nearest) const
{
	int k = std::min(k_nearest, int(all_coils.size()));
	std::vector<int> indices(all_coils.size());
	std::vector<double> cost_coil_to_other(all_coils.size());
	std::iota(indices.begin(), indices.end(), 0); //0,1,2,... not real index
	cost_coil_to_other[0] = 60000;
	
	for (size_t i = 1; i < all_coils.size(); i++) {
		cost_coil_to_other[i] = problem->getCost(coil_index - 1, all_coils[i] - 1);
	}

	std::partial_sort(
		indices.begin(),
		indices.begin() + k,
		indices.end(),
		[&](int i, int j) {
			return cost_coil_to_other[i] < cost_coil_to_other[j];
		}
	);

	return std::vector<size_t>(indices.begin(), indices.begin() + k);
}
