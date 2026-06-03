#pragma once
#include "Node.h"
#include <queue>
#include <bitset>
#include "xgboost/c_api.h"

enum class BoundType {
	FROM_LP, FROM_MIP
};

class BBTree {
public:
	std::priority_queue<Node*, std::vector<Node*>, BBNodeCompare> unexploredNodes;
	Problem* problem{};
	Parameters* para{};

	std::vector<Column>* globalPool = new std::vector<Column>();

	double lb{};
	double ub{};

	int nodesGenerated{};

	Node* nodeAttainUb{};
	BoundType nodeBoundType = BoundType::FROM_LP;

	double gapRoot{};
	double gap{};
	double elapsedTime{};
	int maxDepth{};
	double total_time_on_master{};
	double total_time_on_pricing{};

	//for ng routes
	std::map<std::string, std::vector<int>> poolIndex;
	std::map<std::string, std::vector<std::bitset<Config::MAX_COILS>>> pool_coil_ngsets;

	std::string name;

	BBTree() = default;
	BBTree(const std::string& _name, Problem* _problem, Parameters* _para, std::vector<Column*>& _initColumns);
	~BBTree();
	void explore_tree();

private:
	void printHeader();
	void printRow(const Node* currentNode, double gapNode, double current_time);
	void printResult();

	void try_to_obtain_ub(Node* current_node);
	void update_lb(Node* current_node, unsigned int nodeNumber);

	// create two nodes to unexplored node list
	void branch(Node* current_node);
	bool branch_on_vertex_select(Node* current_node);
	bool branch_on_successive_coils(Node* current_node);
	void branch_strong(Node* current_node, unsigned int nodeNumber);

	//202602: for coding without xgboost
	void branch_strong_learn(BoosterHandle pre_model, Node* current_node, unsigned int nodeNumber);

	void initilize_features(std::map<int, std::map<std::string, double>>& _node_features);

	//the candidate variables' features
	//std::map<int, std::map<std::string, double>> node_variables_features;
	std::vector<std::map<int, std::map<std::string, double>>> all_node_features;
	void saveFatures(std::vector<std::map<int, std::map<std::string, double>>>& all_features);

	//split the sub-problems and get the nearest coils
	void splitsplitPoolIndex();
	void initialize_ng_sets();

	std::vector<size_t> get_k_nearest_neighbors(int coil_index, std::vector<int>& all_coils, int k_nearest) const;

};

