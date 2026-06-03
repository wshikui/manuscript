#pragma once

#include "Column.h"
#include "MasterProblem.h"
#include <optional>
#include "chrono"
#include <set>
#include <bitset>


class Label {
public:
	int curIndex;
	double cost;
	double demand;

	static constexpr double EPS = 0.001;

	Label() {
		cost = 0;
		demand = 0;
		curIndex = 0;
	}

	Label(double _cost, double _demand, int _curIndex) :
		cost{ _cost },
		demand{ _demand },
		curIndex{ _curIndex } {}
};


class ElementaryLabel :public Label {
public:
	std::vector<Material* > materials;

	//ensure elementary
	std::vector<int> compatibleVertex;

	//erase the forbid coils based on rhe compatibleVertex
	std::vector<int> extensionVertex;

	ElementaryLabel(std::vector<Material* > _coils) :
		Label{},
		materials{ _coils } {}

	ElementaryLabel(std::vector<int> _index) :
		Label{},
		extensionVertex{ _index },
		compatibleVertex{ _index } {}

	ElementaryLabel(double _cost, double _demand, int _curIndex, std::vector<Material* > _coils, std::vector<int> _comVertex) :
		Label{ _cost, _demand, _curIndex },
		compatibleVertex{ _comVertex },
		materials{ _coils } {}
};


class ElementaryLabelSR : public Label {
public:
	std::vector<Material* > materials;
	std::vector<int> compatibleVertex;

	//erase the forbid coils based on rhe compatibleVertex
	std::vector<int> extensionVertex;
	std::vector<std::vector<int>> multiS;
	std::vector<int> numVertexinS;

	MasterProblem* rmp;
	int consNum;

	std::vector<int> allIndex;

	//initialize the label_SR : not visit any vertex in S
	ElementaryLabelSR(std::vector<int> _index, std::vector<std::vector<int>> _multiS, MasterProblem* _rmp, int _cons) :
		Label{},
		rmp{ _rmp },
		consNum{ _cons },
		multiS{ _multiS },
		numVertexinS(_multiS.size(), 0),
		allIndex{ _index },
		extensionVertex{ _index },
		compatibleVertex{ _index } {}
};

class LabelSR :public Label {
public:
	//get lower bound of plan
	Parameters* para;

	std::vector<Material* > materials;
	std::vector<int> extensionVertex;
	std::vector<int> compatibleVertex;
	std::bitset<Config::MAX_COILS> ng_memory;

	//the sub-set from master problem
	std::vector<std::vector<int>> multiS;

	//the set to contain the coils belong to sub-sets
	std::vector<std::set<int>> recordSets;
	std::vector<int> numVertexinS;

	MasterProblem* rmp;
	int consNum;

	//1,2,...n, n+1, maybe remove some coils due to branching rule.
	std::vector<int> allIndex;

	//initialize the label_SR : not visit any vertex in S
	LabelSR(Parameters* _para, std::vector<int> _index, std::vector<std::vector<int>> _multiS, MasterProblem* _rmp, int _cons) :
		Label{},
		para{ _para },
		rmp{ _rmp },
		consNum{ _cons },
		multiS{ _multiS },
		numVertexinS(_multiS.size(), 0),
		recordSets(_multiS.size()),
		allIndex{ _index },
		extensionVertex{ _index } {}

	bool is_feasible() const {
		return this->demand >= para->planLower;
	}

	int get_bucket() const {
		return static_cast<int>(this->demand / Config::BUCKET_SIZE);
	}

	bool operator<(const LabelSR& other) const {
		double lhs_cost = this->cost;
		for (int i = 0; i < this->multiS.size(); i++) {
			if (this->numVertexinS[i] > other.numVertexinS[i]) {
				lhs_cost -= this->rmp->getDualVariable(this->consNum + i);
			}
		}
		if (other.cost < lhs_cost - Label::EPS) { return false; }

		//this label visit some coils additionally.
		if ((this->ng_memory & ~other.ng_memory).any()) { return false; }

		//demand 
		bool this_feas = this->is_feasible();
		bool other_feas = other.is_feasible();

		if (this_feas && other_feas) {
			return this->demand <= other.demand + Label::EPS;
		}
		else if (!this_feas && !other_feas) {
			int this_bucket = this->get_bucket();
			int other_bucket = other.get_bucket();
			if (this_bucket == other_bucket) {
				return this->cost <= other.cost + Label::EPS;
			}
			else {
				return false;
			}
		}
		else {
			//in different buckets, different of demand is so large.
			return false;
		}

		return true;
	}
};


class LabelExtender {
public:
	const std::vector<int>& erasedVertex;
	Problem* problem;
	Parameters* parameters;
	MasterProblem* rmp;

	//exclude coils based on branching rule
	std::vector<int> all_index;
	std::map<int, std::vector<int>> forbidIndex;

	//get the true index of nearest set of coils
	std::vector<int> original_index;
	std::vector<std::bitset<Config::MAX_COILS>> coils_ng_sets;

	LabelExtender(const std::vector<int>& _indexes, Problem* _problem, MasterProblem* _rmp, Parameters* _para) :
		erasedVertex{ _indexes },
		problem{ _problem },
		parameters{ _para },
		forbidIndex{ _problem->forbidVertex },
		rmp{ _rmp } {}

	//with ng - nearest set
	LabelExtender(const std::vector<int>& _all_index,
		const std::vector<int>& _indexes, 
		std::map<int, std::vector<int>>& _forbid_coils,
		std::vector<int>& _original_index,
		std::vector<std::bitset<Config::MAX_COILS>>& _coils_ng_sets,
		Problem* _problem, MasterProblem* _rmp, Parameters* _para) :
		all_index{ _all_index },
		erasedVertex{ _indexes },
		original_index{ _original_index },
		coils_ng_sets{ _coils_ng_sets },
		problem{ _problem },
		parameters{ _para },
		forbidIndex{ _forbid_coils },
		rmp{ _rmp } {}


	std::optional<ElementaryLabel> operator()(const ElementaryLabel& label, int targetIndex, int end) const;
	std::optional<ElementaryLabelSR> operator()(const ElementaryLabelSR& label, int targetIndex, int end) const;
	std::optional<LabelSR> operator()(const LabelSR& label, int targetIndex, int end) const;
};


bool operator==(const Label& lhs, const Label& rhs);
bool operator!=(const Label& lhs, const Label& rhs);
bool operator<=(const Label& lhs, const Label& rhs);
bool operator<(const Label& lhs, const Label& rhs);
std::ostream& operator<<(std::ostream& out, const Label& l);

bool operator==(const ElementaryLabel& lhs, const ElementaryLabel& rhs);
bool operator!=(const ElementaryLabel& lhs, const ElementaryLabel& rhs);
bool operator<=(const ElementaryLabel& lhs, const ElementaryLabel& rhs);
bool operator<(const ElementaryLabel& lhs, const ElementaryLabel& rhs);
std::ostream& operator<<(std::ostream& out, const ElementaryLabel& l);


bool operator==(const ElementaryLabelSR& lhs, const ElementaryLabelSR& rhs);
bool operator!=(const ElementaryLabelSR& lhs, const ElementaryLabelSR& rhs);
bool operator<=(const ElementaryLabelSR& lhs, const ElementaryLabelSR& rhs);
bool operator<(const ElementaryLabelSR& lhs, const ElementaryLabelSR& rhs);

bool operator==(const LabelSR& lhs, const LabelSR& rhs);
bool operator!=(const LabelSR& lhs, const LabelSR& rhs);
bool operator<=(const LabelSR& lhs, const LabelSR& rhs);
//bool operator<(const LabelSR& lhs, const LabelSR& rhs);


double getLabelLowerBound(const ElementaryLabel& label, const LabelExtender& extension);

//solve the fathom problem by dynamic programe
double getLabelLowerBoundDy(const ElementaryLabel& label, const LabelExtender& extension);
double getLabelLowerBoundDy(const ElementaryLabelSR& label, const LabelExtender& extension);
double getLabelLowerBoundDy(const LabelSR& label, const LabelExtender& extension);


template<typename Lbl>
class LblContainer {
public:
	Lbl label;
	const LblContainer* preContainer;
	std::optional<int> preIndex;

	mutable bool isDominated;

	// the label is not dominated by default - isDominated - fales
	LblContainer(Lbl label,
		const LblContainer* preContainer,
		const int preIndex) :
		isDominated{ false },
		label{ label },
		preContainer{ preContainer },
		preIndex{ preIndex } {}

	//for the begin label container
	LblContainer(Lbl label) :
		label{ label },
		preContainer{ nullptr },
		preIndex{ std::nullopt },
		isDominated{ false } {}

};

//for printing the coils in a label
template <typename Lbl>
std::string get_coils_path(const LblContainer<Lbl>& c1) {
	std::vector<int> indexes;
	const LblContainer<Lbl>* current = &c1;
	indexes.push_back(current->label.curIndex);

	while (current->preContainer != nullptr) {
		indexes.push_back(current->preIndex.value());
		current = current->preContainer;
	}
	std::string name = "";
	for (auto rit = indexes.rbegin(); rit != indexes.rend(); ++rit) {
		name += std::to_string(*rit);
		name += "_";
	}

	return "label_" + name;
}

//not consider the Vi and from small to big
//TODO no consider the compatible vertex

//template<typename Lbl>
//struct LblContainerComp {
//	bool operator()(const LblContainer<Lbl>& c1, const LblContainer<Lbl>& c2) const {
//		if (c1.label.demand < c2.label.demand) { return true; }
//		return (c1.label.cost < c2.label.cost);
//	}
//};

template<typename Lbl>
struct LblContainerComp {
	bool operator()(const LblContainer<Lbl>& c1, const LblContainer<Lbl>& c2) const {
		if (c1.label.demand < c2.label.demand) { return true; }
		if (c1.label.cost < c2.label.cost) { return true; }
		return c1.label.compatibleVertex != c2.label.compatibleVertex;
	}
};


template<typename Lbl>
class ContainersSet {
	std::set<LblContainer<Lbl>, LblContainerComp<Lbl>> set;
public:
	using iterator = typename std::set<LblContainer<Lbl>, LblContainerComp<Lbl>>::iterator;

	ContainersSet(std::initializer_list<LblContainer<Lbl>> l) : set{ l } {}
	ContainersSet() {}

	friend auto begin(const ContainersSet& s) { return s.set.begin; }
	friend auto end(const ContainersSet& s) { return s.set.end; }

	//is there a non-dominated label
	bool has_undominated_labels() const {
		if (set.empty()) { return false; }
		if (std::any_of(set.begin(), set.end(),
			[](const LblContainer<Lbl>& c) {return !c.isDominated; })) {
			return true;
		}
	}

	auto first_undominated_container() const {
		for (auto it = set.begin(); it != set.end(); ++it) {
			if (!it->isDominated) {
				return it;
			}
		}
		return set.end();
	}

	//get the first label container
	auto first_label_container() const {
		return set.begin();
	}

	const iterator begin() const { return set.begin(); }
	const iterator end() const { return set.end(); }
	iterator begin() { return set.begin(); }
	iterator end() { return set.end(); }

	auto insert(LblContainer<Lbl> l) { return set.insert(l); }
	auto erase(const iterator& l) { return set.erase(l); }
	auto size() const { return set.size(); }
	auto empty() const { return set.empty(); }

	void mark_dominated(iterator l) {
		//assert(l != set.end());
		l->isDominated = true;
	}

};


template<typename Lbl>
class VertexContainersMap {
	std::map<int, ContainersSet<Lbl>> map;
public:
	using iterator = typename std::map<int, ContainersSet<Lbl>>::iterator;

	friend auto begin(const VertexContainersMap& m) { return m.map.begin(); }
	friend auto end(const VertexContainersMap& m) { return m.map.end(); }

	//vertex has labels(no empty)
	bool has_undominated_labels() const {
		if (map.empty()) { return false; }
		for (const auto& vs : map) {
			if (!vs.second.empty()) {
				return true;
			}
		}
		return false;
	}

	bool is_exist_label() const {
		for (const auto& vs : map) {
			if (!vs.second.empty()) {
				return true;
			}
		}
		return false;
	}

	//the first un-dominated container in one of set
	auto first_with_undominated_container() const {
		for (auto it = map.begin(); it != map.end(); ++it) {
			if (std::any_of(it->second.begin(), it->second.end(),
				[](const auto& container) -> bool {return !container.isDominated; })) {
				return it;
			}
		}
		return map.end();
	}

	//get the first no-empty label ContainersSet
	auto first_containerSet() const {
		for (auto it = map.begin(); it != map.end(); ++it) {
			if (!it->second.empty()) {
				return it;
			}
		}
		return map.end();
	}

	const iterator begin() const { return map.begin(); }
	const iterator end() const { return map.end(); }
	iterator begin() { return map.begin(); }
	iterator end() { return map.end(); }
	auto find(const int& index) const { return map.find(index); }
	auto erase(const int& index) { return map.erase(index); }
	auto empty() const { return map.empty(); }
	auto size() const { return map.size(); }
	const ContainersSet<Lbl>& at(const int& index) const { return map.at(index); }
	ContainersSet<Lbl>& at(const int& index) { return map.at(index); }
	ContainersSet<Lbl>& operator[](const int& index) { return map[index]; }
};


template<typename Lbl>
bool operator==(const LblContainer<Lbl>& c1, const LblContainer<Lbl>& c2) {
	if (c1.preIndex != c2.preIndex) { return false; }
	if (c1.preContainer != c2.preContainer) { return false; }
	return c1.label == c2.label;
}


template<typename Lbl, typename LblExt>
class LabellingAlgorithm {
public:
	LabellingAlgorithm() {}
	std::vector<Column> solve(int start_v, int end_v, Lbl start_label, LblExt extension, Problem* problem) const;
	std::vector<Column> solveR(int start_v, int end_v, Lbl start_label, LblExt extension, Problem* problem) const;
	std::vector<Column> solveS(int start_v, int end_v, Lbl start_label, LblExt extension, Problem* problem) const;
};

template<typename Lbl, typename LblExt>
std::vector <Column>
LabellingAlgorithm<Lbl, LblExt>::solve(int start_v, int end_v, Lbl start_label, LblExt extension, Problem* problem) const {
	VertexContainersMap<Lbl> undominated;
	VertexContainersMap<Lbl> unprocessed;

	int fathomNum = 0;

	// In the beginning we only have the starting label, as an unprocessed label at the starting vertex
	unprocessed[start_v] = { LblContainer<Lbl>(start_label) };
	//undominated[start_v] = ContainersSet<Lbl>();
	//for (auto it = start_label.compatibleVertex.begin(); it != start_label.compatibleVertex.end(); it++) {
	//	unprocessed[*it] = ContainersSet<Lbl>();
	//	undominated[*it] = ContainersSet<Lbl>();
	//}

	// While there are unprocessed labels...
	while (unprocessed.has_undominated_labels()) {
		//std::cout << unprocessed.size() << std::endl;

		auto any_set_it = unprocessed.first_with_undominated_container();
		int cur_vertex = any_set_it->first;
		const ContainersSet<Lbl>& containers_at_cur_vertex = any_set_it->second;

		// Get the first unprocessed labelContainer at the selected vertex
		auto any_cnt_it = containers_at_cur_vertex.first_undominated_container();

		// Make a copy of the chosen label, need to modify the preIndex and preContainer in extension
		auto cur_container = LblContainer<Lbl>(*any_cnt_it);

		typename ContainersSet<Lbl>::iterator cur_inserted_it;
		bool cur_inserted;

		// Insert the current label in undominated
		if (undominated.find(cur_vertex) == undominated.end()) { undominated[cur_vertex] = ContainersSet<Lbl>(); }
		std::tie(cur_inserted_it, cur_inserted) = undominated.at(cur_vertex).insert(cur_container);

		unprocessed.at(cur_vertex).erase(any_cnt_it);
		if (unprocessed.at(cur_vertex).empty()) {
			unprocessed.erase(cur_vertex);
		}

		//std::vector<int> connectIndexes = cur_container.label.compatibleVertex; 
		std::vector<int> connectIndexes = cur_container.label.extensionVertex;
		//extension for each compatible vertex
		for (auto oe = connectIndexes.begin(); oe != connectIndexes.end(); ++oe) {
			int dest_vertex = *oe;

			// Call to the extension function
			auto new_label = extension(cur_container.label, dest_vertex, end_v);

			// Extension didn't succeed: skip the rest
			if (new_label == std::nullopt) { continue; }

			//fathom
			double temp_ = getLabelLowerBoundDy(new_label.value(), extension);
			if (new_label.value().cost + temp_ > 0) {
				//auto tempContainer = LblContainer<Lbl>(new_label.value(), &(*cur_inserted_it), cur_vertex);
				fathomNum++;
				continue;
			}
			//if (new_label.value().cost > 0) {
			//	double temp_ = getLabelLowerBoundDy(new_label.value(), extension);
			//	if (new_label.value().cost + temp_ > 0) {
			//		continue;
			//	}
			//}

			// Extension succeeded! Create a container for the new label
			auto new_container = LblContainer<Lbl>(new_label.value(), &(*cur_inserted_it), cur_vertex);
			bool new_container_dominated = false;

			// If there are unprocessed labels at the destination vertex,
			// if any of them dominates the new label, then discard the new label;
			// if the new label dominates any of them, then discard them.
			if (unprocessed.find(dest_vertex) != unprocessed.end()) {
				auto dest_unp_cnt_it = unprocessed.at(dest_vertex).begin();

				while (dest_unp_cnt_it != unprocessed.at(dest_vertex).end()) {
					const LblContainer<Lbl>& dest_contaienr = *dest_unp_cnt_it;

					if (new_container.label < dest_contaienr.label && !(dest_contaienr.label < new_container.label)) {
						unprocessed.at(dest_vertex).erase(dest_unp_cnt_it++);
					}
					else if (dest_contaienr.label < new_container.label && !(new_container.label < dest_contaienr.label)) {
						new_container_dominated = true;
						break;
					}
					else {
						++dest_unp_cnt_it;
					}
				}

				if (unprocessed.at(dest_vertex).empty()) { unprocessed.erase(dest_vertex); }
				if (new_container_dominated) { continue; }
			}

			// If there are undominated labels at the destination vertex,
			// if any of them dominates the new label, then discard the new label;
			// if the new label dominates any of them, then discard them.
			if (undominated.find(dest_vertex) != undominated.end()) {
				auto dest_und_cnt_it = undominated.at(dest_vertex).begin();
				//for each undominated LblContainer at the dest vertex
				while (dest_und_cnt_it != undominated.at(dest_vertex).end()) {
					const LblContainer<Lbl>& dest_container = *dest_und_cnt_it;
					if (!dest_container.isDominated && new_container.label < dest_container.label &&
						!(dest_container.label < new_container.label)) {
						undominated.at(dest_vertex).mark_dominated(dest_und_cnt_it++); //erase the dest_und_cnt_it and the ++
					}
					else if (!dest_container.isDominated && dest_container.label < new_container.label &&
						!(new_container.label < dest_container.label)) {
						new_container_dominated = true;
						break;
					}
					else {
						++dest_und_cnt_it;
					}
				}
				if (undominated.at(dest_vertex).empty()) { undominated.erase(dest_vertex); }
				if (new_container_dominated) { continue; }
			}

			// if we arrived up to here, it means that the new extension label is un-dominated
			// by any existing label at the destination vertex, so we can place it in
			// the set of unprocessed labels at the destination vertex.
			typename ContainersSet<Lbl>::iterator new_insert_it;
			bool new_inserted;
			if (unprocessed.find(dest_vertex) == unprocessed.end()) { unprocessed[dest_vertex] = ContainersSet<Lbl>(); }
			std::tie(new_insert_it, new_inserted) = unprocessed.at(dest_vertex).insert(new_container);
		}
	}

	if (undominated.find(end_v) == undominated.end()) { return std::vector<Column>(); }


	// We now get the undominated labels at the end vertex
	const ContainersSet<Lbl>& pareto_optimal_containers = undominated.at(end_v);
	std::vector<Column> pareto_optimal_solutions;
	pareto_optimal_solutions.reserve(pareto_optimal_containers.size());

	// And, for each of them, we reconstruct the corresponding optimal path
	for (const auto& oc : pareto_optimal_containers) {
		if (oc.isDominated) {
			continue;
		}
		if (oc.label.cost >= 0) {
			continue;
		}

		Column column(problem);
		std::vector<int> coils;
		const LblContainer<Lbl>* current = &oc;
		coils.emplace_back(current->label.curIndex);

		while (current->preContainer != nullptr) {
			coils.emplace_back(current->preIndex.value());
			current = current->preContainer;
		}

		for (auto rit = coils.rbegin() + 1; rit != coils.rend() - 1; ++rit) {
			//std::cout << *rit << " ";
			column.addVertex(*rit);
		}
		column.calculateCost();
		column.setRC(oc.label.cost);

		//std::cout << "cost label = " << oc.label.cost << std::endl;
		//std::cout << "cost column =	" << column.getCost() << std::endl;
		//std::cout << "-------------------" << std::endl;
		//std::cout << std::endl;

		pareto_optimal_solutions.emplace_back(column);

		//column.printRoute();
	}

	//double mass = 0;
	//int index = 0;
	//for (int i = 0; i < pareto_optimal_solutions.size(); i++) {
	//	if (pareto_optimal_solutions[i].getDemand() > mass) {
	//		mass = pareto_optimal_solutions[i].getDemand();
	//		index = i;
	//	}
	//}

	std::cout << "end and fathom label : " << fathomNum << std::endl;

	return pareto_optimal_solutions;
}


template<typename Lbl, typename LblExt>
std::vector <Column>
LabellingAlgorithm<Lbl, LblExt>::solveR(int start_v, int end_v, Lbl start_label, LblExt extension, Problem* problem) const {
	//auto start = std::chrono::high_resolution_clock::now();

	VertexContainersMap<Lbl> undominated;
	VertexContainersMap<Lbl> unprocessed;

	int fathomNum = 0;

	// In the beginning we only have the starting label, as an unprocessed label at the starting vertex
	unprocessed[start_v] = { LblContainer<Lbl>(start_label) };

	// While there are unprocessed labels...
	while (unprocessed.has_undominated_labels()) {

		auto any_set_it = unprocessed.first_with_undominated_container();
		int cur_vertex = any_set_it->first;
		const ContainersSet<Lbl>& containers_at_cur_vertex = any_set_it->second;

		// Get the first unprocessed labelContainer at the selected vertex
		auto any_cnt_it = containers_at_cur_vertex.first_undominated_container();

		//debug
		//if (any_cnt_it->label.curIndex == 50) {
		//	std::cout << "debug: current coil = 50" << std::endl;
		//}

		// Make a copy of the chosen label, need to modify the preIndex and preContainer in extension
		auto cur_container = LblContainer<Lbl>(*any_cnt_it);

		typename ContainersSet<Lbl>::iterator cur_inserted_it;
		bool cur_inserted;

		// Insert the current label in undominated
		if (undominated.find(cur_vertex) == undominated.end()) { undominated[cur_vertex] = ContainersSet<Lbl>(); }
		std::tie(cur_inserted_it, cur_inserted) = undominated.at(cur_vertex).insert(cur_container);

		//Remove the label from unprocessed
		unprocessed.at(cur_vertex).erase(any_cnt_it);
		if (unprocessed.at(cur_vertex).empty()) {
			unprocessed.erase(cur_vertex);
		}

		std::vector<int> connectIndexes = cur_container.label.extensionVertex;		
		//extension for each compatible vertex
		for (auto oe = connectIndexes.begin(); oe != connectIndexes.end(); ++oe) {
			int dest_vertex = *oe;

			//debug
			//if (dest_vertex == 70) {
			//	std::cout << "debug: target coil = 70" << std::endl;
			//}

			// Call to the extension function
			auto new_label = extension(cur_container.label, dest_vertex, end_v);

			// Extension didn't succeed: skip the rest
			if (new_label == std::nullopt) { continue; }

			//fathom
			double temp_ = getLabelLowerBoundDy(new_label.value(), extension);
			if (new_label.value().cost + temp_ > 0) {
				fathomNum++;
				continue;
			}

			// Extension succeeded! Create a container for the new label
			auto new_container = LblContainer<Lbl>(new_label.value(), &(*cur_inserted_it), cur_vertex);
			bool new_container_dominated = false;

			// If there are unprocessed labels at the destination vertex,
			// if any of them dominates the new label, then discard the new label;
			// if the new label dominates any of them, then discard them.
			if (unprocessed.find(dest_vertex) != unprocessed.end()) {
				auto dest_unp_cnt_it = unprocessed.at(dest_vertex).begin();

				while (dest_unp_cnt_it != unprocessed.at(dest_vertex).end()) {
					const LblContainer<Lbl>& dest_contaienr = *dest_unp_cnt_it;

					if (new_container.label < dest_contaienr.label && !(dest_contaienr.label < new_container.label)) {
						unprocessed.at(dest_vertex).erase(dest_unp_cnt_it++);
					}
					else if (dest_contaienr.label < new_container.label && !(new_container.label < dest_contaienr.label)) {
						//debug
						//std::string new_label_str = get_coils_path(new_container);
						//std::string dest_label_str = get_coils_path(dest_contaienr);
						//std::cout << dest_label_str << " dominates " << new_label_str 
						//	<< " with cost " << new_container.label.cost << " < " << dest_contaienr.label.cost
						//	<< std::endl;

						new_container_dominated = true;
						break;
					}
					else {
						++dest_unp_cnt_it;
					}
				}

				if (unprocessed.at(dest_vertex).empty()) { unprocessed.erase(dest_vertex); }
				if (new_container_dominated) { continue; }
			}
			if (undominated.find(dest_vertex) != undominated.end()) {
				auto dest_und_cnt_it = undominated.at(dest_vertex).begin();
				//for each undominated LblContainer at the dest vertex
				while (dest_und_cnt_it != undominated.at(dest_vertex).end()) {
					const LblContainer<Lbl>& dest_container = *dest_und_cnt_it;
					if (!dest_container.isDominated && new_container.label < dest_container.label &&
						!(dest_container.label < new_container.label)) {
						undominated.at(dest_vertex).mark_dominated(dest_und_cnt_it++); //erase the dest_und_cnt_it and the ++
					}
					else if (!dest_container.isDominated && dest_container.label < new_container.label &&
						!(new_container.label < dest_container.label)) {
						//debug
						//std::string new_label_str = get_coils_path(new_container);
						//std::string dest_label_str = get_coils_path(dest_container);
						//std::cout << dest_label_str << " dominates " << new_label_str
						//	<< " with cost " << new_container.label.cost << " < " << dest_container.label.cost
						//	<< std::endl;

						new_container_dominated = true;
						break;
					}
					else {
						++dest_und_cnt_it;
					}
				}
				if (undominated.at(dest_vertex).empty()) { undominated.erase(dest_vertex); }
				if (new_container_dominated) { continue; }
			}
			// if we arrived up to here, it means that the new extension label is un-dominated
			// by any existing label at the destination vertex, so we can place it in
			// the set of unprocessed labels at the destination vertex.
			typename ContainersSet<Lbl>::iterator new_insert_it;
			bool new_inserted;
			if (unprocessed.find(dest_vertex) == unprocessed.end()) { unprocessed[dest_vertex] = ContainersSet<Lbl>(); }
			std::tie(new_insert_it, new_inserted) = unprocessed.at(dest_vertex).insert(new_container);
		}
	}

	//auto end = std::chrono::high_resolution_clock::now();
	//std::cout << "label setting time = " << std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count() << std::endl;

	if (undominated.find(end_v) == undominated.end()) { return std::vector<Column>(); }

	// We now get the undominated labels at the end vertex
	const ContainersSet<Lbl>& pareto_optimal_containers = undominated.at(end_v);
	std::vector<Column> pareto_optimal_solutions;
	pareto_optimal_solutions.reserve(pareto_optimal_containers.size());

	// And, for each of them, we reconstruct the corresponding optimal path
	for (const auto& oc : pareto_optimal_containers) {
		if (oc.isDominated) {
			continue;
		}
		if (oc.label.cost >= 0) {
			continue;
		}

		Column column(problem);
		std::vector<int> coils;
		const LblContainer<Lbl>* current = &oc;
		coils.emplace_back(current->label.curIndex);

		while (current->preContainer != nullptr) {
			coils.emplace_back(current->preIndex.value());
			current = current->preContainer;
		}

		for (auto rit = coils.rbegin() + 1; rit != coils.rend() - 1; ++rit) {
			//std::cout << *rit << " ";
			column.addVertex(*rit);
		}
		column.calculateCost();
		column.setRC(oc.label.cost);

		pareto_optimal_solutions.emplace_back(column);

		//column.printRoute();
	}

	// std::cout << "end and fathom label : " << fathomNum << std::endl;

	return pareto_optimal_solutions;
}

template<typename Lbl, typename LblExt>
std::vector <Column>
LabellingAlgorithm<Lbl, LblExt>::solveS(int start_v, int end_v, Lbl start_label, LblExt extension, Problem* problem) const {
	auto start = std::chrono::high_resolution_clock::now();

	VertexContainersMap<Lbl> unprocessed;
	VertexContainersMap<Lbl> undominated;

	int fathomNum = 0;

	//initialize
	unprocessed[start_v] = { LblContainer<Lbl>(start_label) };
	undominated[start_v] = ContainersSet<Lbl>();
	for (const auto& item : start_label.compatibleVertex) {
		unprocessed[item] = ContainersSet<Lbl>();
		undominated[item] = ContainersSet<Lbl>();
	}

	do {
		//get the first ContainerSet and label(container)
		auto any_set_it = unprocessed.first_containerSet();
		int cur_vertex = any_set_it->first;
		const ContainersSet<Lbl>& containers_at_cur_vertex = any_set_it->second;

		auto any_cnt_it = containers_at_cur_vertex.first_label_container();
		auto cur_container = LblContainer<Lbl>(*any_cnt_it);

		//pop the first label
		unprocessed.at(cur_vertex).erase(any_cnt_it);

		typename ContainersSet<Lbl>::iterator cur_inserted_it;
		bool cur_inserted;

		// Insert the current label in undominated
		std::tie(cur_inserted_it, cur_inserted) = undominated.at(cur_vertex).insert(cur_container);


		bool is_dominated = false;
		auto cur_unp_cnt_it = unprocessed.at(cur_vertex).begin();

		while (cur_unp_cnt_it != unprocessed[cur_vertex].end())
		{
			const LblContainer<Lbl>& other_label_contanier = *cur_unp_cnt_it;
			if (other_label_contanier.label < cur_container.label) {
				is_dominated = true;
				break;
			}
			cur_unp_cnt_it++;
		}
		//no other label can dominate this label
		if (!is_dominated) {
			//extension with feasible vertices
			std::vector<int> connectIndexes = cur_container.label.extensionVertex;
			for (const auto& item : connectIndexes) {
				auto new_label = extension(cur_container.label, item, end_v);
				if (new_label == std::nullopt) {
					continue;
				}
				double temp_ = getLabelLowerBoundDy(new_label.value(), extension);
				if (new_label.value().cost + temp_ > 0) {
					continue;
				}
				auto new_container = LblContainer<Lbl>(new_label.value(), &(*cur_inserted_it), cur_vertex);

				//put the new label container into the unprocessed
				unprocessed.at(item).insert(new_container);
			}
		}

		// if no other label can dominate this label in end vertex, save this label
		if (cur_vertex == end_v) {
			//...
		}
	} while (unprocessed.is_exist_label());

	auto end = std::chrono::high_resolution_clock::now();
	std::cout << "label setting time = " << std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count() << std::endl;

	return std::vector<Column>();
}
