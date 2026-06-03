#include "Labelling.h"

bool operator==(const Label& lhs, const Label& rhs)
{
	return (lhs.demand == rhs.demand) &&
		std::abs(lhs.cost - rhs.cost) < Label::EPS;
}

bool operator!=(const Label& lhs, const Label& rhs)
{
	return !(lhs == rhs);
}

bool operator<=(const Label& lhs, const Label& rhs)
{
	//dominant rule only consider the demand and cost
	if (rhs.cost < lhs.cost - Label::EPS) { return false; }
	if (rhs.demand < lhs.demand) { return false; }
	return true;
}

bool operator<(const Label& lhs, const Label& rhs)
{
	return lhs <= rhs && lhs != rhs;
}

std::ostream& operator<<(std::ostream& out, const Label& l)
{
	out << "cost : " << l.cost << ", demand : " << l.demand << std::endl;
	return out;
}

bool operator==(const ElementaryLabel& lhs, const ElementaryLabel& rhs)
{
	return (lhs.demand == rhs.demand &&
		std::abs(lhs.cost - rhs.cost) < Label::EPS
		);
}

bool operator!=(const ElementaryLabel& lhs, const ElementaryLabel& rhs)
{
	return !(lhs == rhs);
}

bool operator<=(const ElementaryLabel& lhs, const ElementaryLabel& rhs)
{
	// the compatible vector has the can connect coils
	if (rhs.cost < lhs.cost - Label::EPS) { return false; }
	if (rhs.demand < lhs.demand) { return false; }
	if (std::any_of(
		rhs.compatibleVertex.begin(),
		rhs.compatibleVertex.end(),
		[&](const auto& index) {
			return std::find(lhs.compatibleVertex.begin(), lhs.compatibleVertex.end(), index) == lhs.compatibleVertex.end();
		}
	)) {
		return false;
	}
	return true;
}

bool operator<(const ElementaryLabel& lhs, const ElementaryLabel& rhs)
{
	//return lhs <= rhs && lhs != rhs;
	return lhs <= rhs;
}

std::ostream& operator<<(std::ostream& out, const ElementaryLabel& l)
{
	out << "cost : " << l.cost << ", demand : " << l.demand << " " << std::endl;
}


// SR
bool operator==(const ElementaryLabelSR& lhs, const ElementaryLabelSR& rhs)
{
	return (lhs.demand == rhs.demand &&
		lhs.compatibleVertex == rhs.compatibleVertex &&
		std::abs(lhs.cost - rhs.cost) < Label::EPS
		);
}

bool operator!=(const ElementaryLabelSR& lhs, const ElementaryLabelSR& rhs)
{
	return !(lhs == rhs);
}

bool operator<=(const ElementaryLabelSR& lhs, const ElementaryLabelSR& rhs)
{
	double lhs_cost = lhs.cost;
	for (int i = 0; i < lhs.multiS.size(); i++) {
		if (lhs.numVertexinS[i] > rhs.numVertexinS[i]) {
			lhs_cost -= lhs.rmp->getDualVariable(lhs.consNum + i);
		}
	}

	if (rhs.cost < lhs_cost - Label::EPS) { return false; }

	if (rhs.demand < lhs.demand) { return false; }
	if (std::any_of(
		rhs.compatibleVertex.begin(),
		rhs.compatibleVertex.end(),
		[&](const auto& index) {
			return std::find(lhs.compatibleVertex.begin(), lhs.compatibleVertex.end(), index) == lhs.compatibleVertex.end();
		}
	)) {
		return false;
	}
	return true;
}

bool operator<(const ElementaryLabelSR& lhs, const ElementaryLabelSR& rhs)
{
	//double lhs_cost = lhs.cost;
	//for (int i = 0; i < lhs.multiS.size(); i++) {
	//	if (lhs.numVertexinS[i] > rhs.numVertexinS[i]) {
	//		lhs_cost -= lhs.rmp->getDualVariable(lhs.consNum + i);
	//	}
	//}

	//if (rhs.cost < lhs_cost - Label::EPS) { return false; }

	//if (rhs.demand < lhs.demand) { return false; }
	//if (std::any_of(
	//	rhs.compatibleVertex.begin(),
	//	rhs.compatibleVertex.end(),
	//	[&](const auto& index) {
	//		return std::find(lhs.compatibleVertex.begin(), lhs.compatibleVertex.end(), index) == lhs.compatibleVertex.end();
	//	}
	//)) {
	//	return false;
	//}
	//return true;
	return lhs <= rhs && lhs != rhs;
}

bool operator==(const LabelSR& lhs, const LabelSR& rhs)
{
	//ng set
	return (lhs.demand == rhs.demand &&
		lhs.ng_memory == rhs.ng_memory &&
		std::abs(lhs.cost - rhs.cost) < Label::EPS
		);
}

bool operator!=(const LabelSR& lhs, const LabelSR& rhs)
{
	return !(lhs == rhs);
}

bool operator<=(const LabelSR& lhs, const LabelSR& rhs)
{
	double lhs_cost = lhs.cost;
	for (int i = 0; i < lhs.multiS.size(); i++) {
		if (lhs.numVertexinS[i] > rhs.numVertexinS[i]) {
			lhs_cost -= lhs.rmp->getDualVariable(lhs.consNum + i);
		}
	}

	if (rhs.cost < lhs_cost - Label::EPS) { return false; }

	if (rhs.demand < lhs.demand) { return false; }

	if ((lhs.ng_memory & ~rhs.ng_memory).any()) {
		return false;
	}

	return true;
}

//bool operator<(const LabelSR& lhs, const LabelSR& rhs)
//{
//	double lhs_cost = lhs.cost;
//	for (int i = 0; i < lhs.multiS.size(); i++) {
//		if (lhs.numVertexinS[i] > rhs.numVertexinS[i]) {
//			lhs_cost -= lhs.rmp->getDualVariable(lhs.consNum + i);
//		}
//	}
//
//	if (rhs.cost < lhs_cost - Label::EPS) { return false; }
//
//	if (rhs.demand < lhs.demand) { return false; }
//	
//	if ((lhs.ng_memory & ~rhs.ng_memory).any()) {
//		return false;
//	}
//
//	return true;
//}

//bool operator<(const ElementaryLabelSR& lhs, const ElementaryLabelSR& rhs)
//{
//	double lhs_cost = lhs.cost;
//
//	std::vector<int> index_subset;
//
//	std::set<int> set_lhs_compatible(lhs.compatibleVertex.begin(), lhs.compatibleVertex.end());
//
//	std::vector<int> lhs_visited = lhs.allIndex;
//	lhs_visited.erase(std::remove_if(lhs_visited.begin(), lhs_visited.end(),
//		[&set_lhs_compatible](int a) { return set_lhs_compatible.find(a) != set_lhs_compatible.end(); }),
//		lhs_visited.end());
//	std::vector<int> N;
//	N.insert(N.begin(), lhs_visited.begin(), lhs_visited.end());
//	N.insert(N.begin(), rhs.compatibleVertex.begin(), rhs.compatibleVertex.end());
//
//	for (int i = 0; i < lhs.multiS.size(); i++) {
//		if (lhs.numVertexinS[i] > rhs.numVertexinS[i]) {
//			//lhs_cost -= lhs.rmp->getDualVariable(lhs.consNum + i);
//			index_subset.emplace_back(i);
//		}
//	}
//	for (const auto& item : index_subset) {
//		int count = std::count_if(N.begin(), N.end(), [&](int b) {
//			return std::find(lhs.multiS[item].begin(), lhs.multiS[item].end(), b) != lhs.multiS[item].end();
//			});
//		if (count >= 2) {
//			lhs_cost -= lhs.rmp->getDualVariable(lhs.consNum + item);
//		}
//	}
//
//	if (rhs.cost < lhs_cost - Label::EPS) { return false; }
//
//	if (rhs.demand < lhs.demand) { return false; }
//	if (std::any_of(
//		rhs.compatibleVertex.begin(),
//		rhs.compatibleVertex.end(),
//		[&](const auto& index) {
//			return std::find(lhs.compatibleVertex.begin(), lhs.compatibleVertex.end(), index) == lhs.compatibleVertex.end();
//		}
//	)) {
//		return false;
//	}
//	return true;
//}


double getLabelLowerBound(const ElementaryLabel& label, const LabelExtender& extension) {
	double remainWeight = extension.parameters->planUp - label.demand;
	int numCoils = static_cast<int>(label.compatibleVertex.size() - 1);

	IloEnv env;
	IloModel model(env);
	double lowerBound = std::numeric_limits<double>::max();


	try {
		//define the variable
		IloNumVarArray Y(env, numCoils);
		for (int i = 0; i < numCoils; i++) {
			Y[i] = IloNumVar(env, 0, 1, ILOINT);
		}

		//constraint
		IloExpr cons(env);
		for (int i = 0; i < numCoils; i++) {
			int indexRealCoil = label.compatibleVertex[i];
			cons += Y[i] * extension.problem->getWeight(indexRealCoil);
		}
		model.add(cons <= remainWeight);
		cons.end();

		//obj
		IloExpr obj(env);
		for (int l = 0; l < extension.parameters->L; l++) {
			for (int i = 0; i < numCoils; i++) {
				int indexRealCoil = label.compatibleVertex[i];
				if (extension.problem->getFlow(indexRealCoil) == extension.parameters->target_flow[l]) {
					int indexFlow = extension.parameters->N + l;
					obj -= Y[i] * extension.problem->getWeight(indexRealCoil) * extension.rmp->getDualVariable(indexFlow);
				}
			}
		}
		obj -= extension.rmp->getDualVariable(extension.parameters->N + extension.parameters->L + 0);
		obj -= extension.rmp->getDualVariable(extension.parameters->N + extension.parameters->L + 1);
		model.add(IloMinimize(env, obj));
		obj.end();

		IloCplex cplex(model);
		//cplex.setParam(IloCplex::Param::MIP::Display, 0);
		if (!cplex.solve()) {
			env.error() << "Failed!" << std::endl;
		}
		env.out() << "solution status : " << cplex.getStatus() << std::endl;

		//print result 

		for (int i = 0; i < numCoils; i++) {
			if (cplex.getValue(Y[i]) > 0) {
				std::cout << label.compatibleVertex[i] << " weight : " << extension.problem->getWeight(label.compatibleVertex[i])
					<< " flow : " << extension.problem->getFlow(label.compatibleVertex[i])
					<< std::endl;
			}
		}
		lowerBound = cplex.getObjValue();
	}
	catch (IloException& ex) {
		std::cerr << "Error: " << ex << std::endl;
	}
	catch (...) {
		std::cerr << "unknown error" << std::endl;
	}
	env.end();

	//print the coefficient of variable
	for (int i = 0; i < numCoils; i++) {
		int index_ = label.compatibleVertex[i];
		double coef = 0;
		for (int l = 0; l < extension.parameters->L; l++) {
			if (extension.problem->getFlow(index_) == extension.parameters->target_flow[l]) {
				coef += extension.problem->getWeight(index_) * extension.rmp->getDualVariable(extension.parameters->N + l);
			}
		}
		std::cout << "index_origin: " << i << " index: " << index_ << " coef: " << coef << std::endl;
	}
	return lowerBound;
}


double getLabelLowerBoundDy(const ElementaryLabel& label, const LabelExtender& extension)
{
	if (label.curIndex == extension.parameters->N + 1) {
		return 0.0;
	}

	int n = static_cast<int>(label.compatibleVertex.size()) - 1;
	int v = static_cast<int>(std::ceil(extension.parameters->planUp - label.demand));

	//get array of capacity(weight, convert to int) and worth(value)
	int capacity[100]{};
	double worth[100]{};
	double dp[500]{};

	//initialize the coil's capacity and worth(dual variable)
	for (int i = 0; i < n; i++) {
		capacity[i + 1] = std::floor(extension.problem->getWeight(label.compatibleVertex[i]));
		if (extension.problem->getFlow(label.compatibleVertex[i]) == "C512") {
			worth[i + 1] += extension.problem->getWeight(label.compatibleVertex[i]) * extension.rmp->getDualVariable(extension.parameters->N);
		}
		else {
			worth[i + 1] += extension.problem->getWeight(label.compatibleVertex[i]) * extension.rmp->getDualVariable(extension.parameters->N + 1);
		}
	}

	//dynamic programe - 2
	for (int i = 1; i <= n; i++) {
		for (int j = v; j >= 0; j--) {
			if (j - capacity[i] >= 0) {
				dp[j] = std::max(dp[j - capacity[i]] + worth[i], dp[j]);
			}
		}
	}
	//std::cout << "result dp2 = " << dp[v] << std::endl;
	double lowerBound = -dp[v];
	lowerBound -= extension.rmp->getDualVariable(extension.parameters->N + extension.parameters->L + 0);
	lowerBound -= extension.rmp->getDualVariable(extension.parameters->N + extension.parameters->L + 1);
	return lowerBound;
}

double getLabelLowerBoundDy(const ElementaryLabelSR& label, const LabelExtender& extension)
{
	if (label.curIndex == extension.parameters->N + 1) {
		return 0.0;
	}

	int n = static_cast<int>(label.compatibleVertex.size()) - 1;
	int v = static_cast<int>(std::ceil(extension.parameters->planUp - label.demand));

	//get array of capacity(weight, convert to int) and worth(value)
	int capacity[100]{};
	double worth[100]{};
	double dp[500]{};

	//initialize the coil's capacity and worth(dual variable)
	for (int i = 0; i < n; i++) {
		capacity[i + 1] = std::floor(extension.problem->getWeight(label.compatibleVertex[i]));
		if (extension.problem->getFlow(label.compatibleVertex[i]) == "C512") {
			worth[i + 1] += extension.problem->getWeight(label.compatibleVertex[i]) * extension.rmp->getDualVariable(extension.parameters->N);
		}
		else {
			worth[i + 1] += extension.problem->getWeight(label.compatibleVertex[i]) * extension.rmp->getDualVariable(extension.parameters->N + 1);
		}
	}

	//dynamic programe - 2
	for (int i = 1; i <= n; i++) {
		for (int j = v; j >= 0; j--) {
			if (j - capacity[i] >= 0) {
				dp[j] = std::max(dp[j - capacity[i]] + worth[i], dp[j]);
			}
		}
	}
	//std::cout << "result dp2 = " << dp[v] << std::endl;
	double lowerBound = -dp[v];
	lowerBound -= extension.rmp->getDualVariable(extension.parameters->N + extension.parameters->L + 0);
	lowerBound -= extension.rmp->getDualVariable(extension.parameters->N + extension.parameters->L + 1);
	return lowerBound;
}

double getLabelLowerBoundDy(const LabelSR& label, const LabelExtender& extension)
{
	if (label.curIndex == extension.parameters->N + 1) {
		return 0.0;
	}

	int n = static_cast<int>(label.compatibleVertex.size()) - 1;
	int v = static_cast<int>(std::ceil(extension.parameters->planUp - label.demand));

	//get array of capacity(weight, convert to int) and worth(value)
	int capacity[100]{};
	double worth[100]{};
	double dp[500]{};

	//initialize the coil's capacity and worth(dual variable)
	for (int i = 0; i < n; i++) {
		capacity[i + 1] = std::floor(extension.problem->getWeight(label.compatibleVertex[i]));
		if (extension.problem->getFlow(label.compatibleVertex[i]) == "C512") {
			worth[i + 1] += extension.problem->getWeight(label.compatibleVertex[i]) * extension.rmp->getDualVariable(extension.parameters->N);
		}
		else {
			worth[i + 1] += extension.problem->getWeight(label.compatibleVertex[i]) * extension.rmp->getDualVariable(extension.parameters->N + 1);
		}
	}

	//dynamic programe - 2
	for (int i = 1; i <= n; i++) {
		for (int j = v; j >= 0; j--) {
			if (j - capacity[i] >= 0) {
				dp[j] = std::max(dp[j - capacity[i]] + worth[i], dp[j]);
			}
		}
	}
	//std::cout << "result dp2 = " << dp[v] << std::endl;
	double lowerBound = -dp[v];
	lowerBound -= extension.rmp->getDualVariable(extension.parameters->N + extension.parameters->L + 0);
	lowerBound -= extension.rmp->getDualVariable(extension.parameters->N + extension.parameters->L + 1);
	return lowerBound;
}


std::optional<ElementaryLabel> LabelExtender::operator()(const ElementaryLabel& label, int targetIndex, int end) const
{
	//if (std::find(erasedVertex.begin(), erasedVertex.end(), targetIndex) != erasedVertex.end()) {
	//	return std::nullopt;
	//}

	if (std::find(label.compatibleVertex.begin(), label.compatibleVertex.end(), targetIndex) == label.compatibleVertex.end()) {
		return std::nullopt;
	}

	//stop when reach to the end node
	if (label.curIndex == end) {
		return std::nullopt;
	}

	if (label.curIndex == 0 && targetIndex == end) {
		return std::nullopt;
	}

	//a copy
	ElementaryLabel newLabel = label;

	//erase the target vertex in the compatible vertex
	newLabel.compatibleVertex.erase(std::remove(newLabel.compatibleVertex.begin(), newLabel.compatibleVertex.end(), targetIndex),
		newLabel.compatibleVertex.end());

	//TODO erase more incompatible coils
	newLabel.extensionVertex = newLabel.compatibleVertex;

	auto it = forbidIndex.find(targetIndex);
	if (it != forbidIndex.end()) {
		std::vector<int> vec = it->second;
		newLabel.extensionVertex.erase(std::remove_if(newLabel.extensionVertex.begin(), newLabel.extensionVertex.end(), [&vec](int elem) {
			return std::find(vec.begin(), vec.end(), elem) != vec.end();
			}), newLabel.extensionVertex.end());
	}

	//the accumulated resource consumption must not be exceed the capacity
	if (newLabel.demand + problem->getWeight(targetIndex) > parameters->planUp) {
		return std::nullopt;
	}
	newLabel.demand += problem->getWeight(targetIndex);

	//calculate the RC need consider the terminal node and begin node (c_ij = 0)
	if (targetIndex == end) {
		if (newLabel.demand < parameters->planLower || newLabel.demand > parameters->planUp) {
			return std::nullopt;
		}
		newLabel.cost -= rmp->getDualVariable(parameters->N + parameters->L) - rmp->getDualVariable(parameters->N + parameters->L + 1);
		newLabel.cost += parameters->alpha2 * parameters->setUpCost;
		newLabel.curIndex = targetIndex;
		return newLabel;
	}

	if (label.curIndex == 0) {
		newLabel.cost -= rmp->getDualVariable(targetIndex - 1);
		newLabel.cost = newLabel.cost - problem->getWeight(targetIndex) * (problem->getFlow(targetIndex) == "C512" ? 1 : 0) * rmp->getDualVariable(problem->getVertices() + 0);
		newLabel.cost = newLabel.cost - problem->getWeight(targetIndex) * (problem->getFlow(targetIndex) == "C008" ? 1 : 0) * rmp->getDualVariable(problem->getVertices() + 1);
		newLabel.curIndex = targetIndex;
		return newLabel;
	}

	//std::cout << " i and j = " << newLabel.curIndex << " " << targetIndex << std::endl;
	//std::cout << "cost i and j = " << problem->getCost(newLabel.curIndex - 1, targetIndex - 1) << std::endl;
	//std::cout << "cost_" << newLabel.curIndex << targetIndex << " : " << problem->getCost(newLabel.curIndex - 1, targetIndex - 1) << std::endl;
	newLabel.cost = newLabel.cost + parameters->alpha1 * problem->getCost(newLabel.curIndex - 1, targetIndex - 1) - rmp->getDualVariable(targetIndex - 1);

	newLabel.cost = newLabel.cost - problem->getWeight(targetIndex) * (problem->getFlow(targetIndex) == "C512" ? 1 : 0) * rmp->getDualVariable(problem->getVertices() + 0);
	newLabel.cost = newLabel.cost - problem->getWeight(targetIndex) * (problem->getFlow(targetIndex) == "C008" ? 1 : 0) * rmp->getDualVariable(problem->getVertices() + 1);

	newLabel.curIndex = targetIndex;

	return newLabel;
}

std::optional<ElementaryLabelSR> LabelExtender::operator()(const ElementaryLabelSR& label, int targetIndex, int end) const
{
	if (std::find(label.compatibleVertex.begin(), label.compatibleVertex.end(), targetIndex) == label.compatibleVertex.end()) {
		return std::nullopt;
	}

	if (label.curIndex == end) {
		return std::nullopt;
	}

	if (label.curIndex == 0 && targetIndex == end) {
		return std::nullopt;
	}

	//a copy
	ElementaryLabelSR newLabel = label;

	//erase the target vertex in the compatible vertex
	newLabel.compatibleVertex.erase(std::remove(newLabel.compatibleVertex.begin(), newLabel.compatibleVertex.end(), targetIndex),
		newLabel.compatibleVertex.end());

	newLabel.extensionVertex = newLabel.compatibleVertex;

	auto it = forbidIndex.find(targetIndex);
	if (it != forbidIndex.end()) {
		std::vector<int> vec = it->second;
		newLabel.extensionVertex.erase(std::remove_if(newLabel.extensionVertex.begin(), newLabel.extensionVertex.end(), [&vec](int elem) {
			return std::find(vec.begin(), vec.end(), elem) != vec.end();
			}), newLabel.extensionVertex.end());
	}

	//the accumulated resource consumption must not be exceed the capacity
	if (newLabel.demand + problem->getWeight(targetIndex) > parameters->planUp) {
		return std::nullopt;
	}
	newLabel.demand += problem->getWeight(targetIndex);

	//calculate the number of coils involved in the cut
	for (int i = 0; i < newLabel.multiS.size(); i++) {
		if (std::find(newLabel.multiS[i].begin(), newLabel.multiS[i].end(), targetIndex) != newLabel.multiS[i].end()) {
			newLabel.numVertexinS[i] = newLabel.numVertexinS[i] + 1;
		}
	}
	//calculate the penalty of SR
	for (int i = 0; i < newLabel.multiS.size(); i++) {
		if (label.numVertexinS[i] < 2 && newLabel.numVertexinS[i] == 2) {
			newLabel.cost -= rmp->getDualVariable(parameters->N + parameters->L + 2 + i);
		}
	}


	//calculate the RC need consider the terminal node and begin node (c_ij = 0)
	if (targetIndex == end) {
		//For plan lower
		if (newLabel.demand < parameters->planLower || newLabel.demand > parameters->planUp) {
			return std::nullopt;
		}
		// if the target vertex is sink node, return label directly
		newLabel.cost -= rmp->getDualVariable(parameters->N + parameters->L) - rmp->getDualVariable(parameters->N + parameters->L + 1);
		newLabel.cost += parameters->alpha2 * parameters->setUpCost;
		newLabel.curIndex = targetIndex;
		return newLabel;
	}

	if (label.curIndex == 0) {
		newLabel.cost -= rmp->getDualVariable(targetIndex - 1);
		newLabel.cost = newLabel.cost - problem->getWeight(targetIndex) * (problem->getFlow(targetIndex) == "C512" ? 1 : 0) * rmp->getDualVariable(problem->getVertices() + 0);
		newLabel.cost = newLabel.cost - problem->getWeight(targetIndex) * (problem->getFlow(targetIndex) == "C008" ? 1 : 0) * rmp->getDualVariable(problem->getVertices() + 1);
		newLabel.curIndex = targetIndex;
		return newLabel;
	}
	newLabel.cost = newLabel.cost + parameters->alpha1 * problem->getCost(newLabel.curIndex - 1, targetIndex - 1) - rmp->getDualVariable(targetIndex - 1);
	newLabel.cost = newLabel.cost - problem->getWeight(targetIndex) * (problem->getFlow(targetIndex) == "C512" ? 1 : 0) * rmp->getDualVariable(problem->getVertices() + 0);
	newLabel.cost = newLabel.cost - problem->getWeight(targetIndex) * (problem->getFlow(targetIndex) == "C008" ? 1 : 0) * rmp->getDualVariable(problem->getVertices() + 1);

	newLabel.curIndex = targetIndex;
	return newLabel;
}

std::optional<LabelSR> LabelExtender::operator()(const LabelSR& label, int targetIndex, int end) const
{
	if (label.curIndex == end) {
		return std::nullopt;
	}
	if (label.curIndex == 0 && targetIndex == end) {
		return std::nullopt;
	}
	if (label.ng_memory.test(targetIndex)) {
		//the target coil is in the memory of path
		return std::nullopt;
	}
	LabelSR newLabel = label;
	if (targetIndex != end) {
		//auto trans_it = std::find(all_index.begin(), all_index.end(), targetIndex);
		//int trans_index = std::distance(all_index.begin(), trans_it) + 1;

		auto trans_it = std::find(original_index.begin(), original_index.end(), targetIndex);
		int trans_index = std::distance(original_index.begin(), trans_it);
		
		newLabel.ng_memory = (label.ng_memory & coils_ng_sets[trans_index]);
		newLabel.ng_memory.set(targetIndex);
	}

	newLabel.compatibleVertex.erase(std::remove(newLabel.compatibleVertex.begin(), newLabel.compatibleVertex.end(), targetIndex),
		newLabel.compatibleVertex.end());

	newLabel.extensionVertex = all_index;
	auto it = forbidIndex.find(targetIndex);
	if (it != forbidIndex.end()) {
		std::vector<int> vec = it->second;
		newLabel.extensionVertex.erase(std::remove_if(newLabel.extensionVertex.begin(), newLabel.extensionVertex.end(), [&vec](int elem) {
			return std::find(vec.begin(), vec.end(), elem) != vec.end();
			}), newLabel.extensionVertex.end());
	}

	//the accumulated resource consumption must not be exceed the capacity
	if (newLabel.demand + problem->getWeight(targetIndex) > parameters->planUp) {
		return std::nullopt;
	}
	newLabel.demand += problem->getWeight(targetIndex);

	//calculate the number of coils involved in the cut
	for (int i = 0; i < newLabel.multiS.size(); i++) {
		if (std::find(newLabel.multiS[i].begin(), newLabel.multiS[i].end(), targetIndex) != newLabel.multiS[i].end()) {
			newLabel.recordSets[i].insert(targetIndex);
			newLabel.numVertexinS[i] = newLabel.recordSets[i].size();
			//newLabel.numVertexinS[i] = newLabel.numVertexinS[i] + 1;
		}
	}
	//calculate the penalty of SR
	for (int i = 0; i < newLabel.multiS.size(); i++) {
		if (label.numVertexinS[i] < 2 && newLabel.numVertexinS[i] == 2) {
			newLabel.cost -= rmp->getDualVariable(parameters->N + parameters->L + 2 + i);
		}
	}
	//calculate the RC need consider the terminal node and begin node (c_ij = 0)
	if (targetIndex == end) {
		//For plan lower
		if (newLabel.demand < parameters->planLower || newLabel.demand > parameters->planUp) {
			return std::nullopt;
		}
		// if the target vertex is sink node, return label directly
		newLabel.cost -= rmp->getDualVariable(parameters->N + parameters->L) - rmp->getDualVariable(parameters->N + parameters->L + 1);
		newLabel.cost += parameters->alpha2 * parameters->setUpCost;
		newLabel.curIndex = targetIndex;
		return newLabel;
	}

	if (label.curIndex == 0) {
		newLabel.cost -= rmp->getDualVariable(targetIndex - 1);
		newLabel.cost = newLabel.cost - problem->getWeight(targetIndex) * (problem->getFlow(targetIndex) == "C512" ? 1 : 0) * rmp->getDualVariable(problem->getVertices() + 0);
		newLabel.cost = newLabel.cost - problem->getWeight(targetIndex) * (problem->getFlow(targetIndex) == "C008" ? 1 : 0) * rmp->getDualVariable(problem->getVertices() + 1);
		newLabel.curIndex = targetIndex;
		return newLabel;
	}
	newLabel.cost = newLabel.cost + parameters->alpha1 * problem->getCost(newLabel.curIndex - 1, targetIndex - 1) - rmp->getDualVariable(targetIndex - 1);
	newLabel.cost = newLabel.cost - problem->getWeight(targetIndex) * (problem->getFlow(targetIndex) == "C512" ? 1 : 0) * rmp->getDualVariable(problem->getVertices() + 0);
	newLabel.cost = newLabel.cost - problem->getWeight(targetIndex) * (problem->getFlow(targetIndex) == "C008" ? 1 : 0) * rmp->getDualVariable(problem->getVertices() + 1);

	newLabel.curIndex = targetIndex;
	return newLabel;
}

