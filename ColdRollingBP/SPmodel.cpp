#include "SPmodel.h"

SPmodel::SPmodel(const std::vector<std::pair<Column, double>>& columns)
{
	createModel();

	//basic columns
	P = static_cast<int>(columns.size());
	N = columns.front().first.numRealCoils;

	mModel.setName("SPmodel");

	try {
		// create the SP model
		X = IloNumVarArray(mEnv, N);
		for (int i = 0; i < N; i++) {
			X[i] = IloNumVar(mEnv, 0, 1, ILOINT);
		}

		IloExpr con(mEnv);
		for (int i = 0; i < N; i++) {
			con += X[i];
		}
		mModel.add(con == n);
		con.end();

		//for (int p = 0; p < P; p++) {
		//	IloExpr conP(mEnv);
		//	for (int i = 0; i < N; i++) {
		//		conP += columns[p].first.contains(i + 1);
		//	}
		//	mModel.add(IloIfThen(mEnv, conP <= k - 1, conP == 0));
		//	mModel.add(IloIfThen(mEnv, conP >= k, conP == k));
		//	conP.end();
		//}

		IloExpr obj(mEnv);
		for (int p = 0; p < P; p++) {
			for (int i = 0; i < N; i++) {
				obj += (static_cast<double>(1) / k) * columns[p].first.contains(i + 1) * X[i] * columns[p].second;
			}
		}
		obj -= std::floor(n / k);

		mModel.add(IloMaximize(mEnv, obj));
		obj.end();
	}
	catch (IloException& ex)
	{
		cerr << "SP model Error -> " << ex << endl;
	}
	catch (...)
	{
		cerr << "Error when SP model" << endl;
	}

}

bool SPmodel::getS(const std::vector<std::pair<Column, double>>& columns)
{
	bool flag = optimize();
	bool checkFlag = false;

	//the items in S are index of real coils

	if (flag) {
		//maybe exist a SR inequality
		if (cplex.getObjValue() > 0) {
			// std::cout << "maybe exist SR" << std::endl;
			for (int i = 0; i < N; i++) {
				if (cplex.getValue(X[i]) > 0.0001) {
					//std::cout << i + 1 << " -> " << cplex.getValue(X[i]) << std::endl;
					S.emplace_back(i + 1);
				}
			}
			//check the SR
			std::cout << "check the columns in sr separation problem" << std::endl;
			for (int i = 0; i < columns.size(); i++) {

				//columns[i].first.printRoute();

				//std::cout << "------------------" << std::endl;
			}
			double result = 0;
			for (int p = 0; p < columns.size(); p++) {
				int coff = 0;
				int num = 0;
				for (const auto& item : S) {
					if (columns[p].first.contains(item)) {
						num++;
					}
				}
				if (num >= k) {
					coff = 1;
				}
				result += columns[p].second * coff;
			}
			if (result > 1 + EPS) {
				//find a SR inequality and add this inequality to RMP
				checkFlag = true;
				std::cout << "find a SR inequality" << std::endl;
				//std::cout << "sum of columns with SR -> " << std::fixed << std::setprecision(18) << result << std::endl;
				//std::cout << std::setprecision(-1);
				std::cout << "sum of columns with SR -> " << result << std::endl;
				for (const auto& item : S) {
					std::cout << item << " ";
				}
				std::cout << std::endl;
			}
		}
	}

	return checkFlag;
}

bool SPmodel::getSenum(const std::vector<std::pair<Column, double>>& columns)
{
	//find the SR inequality
	bool isFound = false;
	int P = static_cast<int>(columns.size());

	std::set<int> coilIndexes;
	std::vector<std::vector<int>> condidateS;
	for (const auto& item : columns) {
		std::vector<int> indexes = item.first.getRoute();
		for (const auto& index : indexes) {
			coilIndexes.insert(index);
		}
	}

	auto it1 = coilIndexes.begin();
	auto it2 = std::next(it1);
	auto it3 = std::next(it2);

	for (; it1 != std::prev(coilIndexes.end(), 2); ++it1) {
		it2 = std::next(it1);
		for (; it2 != std::prev(coilIndexes.end(), 1); ++it2) {
			it3 = std::next(it2);
			for (; it3 != coilIndexes.end(); ++it3) {
				condidateS.emplace_back(std::initializer_list<int>{ *it1, * it2, * it3 });
			}
		}
	}

	double temp = 0;
	std::vector<int> result;
	for (const auto& subset : condidateS) {
		double obj = 0;
		for (int p = 0; p < P; p++) {
			//change the parameters____
			if (SubProblem::countContainedElements(columns[p].first.getRoute(), subset) >= 2) {
				obj += columns[p].second;
			}
		}
		if (obj > temp) {
			temp = obj;
			result = subset;
		}
	}
	if (temp > 1) {
		S = result;
		//std::cout << "found S by enum -> " << S[0] << " " << S[1] << " " << S[2] << " " << std::endl;
		
		double obj_ = 0;
		for (int p = 0; p < P; p++) {
			if (SubProblem::countContainedElements(columns[p].first.getRoute(), S) >= 2) {
				obj_ += columns[p].second;
			}
		}
		//std::cout << "the sum of SR -> " << obj_ << std::endl;
		isFound = true;
	}

	return isFound;
}

