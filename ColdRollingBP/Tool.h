#pragma once
#include "Attribute.h"
#include <algorithm>
//#include "Algorithm.h"
#include <numeric>
#include <map>
#include <set>
#include <iostream>
#include <sstream>
#include <istream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <random>
#include <cmath>

class Tool
{
public:

	static std::string eNumToString(Attribute attribute)
	{
		// 定义eNumerate和string对应关系的map
		std::map<Attribute, std::string> attributeToString = {
				std::pair<Attribute, std::string>(Attribute::IN_MAT_THICK, "IN_MAT_THICK"),
				std::pair<Attribute, std::string>(Attribute::MAT_NO, "MAT_NO"),
				std::pair<Attribute, std::string>(Attribute::IN_MAT_WIDTH, "IN_MAT_WIDTH"),
				std::pair<Attribute, std::string>(Attribute::OUT_MAT_THICK, "OUT_MAT_THICK"),
				std::pair<Attribute, std::string>(Attribute::OUT_MAT_WIDTH, "OUT_MAT_WIDTH"),
				std::pair<Attribute, std::string>(Attribute::MAT_STATUS, "MAT_STATUS"),
				std::pair<Attribute, std::string>(Attribute::IN_MAT_WT, "IN_MAT_WT"),
				std::pair<Attribute, std::string>(Attribute::CHS_TYPE, "CHS_TYPE"),
				std::pair<Attribute, std::string>(Attribute::RULE_CODE, "RULE_CODE"),
				std::pair<Attribute, std::string>(Attribute::FLOW, "FLOW"),
		};
		return attributeToString[attribute];
	}

	static Attribute stringToeNum(const std::string& label)
	{
		// the map: label -> attribute (eNum)
		std::map<Attribute, std::string> attributeToString = {
				std::pair<Attribute, std::string>(Attribute::IN_MAT_THICK, "IN_MAT_THICK"),
				std::pair<Attribute, std::string>(Attribute::MAT_NO, "MAT_NO"),
				std::pair<Attribute, std::string>(Attribute::IN_MAT_WIDTH, "IN_MAT_WIDTH"),
				std::pair<Attribute, std::string>(Attribute::OUT_MAT_THICK, "OUT_MAT_THICK"),
				std::pair<Attribute, std::string>(Attribute::OUT_MAT_WIDTH, "OUT_MAT_WIDTH"),
				std::pair<Attribute, std::string>(Attribute::MAT_STATUS, "MAT_STATUS"),
				std::pair<Attribute, std::string>(Attribute::IN_MAT_WT, "IN_MAT_WT"),
				std::pair<Attribute, std::string>(Attribute::CHS_TYPE, "CHS_TYPE"),
				std::pair<Attribute, std::string>(Attribute::RULE_CODE, "RULE_CODE"),
				std::pair<Attribute, std::string>(Attribute::FLOW, "FLOW"),
		};

		Attribute attribute = Attribute::MAT_STATUS;
		for (const auto& item : attributeToString)
		{
			if (item.second == label)
			{
				attribute = item.first;
			}
		}
		return attribute;
	}

	static std::vector<std::string> super_split(const std::string& str, const std::string& sep)
	{
		std::vector<std::string> result;

		std::string::size_type pos1, pos2;
		size_t len = str.length();
		pos2 = str.find(sep);
		pos1 = 0;
		while (std::string::npos != pos2)
		{
			result.push_back(str.substr(pos1, pos2 - pos1));

			pos1 = pos2 + sep.size();
			pos2 = str.find(sep, pos1);
		}

		if (pos1 != len)
		{
			result.push_back(str.substr(pos1));
		}

		return result;
	}

	template<typename T>
	static T to_number(const std::string& str)
	{
		std::istringstream iss(str);

		T ret;
		iss >> ret;

		return ret;
	}

	template<typename T>
	static void swap(T& elem1, T& elem2) {
		T tmp = elem1;
		elem1 = elem2;
		elem2 = tmp;
	}

	static void writeToFile(const std::map<int, std::vector<std::string>>& myMap, const std::string& filename) {
		std::ofstream outFile(filename, std::ios::app);
		if (!outFile) {
			std::cerr << "can not open this file: " << filename << std::endl;
			return;
		}
		for (const auto& pair : myMap) {
			int key = pair.first;
			const std::vector<std::string>& values = pair.second;

			if (values.empty()) {
				continue;
			}

			//outFile << "Key: " << key << "\n";

			for (const auto& value : values) {
				outFile << value << "\t";
			}
			outFile << "\n";
		}
		outFile << "---------------------" << "\n";
		outFile.close();
	}


	//template<typename T>
	//static T calculateVariance(const std::vector<T>& vec) {
	//	// Step 1: Calculate mean
	//	T mean = std::accumulate(vec.begin(), vec.end(), 0.0) / vec.size();

	//	// Step 2: Calculate squared differences from the mean
	//	std::vector<T> squared_diffs(vec.size());
	//	std::transform(vec.begin(), vec.end(), squared_diffs.begin(),
	//		[mean](T x) { return std::pow(x - mean, 2); });

	//	// Step 3: Calculate variance (mean of squared differences)
	//	T variance = std::accumulate(squared_diffs.begin(), squared_diffs.end(), 0.0) / vec.size();
	//	return variance;
	//}

	static double calculateVariance(const std::vector<double>& vec) {
		// Step 1: Calculate mean
		double mean = std::accumulate(vec.begin(), vec.end(), 0.0) / vec.size();

		// Step 2: Calculate squared differences from the mean
		std::vector<double> squared_diffs(vec.size());
		std::transform(vec.begin(), vec.end(), squared_diffs.begin(),
			[mean](double x) { return std::pow(x - mean, 2); });

		// Step 3: Calculate variance (mean of squared differences)
		double variance = std::accumulate(squared_diffs.begin(), squared_diffs.end(), 0.0) / vec.size();
		return variance;
	}
};

