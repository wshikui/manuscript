#pragma once
#include "Tool.h"

class Material
{
public:
	int index = 0;

	Material() = default;

	Material(const std::vector<std::string>& names, const std::vector<std::string>& values);

	Material(const Material& material);

	explicit Material(const Material*& material);
	const std::map<Attribute, std::string>& getInfo() const;

	// 返回double属性
	double getValue(Attribute attribute);

	//返回string属性
	std::string getStr(Attribute attribute);

	void setStatus();
	bool getStatus();


private:
	std::map<Attribute, std::string> infos;

	bool isSelected = false;
};

