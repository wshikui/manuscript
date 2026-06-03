//
// Created by wshikui on 2023/9/14.
//
#include "Material.h"

Material::Material(const std::vector<std::string>& names, const std::vector<std::string>& values)
{
	if (names.size() != values.size())
	{
		std::cout << "Construct Error -> input size inequality" << std::endl;
	}
	for (size_t i = 0; i < names.size(); i++)
	{
		Attribute attribute = Tool::stringToeNum(names[i]);
		infos.insert(std::pair<Attribute, std::string>(attribute, values[i]));
	}
}

Material::Material(const Material& material)
{
	//清空已存在内容
	infos.clear();
	std::map<Attribute, std::string> oldInfo = material.getInfo();
	for (const auto& item : oldInfo)
	{
		infos.insert(item);
	}
}

Material::Material(const Material*& material)
{
	infos.clear();
	std::map<Attribute, std::string> oldInfo = material->getInfo();
	for (const auto& item : oldInfo)
	{
		infos.insert(item);
	}
}

const std::map<Attribute, std::string>& Material::getInfo() const
{
	return infos;
}

double Material::getValue(Attribute attribute)
{
	std::string str = infos[attribute];
	auto value = Tool::to_number<double>(str);
	return value;
}

std::string Material::getStr(Attribute attribute)
{
	return infos[attribute];
}

void Material::setStatus()
{
	this->isSelected = true;
}

bool Material::getStatus()
{
	return this->isSelected;
}
