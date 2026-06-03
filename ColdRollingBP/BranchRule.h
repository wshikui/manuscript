#pragma once
#include "Column.h"

class BranchingRule {
public:
	virtual bool is_column_compatible(const Column& column) const = 0;
	virtual bool should_row_be_equality(int vertexIndex) const = 0;
	virtual bool remove_vertex(int vertexIndex) const = 0;
	virtual bool remove_arc(int source, int target) const = 0;
};

class IncludeCoil : public BranchingRule {
public:
	const int index{};
	IncludeCoil(const int _index) : index{ _index } {}

	// override virtual
	bool is_column_compatible(const Column& column) const override;
	bool should_row_be_equality(int vertexIndex) const override;
	bool remove_vertex(int vertexIndex) const override;
	bool remove_arc(int source, int target) const override;
};

class ExcludeCoil : public BranchingRule {
public:
	const int index{};
	ExcludeCoil(const int _index) : index(_index) {};

	bool is_column_compatible(const Column& column) const override;
	bool should_row_be_equality(int vertexIndex) const override;
	bool remove_vertex(int vertexIndex) const override;
	bool remove_arc(int source, int target) const override;
};

class ForceConsecutive : public BranchingRule {
public:
	const std::pair<int, int> arc{};
	ForceConsecutive(const std::pair<int, int> _arc) : arc(_arc) {};

	bool is_column_compatible(const Column& column) const override;
	bool should_row_be_equality(int vertexIndex) const override;
	bool remove_vertex(int vertexIndex) const override;
	bool remove_arc(int source, int target) const override;
};

class ForbidConsecutive : public BranchingRule {
public:
	const std::pair<int, int> arc{};
	ForbidConsecutive(const std::pair<int, int> _arc) : arc(_arc) {};

	bool is_column_compatible(const Column& column) const override;
	bool should_row_be_equality(int vertexIndex) const override;
	bool remove_vertex(int vertexIndex) const override;
	bool remove_arc(int source, int target) const override;
};
