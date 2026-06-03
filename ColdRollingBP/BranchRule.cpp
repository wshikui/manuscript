#include "BranchRule.h"

/* ------------------- Include coil ------------------- */
// the columns are all compatible
bool IncludeCoil::is_column_compatible(const Column& column) const
{
	return true;
}

// Determine whether the constraint corresponding to this node index is an equation
bool IncludeCoil::should_row_be_equality(int vertexIndex) const
{
	return index == vertexIndex;
}

bool IncludeCoil::remove_vertex(int vertexIndex) const
{
	return false;
}

bool IncludeCoil::remove_arc(int source, int target) const
{
	return false;
}

/* ------------------- exclude coil ------------------- */
// if this column contains the coil, return false, this column is not compatible;
// compatible (don't contain) -> true
// incompatible (contain) -> false
bool ExcludeCoil::is_column_compatible(const Column& column) const
{
	return !column.contains(index);
}

// do not change the constraint of RMP
bool ExcludeCoil::should_row_be_equality(int vertexIndex) const
{
	return false;
}

bool ExcludeCoil::remove_vertex(int vertexIndex) const
{
	return index == vertexIndex;
}

bool ExcludeCoil::remove_arc(int source, int target) const
{
	return false;
}

bool ForceConsecutive::is_column_compatible(const Column& column) const
{
	const auto& route = column.getRoute();
	for (size_t i = 0; i < route.size(); i++) {
		//column contains arc.first but not contain arc.second (incompatible)
		if (route[i] == arc.first) {
			if (i + 1 == route.size()) {
				return false;
			}
			if (route[i + 1] != arc.second) {
				return false;
			}
		}

		//column not contains arc.first but contain arc.second (incompatible)
		if (route[i] == arc.second) {
			if (i == 0) {
				return false;
			}
			if (route[i - 1] != arc.first) {
				return false;
			}
		}
	}

	//not contain arc.first or arc.second
	return true;
}

bool ForceConsecutive::should_row_be_equality(int vertexIndex) const
{
	//only restrict the arc selection
	return false;
}

bool ForceConsecutive::remove_vertex(int vertexIndex) const
{
	//branching on arc do not remove any coil
	return false;
}

bool ForceConsecutive::remove_arc(int source, int target) const
{
	if (source == arc.first && target != arc.second) { return true; }
	if (source != arc.first && target == arc.second) { return true; }
	return false;
}

bool ForbidConsecutive::is_column_compatible(const Column& column) const
{
	const auto& route = column.getRoute();
	for (size_t i = 0; i < route.size() - 1; i++) {
		if (route[i] == arc.first && route[i + 1] == arc.second) {
			return false;
		}
	}
	return true;
}

bool ForbidConsecutive::should_row_be_equality(int vertexIndex) const
{
	return false;
}

bool ForbidConsecutive::remove_vertex(int vertexIndex) const
{
	return false;
}

bool ForbidConsecutive::remove_arc(int source, int target) const
{
	if (source == arc.first && target == arc.second) { return true; }

	return false;
}
