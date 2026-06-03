#pragma once
#include "BaseModel.h"
#include "Column.h"
#include "SubProblem.h"

class SPmodel : public BaseModel {

public:
    SPmodel() = default;

    SPmodel(const std::vector<std::pair<Column, double>>& columns);

    bool getS(const std::vector<std::pair<Column, double>>& columns);

    bool getSenum(const std::vector<std::pair<Column, double>>& columns);

    std::vector<int> S{};

private:

    //the number of pattern and coil
    int P = 0;
    int N = 0;
    IloNumVarArray X;

    // |S|
    int n = 3;
    int k = 2;
    //std::vector<int> S{};

    //double EPS = 1e-5;

    double EPS = 0.1;

    //double EPS = 0.2;

    //double EPS = 0.2;
};