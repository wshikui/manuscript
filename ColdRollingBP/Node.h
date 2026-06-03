#pragma once

#include <optional>
#include "ColumnGeneration.h"

class Node {
public:
    Problem* pro{};
    Parameters* para{};

    // the real vertex's index: {1, 2, ...}, is not constraint's index, {0, 1, 2, ...}
    std::vector<int> verticesWithEquality;
    std::vector<int> vertexWithZero; 
    std::map<int, std::vector<int>> forbid_coils;

    std::vector<std::vector<int>> S{};
    std::map<std::string, std::vector<int>> poolIndex;
    std::map<std::string, std::vector<std::bitset<Config::MAX_COILS>>> pool_coil_ngsets;

    std::vector<Column>* globalPool{};
    std::vector<Column> localPool;
    std::shared_ptr<BranchingRule> branchRule{};

    ColumnGeneration* cg{};
    std::vector<std::pair<Column, double>> basicColumns;
    std::vector<std::pair<Column, double>> mipBaseColumns;

    // LP solution (lower bound)
    double solValue{};
    double mipSolValue{};

    std::optional<double> fatherLb;
    int depth{};

    // describe branch rule
    std::string name;
    static constexpr double cplexEpsilon = 1e-5;

    bool cgSolve;

    double totalTimeMP{};
    double totalTimeSp{};
    double totalTime{};

    Node() = default;

    //for root node with ng routes, the subset is empty by default
    Node(
        Problem* _pro,
        Parameters* _para,
        std::vector<int>& _vertices,
        std::vector<Column>* _pool,
        std::vector<Column>& _localPool,
        std::optional<double> _fatherLb,
        std::map<std::string, std::vector<int>> _poolIndex,
        std::map<std::string, std::vector<std::bitset<Config::MAX_COILS>>> _pool_coil_ngsets
    );


    //copy constructor for child node
    Node(const Node& father_node, std::shared_ptr<BranchingRule> branchingRule, std::string _name);
    //Node(const Node& father_node, std::shared_ptr<BranchingRule> branchingRule, std::string _name)

    //copy constructor to save the integer node 
    Node(const Node& father_node);

    ~Node();

    /**
     * @brief run column generation without SR inequalities
     * (not consider the SR)
     * @param node_number
    */
    void solver(unsigned int node_number);

    /**
     * @brief run the column generation with SR inequalities
     * @param node_number
    */
    void solverAddSR(unsigned int node_number);


    //bool solve_lp(std::vector<Column>& pool, std::vector<int>& vertexEquality);


    /**
     * @brief solve the master problem as integer problem
     * @param pool the columns in this node
     * @return
    */
    bool solve_mip(std::vector<Column>& pool);

    //collect feature for branching
    std::tuple<bool, double> get_lower_sr(unsigned int node_number);


    [[nodiscard]] bool isFeasible() const;
    [[nodiscard]] bool has_fractional_solution() const;

    bool has_basic_column_with_cycles() const;
    bool has_integer_coils_with_frac_sol() const;


    
private:
    void determine_equality_constraints();

    void remove_incompatible_columns();

    void no_select_vertex();

    void add_forbid_coils();

    bool hasDuplicates(const std::vector<std::vector<int>>& A);

};

// To create a priority queue
// lower bound priority
class BBNodeCompare {
public:
    bool operator()(const Node* const& n1, const Node* const& n2) const {
        if (!n1->fatherLb) {
            return true;
        }
        if (!n2->fatherLb) {
            return false;
        }
        return (*n1->fatherLb > *n2->fatherLb);
    }
};

