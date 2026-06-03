#include "MathModel.h"
#include "InitColumns.h"
#include "BBTree.h"
#include "Labelling.h"
#include "GRB.h"
#include "gurobi_c++.h"
#include <ilcplex/ilocplex.h>

std::vector<Material*> load_standard_file(std::string path);

void printCoils(const std::vector<Material*>& _coils);
std::string get_path_name(std::string path);
int getPlansUp(const std::vector<Material*>& _coils, Parameters* parameters);

int main() {

    int argc = 10;
    const char* argv[] = {
        "ColdRollingBP.exe",
        "8",
        "50",
        "90",
        "150",
        "150",
        "Coils_32.txt",
        "0",    //cpx
        "0",    //gub
        "1",     //bp
        "0.5",
        "0.5"
    };

    std::cout << "start" << std::endl;
    Inputs input{};

    input.coilInPool = std::stoi(argv[1]);
    input.planLower = std::stoi(argv[2]);
    input.planUp = std::stoi(argv[3]);
    input.target_512 = std::stoi(argv[4]);
    input.target_008 = std::stoi(argv[5]);

    std::string data_path = argv[6];
    int is_solve_cpx_flag = std::stoi(argv[7]);
    int is_solve_gub_flag= std::stoi(argv[8]);
    int is_explore_tree = std::stoi(argv[9]);
    double alpha_1 = std::stod(argv[10]);
    double alpha_2 = std::stod(argv[11]);

    input.alpha_1 = alpha_1;
    input.alpha_2 = alpha_2;

    std::cout << "file name:" << get_path_name(data_path) << std::endl;

    /// get parameters and problem
    auto parameters = new Parameters(input);

    std::vector<Material*> materials = load_standard_file(data_path);

    //calucateRouteCost(materials);
    int plansUP = getPlansUp(materials, parameters);
    parameters->setMaxPlans(plansUP);

    std::cout << "coils in each pool:" << input.coilInPool << std::endl;
    std::cout << "plan capacity:" << input.planLower << " " << input.planUp << std::endl;
    std::cout << "target mass:" << input.target_512 << " " << input.target_008 << std::endl;
    std::cout << "maximum number of plans:" << plansUP << std::endl;

    auto problemUnit = new Problem(materials);

    /// get init columns
    InitColumns initColumns(materials, problemUnit, parameters);
    initColumns.getSolutions();
    std::vector<Column*> initColumnPool = initColumns.getInitPool();
    initColumns.getInitCost();


    /// math model
    if (is_solve_cpx_flag == 1) {
        std::cout << "---------------CPX-------------------" << std::endl;
        //solve the math model
        Parameters para(input);
        para.setMaxPlans(plansUP);

        MathModel model(para, materials);
        model.optimize();
        model.printResult(problemUnit);

    }
    
    if (is_solve_gub_flag == 1) {
        Parameters para(input);
        para.setMaxPlans(plansUP);

        GRBEnv env(true);
        env.start();
        GRBModel model = GRBModel(env);
        SolverModel grb_math_model(para, materials, env, model);
        grb_math_model.optimize(model);
        grb_math_model.print_result(model);
        grb_math_model.print_variables(model, problemUnit, para);
    }


    //-------------- branch and bound test --------------//
    if (is_explore_tree == 1) {
        //run the BP algorithm
        std::string fileName = get_path_name(data_path);
        BBTree tree(fileName, problemUnit, parameters, initColumnPool);
        tree.explore_tree();
    }

    //----------------------closing work-----------------//
    delete parameters;
    for (auto& item : initColumnPool) {
        delete item;
    }
    initColumnPool.clear();
    delete problemUnit;
    for (auto& item : materials) {
        delete item;
    }
    materials.clear();

    return 0;
}


std::vector<Material*> load_standard_file(std::string path)
{
    std::vector<Material*> materials;
    std::ifstream infile;
    infile.open(path.data());

    if (!infile.is_open()) {
        std::cout << "Open materials txt Failed" << std::endl;
    }

    std::string content;
    getline(infile, content);
    getline(infile, content);
    getline(infile, content);
    getline(infile, content);
    getline(infile, content);
    getline(infile, content);
    getline(infile, content);
    //getline(infile, content);

    int _index = 1;
    while (getline(infile, content)) {
        if (content.empty()) { continue; }

        if (content.find("DEPOT_SECTION") != std::string::npos) { break; }

        std::vector<std::string> keys;
        std::vector<std::string> values;

        std::vector<std::string> kvs = Tool::super_split(content, "; ");
        for (auto& item : kvs) {
            std::vector<std::string> kv = Tool::super_split(item, "=>");
            if (kv.size() < 2) {
                kv.emplace_back("");
                std::cout << "Load data: this line has empty attribute" << std::endl;
            }
            keys.push_back(kv.at(0));
            values.push_back(kv.at(1));
        }
        auto* material = new Material(keys, values);
        material->index = _index;
        materials.emplace_back(material);
        _index++;
    }


    return materials;
}


void printCoils(const std::vector<Material*>& _coils) {
    std::cout << "__________coil info__________" << std::endl;

    int width = 10;
    std::cout << std::setw(width) << std::left << "index" << "|"
        << std::setw(width) << std::left << "NO" << "|"
        << std::setw(width) << std::left << "OutWidth" << "|"
        << std::setw(width) << std::left << "WT" << "|"
        << std::setw(width) << std::left << "FLOW" << "|"
        << std::endl;
    for (const auto& item : _coils) {
        std::cout << std::setw(width) << std::left << item->index << "|"
            << std::setw(width) << std::left << item->getStr(Attribute::MAT_NO) << "|"
            << std::setw(width) << std::left << item->getStr(Attribute::OUT_MAT_WIDTH) << "|"
            << std::setw(width) << std::left << item->getValue(Attribute::IN_MAT_WT) << "|"
            << std::setw(width) << std::left << item->getStr(Attribute::FLOW) << "|"
            << std::endl;
    }

    double massAll = 0;
    for (const auto& coil : _coils) {
        massAll += coil->getValue(Attribute::IN_MAT_WT);
    }
    std::cout << "total weight : " << massAll << std::endl;
    std::cout << "average weight : " << massAll / (int)_coils.size() << std::endl;
}

std::string get_path_name(std::string path)
{
    size_t last_backslash = path.rfind('\\');
    if (last_backslash == std::string::npos) {
        size_t last_slash = path.rfind('/');
        std::string _file_name = path.substr(last_slash + 1);
        size_t _dot_pos = _file_name.rfind('.');
        std::string _result = _file_name.substr(0, _dot_pos);
        return _result;
    }
    std::string filename = path.substr(last_backslash + 1);
    size_t dot_pos = filename.rfind('.');
    std::string result = filename.substr(0, dot_pos);

    return result;
}

int getPlansUp(const std::vector<Material*>& _coils, Parameters* parameters) {
    double mass = 0.0;
    for (const auto& _coil : _coils) {
        mass += _coil->getValue(Attribute::IN_MAT_WT);
    }
    int num = static_cast<int>(std::floor(mass / parameters->planLower));
    return num;
}



