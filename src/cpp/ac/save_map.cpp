#include "util/json.hpp"

#include <vector>
#include <fstream>
#include <string>


using nlohmann::json;

void save_map(unsigned int NResidues, std::vector < std::vector <double>> first, std::vector < std::vector <double>> second, std::vector < std::vector <double>> third, std::vector < std::vector <double>> first_d, std::vector < std::vector <double>> second_d, std::vector < std::vector <double>> third_d, std::vector <std::string> names, std::vector <int> real_numbers, std::string name) {
    json save;
    save["NResidues"] = NResidues;
    save["names"] = names;
    save["first"] = first;
    save["second"] = second;
    save["third"] = third;
    save["first_with_dihedrals"] = first_d;
    save["second_with_dihedrals"] = second_d;
    save["third_with_dihedrals"] = third_d;
    save["real_numbers"] = real_numbers;
    
    std::fstream file;
    file.open("output/ac_" + name + ".json", std::ios::trunc | std::ios::out);
    file << std::setw(0) << save;
}
