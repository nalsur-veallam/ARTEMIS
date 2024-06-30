#include "util/json.hpp"

#include <vector>
#include <fstream>
#include <string>


using nlohmann::json;

void save_map(unsigned int NResidues, std::vector <std::vector <double>> map, std::vector <std::string> names, std::vector <int> real_numbers, std::string name) {
    json save;
    save["NResidues"] = NResidues;
    save["names"] = names;
    save["map"] = map;
    save["real_numbers"] = real_numbers;
    
    std::fstream file;
    file.open(name, std::ios::trunc | std::ios::out);
    file << std::setw(0) << save;
}
