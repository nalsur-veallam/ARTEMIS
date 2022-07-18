#include "util/json.hpp"

#include <vector>
#include <fstream>
#include <string>

using nlohmann::json;

std::vector <std::vector <double>> map;
std::vector <std::string> names;
unsigned int NResidues;


void read_map() {
    json save;
    
    std::fstream file;
    file.open("map.json", std::ios::in);
    file >> save;
    
    NResidues = save["NResidues"];
    map = save["map"];
    for (int i = 0; i < save["names"].size(); i++) {
        names.push_back(save["names"][i]);
    }
}

std::vector <std::vector <double>> get_map() {
    return map;
}
        
std::vector <std::string> get_names() {
    return names;
}

unsigned int get_NResidues() {
    return NResidues;
}
