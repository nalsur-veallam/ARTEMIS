#ifndef MY_SAVEMAP_H
#define MY_SAVEMAP_H

#include "util/json.hpp"

#include <vector>
#include <fstream>
#include <string>

using nlohmann::json;

void save_map(unsigned int NResidues, std::vector <std::vector <double>> map, std::vector <std::string> names) {
    json save;
    save["NResidues"] = NResidues;
    
    for (unsigned int resid = 0; resid < NResidues; resid++) {
        save["MIE for the residue number " +  std::to_string(resid+1) + " (" + names[resid] + ") with all residues with a higher number"] = map[resid];
    }
    
    std::fstream file;
    file.open("map.json", std::ios::trunc | std::ios::out);
    file << std::setw(0) << save;
}

#endif
