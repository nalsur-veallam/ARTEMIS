// A program analysing the mutual information between two residues
// Copyright (C) 2020  Markus Fleck
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License version 3 as 
// published by the Free Software Foundation.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.

#include <iostream>
#include <string>

#include "util/Residue_Representation.h"
#include "util/Arg_Parser.h"
#include "../save_map.cpp"

using namespace std;


int main(int argc, char* argv[]){
    
    Arg_Parser arg_parser(argc, argv);
  
    if( !( arg_parser.exists( string("-f") ) && arg_parser.exists( string("-n") ) && (argc==5) ) ){
        cerr<<"USAGE:\n"<<argv[0]<<" -f input.par -n name"<<endl;
        return 1;
    }
    
    if ( strcmp( arg_parser.get_ext( arg_parser.get("-f") ) , "par") ) {
    // check for the extensions of the input file
    cerr<<"USAGE:\n"<<argv[0]<<" -f input.par -n name"<<endl;
    exit(EXIT_FAILURE);
   }
    
    Residue_Representation rep(arg_parser.get("-f"), false, MODE_TOTAL);
    Entropy_Matrix* mat = rep.getEntropy_Matrix();
    
    
    vector <vector <double>> map;
    vector <string> names;
    
    unsigned int NResidues = rep.getNResidues();
    
    for (unsigned int resid1 = 1; resid1 <= NResidues; resid1++) {
        
        vector <double> mie;
        names.push_back(rep.getResidueName(resid1));
        
        for (unsigned int resid2 = 1; resid2 <= NResidues; resid2++) {
            
            vector <vector< int > >dofs1;
            dofs1.push_back(rep.getBonds(resid1));
            dofs1.push_back(rep.getAngles(resid1));
            dofs1.push_back(rep.getDihedrals(resid1));
            
            vector <vector< int > >dofs2;
            dofs2.push_back(rep.getBonds(resid2));
            dofs2.push_back(rep.getAngles(resid2));
            dofs2.push_back(rep.getDihedrals(resid2));
            
            
            double mutual = 0;
            
            for (unsigned char type1 = 0; type1 < 3; type1++){ // for bonds, angles and dihedrals of the first member of the dof pair
                for(unsigned int idx1 = 0; idx1 < dofs1[type1].size(); idx1++){ // and all dofs of the current type of the first member of the dof pair
                    for (unsigned char type2 = 0; type2 < 3; type2++){ // and all "later" types for the second member of the dof pair
                        for(unsigned int idx2 = 0; idx2 < dofs2[type2].size(); idx2++ ){ // and all "later" dofs for the second member of the dof pair
                            
                            mutual += mat->getMutual(type1, type2, dofs1[type1][idx1], dofs2[type2][idx2]);
                                
                        }
                    }
                }
            }
            
            mie.push_back(mutual); 
            
        }
        
        map.push_back(mie);
    }
    
    vector <int> real_numbers;
    for (unsigned int i = 0; i < NResidues; i++) {
        real_numbers.push_back(rep.getResidueNumber(i + 1));
    }
    
    save_map(NResidues, map, names, real_numbers, arg_parser.get("-n"));

    return 0;
}
