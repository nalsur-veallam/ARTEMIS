#include <iostream>
#include <string>
#include <cmath>

#include "util/Residue_Representation.h"
#include "util/Arg_Parser.h"
#include "save_map.cpp"

using namespace std;

#define h 6.62607015e-34

int main(int argc, char* argv[]){
    
    Arg_Parser arg_parser(argc, argv);
  
    if( !( arg_parser.exists( string("-f1") ) && arg_parser.exists( string("-f2") ) && arg_parser.exists( string("-dt1") ) && arg_parser.exists( string("-dt0") )&& arg_parser.exists( string("-dt2") ) && arg_parser.exists( string("-n") ) && (argc==13) ) ){
        cerr<<"USAGE:\n"<<argv[0]<<" -f1 input1.par -f2 input2.par -dt0 dt0 -dt1 dt1 -dt2 dt2 -n name"<<endl;
        return 1;
    }
    
    if ( strcmp( arg_parser.get_ext( arg_parser.get("-f1") ) , "par") ) {
    // check for the extensions of the input file
    cerr<<"USAGE:\n"<<argv[0]<<" -f1 input1.par -f2 input2.par -dt0 dt0 -dt1 dt1 -dt2 dt2 -n name"<<endl;
    exit(EXIT_FAILURE);
    }
   
   if ( strcmp( arg_parser.get_ext( arg_parser.get("-f2") ) , "par") ) {
    // check for the extensions of the input file
    cerr<<"USAGE:\n"<<argv[0]<<" -f1 input1.par -f2 input2.par -dt0 dt0 -dt1 dt1 -dt2 dt2 -n name"<<endl;
    exit(EXIT_FAILURE);
    }
    
    double dt0 = std::stod(arg_parser.get("-dt0"));
    double dt1 = std::stod(arg_parser.get("-dt1"));
    double dt2 = std::stod(arg_parser.get("-dt2"));
    
    Residue_Representation rep1(arg_parser.get("-f1"), false, MODE_TOTAL);
    Entropy_Matrix* mat1 = rep1.getEntropy_Matrix();
    
    Residue_Representation rep2(arg_parser.get("-f2"), false, MODE_TOTAL);
    Entropy_Matrix* mat2 = rep2.getEntropy_Matrix();
    
    vector <vector <double>> map;
    vector <string> names;
    
    unsigned int NResidues = rep1.getNResidues();  // TODO: Errors 
    
    
    for (unsigned int resid1 = 1; resid1 <= NResidues; resid1++) {
        
        vector <double> mie;
        names.push_back(rep1.getResidueName(resid1));
        
        for (unsigned int resid2 = 1; resid2 <= NResidues; resid2++) {
            
            vector <vector< int > >dofs1;
            dofs1.push_back(rep1.getBonds(resid1));
            dofs1.push_back(rep1.getAngles(resid1));
            dofs1.push_back(rep1.getDihedrals(resid1));
            
            vector <vector< int > >dofs2;
            dofs2.push_back(rep1.getBonds(resid2));
            dofs2.push_back(rep1.getAngles(resid2));
            dofs2.push_back(rep1.getDihedrals(resid2));
            
            
            double mutual = 0;
            
            for (unsigned char type1 = 0; type1 < 3; type1++){ // for bonds, angles and dihedrals of the first member of the dof pair
                for(unsigned int idx1 = 0; idx1 < dofs1[type1].size(); idx1++){ // and all dofs of the current type of the first member of the dof pair
                    for (unsigned char type2 = 0; type2 < 3; type2++){ // and all "later" types for the second member of the dof pair
                        for(unsigned int idx2 = 0; idx2 < dofs2[type2].size(); idx2++ ){ // and all "later" dofs for the second member of the dof pair
                            
                            mutual += (dt2-dt0)/(dt2-dt1)*mat1->getMutual(type1, type2, dofs1[type1][idx1], dofs2[type2][idx2]);
                            mutual += (dt1-dt0)/(dt2-dt1)*mat2->getMutual(type1, type2, dofs1[type1][idx1], dofs2[type2][idx2]);
                                
                        }
                    }
                }
            }
            if (mutual < 0) {mutual = 0;}
            mie.push_back(mutual); 
            
        }
        
        map.push_back(mie);
    }
    
    vector <int> real_numbers;
    for (unsigned int i = 0; i < NResidues; i++) {
        real_numbers.push_back(rep1.getResidueNumber(i + 1));
    }
    
    vector <double> entropies;
    for (unsigned int resid = 1; resid <= NResidues; resid++) {
        
        entropies.push_back(0);
        
    }
    
    save_map(NResidues, map, names, real_numbers, arg_parser.get("-n"), entropies);

    return 0;
}
