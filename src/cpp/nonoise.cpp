#include <iostream>
#include <string>
#include <cmath>

#include "util/Residue_Representation.h"
#include "util/Arg_Parser.h"
#include "save_map.cpp"
#include "util/data.h"

using namespace std;

#define h 6.62607015e-34

int main(int argc, char* argv[]){
    
    Arg_Parser arg_parser(argc, argv);
  
    if( !( arg_parser.exists( string("-f1") ) && arg_parser.exists( string("-f2") ) && arg_parser.exists( string("-dt1") ) && arg_parser.exists( string("-dt2") ) && arg_parser.exists( string("-n") ) && (argc==11) ) ){
        cerr<<"USAGE:\n"<<argv[0]<<" -f1 input1.par -f2 input2.par -dt0 dt0 -dt1 dt1 -dt2 dt2 -n name"<<endl;
        return 1;
    }
    
    if ( strcmp( arg_parser.get_ext( arg_parser.get("-f1") ) , "par") ) {
    // check for the extensions of the input file
    cerr<<"USAGE:\n"<<argv[0]<<" -f1 input1.par -f2 input2.par -dt1 dt1 -dt2 dt2 -n name"<<endl;
    exit(EXIT_FAILURE);
    }
   
   if ( strcmp( arg_parser.get_ext( arg_parser.get("-f2") ) , "par") ) {
    // check for the extensions of the input file
    cerr<<"USAGE:\n"<<argv[0]<<" -f1 input1.par -f2 input2.par -dt1 dt1 -dt2 dt2 -n name"<<endl;
    exit(EXIT_FAILURE);
    }
    
    double dt1 = std::stod(arg_parser.get("-dt1"));
    double dt2 = std::stod(arg_parser.get("-dt2"));

    unsigned int it1 = 0;
    unsigned int it2 = 0;

    double C1 = 1;
    double C2 = 1;
    
    if (std::find(data::dts.begin(), data::dts.end(), dt1) != data::dts.end()) {
    
        it1 = std::find(data::dts.begin(), data::dts.end(), dt1) - data::dts.begin();

    } else {
    
        double diff = std::abs(data::dts[it1] - dt1);

        for (size_t it = 0; it < data::dts.size(); it++) {
            if (diff > std::abs(data::dts[it] - dt1) ) {

                it1 = (unsigned int) it;
                diff = std::abs(data::dts[it] - dt1);

            }
        }

        if (dt1 < 1) { C1 = data::approx_linear(dt1, it1);}
        else { C1 = data::approx(dt1, it1);}

    }

    if (std::find(data::dts.begin(), data::dts.end(), dt2) != data::dts.end()) {

        it2 = std::find(data::dts.begin(), data::dts.end(), dt2) - data::dts.begin();

    } else {

        double diff = std::abs(data::dts[it2] - dt2);

        for (size_t it = 0; it < data::dts.size(); it++) {
            if (diff > std::abs(data::dts[it] - dt2) ) {

                it2 = (unsigned int) it;
                diff = std::abs(data::dts[it] - dt2);

            }
        }

        if (dt2 < 1) { C2 = data::approx_linear(dt2, it2);}
        else { C2 = data::approx(dt2, it2);}
    }

    double C = data::noise[it2]/data::noise[it1];

    //std::cout << C1 << " " << C2 << " " << C << std::endl;
    
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

                            mutual -= (C2*C/(C1 - C2*C)) * (mat1->getMutual(type1, type2, dofs1[type1][idx1], dofs2[type2][idx2]));
                            mutual += (C1/(C1 - C2*C)) * (mat2->getMutual(type1, type2, dofs1[type1][idx1], dofs2[type2][idx2]));
                                
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
    
    save_map(NResidues, map, names, real_numbers, arg_parser.get("-n"));

    return 0;
}
