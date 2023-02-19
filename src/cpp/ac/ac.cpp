#define MUTUAL_ARGS nDihedrals,bondsEntropy1D,anglesEntropy1D,dihedralsEntropy1D,bbEntropy,baEntropy,bdEntropy,aaEntropy,adEntropy,ddEntropy


#include <iostream>
#include <fstream>
#include <vector>
#include <string.h>
#include <cstdlib>
#include <algorithm>
#include <numeric>



using namespace std;


#include "util/io.h"
#include "util/util.h"
#include "util/Arg_Parser.h"
#include "save_map.cpp"

struct X { vector <vector <vector <double>>> map; vector <vector <vector <double>>> map_with_dihedrals; unsigned int NResidues = 0; vector <string> names; vector <int> real_numbers;};


struct X f(string file)
{
    
    int double_prec, numFrames, nDihedrals;
    int version, bDens, aDens, dDens, bDens1D, aDens1D, dDens1D;
    vector< vector <int> > dihedrals_top;
    vector <float> masses;
    vector <string> residues;
    vector <int> residueNumbers;
    vector <string> atomNames;
    vector <string> belongsToMolecule;

    double *bondsEntropy1D,*anglesEntropy1D,*dihedralsEntropy1D;
    double *bbEntropy,*aaEntropy,*ddEntropy,*baEntropy,*bdEntropy,*adEntropy;

    vector <vector <vector <double>>> map;
    vector <vector <vector <double>>> map_with_dihedrals;
    vector <unsigned int> bonds;
    vector <unsigned int> angles;
    vector <unsigned int> dihedrals;
    vector <string> names;
    unsigned int NResidues;
    vector <int> real_numbers;

    ifstream infile(file, ios::binary | ios::in);//open the .par file
    if(infile.is_open()) {
        if(read_PAR_header(&infile,&nDihedrals,&double_prec,&numFrames,&dihedrals_top, &masses, &version, &bDens, &aDens, &dDens, &bDens1D, &aDens1D, &dDens1D,&residues,&residueNumbers,&atomNames,&belongsToMolecule)!=0) { //and read the header information
            cerr<<"ERROR READING HEADER OF FILE "<<" !"<<endl;
            struct X tmp;
            return tmp;
        }
        
        string name = residues[0];
        names.push_back(name);
        int j = 0;
        int size = residues.size();
        
        while (j < size - 1) {
            while ((residueNumbers[j+1] == residueNumbers[j]) && (j < size - 1)) {
                j++;
            }
            names.push_back(residues[j+1]);
            name = residues[j+1];
            j++;
        }
        
        names.pop_back();
        real_numbers = residueNumbers;
        int current = 1;
        int total = residueNumbers[0];
        for (int i = 0; i < size; i ++) {
            if (residueNumbers[i] == total)
                residueNumbers[i] = current;
            else {
                total = residueNumbers[i];
                current++;
                residueNumbers[i] = current;
            }
        }
        
        NResidues = names.size();
        
        for (unsigned int i = 0; i < NResidues; i++) {
            vector <vector <double>> line;
            for (unsigned int k = 0; k < NResidues; k++) {
                vector <double> vec;
                line.push_back(vec);
            }
            map.push_back(line);
            map_with_dihedrals.push_back(line);
        }
        
        if (read_PAR_body(&infile,nDihedrals,&bondsEntropy1D, &anglesEntropy1D, &dihedralsEntropy1D, &bbEntropy, &baEntropy, &bdEntropy, &aaEntropy, &adEntropy, &ddEntropy)!=0) {
            cerr<<"ERROR READING FILE "<<" !"<<endl;
            struct X tmp;
            return tmp;
        }

        if (residueNumbers[dihedrals_top[0][0]] == residueNumbers[dihedrals_top[0][1]]) //first bond
            bonds.push_back(residueNumbers[dihedrals_top[0][0]]);
        else
            bonds.push_back(0);
            
        if (residueNumbers[dihedrals_top[0][1]] == residueNumbers[dihedrals_top[0][2]]) //second bond
            bonds.push_back(residueNumbers[dihedrals_top[0][1]]);
        else
            bonds.push_back(0);
        
        if ((residueNumbers[dihedrals_top[0][1]] == residueNumbers[dihedrals_top[0][2]]) && (residueNumbers[dihedrals_top[0][2]] == residueNumbers[dihedrals_top[0][0]])) //first angle
            angles.push_back(residueNumbers[dihedrals_top[0][1]]);
        else
            angles.push_back(0);
        
        for(int i = 0; i < nDihedrals; i++) {// correlate the degrees of freedom with their residues
            
            if (residueNumbers[dihedrals_top[i][2]] == residueNumbers[dihedrals_top[i][3]]) //bonds
                bonds.push_back(residueNumbers[dihedrals_top[i][2]]);
            else
                bonds.push_back(0);
            
            if ((residueNumbers[dihedrals_top[i][1]] == residueNumbers[dihedrals_top[i][2]]) && (residueNumbers[dihedrals_top[i][2]] == residueNumbers[dihedrals_top[i][3]])) //angles
                angles.push_back(residueNumbers[dihedrals_top[i][1]]);
            else
                angles.push_back(0);
            
            if ((residueNumbers[dihedrals_top[i][0]] == residueNumbers[dihedrals_top[i][1]]) && (residueNumbers[dihedrals_top[i][1]] == residueNumbers[dihedrals_top[i][2]]) && (residueNumbers[dihedrals_top[i][2]] == residueNumbers[dihedrals_top[i][3]])) //dihedrals
                dihedrals.push_back(residueNumbers[dihedrals_top[i][1]]);
            else
                dihedrals.push_back(0);
                
        }
        
        
        //then for every bond-bond pair calculate the mutual information
        for(int i = 0; i < nDihedrals + 1; i++) {
            for(int j = i + 1; j < nDihedrals + 2; j++) {
                if (bonds[i] != 0 && bonds[j] != 0) {
                    if (bonds[i] >= bonds[j]) {
                        
                        map[bonds[j] - 1][bonds[i] - 1].push_back(get_mutual(TYPE_BB,i,j,MUTUAL_ARGS));
                
                    }
                    else {
                        
                        map[bonds[i] - 1][bonds[j] - 1].push_back(get_mutual(TYPE_BB,i,j,MUTUAL_ARGS));
                        
                    }
                }
            }
        }
        
        for(int i = 0; i < nDihedrals + 2; i++) { //same for all bond-angle pairs
            for(int j = 0; j < nDihedrals + 1; j++) {
                if (bonds[i] != 0 && angles[j] != 0) {
                    if (bonds[i] >= angles[j]) {
                        
                        map[angles[j] - 1][bonds[i] - 1].push_back(get_mutual(TYPE_BA,i,j,MUTUAL_ARGS));
                
                    }
                    else {
                        
                        map[bonds[i] - 1][angles[j] - 1].push_back(get_mutual(TYPE_BA,i,j,MUTUAL_ARGS));
                        
                    }
                }
            }
        }
        
        for(int i = 0; i < nDihedrals + 2; i++) { //bond-dihedral pairs
            for(int j = 0; j < nDihedrals; j++) {
                if (bonds[i] != 0 && dihedrals[j] != 0) {
                    if (bonds[i] >= dihedrals[j]) {
                        
                        map[dihedrals[j] - 1][bonds[i] - 1].push_back(get_mutual(TYPE_BD,i,j,MUTUAL_ARGS));
                
                    }
                    else {
                        
                        map[bonds[i] - 1][dihedrals[j] - 1].push_back(get_mutual(TYPE_BD,i,j,MUTUAL_ARGS));
                        
                    }
                }
            }
        }
        
        for(int i = 0; i < nDihedrals; i++) { //angle-angle pairs
            for(int j = i + 1; j < nDihedrals + 1; j++) {
                if (angles[i] != 0 && angles[j] != 0) {
                    if (angles[i] >= angles[j]) {
                        
                        map[angles[j] - 1][angles[i] - 1].push_back(get_mutual(TYPE_AA,i,j,MUTUAL_ARGS));
                
                    }
                    else {
                        
                        map[angles[i] - 1][angles[j] - 1].push_back(get_mutual(TYPE_AA,i,j,MUTUAL_ARGS));
                        
                    }
                }
            }
        }
        
        for(int i = 0; i < nDihedrals + 1; i++) { //angle-dihedral pairs
            for(int j = 0; j < nDihedrals; j++) {
                if (angles[i] != 0 && dihedrals[j] != 0) {
                    if (angles[i] >= dihedrals[j]) {
                        
                        map[dihedrals[j] - 1][angles[i] - 1].push_back(get_mutual(TYPE_AD,i,j,MUTUAL_ARGS));
                
                    }
                    else {
                        
                        map[angles[i] - 1][dihedrals[j] - 1].push_back(get_mutual(TYPE_AD,i,j,MUTUAL_ARGS));
                        
                    }
                }
            }
        }
        
        for(int i = 0; i < nDihedrals - 1; i++) { //dihedral-dihedral pairs
            for(int j = i + 1; j < nDihedrals; j++) {
                if (dihedrals[i] != 0 && dihedrals[j] != 0) {
                    if (dihedrals[i] >= dihedrals[j]) {
                        
                        map[dihedrals[j] - 1][dihedrals[i] - 1].push_back(get_mutual(TYPE_DD,i,j,MUTUAL_ARGS));
                        map_with_dihedrals[dihedrals[j] - 1][dihedrals[i] - 1].push_back(get_mutual(TYPE_DD,i,j,MUTUAL_ARGS));
                        
                    }
                    else {
                        
                        map[dihedrals[i] - 1][dihedrals[j] - 1].push_back(get_mutual(TYPE_DD,i,j,MUTUAL_ARGS));
                        map_with_dihedrals[dihedrals[i] - 1][dihedrals[j] - 1].push_back(get_mutual(TYPE_DD,i,j,MUTUAL_ARGS));
                        
                    }
                }
            }
        }        
    }
    else {
        cerr<<"ERROR: COULD NOT OPEN FILE !\n\nUSAGE:\n"<<" input.par"<<endl;
        struct X tmp;
        return tmp;
    }
    infile.close();
    
    delete[]bondsEntropy1D;
    delete[]anglesEntropy1D;
    delete[]dihedralsEntropy1D;
    
    delete[]bbEntropy;
    delete[]aaEntropy;
    delete[]ddEntropy;
    delete[]baEntropy;
    delete[]bdEntropy;
    delete[]adEntropy;
    
    vector <vector <vector <double>>> m(map);
    vector <vector <vector <double>>> mwd(map_with_dihedrals);
    unsigned int NR = 0 + NResidues;
    vector <string> n(names);
    vector <int> rn(real_numbers);
    
    struct X tmp;
    tmp.map = m;
    tmp.map_with_dihedrals = mwd;
    tmp.NResidues = NR;
    tmp.names = n;
    tmp.real_numbers = rn;
    
    return tmp;
}


int main(int argc, char* argv[]) {
    
    Arg_Parser arg_parser(argc, argv);
  
    if (!( arg_parser.exists( string("-f1") ) && arg_parser.exists( string("-f2") ) && arg_parser.exists( string("-n") ) && (argc==7) )) {
        cerr<<"USAGE:\n"<<argv[0]<<" -f1 input1.par -f2 input2.par -n name"<<endl;
        return 1;
    }
    
    if ( strcmp( arg_parser.get_ext( arg_parser.get("-f1") ) , "par") ) {
    // check for the extensions of the input file
    cerr<<"USAGE:\n"<<argv[0]<<" -f1 input1.par -f2 input2.par -n name"<<endl;
    exit(EXIT_FAILURE);
    }
   
   if ( strcmp( arg_parser.get_ext( arg_parser.get("-f2") ) , "par") ) {
    // check for the extensions of the input file
    cerr<<"USAGE:\n"<<argv[0]<<" -f1 input1.par -f2 input2.par -n name"<<endl;
    exit(EXIT_FAILURE);
    }
    X x1 = f(arg_parser.get("-f1"));
    vector <vector <vector <double>>> map1 = x1.map;
    vector <vector <vector <double>>> map1_with_dihedrals = x1.map_with_dihedrals;
    unsigned int NResidues = x1.NResidues;
    vector <string> names = x1.names;
    vector <int> real_numbers = x1.real_numbers;
    int size = NResidues;
    
    X x2 = f(arg_parser.get("-f2"));
    vector <vector <vector <double>>> map2 = x2.map;
    vector <vector <vector <double>>> map2_with_dihedrals = x2.map_with_dihedrals;
    
    vector <vector <double>> first(NResidues, std::vector<double>(NResidues, 0)); // Xa(t)Xb(t+tau)
    vector <vector <double>> second(NResidues, std::vector<double>(NResidues, 0)); // Xa(t)Xb(t)
    vector <vector <double>> third(NResidues, std::vector<double>(NResidues, 0)); //Xa(t)
    vector <vector <double>> first_d(NResidues, std::vector<double>(NResidues, 0)); // Xa(t)Xb(t+tau)
    vector <vector <double>> second_d(NResidues, std::vector<double>(NResidues, 0)); // Xa(t)Xb(t)
    vector <vector <double>> third_d(NResidues, std::vector<double>(NResidues, 0)); //Xa(t)
    
    
    for (int i = 0; i < size; i++) {
        for (int j = i; j < size; j++) {
            for (int k = 0; k < int(map1[i][j].size()); k++) {
                third[i][j] += map1[i][j][k];
                
            }
            third[i][j] = third[i][j]/map1[i][j].size();
            third[j][i] = third[i][j];
            
            for (int k = 0; k < int(map1_with_dihedrals[i][j].size()); k++) {
                third_d[i][j] += map1_with_dihedrals[i][j][k];
                
            }
            third_d[i][j] = third_d[i][j]/map1_with_dihedrals[i][j].size();
            third_d[j][i] = third_d[i][j];
            
            
            double sum = 0;
            for (int k = 0; k < int(map1[i][j].size()); k++) {
                for (int l = 0; l < int(map1[i][j].size()); l++) {
                    sum += map1[i][j][k]*map1[i][j][l];
                }
            }
            sum = sum / map1[i][j].size() / map1[i][j].size();
            
            second[i][j] = sum;
            second[j][i] = second[i][j];
            
            sum = 0;
            for (int k = 0; k < int(map1_with_dihedrals[i][j].size()); k++) {
                for (int l = 0; l < int(map1_with_dihedrals[i][j].size()); l++) {
                    sum += map1_with_dihedrals[i][j][k]*map1_with_dihedrals[i][j][l];
                }
            }
            sum = sum / map1_with_dihedrals[i][j].size() / map1_with_dihedrals[i][j].size();
            
            second_d[i][j] = sum;
            second_d[j][i] = second_d[i][j];
            
            
            sum = 0;
            for (int k = 0; k < int(map1[i][j].size()); k++) {
                for (int l = 0; l < int(map2[i][j].size()); l++) {
                    sum += map1[i][j][k]*map2[i][j][l];
                }
            }
            sum = sum / map1[i][j].size() / map2[i][j].size();
            
            first[i][j] = sum;
            first[j][i] = first[i][j];
            
            sum = 0;
            for (int k = 0; k < int(map1_with_dihedrals[i][j].size()); k++) {
                for (int l = 0; l < int(map2_with_dihedrals[i][j].size()); l++) {
                    sum += map1_with_dihedrals[i][j][k]*map2_with_dihedrals[i][j][l];
                }
            }
            sum = sum / map1_with_dihedrals[i][j].size() / map2_with_dihedrals[i][j].size();
            
            first_d[i][j] = sum;
            first_d[j][i] = first_d[i][j];     
        }
    }

    save_map(NResidues, first, second, third, first_d, second_d, third_d, names, real_numbers, arg_parser.get("-n"));
    return 0;
}
