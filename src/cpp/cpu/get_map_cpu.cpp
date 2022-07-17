//    A program to extract the entropy and mutual information terms from the binary output of the PARENT program
//    Copyright (C) 2015  Markus Fleck (member of the laboratory of Bojan Zagrovic, University of Vienna)
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License  version 3
//    as published by the Free Software Foundation.

//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.

//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.





//    A scientific publication about this program has been released in the Journal of Chemical Theory and Computation:
//		"PARENT: A Parallel Software Suite for the Calculation of Configurational Entropy in Biomolecular Systems"
//		DOI: 10.1021/acs.jctc.5b01217

//   We kindly ask you to include this citation in works that publish
//   results generated using this program or any modifications of it.


#define MUTUAL_ARGS nDihedrals,bondsEntropy1D,anglesEntropy1D,dihedralsEntropy1D,bbEntropy,baEntropy,bdEntropy,aaEntropy,adEntropy,ddEntropy


#include <iostream>
#include <fstream>
#include <vector>
#include <string.h>
#include <cstdlib>
#include <algorithm>


using namespace std;


#include "util/io.h"
#include "util/util.h"
#include "util/Arg_Parser.h"
#include "../save_map.cpp"


int main(int argc, char* argv[]) {
    
    Arg_Parser arg_parser(argc, argv);
  
    if( !( arg_parser.exists( string("-f") ) && (argc==3) ) ){
        cerr<<"USAGE:\n"<<argv[0]<<" -f input.par"<<endl;
        return 1;
    }
    
    if ( strcmp( arg_parser.get_ext( arg_parser.get("-f") ) , "par") ) {
    // check for the extensions of the input file
    cerr<<"USAGE:\n"<<argv[0]<<" -f input.par"<<endl;
    exit(EXIT_FAILURE);
    }

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

    vector <vector <double>> map;
    vector <unsigned int> bonds;
    vector <unsigned int> angles;
    vector <unsigned int> dihedrals;
    vector <string> names;



    ifstream infile(arg_parser.get("-f"), ios::binary | ios::in);//open the .par file
    if(infile.is_open()) {
        if(read_PAR_header(&infile,&nDihedrals,&double_prec,&numFrames,&dihedrals_top, &masses, &version, &bDens, &aDens, &dDens, &bDens1D, &aDens1D, &dDens1D,&residues,&residueNumbers,&atomNames,&belongsToMolecule)!=0) { //and read the header information
            cerr<<"ERROR READING HEADER OF FILE "<<argv[0]<<" !"<<endl;
            return 1;
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
        
        int current = 1;
        int total = 1;
        for (int i = 0; i < size; i ++) {
            if (residueNumbers[i] == total)
                residueNumbers[i] = current;
            else {
                total = residueNumbers[i];
                current++;
                residueNumbers[i] = current;
            }
        }
        
        unsigned int NResidues = names.size();
        
        
        for (unsigned int i = 0; i < NResidues; i++) {
            vector <double> line;
            for (unsigned int k = 0; k < NResidues; k++) {
                line.push_back(0);
            }
            map.push_back(line);
        }

        if (read_PAR_body(&infile,nDihedrals,&bondsEntropy1D, &anglesEntropy1D, &dihedralsEntropy1D, &bbEntropy, &baEntropy, &bdEntropy, &aaEntropy, &adEntropy, &ddEntropy)!=0) {
            cerr<<"ERROR READING FILE "<<argv[0]<<" !"<<endl;
            return 1;
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
                        
                        map[bonds[j] - 1][bonds[i] - 1] += get_mutual(TYPE_BB,i,j,MUTUAL_ARGS);
                
                    }
                    else {
                        
                        map[bonds[i] - 1][bonds[j] - 1] += get_mutual(TYPE_BB,i,j,MUTUAL_ARGS);
                        
                    }
                }
            }
        }
        
        for(int i = 0; i < nDihedrals + 2; i++) { //same for all bond-angle pairs
            for(int j = 0; j < nDihedrals + 1; j++) {
                if (bonds[i] != 0 && angles[j] != 0) {
                    if (bonds[i] >= angles[j]) {
                        
                        map[angles[j] - 1][bonds[i] - 1] += get_mutual(TYPE_BA,i,j,MUTUAL_ARGS);
                
                    }
                    else {
                        
                        map[bonds[i] - 1][angles[j] - 1] += get_mutual(TYPE_BA,i,j,MUTUAL_ARGS);
                        
                    }
                }
            }
        }
        
        for(int i = 0; i < nDihedrals + 2; i++) { //bond-dihedral pairs
            for(int j = 0; j < nDihedrals; j++) {
                if (bonds[i] != 0 && dihedrals[j] != 0) {
                    if (bonds[i] >= dihedrals[j]) {
                        
                        map[dihedrals[j] - 1][bonds[i] - 1] += get_mutual(TYPE_BD,i,j,MUTUAL_ARGS);
                
                    }
                    else {
                        
                        map[bonds[i] - 1][dihedrals[j] - 1] += get_mutual(TYPE_BD,i,j,MUTUAL_ARGS);
                        
                    }
                }
            }
        }
        
        for(int i = 0; i < nDihedrals; i++) { //angle-angle pairs
            for(int j = i + 1; j < nDihedrals + 1; j++) {
                if (angles[i] != 0 && angles[j] != 0) {
                    if (angles[i] >= angles[j]) {
                        
                        map[angles[j] - 1][angles[i] - 1] += get_mutual(TYPE_AA,i,j,MUTUAL_ARGS);
                
                    }
                    else {
                        
                        map[angles[i] - 1][angles[j] - 1] += get_mutual(TYPE_AA,i,j,MUTUAL_ARGS);
                        
                    }
                }
            }
        }
        
        for(int i = 0; i < nDihedrals + 1; i++) { //angle-dihedral pairs
            for(int j = 0; j < nDihedrals; j++) {
                if (angles[i] != 0 && dihedrals[j] != 0) {
                    if (angles[i] >= dihedrals[j]) {
                        
                        map[dihedrals[j] - 1][angles[i] - 1] += get_mutual(TYPE_AD,i,j,MUTUAL_ARGS);
                
                    }
                    else {
                        
                        map[angles[i] - 1][dihedrals[j] - 1] += get_mutual(TYPE_AD,i,j,MUTUAL_ARGS);
                        
                    }
                }
            }
        }
        
        for(int i = 0; i < nDihedrals - 1; i++) { //dihedral-dihedral pairs
            for(int j = i + 1; j < nDihedrals; j++) {
                if (dihedrals[i] != 0 && dihedrals[j] != 0) {
                    if (dihedrals[i] >= dihedrals[j]) {
                        
                        map[dihedrals[j] - 1][dihedrals[i] - 1] += get_mutual(TYPE_DD,i,j,MUTUAL_ARGS);
                        
                    }
                    else {
                        
                        map[dihedrals[i] - 1][dihedrals[j] - 1] += get_mutual(TYPE_DD,i,j,MUTUAL_ARGS);
                        
                    }
                }
            }
        }
        
        for (unsigned int i = 0; i < NResidues; i++) {
            for (unsigned int k = i; k < NResidues; k++) {
                map[k][i] = map[i][k];
            }
        }
        
        save_map(NResidues, map, names);
    }
    else {
        cerr<<"ERROR: COULD NOT OPEN FILE !\n\nUSAGE:\n"<<argv[0]<<" input.par"<<endl;
        return 1;
    }
    infile.close();

    return 0;
}
