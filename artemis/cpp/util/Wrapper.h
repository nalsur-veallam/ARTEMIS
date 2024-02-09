#ifndef WRAPPER_H
#define WRAPPER_H


#define MODE_TOTAL 0
#define MODE_AVER 1
#define MODE_MAX 2

#include <iostream>
#include <string>
#include <cmath>

#include "Residue_Representation.h"
#include "data.h"

#include "Entropy_Matrix.h"
#include  "types.h"
#include "My_Error.cpp"
#include <vector>


using namespace std;

class Wrapper {
public:

    Wrapper(char const * file1, char const * file2, double Cm1, double Cm2) {
        rep1 = new Residue_Representation(file1, false, MODE_TOTAL);
        mat1 = rep1->getEntropy_Matrix();

        rep2 = new Residue_Representation(file2, false, MODE_TOTAL);
        mat2 = rep2->getEntropy_Matrix();

        for (unsigned int i = 0; i < NResidues; i++) {
            real_numbers.push_back(rep1->getResidueNumber(i + 1));
        }

        NResidues = rep1->getNResidues();

        computeMap(Cm1, Cm2);

    };
    virtual ~Wrapper() {};

    void computeMap(double Cm1, double Cm2) {

        for (unsigned int resid1 = 1; resid1 <= NResidues; resid1++) {

            vector <double> mie;
            names.push_back(rep1->getResidueName(resid1));

            for (unsigned int resid2 = 1; resid2 <= NResidues; resid2++) {

                vector <vector< int > >dofs1;
                dofs1.push_back(rep1->getBonds(resid1));
                dofs1.push_back(rep1->getAngles(resid1));
                dofs1.push_back(rep1->getDihedrals(resid1));

                vector <vector< int > >dofs2;
                dofs2.push_back(rep1->getBonds(resid2));
                dofs2.push_back(rep1->getAngles(resid2));
                dofs2.push_back(rep1->getDihedrals(resid2));


                double mutual = 0;

                for (unsigned char type1 = 0; type1 < 3; type1++){ // for bonds, angles and dihedrals of the first member of the dof pair
                    for(unsigned int idx1 = 0; idx1 < dofs1[type1].size(); idx1++){ // and all dofs of the current type of the first member of the dof pair
                        for (unsigned char type2 = 0; type2 < 3; type2++){ // and all "later" types for the second member of the dof pair
                            for(unsigned int idx2 = 0; idx2 < dofs2[type2].size(); idx2++ ){ // and all "later" dofs for the second member of the dof pair

                                mutual += Cm1 * (mat1->getMutual(type1, type2, dofs1[type1][idx1], dofs2[type2][idx2]));
                                mutual +=  Cm2 * (mat2->getMutual(type1, type2, dofs1[type1][idx1], dofs2[type2][idx2]));

                            }
                        }
                    }
                }

                if (mutual < 0) {mutual = 0;}
                mie.push_back(mutual);

            }

            map.push_back(mie);
        }
    }


//    virtual std::vector <std::vector <double>> getMap() = 0;
//    virtual std::vector <std::string> getNames() = 0;
    virtual std::vector <int> getNumbers() = 0;


    vector <int> real_numbers;
    vector <vector <double>> map;
    vector <string> names;


private:
    Entropy_Matrix* mat1;
    Entropy_Matrix* mat2;

    Residue_Representation* rep1;
    Residue_Representation* rep2;

    unsigned int NResidues;

};

class Interface : public Wrapper {
public:
    Interface(char const * file1, char const * file2, double Cm1, double Cm2) : Wrapper(file1, file2, Cm1, Cm2){};


//     std::vector <std::vector <double>> getMap() override {
//       return this->map;
//     }
//
//     std::vector <std::string> getNames() override {
//       return this->names;
//     }

    std::vector <int> getNumbers() override {
      return this->real_numbers;
    }


};

#endif
