#define TYPE_B 0
#define TYPE_A 1
#define TYPE_D 2

#define TYPE_BB 0
#define TYPE_BA 1
#define TYPE_BD 2
#define TYPE_AA 3
#define TYPE_AD 4
#define TYPE_DD 5


#include <iostream>
#include <string>
#include <algorithm>

using namespace std;

double get_mutual(int type,int index1, int index2,int nDihedrals, double* bondsEntropy1D, double* anglesEntropy1D, double* dihedralsEntropy1D, double* bbEntropy, double* baEntropy, double* bdEntropy, double* aaEntropy, double* adEntropy, double* ddEntropy);
char *getCmdOption(char **begin, char **end, const std::string &option);
bool cmdOptionExists(char **begin, char **end, const std::string &option);


