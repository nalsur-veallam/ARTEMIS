#include "util.h"

char *getCmdOption(char **begin, char **end, const string &option) {
  char **itr = find(begin, end, option);
  if (itr != end && ++itr != end) {
    return *itr;
  }
  return 0;
}

bool cmdOptionExists(char **begin, char **end, const string &option) {
  return find(begin, end, option) != end;
}

double get_mutual(int type,int index1, int index2,int nDihedrals, double* bondsEntropy1D, double* anglesEntropy1D, double* dihedralsEntropy1D, double* bbEntropy, double* baEntropy, double* bdEntropy, double* aaEntropy, double* adEntropy, double* ddEntropy) {
    int smaller,bigger,index;
    int nBonds=nDihedrals+2;
    int nAngles=nDihedrals+1;

    if((type==TYPE_BB)&&(index1!=index2)) { // for mutual information between two different bonds
        smaller=index1<index2?index1:index2;
        bigger=index1<index2?index2:index1;
        index=(nBonds-smaller)*(nBonds-smaller-1)/2+smaller-bigger;//the 2D-entropies bonds-bonds (also angles-angles and dihedrals-dihedrals) were stored in reverse order, as documented in "Parent.cpp"
        return bondsEntropy1D[smaller]+bondsEntropy1D[bigger]-bbEntropy[index];
    }
    if(type==TYPE_BA) { // for mutual information between a bond and an angle
        return bondsEntropy1D[index1]+anglesEntropy1D[index2]-baEntropy[index1*nAngles+index2];
    }
    if(type==TYPE_BD) { // for mutual information between a bond and a dihedral
        return bondsEntropy1D[index1]+dihedralsEntropy1D[index2]-bdEntropy[index1*nDihedrals+index2];
    }
    if((type==TYPE_AA)&&(index1!=index2)) { // for mutual information between two different angles
        smaller=index1<index2?index1:index2;
        bigger=index1<index2?index2:index1;
        index=(nAngles-smaller)*(nAngles-smaller-1)/2+smaller-bigger;//the 2D-entropies angles-angles were stored in reverse order, as documented in "Parent.cpp"
        return anglesEntropy1D[smaller]+anglesEntropy1D[bigger]-aaEntropy[index];
    }
    if(type==TYPE_AD) { // for mutual information between a bond and an angle
        return anglesEntropy1D[index1]+dihedralsEntropy1D[index2]-adEntropy[index1*nDihedrals+index2];
    }
    if((type==TYPE_DD)&&(index1!=index2)) { // for mutual information between two different angles
        smaller=index1<index2?index1:index2;
        bigger=index1<index2?index2:index1;
        index=(nDihedrals-smaller)*(nDihedrals-smaller-1)/2+smaller-bigger;//the 2D-entropies dihedrals-dihedrals were stored in reverse order, as documented in "Parent.cpp"
        return dihedralsEntropy1D[smaller]+dihedralsEntropy1D[bigger]-ddEntropy[index];
    }
    cerr<<"WARNING: REQUEST FOR MUTUAL INFORMATION TYPE "<<type<<" WITH INDICES "<<index1<<" AND "<<index2<<"YIELDED NO RESULT."<<endl;
    return 0;
}

