#ifndef RESIDUE_REPRESENTATION_H
#define RESIDUE_REPRESENTATION_H


#define MODE_TOTAL 0
#define MODE_AVER 1
#define MODE_MAX 2



#include "Entropy_Matrix.h"
#include  "types.h"
#include "My_Error.cpp"
#include <vector>


class Residue_Representation {
    public:
      Residue_Representation(char const * infileInput,bool include_full=false, int mode=MODE_TOTAL);
      ~Residue_Representation();
      void calculate_matrix();
    
      std::string getResidueName(unsigned int residueIndex); //gives the name of the residue with index "residueIndex" (indexing starts at 1) 
      int getResidueNumber(unsigned int residueIndex); //gives the number (according to the used topology file) of the residue with index "residueIndex" (indexing starts at 1) 
      std::string getMoleculeName(unsigned int residueIndex); //gives the name of the molecule the residue with index "residueIndex" belongs to (indexing starts at 1) 
    
      unsigned int getNResidues(); //returns the number of residues
      unsigned int getNBonds(unsigned int residueIndex); //returns the number of bonds in the residue with index "residueIndex" (indexing starts at 1) 
      unsigned int getNAngles(unsigned int residueIndex); //returns the number of angles in the residue with index "residueIndex" (indexing starts at 1) 
      unsigned int getNDihedrals(unsigned int residueIndex); //returns the number of dihedrals in the residue with index "residueIndex" (indexing starts at 1)
    
      std::vector <int> getAtoms(unsigned int residueIndex);
      std::vector <int> getBonds(unsigned int residueIndex);
      std::vector <int> getAngles(unsigned int residueIndex);
      std::vector <int> getDihedrals(unsigned int residueIndex);
    
    
      std::string getAtomName(unsigned int atomNumber); // get the name of the atom with the according number (atomnumbers start at 1)

      double getMutual(unsigned int residueIndex1,unsigned int residueIndex2);//returns the total mutual information between the according residues (indexing starts at 1)
      void setMutual(unsigned int residueIndex1,unsigned int residueIndex2, double value);//sets the total mutual information between the according residues (indexing starts at 1)
    
      Entropy_Matrix* getEntropy_Matrix();
    
		
    private:
      int calcBonds;
      int calcAngles;
      int calcDihedrals;
    
      int mode;
    
      bool matrix_calculated = false;
    
      Entropy_Matrix* mat;
      double* mutualArray;
    
      unsigned int nResidues;
		
      std::vector < std::vector <int> > groups;
			
      std::vector <std::string> residueNames;
      std::vector <int> residueNumbers;
      std::vector <std::string> moleculeNames;
    
      std::vector < std::vector <int> > bondIndices;
      std::vector < std::vector <int> > angleIndices;
      std::vector < std::vector <int> > dihedralIndices;
    
      std::vector <int> nBondsVec;
      std::vector <int> nAnglesVec;
      std::vector <int> nDihedralsVec;
};

#endif
