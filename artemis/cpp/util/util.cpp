#include "util.h"
#include <algorithm>
#include <iostream>

using namespace std;


unsigned char get_dof_type_from_id(unsigned int dof_id,
                                   unsigned int n_dihedrals) {
  if (dof_id < n_dihedrals + 2)
    return TYPE_B;
  if (dof_id < 2 * n_dihedrals + 3)
    return TYPE_A;
  return TYPE_D;
}

unsigned int get_min_id_for_type(unsigned char type, unsigned int n_dihedrals) {
  switch (type) {
  case TYPE_B:
    return 0;
  case TYPE_A:
    return n_dihedrals + 2;
  case TYPE_D:
    return 2 * n_dihedrals + 3;
  }
  return 42;
}

unsigned int get_max_id_for_type(unsigned char type, unsigned int n_dihedrals) {
  switch (type) {
  case TYPE_B:
    return n_dihedrals + 1;
  case TYPE_A:
    return 2 * n_dihedrals + 2;
  case TYPE_D:
    return 3 * n_dihedrals + 2;
  }
  return 42;
}

Dof get_dof_from_global_id(unsigned int id, unsigned int n_dihedrals){
    Dof tmp_dof;
    
    if (id < n_dihedrals + 2){
        tmp_dof.type = TYPE_B;
        tmp_dof.id = id;
        return tmp_dof;
    }
    
    if (id < 2 * n_dihedrals + 3){
        tmp_dof.type =  TYPE_A;
        tmp_dof.id = id - n_dihedrals - 2;
        return tmp_dof;
    }
  
    tmp_dof.type = TYPE_D;
    tmp_dof.id = id - 2 * n_dihedrals - 3;
    
    return tmp_dof;
}
