#ifndef UTIL_H
#define UTIL_H

#include "types.h"
#include <string>



#ifdef __CUDACC__
    #include "../util/classes/My_Error.cpp"
    #define gpuErrchk(ans)                                                         \
      { gpuAssert((ans), __FILE__, __LINE__); }
    inline void gpuAssert(cudaError_t return_code, const char *file, int line) {
      if (return_code != cudaSuccess) {
        My_Error my_error(std::string("GPUassert: ") + std::string(cudaGetErrorString(return_code)) + std::string(" ") + std::string(file) + std::string(" ") + std::to_string(line) + std::string("\n"));
        throw my_error;
      }
    }
#endif

struct Dof {
    unsigned char type;
    unsigned int id;

};

unsigned char get_dof_type_from_id(unsigned int dof_id,
                                   unsigned int n_dihedrals);
unsigned int get_min_id_for_type(unsigned char type, unsigned int n_dihedrals);
unsigned int get_max_id_for_type(unsigned char type, unsigned int n_dihedrals);
Dof get_dof_from_global_id(unsigned int id, unsigned int n_diherdals);

#endif
