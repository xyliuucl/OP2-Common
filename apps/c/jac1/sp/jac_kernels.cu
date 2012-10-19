//
// auto-generated by op2.m on 19-Oct-2012 16:21:11
//

// header

#include "op_lib_cpp.h"

#include "op_cuda_rt_support.h"
#include "op_cuda_reduction.h"
// global constants

#ifndef MAX_CONST_SIZE
#define MAX_CONST_SIZE 128
#endif

__constant__ float alpha;

void op_decl_const_alpha(int dim, char const *type, float *dat){
  cutilSafeCall(cudaMemcpyToSymbol(alpha , dat, dim*sizeof(float)));
}

// user kernel files

#include "res_kernel.cu"
#include "update_kernel.cu"
