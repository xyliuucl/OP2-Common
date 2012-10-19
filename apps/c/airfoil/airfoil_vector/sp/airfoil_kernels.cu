//
// auto-generated by op2.m on 19-Oct-2012 16:21:08
//

// header

#include "op_lib_cpp.h"

#include "op_cuda_rt_support.h"
#include "op_cuda_reduction.h"
// global constants

#ifndef MAX_CONST_SIZE
#define MAX_CONST_SIZE 128
#endif

__constant__ float gam;
__constant__ float gm1;
__constant__ float cfl;
__constant__ float eps;
__constant__ float mach;
__constant__ float alpha;
__constant__ float qinf[4];

void op_decl_const_gam(int dim, char const *type, float *dat){
  cutilSafeCall(cudaMemcpyToSymbol(gam , dat, dim*sizeof(float)));
}

void op_decl_const_gm1(int dim, char const *type, float *dat){
  cutilSafeCall(cudaMemcpyToSymbol(gm1 , dat, dim*sizeof(float)));
}

void op_decl_const_cfl(int dim, char const *type, float *dat){
  cutilSafeCall(cudaMemcpyToSymbol(cfl , dat, dim*sizeof(float)));
}

void op_decl_const_eps(int dim, char const *type, float *dat){
  cutilSafeCall(cudaMemcpyToSymbol(eps , dat, dim*sizeof(float)));
}

void op_decl_const_mach(int dim, char const *type, float *dat){
  cutilSafeCall(cudaMemcpyToSymbol(mach , dat, dim*sizeof(float)));
}

void op_decl_const_alpha(int dim, char const *type, float *dat){
  cutilSafeCall(cudaMemcpyToSymbol(alpha , dat, dim*sizeof(float)));
}

void op_decl_const_qinf(int dim, char const *type, float *dat){
  cutilSafeCall(cudaMemcpyToSymbol(qinf , dat, dim*sizeof(float)));
}

// user kernel files

#include "save_soln_kernel.cu"
#include "adt_calc_kernel.cu"
#include "res_calc_kernel.cu"
#include "bres_calc_kernel.cu"
#include "update_kernel.cu"
