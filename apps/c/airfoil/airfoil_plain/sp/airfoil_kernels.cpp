//
// auto-generated by op2.m on 19-Oct-2012 16:21:06
//

// header

#include "op_lib_cpp.h"

// global constants

extern float gam;
extern float gm1;
extern float cfl;
extern float eps;
extern float mach;
extern float alpha;
extern float qinf[4];

void op_decl_const_gam(int dim, char const *type, float *dat){};

void op_decl_const_gm1(int dim, char const *type, float *dat){};

void op_decl_const_cfl(int dim, char const *type, float *dat){};

void op_decl_const_eps(int dim, char const *type, float *dat){};

void op_decl_const_mach(int dim, char const *type, float *dat){};

void op_decl_const_alpha(int dim, char const *type, float *dat){};

void op_decl_const_qinf(int dim, char const *type, float *dat){};

// user kernel files

#include "save_soln_kernel.cpp"
#include "adt_calc_kernel.cpp"
#include "res_calc_kernel.cpp"
#include "bres_calc_kernel.cpp"
#include "update_kernel.cpp"
