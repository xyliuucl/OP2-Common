//
// auto-generated by op2.m on 19-Oct-2012 16:21:08
//

// user function

__device__
#include "adt_calc.h"


// CUDA kernel function

__global__ void op_cuda_adt_calc(
  double *ind_arg0,
  int   *ind_map,
  short *arg_map,
  double *arg4,
  double *arg5,
  int   *ind_arg_sizes,
  int   *ind_arg_offs,
  int    block_offset,
  int   *blkmap,
  int   *offset,
  int   *nelems,
  int   *ncolors,
  int   *colors,
  int   nblocks,
  int   set_size) {

  double *arg0_vec[4];

  __shared__ int   *ind_arg0_map, ind_arg0_size;
  __shared__ double *ind_arg0_s;
  __shared__ int    nelem, offset_b;

  extern __shared__ char shared[];

  if (blockIdx.x+blockIdx.y*gridDim.x >= nblocks) return;
  if (threadIdx.x==0) {

    // get sizes and shift pointers and direct-mapped data

    int blockId = blkmap[blockIdx.x + blockIdx.y*gridDim.x  + block_offset];

    nelem    = nelems[blockId];
    offset_b = offset[blockId];

    ind_arg0_size = ind_arg_sizes[0+blockId*1];

    ind_arg0_map = &ind_map[0*set_size] + ind_arg_offs[0+blockId*1];

    // set shared memory pointers

    int nbytes = 0;
    ind_arg0_s = (double *) &shared[nbytes];
  }

  __syncthreads(); // make sure all of above completed

  // copy indirect datasets into shared memory or zero increment

  for (int n=threadIdx.x; n<ind_arg0_size*2; n+=blockDim.x)
    ind_arg0_s[n] = ind_arg0[n%2+ind_arg0_map[n/2]*2];

  __syncthreads();

  // process set elements

  for (int n=threadIdx.x; n<nelem; n+=blockDim.x) {

      arg0_vec[0] = ind_arg0_s+arg_map[0*set_size+n+offset_b]*2;
      arg0_vec[1] = ind_arg0_s+arg_map[1*set_size+n+offset_b]*2;
      arg0_vec[2] = ind_arg0_s+arg_map[2*set_size+n+offset_b]*2;
      arg0_vec[3] = ind_arg0_s+arg_map[3*set_size+n+offset_b]*2;

      // user-supplied kernel call


      adt_calc(  arg0_vec,
                 arg4+(n+offset_b)*4,
                 arg5+(n+offset_b)*1 );
  }

}


// host stub function

void op_par_loop_adt_calc(char const *name, op_set set,
  op_arg arg0,
  op_arg arg4,
  op_arg arg5 ){


  int    nargs   = 6;
  op_arg args[6];

  arg0.idx = 0;
  args[0] = arg0;
  for (int v = 1; v < 4; v++) {
    args[0 + v] = op_arg_dat(arg0.dat, v, arg0.map, 2, "double", OP_READ);
  }
  args[4] = arg4;
  args[5] = arg5;

  int    ninds   = 1;
  int    inds[6] = {0,0,0,0,-1,-1};

  if (OP_diags>2) {
    printf(" kernel routine with indirection: adt_calc\n");
  }

  // get plan

  #ifdef OP_PART_SIZE_1
    int part_size = OP_PART_SIZE_1;
  #else
    int part_size = OP_part_size;
  #endif

  int set_size = op_mpi_halo_exchanges(set, nargs, args);

  // initialise timers

  double cpu_t1, cpu_t2, wall_t1=0, wall_t2=0;
  op_timing_realloc(1);
  OP_kernels[1].name      = name;
  OP_kernels[1].count    += 1;

  if (set->size >0) {

    op_plan *Plan = op_plan_get(name,set,part_size,nargs,args,ninds,inds);

    op_timers_core(&cpu_t1, &wall_t1);

    // execute plan

    int block_offset = 0;

    for (int col=0; col < Plan->ncolors; col++) {

      if (col==Plan->ncolors_core) op_mpi_wait_all(nargs,args);

    #ifdef OP_BLOCK_SIZE_1
      int nthread = OP_BLOCK_SIZE_1;
    #else
      int nthread = OP_block_size;
    #endif

      dim3 nblocks = dim3(Plan->ncolblk[col] >= (1<<16) ? 65535 : Plan->ncolblk[col],
                      Plan->ncolblk[col] >= (1<<16) ? (Plan->ncolblk[col]-1)/65535+1: 1, 1);
      if (Plan->ncolblk[col] > 0) {
        int nshared = Plan->nsharedCol[col];
        op_cuda_adt_calc<<<nblocks,nthread,nshared>>>(
           (double *)arg0.data_d,
           Plan->ind_map,
           Plan->loc_map,
           (double *)arg4.data_d,
           (double *)arg5.data_d,
           Plan->ind_sizes,
           Plan->ind_offs,
           block_offset,
           Plan->blkmap,
           Plan->offset,
           Plan->nelems,
           Plan->nthrcol,
           Plan->thrcol,
           Plan->ncolblk[col],
           set_size);

        cutilSafeCall(cudaDeviceSynchronize());
        cutilCheckMsg("op_cuda_adt_calc execution failed\n");
      }

      block_offset += Plan->ncolblk[col];
    }

    op_timing_realloc(1);
    OP_kernels[1].transfer  += Plan->transfer;
    OP_kernels[1].transfer2 += Plan->transfer2;

  }


  op_mpi_set_dirtybit(nargs, args);

  // update kernel record

  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[1].time     += wall_t2 - wall_t1;
}
