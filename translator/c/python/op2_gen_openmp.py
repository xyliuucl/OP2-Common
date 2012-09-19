##########################################################################
#
# OpenMP code generator
#
# This routine is called by op2 which parses the input files
#
# It produces a file xxx_kernel.cpp for each kernel,
# plus a master kernel file
#
##########################################################################

import re
import datetime
 
def comm(line):
  global file_text, FORTRAN, CPP  
  
  if len(line) == 0:
    file_text +='\n'
  elif FORTRAN:
    file_text +='!  '+line+'\n'
  elif CPP:
    file_text +='//  '+line+'\n'

def rep(line,m):
  global dims, idxs, typs, indtyps, inddims  
  if m < len(inddims):
    line = re.sub('INDDIM',str(inddims[m]),line)
    line = re.sub('INDTYP',str(indtyps[m]),line)
    
  line = re.sub('INDARG','ind_arg'+str(m),line)
  line = re.sub('DIM',str(dims[m]),line)
  line = re.sub('ARG','arg'+str(m),line)
  line = re.sub('TYP',typs[m],line)
  line = re.sub('IDX',str(int(idxs[m])),line)
  return line

def code(text):
  global file_text, FORTRAN, CPP, g_m
  global depth
  prefix = ' '*depth
  file_text += prefix+rep(text,g_m)+'\n'


def op2_gen_openmp(master, date, consts, kernels):

  global dims, idxs, typs, indtyps, inddims
  global FORTRAN, CPP, g_m, file_text, depth
  
  OP_ID   = 1;  OP_GBL   = 2;  OP_MAP = 3;

  OP_READ = 1;  OP_WRITE = 2;  OP_RW  = 3;
  OP_INC  = 4;  OP_MAX   = 5;  OP_MIN = 6;
  
  accsstring = ['OP_READ','OP_WRITE','OP_RW','OP_INC','OP_MAX','OP_MIN' ]
  
  any_soa = 0
  for nk in range (0,len(kernels)):
    any_soa = any_soa or sum(kernels[nk]['soaflags'])    

##########################################################################
#  create new kernel file
##########################################################################

  for nk in range (0,len(kernels)):
      
    name  = kernels[nk]['name']
    nargs = kernels[nk]['nargs']
    dims  = kernels[nk]['dims']
    maps  = kernels[nk]['maps']
    var   = kernels[nk]['var']
    typs  = kernels[nk]['typs']
    accs  = kernels[nk]['accs']
    idxs  = kernels[nk]['idxs']
    inds  = kernels[nk]['inds']
    soaflags = kernels[nk]['soaflags']

    ninds   = kernels[nk]['ninds']
    inddims = kernels[nk]['inddims']
    indaccs = kernels[nk]['indaccs']
    indtyps = kernels[nk]['indtyps']
    invinds = kernels[nk]['invinds']
    
    vec =  [m for m in range(0,nargs) if int(idxs[m])<0 and maps[m] == OP_MAP]
    
    if len(vec) > 0:
      unique_args = [1];
      vec_counter = 1;
      vectorised = []
      new_dims = []
      new_maps = []
      new_vars = []
      new_typs = []
      new_accs = []
      new_idxs = []
      new_inds = []
      new_soaflags = []
      for m in range(0,nargs):
      	  if int(idxs[m])<0 and maps[m] == OP_MAP:
      	    if m > 0:
      	      unique_args = unique_args + [len(new_dims)+1]
      	    temp = [0]*(-1*int(idxs[m]))
      	    for i in range(0,-1*int(idxs[m])): 
      	      temp[i] = var[m]
      	    new_vars = new_vars+temp
      	    for i in range(0,-1*int(idxs[m])): 
      	      temp[i] = typs[m]
      	    new_typs = new_typs+temp
      	    for i in range(0,-1*int(idxs[m])): 
      	      temp[i] = dims[m]
      	    new_dims = new_dims+temp
      	    new_maps = new_maps+[maps[m]]*int(-1*int(idxs[m]))
      	    new_soaflags = new_soaflags+[0]*int(-1*int(idxs[m]))
      	    new_accs = new_accs+[accs[m]]*int(-1*int(idxs[m]))
            for i in range(0,-1*int(idxs[m])):
            	new_idxs = new_idxs+[i]
            new_inds = new_inds+[inds[m]]*int(-1*int(idxs[m]))  
            vectorised = vectorised + [vec_counter]*int(-1*int(idxs[m]))
            vec_counter = vec_counter + 1;
          else:
            if m > 0:
              unique_args = unique_args + [len(new_dims)+1]
            new_dims = new_dims+[dims[m]]
            new_maps = new_maps+[maps[m]]
            new_accs = new_accs+[int(accs[m])]
            new_soaflags = new_soaflags+[soaflags[m]]
            new_idxs = new_idxs+[int(idxs[m])]
            new_inds = new_inds+[inds[m]]
            new_vars = new_vars+[var[m]]
            new_typs = new_typs+[typs[m]]
            vectorised = vectorised+[0]
      dims = new_dims
      maps = new_maps
      accs = new_accs
      idxs = new_idxs
      inds = new_inds
      var = new_vars
      typs = new_typs
      soaflags = new_soaflags;
      nargs = len(vectorised);
      
      for i in range(1,ninds+1):
      	for index in range(0,len(inds)+1):
      	  if inds[index] == i:
      	    invinds[i-1] = index
      	    break      
    else:
    	vectorised = [0]*nargs
    	unique_args = range(1,nargs+1)

    cumulative_indirect_index = [-1]*nargs;
    j = 0;
    for i in range (0,nargs):
      if maps[i] == OP_MAP:
        cumulative_indirect_index[i] = j
        j = j + 1        
#
# set two logicals
#
    j = 0
    for i in range(0,nargs):
      if maps[i] == OP_MAP and accs[i] == OP_INC:
        j = i
    ind_inc = j > 0
    
    j = 0
    for i in range(0,nargs):
      if maps[i] == OP_GBL and accs[i] <> OP_READ:
        j = i
    reduct = j > 0

##########################################################################
#  start with OpenMP kernel function
##########################################################################
  
    FORTRAN = 0;
    CPP     = 1;
    g_m = 0;
    file_text = ''
    depth = 0
    
    comm('user function')
        
    if FORTRAN:
      code('include '+name+'.inc')
    elif CPP:
      code('#include "'+name+'.h"')
    
    comm('')
    comm('x86 kernel function')
    
    if FORTRAN:
      code('subroutine op_x86_'+name+'(')
    elif CPP:
      code('void op_x86_'+name+'(')
    
    depth = 2
    
    if ninds>0:
      if FORTRAN:
        code('integer(4)  blockIdx,')
      elif CPP:
        code('int  blockIdx,')
      
  
    for g_m in range(0,ninds):
      if FORTRAN:
        code('INDTYP *ind_ARG,')
      elif CPP:
        code('INDTYP *ind_ARG,')
    
    if ninds>0:
      if FORTRAN:
        code('int   *ind_map,')
        code('short *arg_map,')
      elif CPP:
        code('int   *ind_map,')
        code('short *arg_map,')
  
    for g_m in range (0,nargs):
      if maps[g_m]==OP_GBL and accs[g_m] == OP_READ:
        # declared const for performance
        if FORTRAN:
          code('const TYP *ARG,')  
        elif CPP:
          code('const TYP *ARG,')          
      elif maps[g_m]==OP_ID and ninds>0:
        if FORTRAN:
          code('ARG,')
        elif CPP:
          code('TYP  *ARG,')
      elif maps[g_m]==OP_GBL or maps[g_m]==OP_ID:
        if FORTRAN:
          code('ARG,')
        elif CPP:
          code('TYP *ARG,')
    
    if ninds>0:
      if FORTRAN:
        code('ind_arg_sizes,')
        code('ind_arg_offs, ')
        code('block_offset, ')
        code('blkmap,       ')
        code('offset,       ')
        code('nelems,       ')
        code('ncolors,      ')
        code('colors,       ')
        code('set_size) {    ')
        code('')
      elif CPP:
        code('int   *ind_arg_sizes,')
        code('int   *ind_arg_offs, ')
        code('int    block_offset, ')
        code('int   *blkmap,       ')
        code('int   *offset,       ')
        code('int   *nelems,       ')
        code('int   *ncolors,      ')
        code('int   *colors,       ')
        code('int   set_size) {    ')
        code('')          
    else:
      code('int   start,    ')
      code('int   finish ) {')
      
    for g_m in range (0,nargs):
      if maps[g_m]==OP_MAP and accs[g_m]==OP_INC:
        if FORTRAN:
          code('TYP  ARG_l[DIM];')  
        elif CPP:
          code('TYP  ARG_l[DIM];')     
        
    
    code('') 
    
    for m in range (1,ninds+1):
      g_m = m-1
      v = [int(inds[i]==m) for i in range(len(inds))]
      v_i = [vectorised[i] for i in range(len(inds)) if inds[i] == m] 
      if sum(v)>1 and sum(v_i)>0:  
        if indaccs[m-1] == OP_INC:
          ind = int(max([idxs[i] for i in range(len(inds)) if inds[i]==m])) + 1
          code('INDTYP *ARG_vec['+str(ind)+'] = {')
          depth += 2
          for n in range(0,nargs):
            if inds[n] == m:
              g_m = n
              code('ARG_l,')
          depth -= 2
          code('};')
        else:
          ind = int(max([idxs[i] for i in range(len(inds)) if inds[i]==m])) + 1
          code('INDTYP *ARG_vec['+str(ind)+'];')
# 
# lengthy code for general case with indirection
#
    if ninds>0:
      code('')
      for g_m in range (0,ninds):
        code('int  *ind_ARG_map, ind_ARG_size;')
      for g_m in range (0,ninds):
        code('INDTYP *ind_ARG_s;')
       
      code('int    nelem, offset_b;')
      code('')
      code('char shared[128000];')
      code('')
      code('if (0==0) {')
      depth += 2
      code('')
      code('// get sizes and shift pointers and direct-mapped data')
      code('')
      code('int blockId = blkmap[blockIdx + block_offset];')
      code('nelem    = nelems[blockId];')
      code('offset_b = offset[blockId];')
      code('')
       
      for g_m in range (0,ninds):
        code('ind_ARG_size = ind_arg_sizes['+str(g_m)+'+blockId*'+ str(ninds)+'];')
      
      for m in range (1,ninds+1):
        g_m = m-1
        c = [i for i in range(len(inds)) if inds[i]==m]
        code('ind_ARG_map = &ind_map['+str(cumulative_indirect_index[c[0]])+\
        '*set_size] + ind_arg_offs['+str(m-1)+'+blockId*'+str(ninds)+'];')
       
      code('')
      code('// set shared memory pointers')
      code('int nbytes = 0;')
       
      for g_m in range(0,ninds):
        code('ind_ARG_s = (INDTYP *) &shared[nbytes];')
        if g_m < ninds-1:
          code('nbytes += ROUND_UP(ind_ARG_size*sizeof(INDTYP)*INDDIM);')
      depth -= 2
      code('}')
      code('')
      code('// copy indirect datasets into shared memory or zero increment')
      code('')
       
      for g_m in range(0,ninds):
        if indaccs[g_m]==OP_READ or indaccs[g_m]==OP_RW or indaccs[g_m]==OP_INC:
          code('for (int n=0; n<INDARG_size; n++)')
          depth += 2
          code('for (int d=0; d<INDDIM; d++)')
          depth += 2
          if indaccs[g_m]==OP_READ or indaccs[g_m]==OP_RW:
            code('INDARG_s[d+n*INDDIM] = INDARG[d+INDARG_map[n]*INDDIM];')
            code('')
          elif indaccs[g_m]==OP_INC:
            code('INDARG_s[d+n*INDDIM] = ZERO_INDTYP;')
          depth -= 4
       
      code('')
      code('// process set elements')
      code('')
      
      if ind_inc:
        code('for (int n=0; n<nelem; n++) {')
        depth += 2
        code('// initialise local variables            ')
        
        for g_m in range(0,nargs):
          if maps[g_m]==OP_MAP and accs[g_m]==OP_INC:
            code('for (int d=0; d<DIM; d++)')
            depth += 2
            code('ARG_l[d] = ZERO_TYP;')
            depth -= 2
      else:
        code('for (int n=0; n<nelem; n++) {')
        depth += 2
      
#
# simple alternative when no indirection
#
    else:
      code('// process set elements')
      code('for (int n=start; n<finish; n++) {')    
      depth += 2
                        
#
# kernel call#   

    # xxx: array of pointers for non-locals 
    for m in range(1,ninds+1):
      s = [i for i in range(len(inds)) if inds[i]==m]
      if sum(s)>1:
        if indaccs[m-1] <> OP_INC:
          code('')
          ctr = 0
          for n in range(0,nargs):
            if inds[n] == m and vectorised[m]:
              file_text +='    arg'+str(m-1)+'_vec['+str(ctr)+'] = ind_arg'+\
              str(inds[n]-1)+'_s+arg_map['+str(cumulative_indirect_index[n])+\
              '*set_size+n+offset_b]*'+str(dims[n])+';\n'
              ctr = ctr+1
    
    code('')
    code('// user-supplied kernel call')
    
    line = name+'('
    prefix = ' '*len(name)
    a = 0 #only apply indentation if its not the 0th argument
    indent =''
    for m in range (0, nargs):
      if a > 0:
        indent = '       '+' '*len(name)
        
      if maps[m] == OP_GBL:
        line += rep(indent+'ARG,\n',m)
        a = a+1
      elif maps[m]==OP_MAP and  accs[m]==OP_INC and vectorised[m]==0:
        line += rep(indent+'ARG_l,\n',m);
        a = a+1
      elif maps[m]==OP_MAP and vectorised[m]==0:
        line += rep(indent+'ind_arg'+str(inds[m]-1)+'_s+arg_map['+\
        str(cumulative_indirect_index[m])+'*set_size+n+offset_b]*DIM,\n',m)
        a = a+1
      elif maps[m]==OP_MAP and m == 0:
        line += rep(indent+'ARG_vec,'+'\n',inds[m]-1)
        a = a+1
      elif maps[m]==OP_MAP and m>0 and vectorised[m] <> vectorised[m-1]: #xxx:vector
        line += rep(indent+'ARG_vec,'+'\n',inds[m]-1)
        a = a+1
      elif maps[m]==OP_MAP and m>0 and vectorised[m] == vectorised[m-1]:
        line = line
        a = a+1
      elif maps[m]==OP_ID:
        if ninds>0:
          line += rep(indent+'ARG+(n+offset_b)*DIM,'+'\n',m)
          a = a+1
        else:
          line += rep(indent+'ARG+n*DIM,'+'\n',m)
          a = a+1
      else:
        print 'internal error 1 '
    
    code(line[0:-2]+');') #remove final ',' and \n  
    
#
# updating for indirect kernels ...
#    
    if ninds>0:
      if ind_inc:
        code('')
        code('// store local variables            ')
        
        for g_m in range(0,nargs):
          if maps[g_m] == OP_MAP and accs[g_m] == OP_INC:
            code('int ARG_map = arg_map['+ str(cumulative_indirect_index[g_m])+'*set_size+n+offset_b];')
        code('')
        
        for g_m in range(0,nargs):
          if maps[g_m] == OP_MAP and accs[g_m] == OP_INC:
            code('for (int d=0; d<'+str(dims[g_m])+'; d++)')
            depth += 2
            code('ind_arg'+str(inds[g_m]-1)+'_s[d+ARG_map*DIM] += ARG_l[d];')
            depth -= 2
      depth -= 2    
      code('}')
      
      s = [i for i in range(1,ninds+1) if indaccs[i-1]<> OP_READ]
      
      if len(s)>0 and max(s)>0:
        code('')
        code('// apply pointered write/increment')
        
      for g_m in range(0,ninds):
        if indaccs[g_m]==OP_WRITE or indaccs[g_m]==OP_RW or indaccs[g_m]==OP_INC:
          code('for (int n=0; n<INDARG_size; n++)'); depth += 2;
          code('for (int d=0; d<INDDIM; d++)'); depth += 2;
          if indaccs[g_m]==OP_WRITE or indaccs[g_m]==OP_RW:
            code('INDARG[d+INDARG_map[n]*INDDIM] = INDARG_s[d+n*INDDIM];')
          elif indaccs[g_m]==OP_INC:
            code('INDARG[d+INDARG_map[n]*INDDIM] += INDARG_s[d+n*INDDIM];')
          depth -= 4;
#
# ... and direct kernels
#
    else:
      depth -= 2
      code('}')
       
#
# global reduction
#
    depth -= 2
    code('}')    
    code('')
    
##########################################################################
# then C++ stub function
##########################################################################

    code('')
    comm('host stub function          ')
    code('void op_par_loop_'+name+'(char const *name, op_set set,')
    depth += 2
    
    for m in unique_args:
      g_m = m - 1
      if m == unique_args[len(unique_args)-1]:
        code('op_arg ARG){'); 
        code('')
      else:
        code('op_arg ARG,')
    
    for g_m in range (0,nargs):
      if maps[g_m]==OP_GBL and accs[g_m] <> OP_READ:
        code('TYP*ARGh = (TYP *)ARG.data;')
    
    code('int nargs = '+str(nargs)+';')
    code('op_arg args['+str(nargs)+'];')
    
    
    for g_m in range (0,nargs):
      u = [i for i in range(0,len(unique_args)) if unique_args[i]-1 == g_m]
      if len(u) > 0 and vectorised[g_m] > 0:
        code('ARG.idx = 0;')
        code('args['+str(g_m)+'] = ARG;')
        
        v = [int(vectorised[i] == vectorised[g_m]) for i in range(0,len(vectorised))]
        first = [i for i in range(0,len(v)) if v[i] == 1]
        first = first[0]
        code('')
        code('for (int v = 1; v < '+str(sum(v))+'; v++) {'); depth += 2;
        code('args['+str(g_m)+' + v] = op_arg_dat(arg'+str(first)+'.dat, v, arg'+\
        str(first)+'.map, DIM, "TYP", '+accsstring[accs[g_m]-1]+');')
        depth -= 2
        code('}') 
        
      elif vectorised[g_m]>0:
        file_text = file_text
      else:
        code('args['+str(g_m)+'] = ARG;')
        
#
#   indirect bits
#
    if ninds>0:
      code('')
      code('int    ninds   = '+str(ninds)+';')
      file_text += '  int    inds['+str(nargs)+'] = {'
      depth += 2
      for m in range(0,nargs):
        file_text += str(inds[m]-1)+','
      file_text = file_text[:-1]
      
      depth -= 2
      code('};')
      code('')
      
      code('if (OP_diags>2) {');depth += 2;
      code('printf(" kernel routine with indirection: '+name+'\\n");');depth -= 2;
      code('}')
      code('')
      code('// get plan')
      code('#ifdef OP_PART_SIZE_'+ str(nk))
      code('  int part_size = OP_PART_SIZE_'+str(nk)+';')
      code('#else')
      code('  int part_size = OP_part_size;')
      code('#endif')
      code('')
      code('int set_size = op_mpi_halo_exchanges(set, nargs, args);')  
    
#
# direct bit
#   
    else:
      code('')
      code('if (OP_diags>2) {');depth += 2;
      code('printf(" kernel routine w/o indirection:  '+ name + '");')
      depth -= 2;
      code('}')
      code('')
      code('op_mpi_halo_exchanges(set, nargs, args);')

#
# start timing
#
    code('')
    code('// initialise timers')
    code('double cpu_t1, cpu_t2, wall_t1, wall_t2;')
    code('op_timers_core(&cpu_t1, &wall_t1);')
    code('')

#
# set number of threads in x86 execution and create arrays for reduction
#

    if reduct or ninds==0:
      code('// set number of threads')
      code('#ifdef _OPENMP')
      code('  int nthreads = omp_get_max_threads();')
      code('#else')
      code('  int nthreads = 1;')
      code('#endif')
     
    if reduct:
      code('')
      code('// allocate and initialise arrays for global reduction')
      for g_m in range(0,nargs):
        if maps[g_m]==OP_GBL and accs[g_m]<>OP_READ:
          code('TYP ARG_l[DIM+64*64];')
          code('for (int thr=0; thr<nthreads; thr++)'); depth += 2;
          if accs[g_m]==OP_INC:
            code('for (int d=0; d<DIM; d++) ARG_l[d+thr*64]=ZERO_TYP;')
          else:
            code('for (int d=0; d<DIM; d++) ARG_l[d+thr*64]=ARGh[d];')
          depth -= 2  
    
    code('')
    code('if (set->size >0) {'); depth += 2;
    code('')
    
#
# kernel call for indirect version
#
    if ninds>0:
      code('op_plan *Plan = op_plan_get(name,set,part_size,nargs,args,ninds,inds);')
      code('')
      code('// execute plan')
      code('int block_offset = 0;')
      code('for (int col=0; col < Plan->ncolors; col++) {'); depth += 2;
      code('if (col==Plan->ncolors_core) op_mpi_wait_all(nargs, args);')
      code('int nblocks = Plan->ncolblk[col];')
      code('')
      code('#pragma omp parallel for')
      code('for (int blockIdx=0; blockIdx<nblocks; blockIdx++)');depth += 2;
      code('op_x86_'+name+'( blockIdx,')

      for m in range(1,ninds+1):
        g_m = invinds[m-1]
        code('(TYP *)ARG.data,')
    
      code('Plan->ind_map,')
      code('Plan->loc_map,')

      for m in range(0,nargs):
        g_m = m
        if inds[m]==0 and maps[m] == OP_GBL and accs[m] <> OP_READ:
          code('&ARG_l[64*omp_get_thread_num()],')
        if inds[m]==0:
          code('(TYP *)ARG.data,')
      

      code('Plan->ind_sizes,')
      code('Plan->ind_offs,')
      code('block_offset,')
      code('Plan->blkmap,')
      code('Plan->offset,')
      code('Plan->nelems,')
      code('Plan->nthrcol,')
      code('Plan->thrcol,')
      code('set_size);')
      depth -= 2
      code('')
      
      if reduct:
        code('// combine reduction data')
        code('if (col == Plan->ncolors_owned-1) {'); depth += 2;
        for m in range(0,nargs):
          if maps[m] == OP_GBL and accs[m] <> OP_READ:
            code('for (int thr=0; thr<nthreads; thr++)')
            if accs[m]==OP_INC:
              code('for(int d=0; d<DIM; d++) ARGh[d] += ARG_l[d+thr*64];')
            elif accs[m]==OP_MIN:
              code('for(int d=0; d<DIM; d++) ARGh[d]  = MIN(ARGh[d],ARG_l[d+thr*64]);')
            elif  accs(m)==OP_MAX:
              code('for(int d=0; d<DIM; d++) ARGh[d]  = MAX(ARGh[d],ARG_l[d+thr*64]);')
            else:
              error('internal error: invalid reduction option')
      
      code('block_offset += nblocks;'); 
      depth -= 2
      code('}')     
      
#
# kernel call for direct version
#
    else:
      code('// execute plan')
      code('#pragma omp parallel for')
      code('for (int thr=0; thr<nthreads; thr++) {');depth += 2
      code('int start  = (set->size* thr)/nthreads;')
      code('int finish = (set->size*(thr+1))/nthreads;')
      file_text += '      op_x86_'+name+'('

      for g_m in range(0,nargs):
        indent = ''
        if g_m <> 0:
          indent = '        '+' '*len(name)
        else:
          indent = ''
        if maps[g_m]==OP_GBL and accs[g_m] <> OP_READ:
          code(indent+'ARG_l + thr*64,')
        else:
          code(indent+'(TYP *) ARG.data,')
        
      code(' '*len(name)+'        start, finish );')
      depth -= 2
      code('}')
      
    
  
    if ninds>0:
      code('op_timing_realloc('+str(nk)+');')
      code('OP_kernels['+str(nk)+'].transfer  += Plan->transfer; ')
      code('OP_kernels['+str(nk)+'].transfer2 += Plan->transfer2;')
    
    depth -= 2
    code('}')
    code('')    

#
# combine reduction data from multiple OpenMP threads
#
    code('// combine reduction data')
    for g_m in range(0,nargs):
      if maps[g_m]==OP_GBL and accs[g_m]<>OP_READ:
        code('for (int thr=0; thr<nthreads; thr++)')
        if accs[g_m]==OP_INC:
          code('for(int d=0; d<DIM; d++) ARGh[d] += ARG_l[d+thr*64];')
        elif accs[g_m]==OP_MIN:
          code('for(int d=0; d<DIM; d++) ARGh[d]  = MIN(arg'+\
          str(g_m)+'h[d],ARG_l[d+thr*64]);')
        elif accs[g_m]==OP_MAX:
          code('for(int d=0; d<DIM; d++) ARGh[d]  = MAX(ARGh[d],ARG_l[d+thr*64]);')
        else:
          print 'internal error: invalid reduction option'
        code('')
        code('op_mpi_reduce(&ARG,ARGh);')
      
    
    code('op_mpi_set_dirtybit(nargs, args);')
    code('')

#
# update kernel record
#

    code('// update kernel record')
    code('op_timers_core(&cpu_t2, &wall_t2);')
    code('op_timing_realloc('+str(nk)+');')
    code('OP_kernels[' +str(nk)+ '].name      = name;')
    code('OP_kernels[' +str(nk)+ '].count    += 1;')
    code('OP_kernels[' +str(nk)+ '].time     += wall_t2 - wall_t1;')

    if ninds == 0:
      line = 'OP_kernels['+str(nk)+'].transfer += (float)set->size *'

      for g_m in range (0,nargs):
        if maps[g_m]<>OP_GBL:
          if accs[g_m]==OP_READ or accs[g_m]==OP_WRITE:
            code(line+' ARG.size;')
          else:
            code(line+' ARG.size * 2.0f;')
            
    depth -= 2        
    code('}')


##########################################################################
#  output individual kernel file
##########################################################################
    fid = open(name+'_kernel.cpp','w')
    date = datetime.datetime.now()
    fid.write('//\n// auto-generated by op2.py on '+date.strftime("%Y-%m-%d %H:%M")+'\n//\n\n')
    fid.write(file_text)
    fid.close()

# end of main kernel call loop


##########################################################################
#  output one master kernel file
##########################################################################
  
  file_text =''
  comm(' header                 ')
  code('#include "op_lib_cpp.h"       ')
  code('')
  comm('// global constants       ')
              
  for nc in range (0,len(consts)):
    if consts[nc]['dim']==1:
      code('extern '+consts[nc]['type'][1:-1]+' '+consts[nc]['name']+';')
    else:
      if consts[nc]['dim'] > 0:
        num = str(consts[nc]['dim'])
      else:
        num = 'MAX_CONST_SIZE'
    
      code('extern '+consts[nc]['type'][1:-1]+' '+consts[nc]['name']+'['+num+'];')
      
  if any_soa:
    code('')
    code('extern int op2_stride;')
    code('#define OP2_STRIDE(arr, idx) arr[idx]')
    code('')

  comm(' user kernel files')

  for nk in range(0,len(kernels)):
    code('#include "'+kernels[nk]['name']+'_kernel.cpp"')
  master = master.split('.')[0] 
  fid = open(master.split('.')[0]+'_kernels.cpp','w')
  fid.write('//\n// auto-generated by op2.py on '+date.strftime("%Y-%m-%d %H:%M")+'\n//\n\n')
  fid.write(file_text)
  fid.close()
  
              
    
