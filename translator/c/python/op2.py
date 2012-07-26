#!/usr/bin/python
#
# OP2 source code transformation tool
#
# This tool parses the user's original source code to produce
# target-specific code to execute the user's kernel functions.
#
# This prototype is written in Python 
#
# usage: op2('file1','file2',...)
#
# This takes as input
#
# file1.cpp, file2.cpp, ...
#
# and produces as output modified versions
#
# file1_op.cpp, file2_op.cpp, ...
#
# then calls a number of target-specific code generators
# to produce individual kernel files of the form
#
# xxx_kernel.cpp  -- for OpenMP x86 execution
# xxx_kernel.cu   -- for CUDA execution
#
# plus a master kernel file of the form
#
# file1_kernels.cpp  -- for OpenMP x86 execution
# file1_kernels.cu   -- for CUDA execution
#

import sys
import re

#
# declare constants
#

ninit = 0
nexit = 0
npart = 0
nhdf5 = 0

nconsts  = 0
nkernels = 0
consts = []
kernels = []

OP_ID   = 1;  OP_GBL   = 2;  OP_MAP = 3;

OP_READ = 1;  OP_WRITE = 2;  OP_RW  = 3;
OP_INC  = 4;  OP_MAX   = 5;  OP_MIN = 6;

OP_accs_labels = ['OP_READ','OP_WRITE','OP_RW','OP_INC','OP_MAX','OP_MIN' ]



##########################################################################
# Remove comments from text
##########################################################################

def comment_remover(text):
    def replacer(match):
        s = match.group(0)
        if s.startswith('/'):
            return ""
        else:
            return s
    pattern = re.compile(
        r'//.*?$|/\*.*?\*/|\'(?:\\.|[^\\\'])*\'|"(?:\\.|[^\\"])*"',
        re.DOTALL | re.MULTILINE
    )
    return re.sub(pattern, replacer, text)
    
##########################################################################
# parsing for op_init/op_exit/op_partition/op_hdf5 calls
##########################################################################

def op_parse_calls(text):
    
    inits = 0
    search = "op_init"
    found = text.find(search)
    while found > -1:
      found=text.find(search, found+1)
      inits = inits + 1
    
    exits = 0
    search = "op_exit"
    found = text.find(search)
    while found > -1:
      found=text.find(search, found+1)
      exits = exits + 1
      
    parts = 0
    search = "op_partition"
    found = text.find(search)
    while found > -1:
      found=text.find(search, found+1)
      parts = parts + 1
            
    hdf5s = 0
    search = "hdf5"
    found = text.find(search)
    while found > -1:
      found=text.find(search, found+1)
      hdf5 = hdf5s + 1
    
    return (inits, exits, parts, hdf5s)

##########################################################################
# parsing for op_decl_const calls
##########################################################################

def op_decl_const_parse(text):
    consts = []
    num_const = 0
    search = "op_decl_const"
    i = text.find(search)
    while i > -1:
      const_string = text[i+14: text.find(')',i+13)]
      #print 'Args found at index i : '+const_string      
      
      #check for syntax errors
      if len(const_string.split(',')) <> 3:
        print 'Error in op_decl_const : must have three arguments'
        return
      
      temp = {'loc': i,
        'dim':const_string.split(',')[0],
        'type': const_string.split(',')[1],
        'name':const_string.split(',')[2]}
      
      consts.append(temp)
      
      i=text.find(search, i+13)
      num_const = num_const + 1
            
    return (consts)
    
    
    

##########################################################################
# parsing for op_par_loop calls
##########################################################################

def op_par_loop_parse(text):
    loop_args = []
    
    search = "op_par_loop"
    i = text.find(search)
    while i > -1:
      #arg_string = text[i+12:text.find(';',i+11)]
      arg_string = text[text.find('(',i)+1:text.find(';',i+11)]
      #print arg_string

      #parse arguments in par loop
      temp_arg_string = []
      num_args = 0
      
      try:
        #parse each op_arg_dat
        search2 = "op_arg_dat"
        j = arg_string.find(search2)
        while j > -1:
          #dat_args_string =  arg_string[j+11: arg_string.find(')',j+12)]
          dat_args_string = arg_string[arg_string.find('(',j)+1:arg_string.find(')',j+12)]
          #print dat_args_string
          
          #check for syntax errors
          if len(dat_args_string.split(',')) <> 6:
            print 'Error in parsing op_arg_dat('+ dat_args_string +'): must have six arguments'
            return
          temp_arg_string.append('op_arg_dat, '+dat_args_string)
          num_args = num_args + 1        
          j= arg_string.find(search2, j+12) 
                
        #parse each op_arg_gbl
        search3 = "op_arg_gbl"
        k = arg_string.find(search3)
        while k > -1:
          #gbl_args_string =  arg_string[k+11: arg_string.find(')',k+12)]
          gbl_args_string = arg_string[arg_string.find('(',k)+1:arg_string.find(')',k+12)]
          #print gbl_args_string
          
          #check for syntax errors
          if len(gbl_args_string.split(',')) <> 4:
            print 'Error in parsing op_arg_gbl('+ dat_args_string +'): must have four arguments'
            return
          temp_arg_string.append('op_arg_gbl, '+gbl_args_string)
          num_args = num_args + 1
          k= arg_string.find(search3, k+12)        
        
        temp = {'loc':i,
          'name1':arg_string.split(',')[0],
          'name2':arg_string.split(',')[1],
          'set':arg_string.split(',')[2],
          'args':temp_arg_string,
          'nargs':num_args
        }
      
      except (RuntimeError, TypeError, NameError):
        print "error parsing op_par_loop"
      
      loop_args.append(temp)
      i=text.find(search, i+12)
    return (loop_args)



    
###################END OF FUNCTIONS DECLARATIONS #########################




    
##########################################################################
#  loop over all input source files
##########################################################################
for i in range(1,len(sys.argv)):
    print 'processing file '+ str(i) + ' of ' + str(len(sys.argv)) + ' '+ \
    str(sys.argv[i]) 
    
    src_file = str(sys.argv[i])
    f = open(src_file,'r')
    text = f.read()
    
    #remove commnets for parsing
    text = comment_remover(text)
    
    #print str(OP_accs_labels[1])
    
##########################################################################
# check for op_init/op_exit/op_partition/op_hdf5 calls
##########################################################################    
    
    inits, exits, parts, hdf5s = op_parse_calls(text)
    
    if inits+exits+parts+hdf5s > 0:
      print ' '
    if inits > 0:
      print'  contains op_init call\n'
    if exits > 0:
      print'  contains op_exit call\n'
    if parts > 0:
      print'  contains op_partition call\n'
    if hdf5s > 0:
      print'  contains op_hdf5 calls\n'
  
    ninit = ninit + inits;
    nexit = nexit + exits;
    npart = npart + parts;
    nhdf5 = nhdf5 + hdf5s;
    
    
##########################################################################
# parse and process constants
##########################################################################    
    
    const_args = op_decl_const_parse(text)
    #print'  number of constants declared '+str(len(const_args))
    #print'  contains constants '+str(const_args)
    
    #cleanup '&' symbols from name and convert dim to integer
    for i  in range(0,len(const_args)):
      if const_args[i]['name'][0] == '&':
        const_args[i]['name'] = const_args[i]['name'][1:]    
        const_args[i]['dim'] = int(const_args[i]['dim'])
    
    #check for repeats
    nconsts = 0
    for i  in range(0,len(const_args)):
      repeat = 0
      name = const_args[i]['name']
      for c in range(0,nconsts):
        if  const_args[i]['name'] == consts[c]['name']:
          repeat = 1
          if const_args[i]['type'] <> consts[c]['type']:
            print 'type mismatch in repeated op_decl_const'
          if const_args[i]['dim'] <> consts[c]['dim']:
            print 'size mismatch in repeated op_decl_const'
          
      if repeat > 0:
        print 'repeated global constant ' + const_args[i]['name']
      else:
        print 'global constant '+ const_args[i]['name'] + ' of size ' + str(const_args[i]['dim'])
        
    #store away in master list
      if repeat == 0:
        nconsts = nconsts + 1
        temp = {'dim': const_args[i]['dim'],
        'type': const_args[i]['type'],
        'name': const_args[i]['name']}
        consts.append(temp)
        

##########################################################################
# parse and process op_par_loop calls
##########################################################################
   
    loop_args = op_par_loop_parse(text)
    print loop_args
    #for i in range (0, len(loop_args)):
    #  print loop_args[i]
    
    f.close()
