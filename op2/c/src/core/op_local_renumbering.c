#include <time.h>
#include <stdlib.h>
#include <stdio.h>

#include <op_lib_c.h>
#include <op_lib_core.h>
#include <op_util.h>

//metis header
#ifdef HAVE_METIS
#include <metis.h>
#endif

//scotch header
#ifdef HAVE_SCOTCH
#include <metis.h>
#endif

#include <op_mpi_core.h>

int get_lone_elems(op_map primary_map, op_set set)
{
  int* adj_i = (int *)xmalloc(primary_map->to->core_size*sizeof(int ));
  for(int i = 0; i<primary_map->to->core_size; i++)adj_i[i] = 0;

  int lone_elems = 0;
  for(int i = 0; i<primary_map->from->core_size; i++)
  {
    int local_index;
    for(int j=0; j < primary_map->dim; j++) { //for each element pointed at by this entry
      local_index = primary_map->map[i*primary_map->dim+j];
      if(local_index < primary_map->to->core_size)
      {
        for(int k = 0; k<primary_map->dim; k++)
        {
          if(primary_map->map[i*primary_map->dim+k] < primary_map->to->core_size &&
            primary_map->map[i*primary_map->dim+k] != local_index)
            adj_i[local_index]++;
        }
      }
    }
  }

  for(int i = 0; i<primary_map->to->core_size; i++)
    if(adj_i[i] == 0) lone_elems++;
  printf("for set %s lone elements using map %s: %d\n", set->name,
    primary_map->name, lone_elems);

  free(adj_i);

  return lone_elems;
}

int get_lone_elements_from_reverse_map(op_map primary_map, op_set set)
{
  //create temporary reverse map
  int** to = (int **)xmalloc(primary_map->to->size*sizeof(int *));
  int* to_i = (int *)xmalloc(primary_map->to->size*sizeof(int ));
  int* to_cap = (int *)xmalloc(primary_map->to->size*sizeof(int ));

  for(int i = 0; i<primary_map->to->size; i++)to_i[i] = 0;
  for(int i = 0; i<primary_map->to->size; i++)to_cap[i] = primary_map->dim;
  for(int i = 0; i<primary_map->to->size; i++)to[i] = (int *)xmalloc(to_cap[i]*sizeof(int));

  for(int i = 0; i<primary_map->from->size;i++)
  {
    int index;
    for(int j=0; j < primary_map->dim; j++) { //for each element pointed at by this entry
      index = primary_map->map[i*primary_map->dim+j];
      if(index < primary_map->to->size && index >= 0)
      {
        if(to_i[index] >= to_cap[index])
        {
          to_cap[index] = 2* to_cap[index];
          to[index] = (int *)xrealloc(to[index], to_cap[index]*sizeof(int));
        }
        to[index][to_i[index]++] = i;
      }
    }
  }

  for(int i = 0; i<primary_map->to->size; i++)
  {
    //sort
    quickSort(to[i], 0, to_i[i]-1);
    //remove duplicates
    if(to_i[i] > 0)to_i[i] = removeDups(to[i], to_i[i]);
    else to_i[i] = 0;
    //realloc
    to[i] = (int *)xrealloc(to[i],to_i[i]*sizeof(int));
  }

  int* adj_i = (int *)xmalloc(primary_map->from->core_size*sizeof(int ));
  int lone_elems = 0;
  //go through each core element of reverse_map and construct adjacency list
  for(int i = 0; i<primary_map->to->core_size; i++)
  {
    int local_index;
    for(int j=0; j < to_i[i]; j++) { //for each element pointed at by this entry
      local_index = to[i][j];
      if(local_index < primary_map->from->core_size /*&& local_index >=0*/ )
      {
        for(int k = 0; k<to_i[i]; k++)
        {
          if(primary_map->map[i*primary_map->dim+k] < primary_map->to->core_size &&
            primary_map->map[i*primary_map->dim+k] != local_index)
            adj_i[local_index]++;
        }
      }
    }
  }

  for(int i = 0; i<primary_map->from->core_size; i++)
    if(adj_i[i] == 0) lone_elems++;
  printf("for set %s lone elements using map reversed %s: %d\n", set->name,
    primary_map->name, lone_elems);

  free(adj_i);free(to);free(to_i);free(to_cap);

  return lone_elems;
}

void remove_lone_elems(op_map primary_map, op_set set, int** adj, int* adj_i)
{
  //create array with new element indicies when removing lone elements
  int* new_index = (int *)xmalloc(set->core_size*sizeof(int ));

  int count = 0; int rev_count = set->core_size -1;
  for(int i = 0; i < set->core_size; i++)
  {
    quickSort(adj[i], 0, adj_i[i]-1);
    if(adj_i[i]>1)adj_i[i] = removeDups(adj[i], adj_i[i]);
    adj[i] = (int *)xrealloc(adj[i],adj_i[i]*sizeof(int));
    if(adj_i[i] > 1)
      new_index[i] = count++;
    else //add the lone elements to the end
      new_index[i] = rev_count--;
  }

  //renumber the adjacancy list
  for(int i = 0; i<set->core_size; i++)
    for(int j = 0; j<adj_i[i]; j++)
      adj[i][j]= new_index[adj[i][j]];

  //move elements in adj[] to the new locations
  int** adj_new = (int **)xmalloc(set->core_size*sizeof(int *));
  int* adj_i_new = (int *)xmalloc(set->core_size*sizeof(int ));
  for(int i = 0; i<set->core_size; i++)
  {
    adj_new[new_index[i]] = adj[i];
    adj_i_new[new_index[i]] = adj_i[i];
  }

  //look for lone nodes/vertices
  int lone_vertices_index = 0;
  for(int i = 0; i<rev_count+1; i++)
    if(adj_i_new[i] < 2) lone_vertices_index++;
  if(lone_vertices_index > 0)
    printf("Lone vertices FOUND from using map %s for set %s = %d\n",
    primary_map->name, set->name, lone_vertices_index);

  adj = adj_new;
  adj_i = adj_i_new;
}

void create_reordering(op_map primary_map, op_set set, int num_part)
{
  int** adj = (int **)xmalloc(primary_map->to->core_size*sizeof(int *));
  int* adj_i = (int *)xmalloc(primary_map->to->core_size*sizeof(int ));
  int* adj_cap = (int *)xmalloc(primary_map->to->core_size*sizeof(int ));

  for(int i = 0; i<primary_map->to->core_size; i++)adj_i[i] = 0;
  for(int i = 0; i<primary_map->to->core_size; i++)adj_cap[i] = primary_map->dim;
  for(int i = 0; i<primary_map->to->core_size; i++)adj[i] = (int *)xmalloc(adj_cap[i]*sizeof(int));

  //go through each from-element of local primary_map and construct adjacency list
  for(int i = 0; i<primary_map->from->core_size; i++)
  {
    int local_index;
    for(int j=0; j < primary_map->dim; j++) { //for each element pointed at by this entry
      local_index = primary_map->map[i*primary_map->dim+j];

      if(local_index < primary_map->to->core_size)
      {
        for(int k = 0; k<primary_map->dim; k++)
        {
          if(adj_i[local_index] >= adj_cap[local_index])
          {
            adj_cap[local_index] = adj_cap[local_index]*2;
            adj[local_index] = (int *)xrealloc(adj[local_index],
              adj_cap[local_index]*sizeof(int));
          }
          if(primary_map->map[i*primary_map->dim+k] < primary_map->to->core_size)
            adj[local_index][adj_i[local_index]++] =
            primary_map->map[i*primary_map->dim+k];
        }
      }
    }
  }

  //create csr table with *only* the connected set elements
  idxtype *xadj = (idxtype *)xmalloc(sizeof(idxtype)*(rev_count+2));
  int cap = (rev_count+1)*primary_map->dim;
  idxtype *adjncy = (idxtype *)xmalloc(sizeof(idxtype)*cap);
  int count = 0;int prev_count = 0;


  for(int i = 0; i<rev_count+1; i++)
  {
    int l_index = i;
    if(adj_i_new[i] > 1){
      for(int j = 0; j<adj_i_new[i]; j++){
        if(adj_new[i][j] != l_index){
          if(count >= cap){
            cap = cap*2;
            adjncy = (idxtype *)xrealloc(adjncy,sizeof(idxtype)*cap);
          }
          adjncy[count++] = (idxtype)adj_new[i][j];
        }
      }
      if(i != 0){
        xadj[i] = prev_count;
        prev_count = count;
      }
      else{
        xadj[i] = 0;
        prev_count = count;
      }
    }
  }
  xadj[rev_count+1] = count;

  printf("before count for set %s = %d\n",set->name,count);

  remove_lone_elems(primary_map, set, adj, adj_i);


  for(int i = 0; i<rev_count+1; i++)
  {
    int l_index = i;
    if(adj_i_new[i] > 1){
      for(int j = 0; j<adj_i_new[i]; j++){
        if(adj_new[i][j] != l_index){
          if(count >= cap){
            cap = cap*2;
            adjncy = (idxtype *)xrealloc(adjncy,sizeof(idxtype)*cap);
          }
          adjncy[count++] = (idxtype)adj_new[i][j];
        }
      }
      if(i != 0){
        xadj[i] = prev_count;
        prev_count = count;
      }
      else{
        xadj[i] = 0;
        prev_count = count;
      }
    }
  }
  xadj[rev_count+1] = count;

  printf("after count for set %s = %d\n",set->name,count);
  //metis_call(primary_map, set, adj, adj_i, adj_cap, num_part);
}

/*******************************************************************************
 * Routine to do local renumbering using Metis k-way
 *******************************************************************************/

void op_local_renumbering_metiskway(int num_part)
{
  /*need to find the best mapping table to create the adjacancy list for each
  * set. The adjacancy list will be used as a parameter in the metis routine
  * call.
  *
  * The best map will be the one with the least amount of lone verticies
  * Only the core elements will be renumbered
  */

  //for each set find the suitable mapping
  for(int s=0; s<OP_set_index; s++){
    op_set set = OP_set_list[s];

    int best_map = -1; int lone_elems = set->size; int temp = 0;
    //select appropriate mapping table to create adjacancy list for this set
    for(int m = 0; m<OP_map_index; m++){
      op_map map = OP_map_list[m];
      if(compare_sets(map->to,set)==1)
      {
        temp = get_lone_elems(map,set);
        if(lone_elems > temp)
        {
          lone_elems = temp;
          best_map = map->index;
        }
      }
      else if (compare_sets(map->from,set)==1)
      {
        temp = get_lone_elements_from_reverse_map(map,set);
        if(lone_elems > temp)
        {
          lone_elems = temp;
          best_map = map->index;
        }
      }
    }

    if(best_map < 0)
    {
      printf("Set %s not reordered : could not find suitable mapping or reverse mapping\n",
        set->name);
    }
    else
    {
      printf("Best map for Set %s : %s \n", set->name,OP_map_list[best_map]->name);
      if(compare_sets(OP_map_list[best_map]->to,set)==1)
      {
        create_reordering(OP_map_list[best_map], set, num_part);
      }
      else if (compare_sets(OP_map_list[best_map]->from,set)==1)
      {

      }
    }
  }
}
