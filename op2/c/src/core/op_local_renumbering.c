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

//set up a list of arrays to hold the new reordered indicies of a set
int** reordered_index;
//array to hold the flag that says if a set is reordered
int *reordered;

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

int get_lone_elements_by_reverse_map(op_map primary_map, op_set set)
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

  for(int i = 0; i<primary_map->to->size; i++)free(to[i]);
  free(adj_i);free(to);free(to_i);free(to_cap);

  return lone_elems;
}

int remove_lone_elems(op_set set, int** adj, int* adj_i, int* new_index)
{
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

  for(int i = 0; i<set->core_size; i++)
  {
    adj[i] = adj_new[i];
    adj_i[i] = adj_i_new[i];
  }
  free(adj_i_new);free(adj_new);
  return rev_count;
}

int* metis_call(op_map primary_map, op_set set, int elem_count,
  int** adj, int* adj_i, int num_part )
{
  //create csr table with *only* the connected set elements
  idxtype *xadj = (idxtype *)xmalloc(sizeof(idxtype)*(elem_count+2));
  int cap = (elem_count+1)*primary_map->dim;
  idxtype *adjncy = (idxtype *)xmalloc(sizeof(idxtype)*cap);
  int count = 0;int prev_count = 0;

  for(int i = 0; i<elem_count+1; i++)
  {
    int l_index = i;
    if(adj_i[i] > 1){
      for(int j = 0; j<adj_i[i]; j++){
        if(adj[i][j] != l_index){
          if(count >= cap){
            cap = cap*2;
            adjncy = (idxtype *)xrealloc(adjncy,sizeof(idxtype)*cap);
          }
          adjncy[count++] = (idxtype)adj[i][j];
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
  xadj[elem_count+1] = count;


  idxtype ncon = 1;
  int edgecut = 0;
  int part = num_part;

  idxtype wgtflag = 0;
  int options[5] = {0};
  idxtype numflag = 0;
  int nn = elem_count+1;

  //free memory before calling METIS
  for(int i = 0; i<set->core_size; i++)free(adj[i]);
  free(adj_i);free(adj);

  //call partitioner with the number of mini-partitions required and get the partition info array
  if(nn > part)
  {
    adjncy = (idxtype *)xrealloc(adjncy,sizeof(idxtype)*count);

    idxtype *p = (idxtype *)xmalloc(sizeof(idxtype)*elem_count+1);
    for(int i = 0; i < elem_count+1; i++){ p[i] = -99; }

    METIS_PartGraphKway(&nn, xadj, adjncy, NULL, NULL, &wgtflag, &numflag,
      &part, options, &edgecut, p);

    //free memory used for METIS
    free(xadj); free(adjncy);

    return (int*)p;
  }

  //free memory used for METIS
  free(xadj); free(adjncy);

  return NULL;
}

void create_map(op_map primary_map, op_set set, int flag, int** map, int* map_i, int* cap)
{
  if(flag == 0) //return same map
  {
    for(int i = 0; i<primary_map->from->size;i++)
    {
      int index;
      for(int j=0; j < primary_map->dim; j++) { //for each element pointed at by this entry
        index = primary_map->map[i*primary_map->dim+j];
        if(map_i[i] >= cap[i])
        {
          cap[i] = 2* cap[i];
          map[i] = (int *)xrealloc(map[i], cap[i]*sizeof(int));
        }
        map[i][map_i[i]++] = index;
      }
    }
  }
  else //return reverse map
  {
    for(int i = 0; i<primary_map->from->size;i++)
    {
      int index;
      for(int j=0; j < primary_map->dim; j++) { //for each element pointed at by this entry
        index = primary_map->map[i*primary_map->dim+j];
        if(index < primary_map->to->size && index >= 0)
        {
          if(map_i[index] >= cap[index])
          {
            cap[index] = 2* cap[index];
            map[index] = (int *)xrealloc(map[index], cap[index]*sizeof(int));
          }
          map[index][map_i[index]++] = i;
        }
      }
    }

    for(int i = 0; i<primary_map->to->size; i++)
    {
      //sort
      quickSort(map[i], 0, map_i[i]-1);
      //remove duplicates
      if(map_i[i] > 0)map_i[i] = removeDups(map[i], map_i[i]);
      else map_i[i] = 0;
      //realloc
      map[i] = (int *)xrealloc(map[i],map_i[i]*sizeof(int));
    }
  }
}


void reorder(op_map primary_map, op_set set, int num_part, int flag)
{
  int** map = NULL;
  int* map_i = NULL;
  int* cap = NULL;

  int count1 = 0;
  int count2 = 0;
  if(flag == 0)
  {
    count1 = primary_map->from->core_size;
    count2 = primary_map->to->core_size;

    map = (int **)xmalloc(primary_map->from->size*sizeof(int *));
    map_i = (int *)xmalloc(primary_map->from->size*sizeof(int ));
    cap = (int *)xmalloc(primary_map->from->size*sizeof(int ));

    for(int i = 0; i<primary_map->from->size; i++)map_i[i] = 0;
    for(int i = 0; i<primary_map->from->size; i++)cap[i] = primary_map->dim;
    for(int i = 0; i<primary_map->from->size; i++)map[i] = (int *)xmalloc(cap[i]*sizeof(int));
  }
  else
  {
    count1 = primary_map->to->core_size;
    count2 = primary_map->from->core_size;

    map = (int **)xmalloc(primary_map->to->size*sizeof(int *));
    map_i = (int *)xmalloc(primary_map->to->size*sizeof(int ));
    cap = (int *)xmalloc(primary_map->to->size*sizeof(int ));

    for(int i = 0; i<primary_map->to->size; i++)map_i[i] = 0;
    for(int i = 0; i<primary_map->to->size; i++)cap[i] = primary_map->dim;
    for(int i = 0; i<primary_map->to->size; i++)map[i] = (int *)xmalloc(cap[i]*sizeof(int));
  }
  create_map(primary_map, set, flag, map, map_i, cap);

  int** adj = (int **)xmalloc(count2*sizeof(int *));
  int* adj_i = (int *)xmalloc(count2*sizeof(int ));
  int* adj_cap = (int *)xmalloc(count2*sizeof(int ));

  for(int i = 0; i<count2; i++)adj_i[i] = 0;
  for(int i = 0; i<count2; i++)adj_cap[i] = primary_map->dim;
  for(int i = 0; i<count2; i++)adj[i] = (int *)xmalloc(adj_cap[i]*sizeof(int));

  //go through each from-element of local primary_map and construct adjacency list
  for(int i = 0; i<count1; i++)
  {
    int local_index;
    for(int j=0; j < map_i[i]; j++) { //for each element pointed at by this entry
      local_index = map[i][j];

      if(local_index < count2 && local_index >=0)
      {
        for(int k = 0; k<map_i[i]; k++)
        {
          if(adj_i[local_index] >= adj_cap[local_index])
          {
            adj_cap[local_index] = adj_cap[local_index]*2;
            adj[local_index] = (int *)xrealloc(adj[local_index],
                adj_cap[local_index]*sizeof(int));
          }
          if(map[i][k] < count2)
            adj[local_index][adj_i[local_index]++] = map[i][k];
        }
      }
    }
  }
  free(adj_cap);
  for(int i = 0; i<count1; i++)free(map[i]);
  free(cap);free(map_i);

  //create array with new element indicies when removing lone elements
  int* new_index = (int *)xmalloc(set->core_size*sizeof(int ));
  int elem_count = remove_lone_elems(set, adj, adj_i, new_index);

  //look for lone nodes/vertices
  int lone_vertices_index = 0;
  for(int i = 0; i<elem_count+1; i++)
    if(adj_i[i] < 2) lone_vertices_index++;
  if(lone_vertices_index > 0)
    printf("Lone vertices FOUND from using map %s for set %s = %d\n",
        primary_map->name, set->name, lone_vertices_index);

  int* p = metis_call(primary_map, set, elem_count, adj, adj_i, num_part);

  if(p != NULL)
  {
    //create index array to hold new element indecies
    int* new_index_2 = (int *)xmalloc(set->core_size*sizeof(int ));
    int* temp = (int *)xmalloc(set->core_size*sizeof(int ));
    for(int i = 0; i<set->core_size; i++)new_index_2[i] = i;

    //sort new_index_2 according to parition information in p[]
    quickSort_2(p, new_index_2, 0, elem_count);

    //reveres new_index_2
    for(int i = 0; i<set->core_size; i++)
      temp[new_index_2[i]] = i;
    free(new_index_2);
    new_index_2 = temp;

    //update new_index using new_index_2
    for(int i = 0; i<set->core_size; i++)
      new_index[i] = new_index_2[new_index[i]];
    free(new_index_2);
    free(p);

    //store reorder information in reordered_index array
    reordered_index[set->index] = new_index;
    reordered[set->index] = 1;
  }
}

void migrate_dats_and_maps()
{
  for(int s=0; s<OP_set_index; s++) { //for each set
    op_set set=OP_set_list[s];

    if(reordered[set->index] == 1)
    {
      //using the possition of the reordered_index array move the data on all
      //sets to the correct possition
      for(int d=0; d<OP_dat_index; d++) { //for each dat
        op_dat dat=OP_dat_list[d];

        if(compare_sets(set,dat->set)==1)//if this data array is defined on this set
        {
          char* new_dat = (char* )xmalloc(set->core_size*dat->size);
          for(int i = 0; i<set->core_size; i++)
          {
            memcpy(&new_dat[reordered_index[set->index][i]*dat->size],
                &dat->data[i*dat->size],dat->size);
          }
          memcpy(&dat->data[0],&new_dat[0], set->core_size*dat->size);
          free(new_dat);
        }
      }

      //using the possition of the reordered_index array move the maps on all
      //set to the correct possition
      for(int m=0; m<OP_map_index; m++) { //for each map
        op_map map=OP_map_list[m];

        if(compare_sets(set,map->from)==1)//if this data array is defined on this set
        {
          int* new_map = (int* )xmalloc(set->core_size*map->dim*sizeof(int));
          for(int i = 0; i<set->core_size; i++)
          {
            memcpy(&new_map[reordered_index[set->index][i]*map->dim],
                &map->map[i*map->dim],map->dim*sizeof(int));
          }
          memcpy(&map->map[0],&new_map[0], set->core_size*map->dim*sizeof(int));
          free(new_map);
        }
      }

      //go thorugh each mapping table - both core and non-core and exec-halo
      //renumber the to set with the new possitions of the elements of the to set
      for(int m=0; m<OP_map_index; m++) { //for each set
        op_map map=OP_map_list[m];

        if(compare_sets(set,map->to)==1)//if this data array is defined on this set
        {
          int map_size = map->from->size+OP_import_exec_list[map->from->index]->size;
          for(int i = 0; i<map_size; i++)
          {
            for(int k = 0; k<map->dim; k++)
            {
              int local_index = map->map[i*map->dim+k];
              if(local_index < map->to->core_size)
              {
                OP_map_list[map->index]->map[i*map->dim+k] =
                  reordered_index[map->to->index][map->map[i*map->dim+k]];
              }
            }
          }
        }
      }

      //go through the export halo lists of this set and update the list with the new indexes
      for(int i = 0; i< OP_export_exec_list[set->index]->size; i++)
      {
        int local_index = OP_export_exec_list[set->index]->list[i];
        if(local_index < set->core_size)
          OP_export_exec_list[set->index]->list[i] =
            reordered_index[set->index][local_index];
      }
      for(int i = 0; i< OP_export_nonexec_list[set->index]->size; i++)
      {
        int local_index = OP_export_nonexec_list[set->index]->list[i];
        if(local_index < set->core_size)
          OP_export_nonexec_list[set->index]->list[i] =
            reordered_index[set->index][local_index];
      }
      printf("Reordering and data migration done for set %s\n",set->name);
    }
  }
}

/*******************************************************************************
 * Routine to do local renumbering using a specified partitioning library
 *******************************************************************************/

void op_local_renumbering(int num_part /* const char *library*/)
{
  /*need to find the best mapping table to create the adjacancy list for each
   * set. The adjacancy list will be used as a parameter in the metis routine
   * call.
   *
   * The best map will be the one with the least amount of lone verticies
   * Only the core elements will be renumbered
   */

  //set up a list of arrays to hold the new reordered indicies of a set
  reordered_index = (int **)xmalloc(OP_set_index*sizeof(int *));

  //array to hold the flag that says if a set is reordered
  reordered = (int *)xmalloc(OP_set_index*sizeof(int));
  for(int i = 0; i<OP_set_index; i++) reordered[i] = 0;

  //for each set find the suitable mapping
  for(int s=0; s<OP_set_index; s++){
    op_set set = OP_set_list[s];

    int best_map = -1; int lone_elems = set->size; int temp = 0;
    //select appropriate mapping table to create adjacancy list for this set

    // - maps that allow to partition the "from set" is prefered
    for(int m = 0; m<OP_map_index; m++){
      op_map map = OP_map_list[m];
      if (compare_sets(map->from,set)==1)
      {
        temp = get_lone_elements_by_reverse_map(map,set);
        if(lone_elems > temp)
        {
          lone_elems = temp;
          best_map = map->index;
        }
      }
    }

    if(best_map < 0)
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
      }

    if(best_map < 0)
    {
      printf("Set %s not reordered : could not find suitable mapping or reverse mapping\n",
        set->name);
    }
    else
    {
      printf("Best map for Set %s : %s \n", set->name,OP_map_list[best_map]->name);
      if (compare_sets(OP_map_list[best_map]->from,set)==1)
      {
        //create_reverse_reordering(OP_map_list[best_map], set, num_part);
        reorder(OP_map_list[best_map], set, num_part, 1);
      }
      else if(compare_sets(OP_map_list[best_map]->to,set)==1)
      {
        //create_reordering(OP_map_list[best_map], set, num_part);
        reorder(OP_map_list[best_map], set, num_part, 0);
      }
    }
  }

  migrate_dats_and_maps();

  //cleanup
  for(int i = 0; i<OP_set_index;i++)if(reordered[i]>0)free(reordered_index[i]);
  free(reordered_index);free(reordered);
}
