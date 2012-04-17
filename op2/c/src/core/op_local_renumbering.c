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

#include <op_mpi_core.h>

//Need OP_reverse_map_list[OP_map_index] to hold the reverse mappings
int OP_reverse_map_index = 0;
op_map * OP_reverse_map_list;

void reverse_map(op_map map)
{
  int cap = map->dim;
  int** rev_map = (int **)xmalloc(map->from->size*sizeof(int));
}
  //for each map create reverse map
  /*for(int m = 0; m<OP_map_index; m++)
  {
    op_map map = OP_map_list[m];

    op_map r_map = ( op_map ) malloc ( sizeof ( op_map_core ) );
    map->index = OP_reverse_map_index;
    r_map->from = map->to;
    r_map->to = map->from;
    r_map->dim = dim;
    r_map->map = imap;
    r_map->name = "r_"+map->name;

    OP_reverse_map_list[OP_reverse_map_index++] = r_map;
  }*/

/*******************************************************************************
 * Routine to do local renumbering using Metis k-way
 *******************************************************************************/

void op_local_renumbering_metiskway(int num_part)
{
  //set up a list os arrays to hold the new reordered indicies of a set
  int** reordered_index = (int **)xmalloc(OP_set_index*sizeof(int *));

  //array to hold the flag that says if a set is reordered
  int *reordered = (int *)xmalloc(OP_set_index*sizeof(int));
  for(int i = 0; i<OP_set_index; i++) reordered[i] = 0;

  //for each set need to do the local renumbering for the core elements
  for(int s=0; s<OP_set_index; s++){
    op_set set = OP_set_list[s];
    //need int array to hold the candidate mapping indexes
    int candidate_maps[OP_map_index*2];
    int found = 0;
    for(int i = 0; i< OP_map_index*2; i++)candidate_maps[i]=-1;

    //select appropriate mapping table or reverse mapping table to create adjacancy list for this set
    for(int m = 0; m<OP_map_index; m++){
      op_map map = OP_map_list[m];

      if(compare_sets(map->to,set)==1) //partition to set
      {
        candidate_maps[found++] = map->index;//select map->index as a a candidate mapping table for this set
      }
      /*else if(compare_sets(map->from,set)==1) //partition from set
      {
        //create reverse mapping table and add it to OP_reverse_map_list[OP_reverse_map_index]
        //add map->index as a a candidate mapping table for this set
      }*/
    }

    if(candidate_maps[0] != -1)
    {
      op_map primary_map = OP_map_list[candidate_maps[0]];

      int** adj = (int **)xmalloc(primary_map->to->core_size*sizeof(int *));
      int* adj_i = (int *)xmalloc(primary_map->to->core_size*sizeof(int ));
      int* adj_cap = (int *)xmalloc(primary_map->to->core_size*sizeof(int ));

      for(int i = 0; i<primary_map->to->core_size; i++)adj_i[i] = 0;
      for(int i = 0; i<primary_map->to->core_size; i++)adj_cap[i] = primary_map->dim;
      for(int i = 0; i<primary_map->to->core_size; i++)adj[i] = (int *)xmalloc(adj_cap[i]*sizeof(int));

      for(int i = 0; i<primary_map->to->core_size; i++)
        for(int j = 0; j<adj_cap[i]; j++)
          adj[i][j] = 0;

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

      //look for lone nodes/vertices - BEFORE THEY ARE REMOVED
      /*int lone_vertices_index = 0;
      for(int i = 0; i<rev_count+1; i++)
        if(adj_i[i] < 2) lone_vertices_index++;
      printf("Lone vertices FOUND from using map %s for set %s = %d\n",
          primary_map->name, set->name, lone_vertices_index);*/

      //create array with new element indicies when removing lone elements
      int* new_index = (int *)xmalloc(primary_map->to->core_size*sizeof(int ));

      int count = 0; int rev_count = primary_map->to->core_size -1;
      for(int i = 0; i < primary_map->to->core_size; i++)
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
      for(int i = 0; i<primary_map->to->core_size; i++)
        for(int j = 0; j<adj_i[i]; j++)
          adj[i][j]= new_index[adj[i][j]];

      //move elements in adj[] to the new locations
      int** adj_new = (int **)xmalloc(primary_map->to->core_size*sizeof(int *));
      int* adj_i_new = (int *)xmalloc(primary_map->to->core_size*sizeof(int ));
      for(int i = 0; i<primary_map->to->core_size; i++)
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

      //create csr table with *only* the connected set elements
      idxtype *xadj = (idxtype *)xmalloc(sizeof(idxtype)*(rev_count+2));
      int cap = (rev_count+1)*primary_map->dim;
      idxtype *adjncy = (idxtype *)xmalloc(sizeof(idxtype)*cap);

      count = 0;int prev_count = 0;
      for(int i = 0; i<rev_count+1; i++)
      {
        int l_index = i;
        if(adj_i_new[i] > 1)
        {
          for(int j = 0; j<adj_i_new[i]; j++)
          {
            if(adj_new[i][j] != l_index)
            {
              if(count >= cap)
              {
                cap = cap*2;
                adjncy = (idxtype *)xrealloc(adjncy,sizeof(idxtype)*cap);
              }
              adjncy[count++] = (idxtype)adj_new[i][j];
            }
          }
          if(i != 0)
          {
            xadj[i] = prev_count;
            prev_count = count;
          }
          else
          {
            xadj[i] = 0;
            prev_count = count;
          }
        }
      }
      xadj[rev_count+1] = count;

      adjncy = (idxtype *)xrealloc(adjncy,sizeof(idxtype)*count);

      idxtype *p = (idxtype *)xmalloc(sizeof(idxtype)*rev_count+1);
      for(int i = 0; i < rev_count+1; i++){ p[i] = -99; }

      idxtype ncon = 1;
      int edgecut = 0;
      int part = num_part;
      idxtype wgtflag = 0;
      int options[5] = {0};
      idxtype numflag = 0;
      int nn = rev_count+1;

      //printf("for set %s set_size %d, set_core_size %d, partitioned size %d\n",
      //  set->name, set->size, set->core_size, nn);

      //free memory before calling METIS
      for(int i = 0; i<primary_map->to->core_size; i++)free(adj[i]);
      free(adj_i);free(adj_i_new);free(adj_cap);free(adj);free(adj_new);

      //call partitioner with the number of mini-partitions required and get the partition info array
      METIS_PartGraphKway(&nn, xadj, adjncy, NULL, NULL, &wgtflag, &numflag, &part, options, &edgecut, p);

      //free memory used for METIS
      free(xadj); free(adjncy);

      //create index array to hold new element indecies
      int* new_index_2 = (int *)xmalloc(primary_map->to->core_size*sizeof(int ));
      int* temp = (int *)xmalloc(primary_map->to->core_size*sizeof(int ));
      for(int i = 0; i<primary_map->to->core_size; i++)new_index_2[i] = i;

      //sort new_index_2 according to parition information in p[]
      quickSort_2(p, new_index_2, 0, nn-1);

      //reveres new_index_2
      for(int i = 0; i<primary_map->to->core_size; i++)
        temp[new_index_2[i]] = i;
      free(new_index_2);
      new_index_2 = temp;

      //update new_index using new_index_2
      for(int i = 0; i<primary_map->to->core_size; i++)
        new_index[i] = new_index_2[new_index[i]];
      free(new_index_2);
      free(p);

      //store reorder information in reordered_index array
      reordered_index[set->index] = new_index;
      reordered[set->index] = 1;
    }
    else //cound not find a suitable mapping table to do the reordering
    {
      printf("Set %s not reordered : could not find suitable mapping or reverse mapping\n",
        set->name);
    }
  }

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

  //update the original global index array (if it exists)
  //with the new possitions so that partition reversal can be done at the end
  if(OP_part_index == OP_set_index)
  {
    for(int s=0; s<OP_set_index; s++) { //for each set
      op_set set=OP_set_list[s];
      if(reordered[set->index] > 0)
      {
        int *temp= (int *)xmalloc(sizeof(int)*set->size);

        for(int i = 0; i<set->size; i++)
        {
          if(i<set->core_size)
            //printf("%d - %d\n",i, reordered_index[set->index][i]);
            temp[reordered_index[set->index][i]] = OP_part_list[set->index]->g_index[i];
        }
        free(OP_part_list[set->index]->g_index);
        OP_part_list[set->index]->g_index = temp;
      }
    }
  }
  else //If not create global index array as this is going to be for single node backend
  {
    //TODO
  }

  //cleanup
  for(int i = 0; i<OP_set_index;i++)if(reordered[i]>0)free(reordered_index[i]);
  free(reordered_index);free(reordered);
}
