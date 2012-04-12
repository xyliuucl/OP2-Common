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

/*******************************************************************************
 * Routine to do local renumbering using Metis k-way
 *******************************************************************************/

void op_local_renumbering_metiskway()
{
  //Need OP_reverse_map_list[OP_map_index] to hold the reverse mappings
  op_map * OP_reverse_map_list;
  int OP_reverse_map_index = 0;

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

  //for each set need to do the local renumbering for the core elements
  for(int s=0; s<OP_set_index; s++) {
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

    //use each candidate mapping found and
    //create list with each index holding an array with the neighbors of this set element
    //select the list that gives a minimum lone nodes/vertices, max dimension of the map used in that order
    /*for (int i = 0; i< found; i++)
      {
      op_map map = OP_map_list[candidate_maps[i]];
      int lone_vertices = 0;
      int dim = map->dim;
      if(compare_sets(map->to,set)==1) //partition to set
      {

      }
      }*/

    if(candidate_maps[0] != -1)
    {
      op_map primary_map = OP_map_list[candidate_maps[0]];

      int** adj = (int **)xmalloc(primary_map->to->core_size*sizeof(int *));
      int* adj_i = (int *)xmalloc(primary_map->to->core_size*sizeof(int ));
      int* adj_cap = (int *)xmalloc(primary_map->to->core_size*sizeof(int ));

      for(int i = 0; i<primary_map->to->core_size; i++)adj_i[i] = 0;
      for(int i = 0; i<primary_map->to->core_size; i++)adj_cap[i] = primary_map->dim;
      for(int i = 0; i<primary_map->to->core_size; i++)adj[i] = (int *)xmalloc(adj_cap[i]*sizeof(int));

      //go through each from-element of local primary_map and construct adjacency list
      for(int i = 0; i<primary_map->from->core_size; i++)
      {
        int part, local_index;
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
              adj[local_index][adj_i[local_index]++] =
                primary_map->map[i*primary_map->dim+k];
            }
          }
        }
      }

      //look for lone nodes/vertices
      int lone_vertices_index = 0;
      int* lone_vertices = (int *)xmalloc(primary_map->to->core_size*sizeof(int ));
      for(int i = 0; i<primary_map->to->core_size; i++)
      {
        if(adj_i[i] == 0)
        {
          lone_vertices[lone_vertices_index++] = i;
        }
      }
      printf("lone vertices from using map %s for set %s = %d\n",
          primary_map->name, set->name, lone_vertices_index);

      //create array that holds the original local-index of each set element
      int* local_index = (int *)xmalloc(primary_map->to->core_size*sizeof(int ));
      for(int i = 0; i<primary_map->to->core_size; i++)local_index[i] = i;

      //move non-connected set elements to the end of the list -
      //update local-index array to point to new possitions of the elements
      int count = 0;
      for(int i = 0; i<primary_map->to->core_size; i++)
      {
        if(adj_i[i] == 0)
        {
          count++;
          local_index[count] = local_index[i];
        }
      }
      for(int i = 0; i< lone_vertices_index; i++) local_index[count] = lone_vertices[i];

      //create csr table with *only* the connected set elements
      idxtype *xadj = (idxtype *)xmalloc(sizeof(idxtype)*(primary_map->to->core_size+1));
      int cap = (primary_map->to->size)*primary_map->dim;
      idxtype *adjncy = (idxtype *)xmalloc(sizeof(idxtype)*cap);

      count = 0;int prev_count = 0;
      for(int i = 0; i<primary_map->to->core_size; i++)
      {
        int l_index = i;
        quickSort(adj[i], 0, adj_i[i]-1);
        adj_i[i] = removeDups(adj[i], adj_i[i]);

        if(adj_i[i] != 0)
        {
          adj[i] = (int *)xrealloc(adj[i],adj_i[i]*sizeof(int));
          for(int j = 0; j<adj_i[i]; j++)
          {
            if(adj[i][j] != l_index)
            {
              if(count >= cap)
              {
                cap = cap*2;
                adjncy = (idxtype *)xrealloc(adjncy,sizeof(idxtype)*cap);
              }
              adjncy[count++] = (idxtype)adj[i][j];
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
      xadj[primary_map->to->size] = count;

      for(int i = 0; i<primary_map->to->core_size; i++)free(adj[i]);
      free(adj_i);free(adj_cap);free(adj);


      idxtype *partition = (idxtype *)xmalloc(sizeof(idxtype)*primary_map->to->core_size);
      for(int i = 0; i < primary_map->to->core_size; i++){ partition[i] = -99; }

      idxtype ncon = 1;
      int edgecut = 0;
      int part = 1024;

      METIS_PartGraphKway(&count, xadj, adjncy, NULL, NULL, NULL, NULL, &part, NULL, &edgecut, partition);

      //call partitioner with the number of mini-partitions required and get the partition info array

      //move elements in local-index array to the correct mini-partitions and sort them within each mini-partition
      //according to original index

      //add the lone-nodes to the end of the local-index array (perhaps already there) and sort them also

      //using the possition of the local-index array move the data and maps on this set to the correct possition
    }
  }

  //go thorugh each mapping table - both core and non-core and exec-halo elements
  //renumber the to set with the new possitions of the elements

  //update the original global index array (if it exists)
  //with the new possitions so that partition revers can be done at the end

  //If not create global index array as this is going to be for single node backend
}
