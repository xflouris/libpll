/*
    Copyright (C) 2016 Tomas Flouri

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as
    published by the Free Software Foundation, either version 3 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    Contact: Tomas Flouri <Tomas.Flouri@h-its.org>,
    Exelixis Lab, Heidelberg Instutute for Theoretical Studies
    Schloss-Wolfsbrunnenweg 35, D-69118 Heidelberg, Germany
*/

#include "pll.h"

/* simulate exactly the non-reentrant glibc srandom() function */
#define RAND_STATE_SIZE 128

typedef struct
{
  int clv_valid;
} node_info_t;

static pll_unode_t ** travbuffer;
static pll_pars_buildop_t * parsops;

static char * xstrdup(const char * s)
{
  size_t len = strlen(s);
  char * p = (char *)malloc(len+1);
  if (!p)
  {
    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf(pll_errmsg, 200, "Memory allocation failed");
    return NULL;
  }
  return strcpy(p,s);
}

/* Fisher-Yates shuffle */
static unsigned int * create_shuffled(unsigned int n, unsigned int seed)
{
  unsigned int i,j;
  char * statebuf;
  struct pll_random_data * buf;
  
  unsigned int * x = (unsigned int *)malloc(n*sizeof(unsigned int));
  if (!x)
  {
    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf(pll_errmsg, 200, "Unable to allocate enough memory.");
    return NULL;
  }

  for (i=0; i<n; ++i)
    x[i] = i;

  /* if seed == 0 then do not shuffle! */
  if (!seed)
    return x;

  /* init re-entrant randomizer */
  buf = (struct pll_random_data *)calloc(1, sizeof(struct pll_random_data));
  statebuf = (char *)calloc(RAND_STATE_SIZE,sizeof(char));

  pll_initstate_r(seed,statebuf,RAND_STATE_SIZE,buf);
  pll_srandom_r(seed,buf);

  /* perform Fisher-Yates shuffle */
  if (n > 1)
  {
    i = n - 1;
    while (1)
    {
      int rint;
      pll_random_r(buf,&rint);
      double r = ((double)rint / RAND_MAX);
      j = (unsigned int)(r * (i+1));

      PLL_SWAP(x[i],x[j]);

      if (i == 0) break;
      --i;
    }
  }

  /* dealloc and return shuffled array */
  free(statebuf);
  free(buf);
  return x;
}

static void dealloc_data_onenode(pll_unode_t * node)
{
  if (node->data)
  {
    free(node->data);
    node->data = NULL;
  }
}

static void dealloc_data(pll_unode_t * node)
{
  dealloc_data_onenode(node);
  dealloc_data_onenode(node->next);
  dealloc_data_onenode(node->next->next);
}

/* a callback function for performing a partial traversal */
static int cb_partial_traversal(pll_unode_t * node)
{
  node_info_t * node_info;

  /* if we don't want tips in the traversal we must return 0 here. For now,
     allow tips */
  if (!node->next) return 1;

  /* get the data element from the node and check if the CLV vector is
     oriented in the direction that we want to traverse. If the data
     element is not yet allocated then we allocate it, set the direction
     and instruct the traversal routine to place the node in the traversal array
     by returning 1 */
  node_info = (node_info_t *)(node->data);
  if (!node_info)
  {
    /* allocate data element. TODO: Check whether allocation was successful */
    node->data             = (node_info_t *)calloc(1,sizeof(node_info_t));
    node->next->data       = (node_info_t *)calloc(1,sizeof(node_info_t));
    node->next->next->data = (node_info_t *)calloc(1,sizeof(node_info_t));

    /* set orientation on selected direction and traverse the subtree */
    node_info = node->data;
    node_info->clv_valid = 1;
    return 1;
  }

  /* if the data element was already there and the CLV on this direction is
     set, i.e. the CLV is valid, we instruct the traversal routine not to
     traverse the subtree rooted in this node/direction by returning 0 */
  if (node_info->clv_valid) return 0;

  /* otherwise, set orientation on selected direction */
  node_info->clv_valid = 1;

  /* reset orientation on the other two directions and return 1 (i.e. traverse
     the subtree */
  node_info = node->next->data;
  node_info->clv_valid = 0;
  node_info = node->next->next->data;
  node_info->clv_valid = 0;

  return 1;
}

static pll_unode_t * utree_inner_create(unsigned int i)
{
  pll_unode_t * node = (pll_unode_t *)calloc(1,sizeof(pll_unode_t));
  if (!node)
    return NULL;
  
  node->next = (pll_unode_t *)calloc(1,sizeof(pll_unode_t));
  if (!node->next)
  {
    free(node);
    return NULL;
  }
  node->next->next = (pll_unode_t *)calloc(1,sizeof(pll_unode_t));
  if (!node->next->next)
  {
    free(node->next);
    free(node);
    return NULL;
  }

  node->next->next->next = node;

  node->clv_index = i;
  node->next->clv_index = i;
  node->next->next->clv_index = i;

  return node;
}

static pll_unode_t * utree_tip_create(unsigned int i)
{
  pll_unode_t * node = (pll_unode_t *)calloc(1,sizeof(pll_unode_t));
  node->next = NULL;
  node->clv_index = i;

  return node;
}

static void utree_link(pll_unode_t * a, pll_unode_t * b)
{
  /*

    *               *               *                * 
     \             /                 \              /
      *---*   *---*        -->        *---*-----*--*
     /    a   b    \                 /    a     b   \
    *               *               *                *

  */

  a->back = b;
  b->back = a;
}

static void utree_edgesplit(pll_unode_t * a, pll_unode_t * b, pll_unode_t * c)
{
  /*
                *                                      *
                |                                      |
                *                                      *
               / \                                    / \
            b *   * c                              b *   * c
                                                    /     \
    *                      *      -->      *       /       \      * 
     \                    /                 \     /         \    /
      *---*----------*---*                   *---*           *--*
     /    a          d    \                 /    a           d   \
    *                      *               *                      *

  */

  /* link d<->c */
  utree_link(a->back,c);

  /* link a<->b */
  utree_link(a,b);
}

static unsigned int utree_iterate(pll_parsimony_t ** list,
                                  pll_unode_t ** edge_list,
                                  pll_unode_t * inner_node,
                                  pll_unode_t * tip_node,
                                  unsigned int edge_count,
                                  unsigned int partition_count)
{
  unsigned int i,j;
  unsigned int min_cost = 0;
  unsigned int best_index = 0;
  unsigned int cost;
  unsigned int ops_count;
  unsigned int traversal_size;

  /* set min cost to maximum possible value */
  min_cost = ~0u;

  /* find first empty slot in edge_list */
  pll_unode_t ** empty_slot = edge_list + edge_count;

  for (i = 0; i < edge_count; ++i)
  {
    /* make the split */
    pll_unode_t * d = edge_list[i]->back;
    utree_edgesplit(edge_list[i], inner_node, inner_node->next); 
    utree_link(inner_node->next->next, tip_node);

    /* add the two new edges to the end of the list */
    empty_slot[0] = inner_node->next;
    empty_slot[1] = inner_node->next->next;

    /* make a partial traversal */
    if (!pll_utree_traverse(tip_node->back,
                            PLL_TREE_TRAVERSE_POSTORDER,
                            cb_partial_traversal,
                            travbuffer,
                            &traversal_size))
      assert(0);

    /* create parsimony operations */
    pll_utree_create_pars_buildops(travbuffer,
                                   traversal_size,
                                   parsops,
                                   &ops_count);

    /* compute the costs for each parsimony partition */
    cost = 0;
    for (j = 0; j < partition_count; ++j)
    {
      /* update parsimony vectors */
      pll_fastparsimony_update_vectors(list[j], parsops, ops_count);

      /* get parsimony score */
      cost += pll_fastparsimony_edge_score(list[j],
                                           tip_node->clv_index,
                                           tip_node->back->clv_index);
    }

    /* if current cost is smaller than minimum cost save topology index */
    if (cost < min_cost)
    {
      min_cost = cost;
      best_index = i;
    }

    /* reset direction for the newly placed inner node */
    node_info_t * node_info = (node_info_t *)(inner_node->next->next->data);
    node_info->clv_valid = 0;

    /* restore tree to its state before placing the tip (and inner) node */
    utree_link(edge_list[i], d);
    inner_node->back = NULL;
    inner_node->next->back = NULL;
    inner_node->next->next->back = NULL;
    tip_node->back = NULL;
  }

  /* perform the placement yielding the lowest cost */
  utree_edgesplit(edge_list[best_index], inner_node, inner_node->next); 
  utree_link(inner_node->next->next, tip_node);

  return min_cost;
}

static void invalidate_node(pll_unode_t * node)
{
  node_info_t * info;

  info = (node_info_t *)(node->data);
  info->clv_valid = 0;
  info = (node_info_t *)(node->next->data);
  info->clv_valid = 0;
  info = (node_info_t *)(node->next->next->data);
  info->clv_valid = 0;
}

PLL_EXPORT pll_utree_t * pll_fastparsimony_stepwise(pll_parsimony_t ** list,
                                                    char * const * labels,
                                                    unsigned int * cost,
                                                    unsigned int count,
                                                    unsigned int seed)
{
  unsigned int i,j;

  unsigned int tips_count = list[0]->tips;
  unsigned int inner_nodes = list[0]->inner_nodes;

  if (tips_count < 3)
  {
    pll_errno = PLL_ERROR_STEPWISE_TIPS;
    snprintf(pll_errmsg, 200,
             "Stepwise parsimony requires at least three tips.");
    return NULL;
  }

  //if (tips_count != inner_nodes + 2)
  if (inner_nodes < tips_count-2)
  {
    pll_errno = PLL_ERROR_STEPWISE_UNSUPPORTED;
    snprintf(pll_errmsg, 200,
             "Stepwise parsimony currently supports only unrooted trees.");
    return NULL;
  }

  *cost = ~0u;


  pll_unode_t * root;

  /* check that all parsimony structures have the same number of tips and
     inner nodes */

  for (i = 1; i < count; ++i)
  {
    if ((list[i]->tips != tips_count) ||
        (list[i]->inner_nodes != inner_nodes))
    {
      pll_errno = PLL_ERROR_STEPWISE_STRUCT;
      snprintf(pll_errmsg, 200,
               "Parsimony structures tips/inner nodes not equal.");
      return NULL;
    }
  }
    

  /* 1. Make all allocations at the beginning and check everything was
        allocated, otherwise return an error */

  travbuffer = (pll_unode_t **)malloc((2*tips_count-2) * sizeof(pll_unode_t *));

  root = utree_inner_create(2*tips_count-3);

  /* allocate parsimony operations container */
  parsops = (pll_pars_buildop_t *)malloc((tips_count-2)*
                                         sizeof(pll_pars_buildop_t));

  /* create tip node list with a terminating NULL element */
  pll_unode_t ** tip_node_list = (pll_unode_t **)calloc(tips_count+1,
                                                        sizeof(pll_unode_t *));

  /* create inner node list for (tips_count - 3) inner nodes (root was already
     created, and leave the last slot NULL for termination */
  pll_unode_t ** inner_node_list = (pll_unode_t **)calloc(tips_count - 2,
                                                          sizeof(pll_unode_t *));

  if (!inner_node_list || !parsops || !tip_node_list || !root || !travbuffer)
  {
    pll_utree_graph_destroy(root,NULL);
    free(parsops);
    free(inner_node_list);
    free(tip_node_list);
    free(travbuffer);

    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf(pll_errmsg, 200, "Unable to allocate enough memory.");
    return NULL;
  }

  /* allocate all inner nodes */
  for (i=0; i<tips_count-3; ++i)
  {
    inner_node_list[i] = utree_inner_create(i+tips_count);
    if (!inner_node_list[i])
    {
      pll_utree_graph_destroy(root,NULL);
      free(parsops);
      free(tip_node_list);
      free(travbuffer);
      for (j = 0; j < i; ++j)
        pll_utree_graph_destroy(inner_node_list[j],NULL);
      free(inner_node_list);

      pll_errno = PLL_ERROR_MEM_ALLOC;
      snprintf(pll_errmsg, 200, "Unable to allocate enough memory.");
      return NULL;
    }
  }

  /* shuffle the order of iterating tip sequences */
  unsigned int * order = create_shuffled(tips_count,seed);
  if (!order) return NULL;

  /* allocate all tips */
  for (i=0; i<tips_count; ++i)
  {
    unsigned int index = order[i];
    tip_node_list[i] = utree_tip_create(index);
    if (tip_node_list[i])
      tip_node_list[i]->label = xstrdup(labels[index]);

    if (!tip_node_list[i] || !tip_node_list[i]->label)
    {
      free(tip_node_list[i]); 

      pll_utree_graph_destroy(root,NULL);
      free(parsops);
      free(inner_node_list);
      free(travbuffer);
      for (j = 0; j < i; ++j)
        pll_utree_graph_destroy(tip_node_list[j],NULL);
      free(tip_node_list);

      pll_errno = PLL_ERROR_MEM_ALLOC;
      snprintf(pll_errmsg, 200, "Unable to allocate enough memory.");
      return NULL;

    }
  }
  free(order);

  /* 2. Create the following topology with three leaves
          
            *
           /
      *---*
           \
            *
  
  */

  /* place first three tips */
  utree_link(root, tip_node_list[0]);
  utree_link(root->next, tip_node_list[1]);
  utree_link(root->next->next, tip_node_list[2]);

  /* available placements */
  pll_unode_t ** edge_list = (pll_unode_t **)calloc(2*tips_count-3,
                                                    sizeof(pll_unode_t *));
  edge_list[0] = root;
  edge_list[1] = root->next;
  edge_list[2] = root->next->next;

  /* 3. The stepwise parsimony. Current topology is the tree with three leaves,
        and repeat the following steps for each remaining tip u:
         (i) compute the parsimony score of all possible topologies by placing
             u at every possible edge of the current tree topology.
        (ii) set current toplogy as the tree with the smallest parsimony score
  */

  if (tips_count > 3)
  {
    unsigned int edge_count = 3;
    
    for (i = 3; i < tips_count; ++i)
    {
      /* printf("%d -- adding %s\n", i, tip_node_list[i]->label); */
      *cost = utree_iterate(list,
                            edge_list,
                            inner_node_list[i-3],
                            tip_node_list[i],
                            edge_count,
                            count);

      /* reset traversal such that parsimony vectors are re-computed */
      for (j = 0; j < i-2; ++j)
        invalidate_node(inner_node_list[j]);
      invalidate_node(root);

      /* after adding a leaf, we have two new edges */
      edge_count += 2;
    }
  }
  else
  {
    *cost = 0;
    for (i = 0; i < count; ++i)
      *cost += list[i]->const_cost;
  }

  /* delete data elements */
  for (i = 0; i < tips_count-3; ++i)
    dealloc_data(inner_node_list[i]);
  dealloc_data(root);

  /* deallocate auxiliary arrays */
  free(inner_node_list);
  free(tip_node_list);
  free(edge_list);
  free(travbuffer);
  free(parsops);

  /* wrap tree */
  pll_utree_t * tree = pll_utree_wraptree(root,tips_count);

  return tree;
}
