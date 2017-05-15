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

#define PLL_BITVECTOR_SIZE 32

static int alloc_pars_structs(pll_parsimony_t * parsimony,
                              unsigned int bitvectors)
{
  unsigned int i,j;

  /* TODO: Test this for compatibility with rooted and unrooted trees */
  unsigned int nodes_count = parsimony->tips + parsimony->inner_nodes;

  parsimony->node_cost = (unsigned int *)calloc(nodes_count,
                                                sizeof(unsigned int));
  if (!parsimony->node_cost)
  {
    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf(pll_errmsg, 200,
             "Cannot allocate parsimony cost array.");
    return PLL_FAILURE;
  }

  /* allocate parsimony vector container */
  parsimony->packedvector = (unsigned int **)malloc(nodes_count *
                                                    sizeof(unsigned int *));
  if (!parsimony->packedvector)
  {
    free(parsimony->node_cost);
    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf (pll_errmsg, 200,
              "Cannot allocate parsimony vector container.");
    return PLL_FAILURE;
  }

  /* allocate individual vectors */
  unsigned int ** vector = parsimony->packedvector;
  for (i = 0; i < nodes_count; ++i)
  {
    vector[i] = (unsigned int *)pll_aligned_alloc(parsimony->states *
                                       bitvectors *
                                       sizeof(unsigned int),
                                       parsimony->alignment);
    if (!vector[i])
    {
      free(parsimony->node_cost);
      pll_errno = PLL_ERROR_MEM_ALLOC;
      snprintf (pll_errmsg, 200,
                "Cannot allocate parsimony vector.");
      
      /* free all allocated vectors */
      for (j = 0; j < i; ++j)
        pll_aligned_free(vector[i]);
      free(vector);

      return PLL_FAILURE;
    }
  }
  return PLL_SUCCESS;
}

static int check_informative_extended(const pll_partition_t * partition,
                                      unsigned int index,
                                      unsigned int * singleton)
{
  unsigned int c;
  int count = 0;
  unsigned int * map;
  unsigned int i,j;
  unsigned int range = (1 << partition->states);

  /* TODO: Move allocation outside of function such that it happens only once */
  //map = (int *)malloc(1024*1024*sizeof(int));
  map = (unsigned int *)malloc(range*sizeof(unsigned int));

  memset(map,0,range*sizeof(unsigned int));

  for (i = 0; i < partition->tips; ++i)
  {
    c = 0;

    double * clv = partition->clv[i] +
                   index*partition->states_padded*partition->rate_cats;

    for (j = 0; j < partition->states; ++j)
       c = (c << 1) | (unsigned int)(clv[j]);

    map[c]++;
  }

  /* TODO: change to only required range */
  for (i=0; i< range; ++i)
    if (map[i] > 1)
      count++;
    else if (map[i] == 1)
      (*singleton)++;

  free(map);

  if (count <= 1)
    return 0;

  return 1;
}

static int check_informative(const pll_partition_t * partition,
                             unsigned int index,
                             unsigned int * singleton)
{
  unsigned int i,j;
  unsigned int map[256];
  int count = 0;
  unsigned int c;

  /* in case the site is non-informative, count the number of states appearing
     only once, and which is equal to the number of mutations */
  *singleton = 0;

  for (i = 0; i < 256; ++i)
    map[i] = 0;

  /* if tips states are presented by characters */
  if (partition->attributes & PLL_ATTRIB_PATTERN_TIP)
  {
    for (i = 0; i < partition->tips; ++i)
    {
      c = partition->tipchars[i][index];
      map[c]++;
    }
  }
  else
  {
    /* otherwise, tips are represented by conditional probabilities */

    /* TODO: Currently this is extremely time-consuming as we have to allocate
    an array of 2^states elements each time we call check_informative_extended
    (i.e. for every character). For more than 20 states, this process becomes
    extremely slow, and we need to think of a better way (or completely disallow
    the non pattern tip case */
    assert(partition->states <= 20);
    if (partition->states > 8)
      return check_informative_extended(partition, index, singleton);

    for (i = 0; i < partition->tips; ++i)
    {
      c = 0;

      double * clv = partition->clv[i] +
                     index*partition->states_padded*partition->rate_cats;

      for (j = 0; j < partition->states; ++j)
         c = (c << 1) | (unsigned int)(clv[j]);

      map[c]++;
    }
  }


  /* TODO: change to only required range */
  for (i=0; i<256; ++i)
    if (map[i] > 1)
      count++;
    else if (map[i] == 1)
      (*singleton)++;

  if (count <= 1)
    return 0;

  return 1;
}

static int fill_parsimony_vectors(const pll_partition_t * partition,
                                  pll_parsimony_t * parsimony)
{
  unsigned int c;
  unsigned int i,j,k;
  unsigned int bitcount = 0;
  unsigned int bitvectors;

  unsigned int states = parsimony->states;

  /*
     Example: assume the following tree and encoding

       /\        A = 0001    seq1 = TTAACT
      /  \       C = 0010    seq2 = TTCAGG
     /\  /\      G = 0100    seq3 = CCGACT
    1  23  4     T = 1000    seq4 = CCTAGG
                      informative = 110011 

    Create 32-bit state vectors for each tip:

    node1:
      vecA    00110011111111111111111111111111
      vecC    00001011111111111111111111111111
      vecG    00000011111111111111111111111111
      vecT    11000111111111111111111111111111

    node2:
      vecA    00010011111111111111111111111111
      vecC    00100011111111111111111111111111
      vecG    00001111111111111111111111111111
      vecT    11000011111111111111111111111111

    node3:
      vecA    00010011111111111111111111111111
      vecC    11001011111111111111111111111111
      vecG    00100011111111111111111111111111
      vecT    00000111111111111111111111111111

    node4:
      vecA    00010011111111111111111111111111
      vecC    11000011111111111111111111111111
      vecG    00001111111111111111111111111111
      vecT    00100011111111111111111111111111
  */

  /* compute total number of bits required */
  for (i=0; i<parsimony->sites; ++i)
    if (parsimony->informative[i])
      bitcount += partition->pattern_weights[i];

  /* number of 32-bit bit-vectors required */
  bitvectors = (bitcount / PLL_BITVECTOR_SIZE) +
               (bitcount % PLL_BITVECTOR_SIZE != 0);
  
#ifdef HAVE_SSE3
  if (parsimony->attributes & PLL_ATTRIB_ARCH_SSE && PLL_STAT(sse3_present))
    bitvectors = (bitvectors+3) & 0xFFFFFFFC;
#endif

#ifdef HAVE_AVX
  if (parsimony->attributes & PLL_ATTRIB_ARCH_AVX && PLL_STAT(avx_present))
    bitvectors = (bitvectors+7) & 0xFFFFFFF8;
#endif

#ifdef HAVE_AVX2
  if (parsimony->attributes & PLL_ATTRIB_ARCH_AVX2 && PLL_STAT(avx2_present))
    bitvectors = (bitvectors+7) & 0xFFFFFFF8;
#endif
  
  /* allocate necessary data structures */
  if (!alloc_pars_structs(parsimony, bitvectors))
    return PLL_FAILURE;
  
  unsigned int ** statevec = (unsigned int **)malloc(states *
                                                     sizeof(unsigned int *));
  unsigned int * val = (unsigned int *)calloc(states,sizeof(unsigned int));
  if (!val || !statevec)
  {
    if (val) free(val);
    if (statevec) free(statevec);

    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf (pll_errmsg, 200,
              "Cannot allocate bitvector data.");
    return PLL_FAILURE;
  }

  /* TODO: Separate tip pattern and clv case */
  for (i = 0; i < parsimony->tips; ++i)
  {
    for (k = 0; k < parsimony->states; ++k) val[k] = 0;
    bitcount = 0;

    for (j = 0; j < parsimony->states; ++j)
      statevec[j] = parsimony->packedvector[i]+bitvectors*j;

    unsigned int vec_index = 0;
    for (j = 0; j < parsimony->sites; ++j)
    {
      if (parsimony->informative[j])
      {
        /* TODO: check for tip map */
        unsigned int m;
        for (m = 0; m < partition->pattern_weights[j]; ++m)
        {
          if (parsimony->attributes & PLL_ATTRIB_PATTERN_TIP)
          {
            c = partition->tipchars[i][j];
            if (states != 4) c = partition->tipmap[c];
            for (k = 0; k < parsimony->states; ++k, c >>= 1)
              if (c & 1) 
                val[k] |= (1 << bitcount);
          }
          else
          {
            double * clv = partition->clv[i] + 
                           j*partition->states_padded * partition->rate_cats;

            for (k = 0; k < states; ++k)
              if ((int)(clv[k]))
              {
                val[k] |= (1 << bitcount);
              }
          }

          bitcount++;

          if (bitcount == PLL_BITVECTOR_SIZE)
          {
            for (k = 0; k < states; ++k)
            {
              statevec[k][vec_index] = val[k];
              val[k] = 0;
            }

            vec_index++;
            bitcount = 0;
          }
        }

      }
    }

    /* fill up the remaining bit entries in the current vector with ones */
    if (bitcount && (bitcount != PLL_BITVECTOR_SIZE))
    {
      for (; bitcount < PLL_BITVECTOR_SIZE; ++bitcount)
        for (k = 0; k < states; ++k)
          val[k] |= (1 << bitcount);
      
      for (k = 0; k < states; ++k)
        statevec[k][vec_index] = val[k];

      vec_index++;
    }

    /* fill up the remaining bit vectors due to padding with ones */
    for (; vec_index < bitvectors; ++vec_index)
      for (k = 0; k < states; ++k)
        statevec[k][vec_index] = ~0u;
  }
  parsimony->packedvector_count = bitvectors;

  free(statevec);
  free(val);

  return PLL_SUCCESS;
}

static int pll_set_informative(const pll_partition_t * partition,
                               pll_parsimony_t * parsimony)
{
  unsigned int i;
  unsigned int singletons = 0;
  unsigned int count = 0;

  /* allocate array for indicating whether a site is informative or not */
  parsimony->informative = (int *)malloc(parsimony->sites * sizeof(int));
  if (!parsimony->informative)
  {
    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf (pll_errmsg, 200,
              "Cannot allocate informative array.");
    return PLL_FAILURE;
  }

  /* identify and mark informative sites */
  for (i = 0; i < parsimony->sites; ++i)
  {
    if (check_informative(partition,i,&singletons))
      parsimony->informative[i] = 1;
    else
    {
      parsimony->informative[i] = 0;
      count++;
      parsimony->const_cost += singletons *
                               partition->pattern_weights[i];
    }
  }

  parsimony->informative_count = parsimony->sites - count;

  return PLL_SUCCESS;
}

PLL_EXPORT unsigned int pll_fastparsimony_edge_score_4x4(const pll_parsimony_t * parsimony,
                                                         unsigned int node1_score_index,
                                                         unsigned int node2_score_index)
{
  unsigned int i;

  unsigned int * node1[4];
  unsigned int * node2[4];

  unsigned int ** vector = parsimony->packedvector;
  unsigned int vector_count = parsimony->packedvector_count;

  unsigned int score = 0;

  /* point to the parsimony vectors for each node and for each state */
  for (i = 0; i < 4; ++i)
  {
    node1[i] = vector[node1_score_index] + i*vector_count;
    node2[i] = vector[node2_score_index] + i*vector_count;
  }

  unsigned int xmm0,xmm1,xmm2,xmm3,xmm4,xmm5,xmm6;

  /* set all bits to one */
  xmm1 = ~0u;

  for (i = 0; i < parsimony->packedvector_count; ++i)
  {
    /* compute AND bit vectors for state 0 */
    xmm2 = node1[0][i] & node2[0][i];

    /* compute AND bit vectors for state 1 */
    xmm3 = node1[1][i] & node2[1][i];

    /* compute AND bit vectors for state 2 */
    xmm4 = node1[2][i] & node2[2][i];

    /* compute AND bit vectors for state 3 */
    xmm5 = node1[3][i] & node2[3][i];

    /* OR the ANDs of states 0,1,2,3 */
    xmm6 = xmm2 | xmm3 | xmm4 | xmm5;

    xmm0 = ~xmm6 & xmm1;

    score += (unsigned int)__builtin_popcount(xmm0);
  }
  unsigned int score1 = parsimony->node_cost[node1_score_index];
  unsigned int score2 = parsimony->node_cost[node2_score_index];

  return score+score1+score2+parsimony->const_cost;
}

PLL_EXPORT void pll_fastparsimony_update_vector_4x4(pll_parsimony_t * parsimony,
                                                    const pll_pars_buildop_t * op)
{
  unsigned int i;
  unsigned int * parent[4];
  unsigned int * child1[4];
  unsigned int * child2[4];

  unsigned int ** vector = parsimony->packedvector;
  unsigned int vector_count = parsimony->packedvector_count;

  unsigned int score = 0;

  /* point to the parsimony vectors for each node and for each state */
  for (i = 0; i < 4; ++i)
  {
    parent[i] = vector[op->parent_score_index] + i*vector_count;
    child1[i] = vector[op->child1_score_index] + i*vector_count;
    child2[i] = vector[op->child2_score_index] + i*vector_count;
  }

  unsigned int xmm0,xmm1,xmm2,xmm3,xmm4,xmm5,xmm6,xmm7,xmm8,xmm9;

  /* set all bits to one */
  xmm9 = ~0u;

  for (i = 0; i < parsimony->packedvector_count; ++i)
  {
    /* load, and, or bit vectors for state 0 */

    xmm0 = child1[0][i] & child2[0][i];
    xmm1 = child1[0][i] | child2[0][i];

    /* load, and, or bit vectors for state 1 */

    xmm2 = child1[1][i] & child2[1][i];
    xmm3 = child1[1][i] | child2[1][i];

    /* load, and, or bit vectors for state 2 */

    xmm4 = child1[2][i] & child2[2][i];
    xmm5 = child1[2][i] | child2[2][i];

    /* load, and, or bit vectors for state 3 */

    xmm6 = child1[3][i] & child2[3][i];
    xmm7 = child1[3][i] | child2[3][i];
    
    /* OR the ANDs of states 0,1,2,3 */

    xmm8 = xmm0 | xmm2 | xmm4 | xmm6;

    /* store them */
    parent[0][i] = xmm0 | (~xmm8 & xmm1);
    parent[1][i] = xmm2 | (~xmm8 & xmm3);
    parent[2][i] = xmm4 | (~xmm8 & xmm5);
    parent[3][i] = xmm6 | (~xmm8 & xmm7);

    score += (unsigned int)__builtin_popcount(~xmm8 & xmm9);
  }
  unsigned int score1 = parsimony->node_cost[op->child1_score_index];
  unsigned int score2 = parsimony->node_cost[op->child2_score_index];
  parsimony->node_cost[op->parent_score_index] = score+score1+score2;
}

PLL_EXPORT pll_parsimony_t * pll_fastparsimony_init(const pll_partition_t * partition)
{
  pll_parsimony_t * parsimony;

  /* TODO: Currently only upto 20 states are supported with non pattern-tip
     compression */
  if (partition->states > 20 && 
      (partition->attributes & PLL_ATTRIB_PATTERN_TIP) == 0)
  {
    pll_errno = PLL_ERROR_STEPWISE_UNSUPPORTED;
    snprintf(pll_errmsg,
             200,
             "Use PLL_ATTRIB_PATTERN_TIP for more than 20 states.");
    return NULL;
  }
  
  parsimony = (pll_parsimony_t *)calloc(1,sizeof(pll_parsimony_t));

  parsimony->tips = partition->tips;
  parsimony->inner_nodes = partition->tips-1;
  parsimony->sites = partition->sites;
  parsimony->attributes = partition->attributes;
  parsimony->states = partition->states;
  parsimony->alignment = partition->alignment;

  if (!pll_set_informative(partition,parsimony))
    return NULL;

  if (!fill_parsimony_vectors(partition,parsimony))
    return NULL;

  return parsimony;
}

PLL_EXPORT void pll_fastparsimony_update_vector(pll_parsimony_t * parsimony,
                                                const pll_pars_buildop_t * op)
{
  unsigned int i,j;
  unsigned int states = parsimony->states;
  unsigned int vector_count = parsimony->packedvector_count;
  unsigned int ** vector = parsimony->packedvector;

  unsigned int * child1;
  unsigned int * child2;
  unsigned int * parent;

  unsigned int score = 0;

  /* set all bits to one */
  unsigned int vones = ~0u;

  for (i = 0; i < parsimony->packedvector_count; ++i)
  {

    /* OR the ANDs of all states */
    unsigned int orvand = 0;
    child1 = vector[op->child1_score_index];
    child2 = vector[op->child2_score_index];
    for (j = 0; j < states; ++j)
    {
      orvand |= (child1[i] & child2[i]);

      child1 += vector_count;
      child2 += vector_count;
    }

    /* store vectors at parent */
    child1 = vector[op->child1_score_index];
    child2 = vector[op->child2_score_index];
    parent = vector[op->parent_score_index];
    for (j = 0; j < states; ++j)
    {
      parent[i] = (child1[i] & child2[i]) | 
                     (~orvand & (child1[i] | child2[i]));

      child1 += vector_count;
      child2 += vector_count;
      parent += vector_count;

    }

    score += (unsigned int)__builtin_popcount(~orvand & vones);
  }
  unsigned int score1 = parsimony->node_cost[op->child1_score_index];
  unsigned int score2 = parsimony->node_cost[op->child2_score_index];
  parsimony->node_cost[op->parent_score_index] = score+score1+score2;
}

static unsigned int fastparsimony_edge_score(const pll_parsimony_t * parsimony,
                                             unsigned int node1_score_index,
                                             unsigned int node2_score_index)
{
  unsigned int i,j;
  unsigned int states = parsimony->states;
  unsigned int vector_count = parsimony->packedvector_count;
  unsigned int ** vector = parsimony->packedvector;

  unsigned int * node1;
  unsigned int * node2;

  unsigned int score = 0;

  /* set all bits to one */
  unsigned int vones = ~0u;

  for (i = 0; i < parsimony->packedvector_count; ++i)
  {
    /* OR the ANDs of all states */
    unsigned int orvand = 0;
    node1 = vector[node1_score_index];
    node2 = vector[node2_score_index];
    for (j = 0; j < states; ++j)
    {
      orvand |= (node1[i] & node2[i]);
      
      node1 += vector_count;
      node2 += vector_count;
    }

    score += (unsigned int )__builtin_popcount(~orvand & vones);
  }
  unsigned int score1 = parsimony->node_cost[node1_score_index];
  unsigned int score2 = parsimony->node_cost[node2_score_index];

  return score+score1+score2+parsimony->const_cost;
}

static void fastparsimony_update_vectors_4x4(pll_parsimony_t * parsimony,
                                             const pll_pars_buildop_t * ops,
                                             unsigned int count)
{
  unsigned int i;
  const pll_pars_buildop_t * op;

  for (i = 0; i < count; ++i)
  {
    op = &(ops[i]);
#ifdef HAVE_SSE3
    if (parsimony->attributes & PLL_ATTRIB_ARCH_SSE && PLL_STAT(sse3_present))
      pll_fastparsimony_update_vector_4x4_sse(parsimony,op);
    else
#endif
#ifdef HAVE_AVX
    if (parsimony->attributes & PLL_ATTRIB_ARCH_AVX && PLL_STAT(avx_present))
      pll_fastparsimony_update_vector_4x4_avx(parsimony,op);
    else
#endif
#ifdef HAVE_AVX2
    if (parsimony->attributes & PLL_ATTRIB_ARCH_AVX2 && PLL_STAT(avx2_present))
      pll_fastparsimony_update_vector_4x4_avx2(parsimony,op);
    else
#endif
      pll_fastparsimony_update_vector_4x4(parsimony,op);
  }
}

static int fastparsimony_update_vectors(pll_parsimony_t * parsimony,
                                        const pll_pars_buildop_t * ops,
                                        unsigned int count)
{
  unsigned int i;
  const pll_pars_buildop_t * op;

  for (i = 0; i < count; ++i)
  {
    op = &(ops[i]);
#ifdef HAVE_SSE3
    if (parsimony->attributes & PLL_ATTRIB_ARCH_SSE && PLL_STAT(sse3_present))
      pll_fastparsimony_update_vector_sse(parsimony,op);
    else
#endif
#ifdef HAVE_AVX
    if (parsimony->attributes & PLL_ATTRIB_ARCH_AVX && PLL_STAT(avx_present))
      pll_fastparsimony_update_vector_avx(parsimony,op);
    else
#endif
#ifdef HAVE_AVX2
    if (parsimony->attributes & PLL_ATTRIB_ARCH_AVX2 && PLL_STAT(avx2_present))
      pll_fastparsimony_update_vector_avx2(parsimony,op);
    else
#endif
      pll_fastparsimony_update_vector(parsimony,op);
  }
  return PLL_SUCCESS;
}

PLL_EXPORT void pll_fastparsimony_update_vectors(pll_parsimony_t * parsimony,
                                                 const pll_pars_buildop_t * ops,
                                                 unsigned int count)
{
  if (parsimony->states == 4)
    fastparsimony_update_vectors_4x4(parsimony,ops,count);
  else
    fastparsimony_update_vectors(parsimony,ops,count);
}

PLL_EXPORT unsigned int pll_fastparsimony_edge_score(const pll_parsimony_t * parsimony,
                                                     unsigned int node1_score_index,
                                                     unsigned int node2_score_index)
{
  if (parsimony->states == 4)
  {
#ifdef HAVE_SSE3
    if (parsimony->attributes & PLL_ATTRIB_ARCH_SSE && PLL_STAT(sse3_present))
      return pll_fastparsimony_edge_score_4x4_sse(parsimony,
                                                  node1_score_index,
                                                  node2_score_index);
#endif
#ifdef HAVE_AVX
    if (parsimony->attributes & PLL_ATTRIB_ARCH_AVX && PLL_STAT(avx_present))
      return pll_fastparsimony_edge_score_4x4_avx(parsimony,
                                                  node1_score_index,
                                                  node2_score_index);
#endif
#ifdef HAVE_AVX2
    if (parsimony->attributes & PLL_ATTRIB_ARCH_AVX2 && PLL_STAT(avx2_present))
      return pll_fastparsimony_edge_score_4x4_avx2(parsimony,
                                                   node1_score_index,
                                                   node2_score_index);
#endif
    return pll_fastparsimony_edge_score_4x4(parsimony,
                                            node1_score_index,
                                            node2_score_index);
  }

#ifdef HAVE_SSE3
  if (parsimony->attributes & PLL_ATTRIB_ARCH_SSE && PLL_STAT(sse3_present))
    return pll_fastparsimony_edge_score_sse(parsimony,
                                            node1_score_index,
                                            node2_score_index);
  else
#endif
#ifdef HAVE_AVX
  if (parsimony->attributes & PLL_ATTRIB_ARCH_AVX && PLL_STAT(avx_present))
    return pll_fastparsimony_edge_score_avx(parsimony,
                                            node1_score_index,
                                            node2_score_index);
  else
#endif
#ifdef HAVE_AVX2
  if (parsimony->attributes & PLL_ATTRIB_ARCH_AVX2 && PLL_STAT(avx2_present))
    return pll_fastparsimony_edge_score_avx2(parsimony,
                                             node1_score_index,
                                             node2_score_index);
  else
#endif
  return fastparsimony_edge_score(parsimony,
                                  node1_score_index,
                                  node2_score_index);

}


PLL_EXPORT unsigned int pll_fastparsimony_root_score(const pll_parsimony_t * parsimony,
                                                     unsigned int root_index)
{
  return parsimony->node_cost[root_index] +
         parsimony->const_cost;
}
