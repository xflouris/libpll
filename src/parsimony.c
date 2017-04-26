/*
    Copyright (C) 2015 Tomas Flouri

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

PLL_EXPORT int pll_set_parsimony_sequence(pll_parsimony_t * pars,
                                          unsigned int tip_index,
                                          const unsigned int * map,
                                          const char * sequence)
{
  unsigned int c;
  unsigned int i,j;

  unsigned int states = pars->states;
  double * tipstate = pars->sbuffer[tip_index];

  double * score_matrix = pars->score_matrix;

  /* set infinity as the highest score in the matrix plus one */
  double inf = score_matrix[0];
  for (i = 1; i < states*states; ++i)
    if (score_matrix[i] > inf)
      inf = score_matrix[i];
  inf++;

  for (i = 0; i < pars->sites; ++i)
  {
    if ((c = map[(int)sequence[i]]) == 0)
    {
      pll_errno = PLL_ERROR_TIPDATA_ILLEGALSTATE;
      snprintf(pll_errmsg, 200, "Illegal state code in tip \"%c\"", sequence[i]);
      printf ("%s\n", pll_errmsg);
      return PLL_FAILURE;
    }

    for (j = 0; j < states; ++j)
    {
      if (c & 1)
        tipstate[j] = 0;
      else
        tipstate[j] = inf;
      c >>= 1;
    }

    tipstate += states;
  }

  return PLL_SUCCESS;
}

PLL_EXPORT void pll_parsimony_destroy(pll_parsimony_t * parsimony)
{
  unsigned int i;
  unsigned int nodes_count = parsimony->tips + parsimony->inner_nodes;

  if (!parsimony) return;

  /* deallocate fast parsimony structures */
  if (parsimony->packedvector)
  {
    for (i=0; i < nodes_count; ++i)
      pll_aligned_free(parsimony->packedvector[i]);
    free(parsimony->packedvector);
  }

  if (parsimony->node_cost)
    free(parsimony->node_cost);

  if (parsimony->informative)
    free(parsimony->informative);

  /* if available, deallocate structures for weighted parsimony */

  /* score buffers */
  if (parsimony->sbuffer)
  {
    for (i = 0; i < parsimony->score_buffers + parsimony->tips; ++i)
      free(parsimony->sbuffer[i]);
    free(parsimony->sbuffer);
  }

  /* ancestral state buffers */
  if (parsimony->anc_states)
  {
    for (i=parsimony->tips; i<parsimony->ancestral_buffers+parsimony->tips; ++i)
      free(parsimony->anc_states[i]);
    free(parsimony->anc_states);
  }

  /* scoring matrix */
  if (parsimony->score_matrix) free(parsimony->score_matrix);

  free(parsimony);
}

PLL_EXPORT pll_parsimony_t * pll_parsimony_create(unsigned int tips,
                                                  unsigned int states,
                                                  unsigned int sites,
                                                  const double * score_matrix,
                                                  unsigned int score_buffers,
                                                  unsigned int ancestral_buffers)
{
  unsigned int i;

  /* create parsimony instance */
  pll_parsimony_t * pars = (pll_parsimony_t *)calloc(1,sizeof(pll_parsimony_t));
  if (!pars)
  {
    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf(pll_errmsg, 200, "Unable to allocate enough memory.");
    return PLL_FAILURE;
  }

  /* store passed parameters */
  pars->tips = tips;
  pars->states = states;
  pars->sites = sites;
  pars->score_buffers = score_buffers;
  pars->ancestral_buffers = ancestral_buffers;

  /* create and copy a states*states scoring matrix */
  pars->score_matrix = (double *)calloc(states*states, sizeof(double));
  if (!pars->score_matrix)
  {
    pll_parsimony_destroy(pars);
    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf(pll_errmsg, 200,
             "Unable to allocate enough memory for scoring matrix.");
    return NULL;
  }
  memcpy(pars->score_matrix, score_matrix, states*states*sizeof(double));

  /* create requested score buffers */
  pars->sbuffer = (double **)calloc(score_buffers+tips, sizeof(double *));
  if (!pars->sbuffer)
  {
    pll_parsimony_destroy(pars);
    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf(pll_errmsg, 200,
             "Unable to allocate enough memory for score buffers.");
    return NULL;
  }
  for (i=0; i < score_buffers+tips; ++i)
  {
    pars->sbuffer[i] = (double *)calloc(sites*states, sizeof(double *));
    if (!pars->sbuffer[i])
    {
      pll_parsimony_destroy(pars);
      pll_errno = PLL_ERROR_MEM_ALLOC;
      snprintf(pll_errmsg, 200,
               "Unable to allocate enough memory for score buffers.");
      return NULL;
    }
  }

  /* create ancestral buffers */
  pars->anc_states = (unsigned int **)calloc(tips+ancestral_buffers,
                                             sizeof(unsigned int *));
  if (!pars->anc_states)
  {
    pll_parsimony_destroy(pars);
    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf(pll_errmsg, 200,
             "Unable to allocate enough memory for score buffers.");
    return NULL;
  }
  for (i=tips; i < ancestral_buffers+tips; ++i)
  {
    pars->anc_states[i] = (unsigned int *)calloc(sites, sizeof(unsigned int));
    if (!pars->anc_states[i])
    {
      pll_parsimony_destroy(pars);
      pll_errno = PLL_ERROR_MEM_ALLOC;
      snprintf(pll_errmsg, 200,
               "Unable to allocate enough memory for score buffers.");
      return NULL;
    }
  }

  return pars;
}
                                                  
PLL_EXPORT double pll_parsimony_build(pll_parsimony_t * pars,
                                      const pll_pars_buildop_t * operations,
                                      unsigned int count)
{
  unsigned int i,j,k,n;
  const pll_pars_buildop_t * op;

  unsigned int sites = pars->sites;
  unsigned int states = pars->states;
  double minimum;

  double * score_buffer;
  double * child1_score_buffer;
  double * child2_score_buffer;

  double * score_matrix = pars->score_matrix;


  /* Implementation of the 'minimum mutation trees' algorithm by David Sankoff.
     For more information see:
     
  Sankoff D: Minimal Mutation Trees of Sequences. SIAM Appl Math 28:35-32, 1975
  doi: 10.1137/0128004

  Sankoff D, Rousseau P: Locating the Vertices of a Steiner Tree in Arbitrary
  Metric Space. Math Program 1975, 9:240-246.
  doi: 10.1007/BF01681346

  */

  /* iterate through inner nodes (post-order traversal) */
  for (i = 0; i < count; ++i)
  {
    op = &(operations[i]);

    /* get parent score buffer */
    score_buffer = pars->sbuffer[op->parent_score_index];

    /* get child1 score buffer if it's not a tip */
    child1_score_buffer = pars->sbuffer[op->child1_score_index];

    /* get child2 score buffer if it's not a tip */
    child2_score_buffer = pars->sbuffer[op->child2_score_index];

    /* iterate through sites */
    for (j = 0; j < sites; ++j)
    {

      for (n = 0; n < states; ++n)
      {
        /* process child 1 */

        minimum = child1_score_buffer[0] + score_matrix[n];
        for (k = 1; k < states; ++k)
          minimum = fmin(child1_score_buffer[k] + score_matrix[k*states+n],
                         minimum);

        score_buffer[n] = minimum;

        /* process child 2 */

        minimum = child2_score_buffer[0] + score_matrix[n];
        for (k = 1; k < states; ++k)
          minimum = fmin(child2_score_buffer[k] + score_matrix[k*states+n],
                         minimum);

        score_buffer[n] += minimum;
      }

      score_buffer += states;
      child1_score_buffer += states;
      child2_score_buffer += states;
    }
  }

  /* get the sum of minimum scores for each site of the last node in the
     postorder traversal */
  op = &(operations[count-1]);

  return pll_parsimony_score(pars,op->parent_score_index);
}

PLL_EXPORT double pll_parsimony_score(pll_parsimony_t * pars,
                                      unsigned int score_buffer_index)
{
  unsigned int i,j,k;
  unsigned int states = pars->states;
  unsigned int sites = pars->sites;
  double sum = 0;
  double minimum;

  double * score_buffer = pars->sbuffer[score_buffer_index];

  for (k = 0, i = 0; i < sites; ++i)
  {
    minimum = score_buffer[k++];
    for (j = 1; j < states; ++j) 
      minimum = fmin(score_buffer[k++],minimum);
    
    sum += minimum;
  }

  return sum;
}

PLL_EXPORT void pll_parsimony_reconstruct(pll_parsimony_t * pars,
                                          const unsigned int * map,
                                          const pll_pars_recop_t * operations,
                                          unsigned int count)
{
  unsigned int i,j,n;
  unsigned int revmap[256];

  double * score_buffer;
  unsigned int * ancestral_buffer;
  double * parent_score_buffer;
  unsigned int * parent_ancestral_buffer;
  unsigned int minindex;

  unsigned int states=pars->states;

  const pll_pars_recop_t * op;

  for (i = 0; i < 256; ++i) revmap[i] = 0;
  for (i = 0; i < 256; ++i)
  {
    if (__builtin_popcount(map[i]) == 1)
    {
      revmap[__builtin_ctz(map[i])] = i;
    }
  }

  /* start from root of given subtree */
  op = &(operations[0]);
  score_buffer = pars->sbuffer[op->node_score_index];
  ancestral_buffer = pars->anc_states[op->node_ancestral_index];
  for (n = 0; n < pars->sites; ++n)
  {
    minindex= 0;
    for (i = 1; i < pars->states; ++i)
    {
      if (score_buffer[n*states+i] < score_buffer[n*states+minindex])
        minindex= i;
    }
    ancestral_buffer[n] = revmap[minindex];
  }

  /* continue with preorder traversal */
  for (i = 1; i < count; ++i)
  {
    op = &(operations[i]);

    /* get node score and ancestral buffer */
    parent_score_buffer = pars->sbuffer[op->parent_score_index];
    parent_ancestral_buffer = pars->anc_states[op->parent_ancestral_index];

    /* get node score and ancestral buffer */
    score_buffer = pars->sbuffer[op->node_score_index];
    ancestral_buffer = pars->anc_states[op->node_ancestral_index];

    for (n = 0; n < pars->sites; ++n)
    {
      /* compute index with minimum value */
      minindex = 0;
      for (j = 1; j < pars->states; ++j)
      {
        if (score_buffer[n*states+j] < score_buffer[n*states+minindex])
          minindex = j;
      }

      double parent_val = parent_score_buffer[n*states + 
                            __builtin_ctz(map[parent_ancestral_buffer[n]])];

      if (score_buffer[n*states+minindex] + 1 > parent_val)
        ancestral_buffer[n] = parent_ancestral_buffer[n];
      else
        ancestral_buffer[n] = revmap[minindex];
    }
  }
}
