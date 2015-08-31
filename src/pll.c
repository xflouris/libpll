/*
    Copyright (C) 2015 Tomas Flouri, Diego Darriba

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

int pll_errno;
char pll_errmsg[200] = {0};

static void dealloc_partition_data(pll_partition_t * partition);

static void dealloc_partition_data(pll_partition_t * partition)
{
  unsigned int i;

  if (!partition) return;

  free(partition->rates);
  free(partition->eigen_decomp_valid);
  if (partition->prop_invar)
    free(partition->prop_invar);
  if (partition->invariant)
    free(partition->invariant);
  if (!partition->pattern_weights)
    free(partition->pattern_weights);

  if (partition->scale_buffer)
    for (i = 0; i < partition->scale_buffers; ++i)
      free(partition->scale_buffer[i]);
  free(partition->scale_buffer);

  if (partition->clv)
    for (i = 0; i < partition->clv_buffers + partition->tips; ++i)
      pll_aligned_free(partition->clv[i]);
  free(partition->clv);

  if (partition->pmatrix)
    for (i = 0; i < partition->prob_matrices; ++i)
      pll_aligned_free(partition->pmatrix[i]);
  free(partition->pmatrix);

  if (partition->subst_params)
    for (i = 0; i < partition->rate_matrices; ++i)
      pll_aligned_free(partition->subst_params[i]);
  free(partition->subst_params);

  if (partition->eigenvecs)
    for (i = 0; i < partition->rate_matrices; ++i)
      pll_aligned_free(partition->eigenvecs[i]);
  free(partition->eigenvecs);

  if (partition->inv_eigenvecs)
    for (i = 0; i < partition->rate_matrices; ++i)
      pll_aligned_free(partition->inv_eigenvecs[i]);
  free(partition->inv_eigenvecs);

  if (partition->eigenvals)
    for (i = 0; i < partition->rate_matrices; ++i)
      pll_aligned_free(partition->eigenvals[i]);
  free(partition->eigenvals);

  if (partition->frequencies)
    for (i = 0; i < partition->rate_matrices; ++i)
      pll_aligned_free(partition->frequencies[i]);
  free(partition->frequencies);

  if (partition->pattern_weights)
    free(partition->pattern_weights);

  free(partition);
}

PLL_EXPORT void * pll_aligned_alloc(size_t size, size_t alignment)
{
  void * mem;

#if (defined(__WIN32__) || defined(__WIN64__))
  mem = _aligned_malloc(size, alignment);
#else
  if (posix_memalign(&mem, alignment, size)) 
    mem = NULL;
#endif
  
  return mem;
}

PLL_EXPORT void pll_aligned_free(void * ptr)
{
#if (defined(__WIN32__) || defined(__WIN64__))
  _aligned_free(ptr);
#else
  free(ptr);
#endif
}

PLL_EXPORT pll_partition_t * pll_partition_create(unsigned int tips,
                                                  unsigned int clv_buffers,
                                                  unsigned int states,
                                                  unsigned int sites,
                                                  unsigned int rate_matrices,
                                                  unsigned int prob_matrices,
                                                  unsigned int rate_cats,
                                                  unsigned int scale_buffers,
                                                  unsigned int attributes)
{
  unsigned int i;

  /* make sure that multiple ARCH were not specified */
  if (__builtin_popcount(attributes & PLL_ATTRIB_ARCH_MASK) > 1)
  {
    pll_errno = PLL_ERROR_MULTIPLE_ARCH;
    return PLL_FAILURE;
  }

  /* allocate partition */
  pll_partition_t * partition = (pll_partition_t *)malloc(sizeof(pll_partition_t));
  if (!partition) return PLL_FAILURE;

  /* extract architecture and set vectorization parameters */
  partition->alignment = PLL_ALIGNMENT_CPU;
  partition->attributes = attributes;
  partition->states_padded = states;
#ifdef HAVE_SSE
  if (attributes | PLL_ATTRIB_ARCH_SSE)
  {
    partition->alignment = PLL_ALIGNMENT_SSE;
    partition->states_padded = (states+1) & 0xFFFFFFFE;
  }
#endif
#ifdef HAVE_AVX
  if (attributes | PLL_ATTRIB_ARCH_AVX)
  {
    partition->alignment = PLL_ALIGNMENT_AVX;
    partition->states_padded = (states+3) & 0xFFFFFFFC;
  }
#endif
  unsigned int states_padded = partition->states_padded;

  /* initialize properties */

  partition->tips = tips;
  partition->clv_buffers = clv_buffers;
  partition->states = states;
  partition->sites = sites;
  partition->rate_matrices = rate_matrices;
  partition->prob_matrices = prob_matrices;
  partition->rate_cats = rate_cats;
  partition->scale_buffers = scale_buffers;

  partition->prop_invar = NULL;
  partition->invariant = NULL;
  partition->pattern_weights = NULL;

  partition->eigenvecs = NULL;
  partition->inv_eigenvecs = NULL;
  partition->eigenvals = NULL;

  partition->rates = NULL;
  partition->subst_params = NULL;
  partition->scale_buffer = NULL;
  partition->frequencies = NULL;
  partition->eigen_decomp_valid = 0;

  /* allocate structures */

  /* eigen_decomp_valid */
  partition->eigen_decomp_valid = (int *)calloc(partition->rate_matrices,
                                                sizeof(int));
  if (!partition->eigen_decomp_valid)
  {
    dealloc_partition_data(partition);
    return PLL_FAILURE;
  }
  /* clv */
  partition->clv = (double **)calloc(partition->tips + partition->clv_buffers, 
                                     sizeof(double *));
  if (!partition->clv)
  {
    dealloc_partition_data(partition);
    return PLL_FAILURE;
  }
  for (i = 0; i < partition->tips + partition->clv_buffers; ++i)
  {
    partition->clv[i] = pll_aligned_alloc(partition->sites * states_padded *
                                             rate_cats * sizeof(double),
                                          partition->alignment);
    if (!partition->clv[i])
    {
      dealloc_partition_data(partition);
      return PLL_FAILURE;
    }
  }

  /* pmatrix */
  partition->pmatrix = (double **)calloc(partition->prob_matrices,
                                         sizeof(double *));
  if (!partition->pmatrix)
  {
    dealloc_partition_data(partition);
    return PLL_FAILURE;
  }
  for (i = 0; i < partition->prob_matrices; ++i)
  {
    partition->pmatrix[i] = pll_aligned_alloc(states_padded * states_padded *
                                              rate_cats * sizeof(double),
                                              partition->alignment);
    if (!partition->pmatrix[i])
    {
      dealloc_partition_data(partition);
      return PLL_FAILURE;
    }
    /* TODO: don't forget to add code for SSE/AVX */
    memset(partition->pmatrix[i],
           0,
           states_padded*states_padded*rate_cats*sizeof(double));
  }

  /* eigenvecs */
  partition->eigenvecs = (double **)calloc(partition->rate_matrices,
                                           sizeof(double *));
  if (!partition->eigenvecs)
  {
    dealloc_partition_data(partition);
    return PLL_FAILURE;
  }
  for (i = 0; i < partition->rate_matrices; ++i)
  {
    partition->eigenvecs[i] = pll_aligned_alloc(states_padded * states_padded *
                                                  sizeof(double),
                                                partition->alignment);
    if (!partition->eigenvecs[i])
    {
      dealloc_partition_data(partition);
      return PLL_FAILURE;
    }
    /* TODO: don't forget to add code for SSE/AVX */
  }

  /* inv_eigenvecs */
  partition->inv_eigenvecs = (double **)calloc(partition->rate_matrices,
                                               sizeof(double *));
  if (!partition->inv_eigenvecs)
  {
    dealloc_partition_data(partition);
    return PLL_FAILURE;
  }
  for (i = 0; i < partition->rate_matrices; ++i)
  {
    partition->inv_eigenvecs[i] = pll_aligned_alloc(states_padded*states_padded*
                                                       sizeof(double),
                                                    partition->alignment);
    if (!partition->inv_eigenvecs[i])
    {
      dealloc_partition_data(partition);
      return PLL_FAILURE;
    }
    /* TODO: don't forget to add code for SSE/AVX */
  }

  /* eigenvals */
  partition->eigenvals = (double **)calloc(partition->rate_matrices,
                                               sizeof(double *));
  if (!partition->eigenvals)
  {
    dealloc_partition_data(partition);
    return PLL_FAILURE;
  }
  for (i = 0; i < partition->rate_matrices; ++i)
  {
    partition->eigenvals[i] = pll_aligned_alloc(states_padded*sizeof(double),
                                                partition->alignment);
    if (!partition->eigenvals[i])
    {
      dealloc_partition_data(partition);
      return PLL_FAILURE;
    }
    /* TODO: don't forget to add code for SSE/AVX */
  }

  /* subst_params */
  partition->subst_params = (double **)calloc(partition->rate_matrices,
                                              sizeof(double *));
  if (!partition->subst_params)
  {
    dealloc_partition_data(partition);
    return PLL_FAILURE;
  }
  for (i = 0; i < partition->rate_matrices; ++i)
  {
    partition->subst_params[i] = pll_aligned_alloc(((states*states-states)/2) *
                                                   sizeof(double),
                                                   partition->alignment);
    if (!partition->subst_params[i])
    {
      dealloc_partition_data(partition);
      return PLL_FAILURE;
    }
    /* TODO: don't forget to add code for SSE/AVX */
  }

  /* frequencies */
  partition->frequencies = (double **)calloc(partition->rate_matrices, 
                                             sizeof(double *));
  if (!partition->frequencies)
  {
    dealloc_partition_data(partition);
    return PLL_FAILURE;
  }
  for (i = 0; i < partition->rate_matrices; ++i)
  {
    partition->frequencies[i] = pll_aligned_alloc(states_padded*sizeof(double),
                                                  partition->alignment);
    if (!partition->frequencies[i])
    {
      dealloc_partition_data(partition);
      return PLL_FAILURE;
    }
    /* TODO: don't forget to add code for SSE/AVX */
    memset(partition->frequencies[i],
           0,
           states_padded*sizeof(double));
  }

  /* rates */
  partition->rates = (double *)calloc(partition->rate_cats,sizeof(double));
  if (!partition->rates)
  {
    dealloc_partition_data(partition);
    return PLL_FAILURE;
  }

  /* proportion of invariant sites */
  partition->prop_invar = (double *)calloc(partition->rate_matrices,
                                             sizeof(double));
  if (!partition->prop_invar)
  {
    dealloc_partition_data(partition);
    return PLL_FAILURE;
  }

  /* site weights */
  partition->pattern_weights = (unsigned int *)malloc(partition->sites *
                                                      sizeof(unsigned int));
  if (!partition->pattern_weights)
  {
    dealloc_partition_data(partition);
    return PLL_FAILURE;
  }
  for (i = 0; i < partition->sites; ++i) partition->pattern_weights[i] = 1;
  
  /* scale_buffer */
  partition->scale_buffer = (unsigned int **)calloc(partition->scale_buffers,
                                                    sizeof(unsigned int *));
  if (!partition->scale_buffer)
  {
    dealloc_partition_data(partition);
    return PLL_FAILURE;
  }
  for (i = 0; i < partition->scale_buffers; ++i)
  {
    partition->scale_buffer[i] = (unsigned int *)calloc(partition->sites,
                                                        sizeof(unsigned int));
    if (!partition->scale_buffer[i])
    {
      dealloc_partition_data(partition);
      return PLL_FAILURE;
    }
  }

  return partition;
}

PLL_EXPORT void pll_partition_destroy(pll_partition_t * partition)
{
  dealloc_partition_data(partition);
}

PLL_EXPORT int pll_set_tip_states(pll_partition_t * partition, 
                                  unsigned int tip_index,
                                  const unsigned int * map,
                                  const char * sequence)
{
  unsigned int c;
  unsigned int i,j;
  double * tipclv = partition->clv[tip_index];

  /* iterate through sites */
  for (i = 0; i < partition->sites; ++i)
  {
    if ((c = map[(int)sequence[i]]) == 0)
    {
      pll_errno = PLL_ERROR_TIP_DATA_ILLEGAL_STATE;
      snprintf(pll_errmsg, 200, "Illegal state code in tip \"%c\"", sequence[i]);
      return PLL_FAILURE;
    }

    /* decompose basecall into the encoded residues and set the appropriate
       positions in the tip vector */
    for (j = 0; j < partition->states; ++j)
    {
      tipclv[j] = c & 1;
      c >>= 1;
    }

    /* fill in the entries for the other gamma values */
    tipclv += partition->states;
    for (j = 0; j < partition->rate_cats - 1; ++j)
    {
      memcpy(tipclv, tipclv - partition->states,
             partition->states * sizeof(double));
      tipclv += partition->states;
    }
  }
  return PLL_SUCCESS;
}

PLL_EXPORT void pll_set_tip_clv(pll_partition_t * partition,
                                unsigned int tip_index,
                                const double * clv)
{
  unsigned int i,j;
  double * tipclv = partition->clv[tip_index];

  for (i = 0; i < partition->sites; ++i)
  {
    for (j = 0; j < partition->rate_cats; ++j)
    {
      memcpy(tipclv, clv, partition->states*sizeof(double));
      tipclv += partition->states;
    }
    clv += partition->states;
  }
}

PLL_EXPORT void pll_set_pattern_weights(pll_partition_t * partition,
                                        const unsigned int * pattern_weights)
{
  memcpy(partition->pattern_weights, 
         pattern_weights,
         sizeof(unsigned int)*partition->sites);
}
