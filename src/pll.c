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

__thread int pll_errno;
__thread char pll_errmsg[200] = {0};

static void dealloc_partition_data(pll_partition_t * partition);

static void dealloc_partition_data(pll_partition_t * partition)
{
  unsigned int i;

  if (!partition) return;

  free(partition->rates);
  free(partition->rate_weights);
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

  if (partition->tipchars)
    for (i = 0; i < partition->tips; ++i)
      pll_aligned_free(partition->tipchars[i]);
  free(partition->tipchars);

  if (partition->ttlookup)
    pll_aligned_free(partition->ttlookup);

  if (partition->charmap)
    free(partition->charmap);

  if (partition->tipmap)
    free(partition->tipmap);

  if (partition->clv)
  {
    int start = (partition->attributes & PLL_ATTRIB_PATTERN_TIP) ?
                    partition->tips : 0;
    for (i = start; i < partition->clv_buffers + partition->tips; ++i)
      pll_aligned_free(partition->clv[i]);
  }
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

/* create a bijective mapping from states to the range <1,maxstates> where
   maxstates is the maximum number of states (including ambiguities). It is 
   neeed to index the precomputed conditional likelihoods for each pair of
   states. The sequences are then encoded using this charmap, and we store
   the precomputated CLV for a charmapped pair i and j, at index:
   
   (i << ceil(log(maxstate)) + j) << log(states) << log(rates) */
static int create_charmap(pll_partition_t * partition)
{
  unsigned int i,j,m = 0;
  char k = 0;
  unsigned int map[PLL_ASCII_SIZE];
  
  memcpy(map, partition->map, PLL_ASCII_SIZE * sizeof(unsigned int));

  if (!(partition->charmap = (char *)calloc(PLL_ASCII_SIZE,sizeof(char))))
  {
    pll_errno = PLL_ERROR_INIT_CHARMAP;
    snprintf (pll_errmsg, 200,
              "Cannot allocate charmap for tip-tip precomputation.");
    return PLL_FAILURE;
  }

  if (!(partition->tipmap= (unsigned int *)calloc(PLL_ASCII_SIZE,
                                                   sizeof(unsigned int))))
  {
    pll_errno = PLL_ERROR_INIT_CHARMAP;
    snprintf (pll_errmsg, 200,
              "Cannot allocate tipmap for tip-tip precomputation.");
    return PLL_FAILURE;
  }


  for (i = 0; i < PLL_ASCII_SIZE; ++i)
  {
    if (map[i])
    {
      if (map[i] > m) m = map[i];

      partition->charmap[i] = k;
      partition->tipmap[(int)k] = map[i];
      for (j = i+1; j < PLL_ASCII_SIZE; ++j)
      {
        if (map[i] == map[j])
        {
          partition->charmap[j] = k;
          map[j] = 0;
        }
      }
      ++k;
    }
  }

  if (partition->states == 4) k = m+1;

  /* set maximum number of states (including ambiguities), its logarithm,
     and the logarithm of states */
  partition->maxstates = (unsigned int)k;
  partition->log2_maxstates = (unsigned int)ceil(log2(partition->maxstates));
  partition->log2_states = (unsigned int)ceil(log2(partition->states));
  partition->log2_rates = (unsigned int)ceil(log2(partition->rate_cats));

  /* allocate space for the precomputed tip-tip likelihood vector */
  unsigned int addressbits = partition->log2_maxstates + 
                             partition->log2_maxstates +
                             partition->log2_states +
                             partition->log2_rates;

  /* dedicated 4x4 function */
  if (partition->states == 4 && (partition->attributes & PLL_ATTRIB_ARCH_AVX))
  {
    partition->ttlookup = pll_aligned_alloc(1024 * partition->rate_cats *
                                            sizeof(double),
                                            partition->alignment);
    if (!partition->ttlookup)
    {
      pll_errno = PLL_ERROR_INIT_CHARMAP;
      snprintf (pll_errmsg, 200,
                "Cannot allocate space for storing precomputed tip-tip CLVs.");
      return PLL_FAILURE;
    }
  }
  else
  {
    partition->ttlookup = pll_aligned_alloc((1 << addressbits) * 
                                            sizeof(double),
                                            partition->alignment);
    if (!partition->ttlookup)
    {
      pll_errno = PLL_ERROR_INIT_CHARMAP;
      snprintf (pll_errmsg, 200,
                "Cannot allocate space for storing precomputed tip-tip CLVs.");
      return PLL_FAILURE;
    }
  }

  /* allocate tip character arrays */
  partition->tipchars= (char **)calloc(partition->tips, sizeof(char *));
  if (!partition->tipchars)
  {
    pll_errno = PLL_ERROR_INIT_CHARMAP;
    snprintf (pll_errmsg, 200,
              "Cannot allocate space for storing tip characters.");
    return PLL_FAILURE;
  }

  for (i = 0; i < partition->tips; ++i)
  {
    partition->tipchars[i] = (char *)malloc(partition->sites * sizeof(char));
    if (!partition->tipchars[i])
    {
      pll_errno = PLL_ERROR_INIT_CHARMAP;
      snprintf (pll_errmsg, 200,
                "Cannot allocate space for storing tip characters.");
      return PLL_FAILURE;
    }
  }

  return PLL_SUCCESS;
}

PLL_EXPORT pll_partition_t * pll_partition_create(unsigned int tips,
                                                  unsigned int clv_buffers,
                                                  unsigned int states,
                                                  unsigned int sites,
                                                  unsigned int rate_matrices,
                                                  unsigned int prob_matrices,
                                                  unsigned int rate_cats,
                                                  unsigned int scale_buffers,
                                                  const unsigned int * map,
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
  if (attributes & PLL_ATTRIB_ARCH_SSE)
  {
    partition->alignment = PLL_ALIGNMENT_SSE;
    partition->states_padded = (states+1) & 0xFFFFFFFE;
  }
#endif
#ifdef HAVE_AVX
  if (attributes & PLL_ATTRIB_ARCH_AVX)
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
  partition->rate_weights = NULL;
  partition->subst_params = NULL;
  partition->scale_buffer = NULL;
  partition->frequencies = NULL;
  partition->eigen_decomp_valid = 0;

  partition->ttlookup = NULL;
  partition->tipchars = NULL;
  partition->charmap = NULL;
  partition->tipmap = NULL;

  partition->map = map;

  /* create character map for tip-tip precomputations */
  if (partition->attributes & PLL_ATTRIB_PATTERN_TIP)
  {
    if (!create_charmap(partition))
    {
      dealloc_partition_data(partition);
      return PLL_FAILURE;
    }
  }

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

  /* if tip pattern precomputation is enabled, then do not allocate CLV space
     for the tip nodes */
  int start = (partition->attributes & PLL_ATTRIB_PATTERN_TIP) ?
                  partition->tips : 0;

  for (i = start; i < partition->tips + partition->clv_buffers; ++i)
  {
    partition->clv[i] = pll_aligned_alloc(partition->sites * states_padded *
                                             rate_cats * sizeof(double),
                                          partition->alignment);
    if (!partition->clv[i])
    {
      dealloc_partition_data(partition);
      return PLL_FAILURE;
    }
    memset(partition->clv[i],
           0,
           partition->sites*states_padded*rate_cats*sizeof(double));
  }

  /* pmatrix */
  partition->pmatrix = (double **)calloc(partition->prob_matrices,
                                         sizeof(double *));
  if (!partition->pmatrix)
  {
    dealloc_partition_data(partition);
    return PLL_FAILURE;
  }

  size_t displacement = (states_padded - states) * (states_padded) * sizeof(double);
  for (i = 0; i < partition->prob_matrices; ++i)
  {
    partition->pmatrix[i] = pll_aligned_alloc(states * states_padded *
                                              rate_cats * sizeof(double) +
                                              displacement,
                                              partition->alignment);
    if (!partition->pmatrix[i])
    {
      dealloc_partition_data(partition);
      return PLL_FAILURE;
    }
    /* TODO: don't forget to add code for SSE/AVX */
    memset(partition->pmatrix[i],
           0,
           states*states_padded*rate_cats*sizeof(double));
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
    partition->eigenvecs[i] = pll_aligned_alloc(states * states_padded *
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
    partition->inv_eigenvecs[i] = pll_aligned_alloc(states*states_padded*
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

  /* rate weights */
  partition->rate_weights = (double *)calloc(partition->rate_cats,sizeof(double));
  if (partition->rate_weights)
  {
    /* initialize to 1/n_rates */
    for (i = 0; i < partition->rate_cats; ++i)
      partition->rate_weights[i] = 1.0/partition->rate_cats;
  }
  else
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

static int set_tipchars_4x4(pll_partition_t * partition,
                            unsigned int tip_index,
                            const unsigned int * map,
                            const char * sequence)
{
  unsigned int c;
  unsigned int i;

  /* iterate through sites */
  for (i = 0; i < partition->sites; ++i)
  {
    if ((c = map[(int)sequence[i]]) == 0)
    {
      pll_errno = PLL_ERROR_TIP_DATA_ILLEGAL_STATE;
      snprintf(pll_errmsg, 200, "Illegal state code in tip \"%c\"", sequence[i]);
      return PLL_FAILURE;
    }

    /* store states as the remapped characters from charmap */
    partition->tipchars[tip_index][i] = (char)c;
  }

  /* tipmap is never used in the 4x4 case */
  memcpy(partition->tipmap, map, PLL_ASCII_SIZE*sizeof(unsigned int));

  return PLL_SUCCESS;
}

static int set_tipchars(pll_partition_t * partition,
                        unsigned int tip_index,
                        const unsigned int * map,
                        const char * sequence)
{
  unsigned int c;
  unsigned int i;

  /* iterate through sites */
  for (i = 0; i < partition->sites; ++i)
  {
    if ((c = map[(int)sequence[i]]) == 0)
    {
      pll_errno = PLL_ERROR_TIP_DATA_ILLEGAL_STATE;
      snprintf(pll_errmsg, 200, "Illegal state code in tip \"%c\"", sequence[i]);
      return PLL_FAILURE;
    }

    /* store states as the remapped characters from charmap */
    partition->tipchars[tip_index][i] = partition->charmap[(int)(sequence[i])];
  }
  return PLL_SUCCESS;
}

static int set_tipclv(pll_partition_t * partition,
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
    tipclv += partition->states_padded;
    for (j = 0; j < partition->rate_cats - 1; ++j)
    {
      memcpy(tipclv, tipclv - partition->states_padded,
             partition->states * sizeof(double));
      tipclv += partition->states_padded;
    }
  }
  return PLL_SUCCESS;
}

PLL_EXPORT int pll_set_tip_states(pll_partition_t * partition, 
                                  unsigned int tip_index,
                                  const unsigned int * map,
                                  const char * sequence)
{
  int rc;

  if (partition->attributes & PLL_ATTRIB_PATTERN_TIP)
  {
    if (partition->states == 4)
      rc = set_tipchars_4x4(partition, tip_index, map, sequence);
    else
      rc = set_tipchars(partition, tip_index, map, sequence);
  }
  else
    rc = set_tipclv(partition, tip_index, map, sequence);

  return rc;
}

//TODO: <DOC> We should account for padding before calling this function
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
      tipclv += partition->states_padded;
    }
    clv += partition->states_padded;
  }
}

PLL_EXPORT void pll_set_pattern_weights(pll_partition_t * partition,
                                        const unsigned int * pattern_weights)
{
  memcpy(partition->pattern_weights, 
         pattern_weights,
         sizeof(unsigned int)*partition->sites);
}
