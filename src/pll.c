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

static pll_dlist_t * partition_list = NULL;
static void dealloc_partition_data(pll_partition_t * partition);



static void dealloc_partition_data(pll_partition_t * partition)
{
  int i;

  if (!partition) return;

  free(partition->rates);
  free(partition->scale_buffer);
  free(partition->eigen_decomp_valid);
  free(partition->invariant);

/*
  if (partition->tip_clv)
    for (i = 0; i < partition->tips; ++i) 
      free(partition->tip_clv[i]);
  free(partition->tip_clv);

  if (partition->inner_clv)
    for (i = 0; i < partition->clv_buffers; ++i) 
      free(partition->inner_clv[i]);
  free(partition->inner_clv);
*/
  if (partition->clv)
    for (i = 0; i < partition->clv_buffers + partition->tips; ++i)
      free(partition->clv[i]);
  free(partition->clv);

  if (partition->pmatrix)
    for (i = 0; i < partition->prob_matrices; ++i)
      free(partition->pmatrix[i]);
  free(partition->pmatrix);

  if (partition->subst_params)
    for (i = 0; i < partition->rate_matrices; ++i)
      free(partition->subst_params[i]);
  free(partition->subst_params);

  if (partition->eigenvecs)
    for (i = 0; i < partition->rate_matrices; ++i)
      free(partition->eigenvecs[i]);
  free(partition->eigenvecs);

  if (partition->inv_eigenvecs)
    for (i = 0; i < partition->rate_matrices; ++i)
      free(partition->inv_eigenvecs[i]);
  free(partition->inv_eigenvecs);

  if (partition->eigenvals)
    for (i = 0; i < partition->rate_matrices; ++i)
      free(partition->eigenvals[i]);
  free(partition->eigenvals);

  if (partition->frequencies)
    for (i = 0; i < partition->rate_matrices; ++i)
      free(partition->frequencies[i]);
  free(partition->frequencies);

  free(partition);
}

PLL_EXPORT pll_partition_t * pll_create_partition(int tips,
                                                  int clv_buffers,
                                                  int states,
                                                  int sites,
                                                  int rate_matrices,
                                                  int prob_matrices,
                                                  int rate_cats,
                                                  int scale_buffers,
                                                  int attributes)
{
  int i;

  pll_partition_t * partition = (pll_partition_t *)malloc(sizeof(pll_partition_t));
  if (!partition) return PLL_FAILURE;

  int alignment = PLL_ALIGNMENT_CPU;
#ifdef HAVE_SSE
  if (attributes | PLL_ATTRIBUTE_SSE) alignment = PLL_ALIGNMENT_SSE;
#endif
#ifdef HAVE_AVX
  if (attributes | PLL_ATTRIBUTE_AVX) alignment = PLL_ALIGNMENT_AVX;
#endif

  /* initialize properties */

  partition->tips = tips;
  partition->clv_buffers = clv_buffers;
  partition->states = states;
  partition->sites = sites;
  partition->rate_matrices = rate_matrices;
  partition->prob_matrices = prob_matrices;
  partition->rate_cats = rate_cats;
  partition->scale_buffers = scale_buffers;

  partition->prop_invar = 0.0;
  partition->invariant = NULL;

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
    if (posix_memalign((void **)&(partition->clv[i]),
                       alignment, 
                       partition->sites * states * rate_cats * sizeof(double)))
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
    /* TODO: don't forget to add code for SSE/AVX */
    if (posix_memalign((void **)&(partition->pmatrix[i]), alignment,
                                   states*states*rate_cats*sizeof(double)))
    {
      dealloc_partition_data(partition);
      return PLL_FAILURE;
    }
  }

  /* eigenvecs */
  partition->eigenvecs = (double **)calloc(partition->rate_matrices,
                                           sizeof(double *));
  if (!partition->pmatrix)
  {
    dealloc_partition_data(partition);
    return PLL_FAILURE;
  }
  for (i = 0; i < partition->rate_matrices; ++i)
  {
    /* TODO: don't forget to add code for SSE/AVX */
    if (posix_memalign((void **)&(partition->eigenvecs[i]), alignment, 
                                   states*states*sizeof(double)))
    {
      dealloc_partition_data(partition);
      return PLL_FAILURE;
    }
  }

  /* inv_eigenvecs */
  partition->inv_eigenvecs = (double **)calloc(partition->rate_matrices,
                                               sizeof(double *));
  if (!partition->pmatrix)
  {
    dealloc_partition_data(partition);
    return PLL_FAILURE;
  }
  for (i = 0; i < partition->rate_matrices; ++i)
  {
    /* TODO: don't forget to add code for SSE/AVX */
    if (posix_memalign((void **)&(partition->inv_eigenvecs[i]), alignment, 
                                   states*states*sizeof(double)))
    {
      dealloc_partition_data(partition);
      return PLL_FAILURE;
    }
  }

  /* eigenvals */
  partition->eigenvals = (double **)calloc(partition->rate_matrices,
                                               sizeof(double *));
  if (!partition->pmatrix)
  {
    dealloc_partition_data(partition);
    return PLL_FAILURE;
  }
  for (i = 0; i < partition->rate_matrices; ++i)
  {
    /* TODO: don't forget to add code for SSE/AVX */
    if (posix_memalign((void **)&(partition->eigenvals[i]), alignment, 
                                   states*sizeof(double)))
    {
      dealloc_partition_data(partition);
      return PLL_FAILURE;
    }
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
    /* TODO: don't forget to add code for SSE/AVX */
    if (posix_memalign((void **)&(partition->subst_params[i]), alignment, 
                                ((states*states - states)/2)*sizeof(double)))
    {
      dealloc_partition_data(partition);
      return PLL_FAILURE;
    }
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
    /* TODO: don't forget to add code for SSE/AVX */
    if (posix_memalign((void **)&(partition->frequencies[i]), 
                                  alignment, 
                                  states*sizeof(double)))
    {
      dealloc_partition_data(partition);
      return PLL_FAILURE;
    }
  }

  /* rates */
  partition->rates = (double *)calloc(partition->rate_cats,sizeof(double));
  if (!partition->rates)
  {
    dealloc_partition_data(partition);
    return PLL_FAILURE;
  }

  if (!pll_dlist_append(&partition_list, (void *)partition))
  {
    dealloc_partition_data(partition);
    return PLL_FAILURE;
  }

  return partition;
}

PLL_EXPORT int pll_destroy_partition(pll_partition_t * partition)
{
  int rc;
  rc = pll_dlist_remove(&partition_list, (void *)partition);
  dealloc_partition_data(partition);
  return rc;
}

PLL_EXPORT int pll_set_tip_states(pll_partition_t * partition, 
                                  int tip_index,
                                  const unsigned int * map,
                                  const char * sequence)
{
  unsigned int c;
  int i,j;
  double * tipclv = partition->clv[tip_index];

  /* iterate through sites */
  for (i = 0; i < partition->sites; ++i)
  {
    if ((c = map[(int)sequence[i]]) == 0)
      return PLL_FAILURE;

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
                                int tip_index,
                                const double * clv)
{
  int i,j;
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
