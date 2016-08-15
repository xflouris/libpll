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

static int sumtable_tipinner(pll_partition_t * partition,
                                unsigned int parent_clv_index,
                                unsigned int child_clv_index,
                                const unsigned int * params_indices,
                                double *sumtable)
{
  unsigned int i, retval;
  unsigned int tip_clv_index;
  unsigned int inner_clv_index;
  unsigned int sites = partition->sites;
  double ** eigenvecs = (double **)malloc(partition->rate_cats *
                                          sizeof(double *));
  double ** inv_eigenvecs = (double **)malloc(partition->rate_cats *
                                              sizeof(double *));
  double ** freqs = (double **)malloc(partition->rate_cats *
                                      sizeof(double *));
  if (!eigenvecs || !inv_eigenvecs || !freqs)
  {
    if (eigenvecs) free(eigenvecs);
    if (inv_eigenvecs) free(inv_eigenvecs);
    if (freqs) free(freqs);

    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf(pll_errmsg, 200, "Unable to allocate enough memory.");
    return PLL_FAILURE;
  }

  /* ascertaiment bias correction */
  if (partition->asc_bias_alloc)
    sites += partition->states;

  for (i = 0; i < partition->rate_cats; ++i)
  {
    eigenvecs[i] = partition->eigenvecs[params_indices[i]];
    inv_eigenvecs[i] = partition->inv_eigenvecs[params_indices[i]];
    freqs[i] = partition->frequencies[params_indices[i]];
  }

  /* find which of the two child nodes is the tip */
  if (parent_clv_index < partition->tips)
  {
    tip_clv_index = parent_clv_index;
    inner_clv_index = child_clv_index;
  }
  else
  {
    tip_clv_index = child_clv_index;
    inner_clv_index = parent_clv_index;
  }

  retval = pll_core_update_sumtable_ti(partition->states,
                                       sites,
                                       partition->rate_cats,
                                       partition->clv[inner_clv_index],
                                       partition->tipchars[tip_clv_index],
                                       eigenvecs,
                                       inv_eigenvecs,
                                       freqs,
                                       partition->tipmap,
                                       sumtable,
                                       partition->attributes);

  free(freqs);
  free(eigenvecs);
  free(inv_eigenvecs);

  return retval;
}

static int sumtable_innerinner(pll_partition_t * partition,
                                unsigned int parent_clv_index,
                                unsigned int child_clv_index,
                                const unsigned int * params_indices,
                                double *sumtable)
{
  unsigned int i, retval;
  unsigned int sites = partition->sites;

  double ** eigenvecs = (double **)malloc(partition->rate_cats *
                                          sizeof(double *));
  double ** inv_eigenvecs = (double **)malloc(partition->rate_cats *
                                              sizeof(double *));
  double ** freqs = (double **)malloc(partition->rate_cats *
                                      sizeof(double *));
  if (!eigenvecs || !inv_eigenvecs || !freqs)
  {
    if (eigenvecs) free(eigenvecs);
    if (inv_eigenvecs) free(inv_eigenvecs);
    if (freqs) free(freqs);

    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf(pll_errmsg, 200, "Unable to allocate enough memory.");
    return PLL_FAILURE;
  }

  /* ascertaiment bias correction */
  if (partition->asc_bias_alloc)
    sites += partition->states;

  for (i = 0; i < partition->rate_cats; ++i)
  {
    eigenvecs[i] = partition->eigenvecs[params_indices[i]];
    inv_eigenvecs[i] = partition->inv_eigenvecs[params_indices[i]];
    freqs[i] = partition->frequencies[params_indices[i]];
  }

  retval =
  pll_core_update_sumtable_ii(partition->states,
                              sites,
                              partition->rate_cats,
                              partition->clv[parent_clv_index],
                              partition->clv[child_clv_index],
                              eigenvecs,
                              inv_eigenvecs,
                              freqs,
                              sumtable,
                              partition->attributes);

  free(freqs);
  free(eigenvecs);
  free(inv_eigenvecs);

  return retval;
}

/* computes the table containing the constant parts of the likelihood function
 * partial derivatives on the branch lengths.
 * sumtable: [output] must be allocated for storing (rates x states_padded) values */
PLL_EXPORT int pll_update_sumtable(pll_partition_t * partition,
                                      unsigned int parent_clv_index,
                                      unsigned int child_clv_index,
                                      const unsigned int * params_indices,
                                      double *sumtable)
{
  int retval;

  if (partition->attributes & PLL_ATTRIB_PATTERN_TIP)
      {
        if ((parent_clv_index < partition->tips) &&
            (child_clv_index < partition->tips))
        {
          /* tip-tip case */
          assert(0);
        }
        else if ((parent_clv_index < partition->tips) ||
                 (child_clv_index < partition->tips))
        {
          /* tip-inner */
          retval = sumtable_tipinner(partition,
                                     parent_clv_index,
                                     child_clv_index,
                                     params_indices,
                                     sumtable);
        }
        else
        {
          /* inner-inner */
          retval = sumtable_innerinner(partition,
                                       parent_clv_index,
                                       child_clv_index,
                                       params_indices,
                                       sumtable);
        }
      }
      else
      {
        /* inner-inner */
        retval = sumtable_innerinner(partition,
                                     parent_clv_index,
                                     child_clv_index,
                                     params_indices,
                                     sumtable);
      }

  return retval;
}

/* Computes partial derivatives on the branch lengths.
 * branch_length: [input] value where the derivative is computed
 * sumtable: [input] must be computed at the edge where the derivatives will
 *                   be computed
 * d_f:  [output] first derivative
 * dd_f: [output] second derivative
 */
PLL_EXPORT int pll_compute_likelihood_derivatives(pll_partition_t * partition,
                                                  int parent_scaler_index,
                                                  int child_scaler_index,
                                                  double branch_length,
                                                  const unsigned int * params_indices,
                                                  const double * sumtable,
                                                  double * d_f,
                                                  double * dd_f)
{
  unsigned int * parent_scaler;
  unsigned int * child_scaler;
  unsigned int i;
  unsigned int rate_cats = partition->rate_cats;

  double ** eigenvals = (double **) malloc(rate_cats * sizeof(double *));
  double ** freqs     = (double **) malloc(rate_cats * sizeof(double *));
  double * prop_invar = (double *)  malloc(rate_cats * sizeof(double));
  if (!eigenvals || !prop_invar || !freqs)
  {
    if (eigenvals) free(eigenvals);
    if (prop_invar) free(prop_invar);
    if (freqs) free(freqs);

    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf(pll_errmsg, 200, "Unable to allocate enough memory.");
    return PLL_FAILURE;
  }

  for (i=0; i<rate_cats; ++i)
  {
    eigenvals[i]  = partition->eigenvals[params_indices[i]];
    freqs[i]      = partition->frequencies[params_indices[i]];
    prop_invar[i] = partition->prop_invar[params_indices[i]];
  }

  /* get parent scaler */
  if (parent_scaler_index == PLL_SCALE_BUFFER_NONE)
    parent_scaler = NULL;
  else
    parent_scaler = partition->scale_buffer[parent_scaler_index];

  if (child_scaler_index == PLL_SCALE_BUFFER_NONE)
    child_scaler = NULL;
  else
    child_scaler = partition->scale_buffer[child_scaler_index];

  int retval = pll_core_likelihood_derivatives(partition->states,
                                               partition->sites,
                                               partition->rate_cats,
                                               partition->rate_weights,
                                               parent_scaler,
                                               child_scaler,
                                               partition->invariant,
                                               partition->pattern_weights,
                                               branch_length,
                                               prop_invar,
                                               freqs,
                                               partition->rates,
                                               eigenvals,
                                               sumtable,
                                               d_f,
                                               dd_f,
                                               partition->attributes);

  free (freqs);
  free (prop_invar);
  free (eigenvals);

  return retval;
}
