/*
    Copyright (C) 2015 Tomas Flouri, Diego Darriba, Alexandros Stamatakis

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

static int mytqli(double *d, double *e, const unsigned int n, double **z)
{
  unsigned int     m, l, iter, i, k;
  double  s, r, p, g, f, dd, c, b;

  for (i = 2; i <= n; i++)
    e[i - 2] = e[i - 1];

  e[n - 1] = 0.0;

  for (l = 1; l <= n; l++)
    {
      iter = 0;
      do
        {
          for (m = l; m <= n - 1; m++)
            {
              dd = fabs(d[m - 1]) + fabs(d[m]);
              if (fabs(e[m - 1]) + dd == dd)
                break;
            }
          if (m != l)
           {
             assert(iter < 30);

             g = (d[l] - d[l - 1]) / (2.0 * e[l - 1]);
             r = sqrt((g * g) + 1.0);
             g = d[m - 1] - d[l - 1] + e[l - 1] / (g + ((g < 0)?-fabs(r):fabs(r)));/*MYSIGN(r, g));*/
             s = c = 1.0;
             p = 0.0;

             for (i = m - 1; i >= l; i--)
               {
                 f = s * e[i - 1];
                 b = c * e[i - 1];
                 if (fabs(f) >= fabs(g))
                   {
                     c = g / f;
                     r = sqrt((c * c) + 1.0);
                     e[i] = f * r;
                     c *= (s = 1.0 / r);
                   }
                 else
                   {
                     s = f / g;
                     r = sqrt((s * s) + 1.0);
                     e[i] = g * r;
                     s *= (c = 1.0 / r);
                   }
                 g = d[i] - p;
                 r = (d[i - 1] - g) * s + 2.0 * c * b;
                 p = s * r;
                 d[i] = g + p;
                 g = c * r - b;
                 for (k = 1; k <= n; k++)
                   {
                     f = z[i][k-1];
                     z[i][k-1] = s * z[i - 1][k - 1] + c * f;
                     z[i - 1][k - 1] = c * z[i - 1][k - 1] - s * f;
                   }
               }

             d[l - 1] = d[l - 1] - p;
             e[l - 1] = g;
             e[m - 1] = 0.0;
           }
        }
      while (m != l);
    }



    return (1);
 }


static void mytred2(double **a, const unsigned int n, double *d, double *e)
{
  unsigned int     l, k, j, i;
  double  scale, hh, h, g, f;

  for (i = n; i > 1; i--)
    {
      l = i - 1;
      h = 0.0;
      scale = 0.0;

      if (l > 1)
        {
          for (k = 1; k <= l; k++)
            scale += fabs(a[k - 1][i - 1]);
          if (scale == 0.0)
            e[i - 1] = a[l - 1][i - 1];
          else
            {
              for (k = 1; k <= l; k++)
                {
                  a[k - 1][i - 1] /= scale;
                  h += a[k - 1][i - 1] * a[k - 1][i - 1];
                }
              f = a[l - 1][i - 1];
              g = ((f > 0) ? -sqrt(h) : sqrt(h)); /* diff */
              e[i - 1] = scale * g;
              h -= f * g;
              a[l - 1][i - 1] = f - g;
              f = 0.0;
              for (j = 1; j <= l; j++)
                {
                  a[i - 1][j - 1] = a[j - 1][i - 1] / h;
                  g = 0.0;
                  for (k = 1; k <= j; k++)
                    g += a[k - 1][j - 1] * a[k - 1][i - 1];
                  for (k = j + 1; k <= l; k++)
                    g += a[j - 1][k - 1] * a[k - 1][i - 1];
                  e[j - 1] = g / h;
                  f += e[j - 1] * a[j - 1][i - 1];
                }
              hh = f / (h + h);
              for (j = 1; j <= l; j++)
                {
                  f = a[j - 1][i - 1];
                  g = e[j - 1] - hh * f;
                  e[j - 1] = g;
                  for (k = 1; k <= j; k++)
                    a[k - 1][j - 1] -= (f * e[k - 1] + g * a[k - 1][i - 1]);
                }
            }
        }
      else
        e[i - 1] = a[l - 1][i - 1];
      d[i - 1] = h;
    }
  d[0] = 0.0;
  e[0] = 0.0;

  for (i = 1; i <= n; i++)
    {
      l = i - 1;
      if (d[i - 1] != 0.0)
        {
          for (j = 1; j <= l; j++)
            {
                g = 0.0;
                for (k = 1; k <= l; k++)
                  g += a[k - 1][i - 1] * a[j - 1][k - 1];
                for(k = 1; k <= l; k++)
                  a[j - 1][k - 1] -= g * a[i - 1][k - 1];
            }
        }
      d[i - 1] = a[i - 1][i - 1];
      a[i - 1][i - 1] = 1.0;
      for (j = 1; j <= l; j++)
        a[i - 1][j - 1] = a[j - 1][i - 1] = 0.0;
    }
}

/* TODO: Add code for SSE/AVX. Perhaps allocate qmatrix in one chunk to avoid the
complex checking when to dealloc */
static double ** create_ratematrix(double * params,
                                   double * frequencies,
                                   unsigned int states)
{
  unsigned int i,j,k,success;

  double ** qmatrix;

  /* normalize substitution parameters */
  unsigned int params_count = (states*states - states) / 2;
  double * params_normalized = (double *)malloc(sizeof(double) * params_count);
  if (!params_normalized)
    return NULL;

  memcpy(params_normalized,params,params_count*sizeof(double));

  if (params_normalized[params_count - 1] > 0.0)
    for (i = 0; i < params_count; ++i)
      params_normalized[i] /= params_normalized[params_count - 1];

  /* allocate qmatrix */
  qmatrix = (double **)malloc(states*sizeof(double *));
  if (!qmatrix)
  {
    free(params_normalized);
    return NULL;
  }

  success = 1;
  for (i = 0; i < states; ++i)
    if (!(qmatrix[i] = (double *)malloc(states*sizeof(double)))) success=0;

  if (!success)
  {
    for(i = 0; i < states; ++i) free(qmatrix[i]);
    free(qmatrix);
    free(params_normalized);
    return NULL;
  }

  /* construct a matrix equal to sqrt(pi) * Q sqrt(pi)^-1 in order to ensure
     it is symmetric */

  for (i = 0; i < states; ++i) qmatrix[i][i] = 0;

  k = 0;
  for (i = 0; i < states; ++i)
  {
    for (j = i+1; j < states; ++j)
    {
      double factor = params_normalized[k++];
      qmatrix[i][j] = qmatrix[j][i] = factor * sqrt(frequencies[i] * frequencies[j]);
      qmatrix[i][i] -= factor * frequencies[j];
      qmatrix[j][j] -= factor * frequencies[i];
    }
  }


  double mean = 0;
  for (i = 0; i < states; ++i) mean += frequencies[i] * (-qmatrix[i][i]);
  for (i = 0; i < states; ++i)
    for (j = 0; j < states; ++j)
      qmatrix[i][j] /= mean;

  free(params_normalized);

  return qmatrix;
}

PLL_EXPORT int pll_update_eigen(pll_partition_t * partition,
                                unsigned int params_index)
{
  unsigned int i,j,k;
  double *e, *d;
  double ** a;

  double * eigenvecs = partition->eigenvecs[params_index];
  double * inv_eigenvecs = partition->inv_eigenvecs[params_index];
  double * eigenvals = partition->eigenvals[params_index];
  double * freqs = partition->frequencies[params_index];
  double * subst_params = partition->subst_params[params_index];

  unsigned int states = partition->states;
  unsigned int states_padded = partition->states_padded;

  a = create_ratematrix(subst_params,
                        freqs,
                        states);
  if (!a)
  {
    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf(pll_errmsg, 200, "Unable to allocate enough memory.");
    return PLL_FAILURE;
  }

  d = (double *)malloc(states*sizeof(double));
  e = (double *)malloc(states*sizeof(double));
  if (!d || !e)
  {
    if (d) free(d);
    if (e) free(e);
    for(i = 0; i < states; ++i) free(a[i]);
    free(a);
    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf(pll_errmsg, 200, "Unable to allocate enough memory.");
    return PLL_FAILURE;
  }

  mytred2(a, states, d, e);
  mytqli(d, e, states, a);

  /* store eigen vectors */
  for (i = 0; i < states; ++i)
    memcpy(eigenvecs + i*states_padded, a[i], states*sizeof(double));

  /* store eigen values */
  memcpy(eigenvals, d, states*sizeof(double));

  /* store inverse eigen vectors */
  for (k = 0, i = 0; i < states; ++i)
  {
    for (j = i; j < states_padded*states; j += states_padded)
      inv_eigenvecs[k++] = eigenvecs[j];

    /* account for padding */
    k += states_padded - states;
  }

  assert(k == states_padded*states);

  /* multiply the inverse eigen vectors from the left with sqrt(pi)^-1 */
  for (i = 0; i < states; ++i)
    for (j = 0; j < states; ++j)
      inv_eigenvecs[i*states_padded+ j] /= sqrt(freqs[i]);

  /* multiply the eigen vectors from the right with sqrt(pi) */
  for (i = 0; i < states; ++i)
    for (j = 0; j < states; ++j)
      eigenvecs[i*states_padded+j] *= sqrt(freqs[j]);

  partition->eigen_decomp_valid[params_index] = 1;

  free(d);
  free(e);
  for (i = 0; i < states; ++i)
    free(a[i]);
  free(a);

  return PLL_SUCCESS;
}

PLL_EXPORT int pll_update_prob_matrices(pll_partition_t * partition,
                                        const unsigned int * params_indices,
                                        const unsigned int * matrix_indices,
                                        const double * branch_lengths,
                                        unsigned int count)
{
  unsigned int n;

  /* check whether we have cached an eigen decomposition. If not, compute it */
  for (n = 0; n < partition->rate_cats; ++n)
  {
    if (!partition->eigen_decomp_valid[params_indices[n]])
    {
      if (!pll_update_eigen(partition, params_indices[n]))
        return PLL_FAILURE;
    }
  }

  return pll_core_update_pmatrix(partition->pmatrix,
                                 partition->states,
                                 partition->rate_cats,
                                 partition->rates,
                                 branch_lengths,
                                 matrix_indices,
                                 params_indices,
                                 partition->prop_invar,
                                 partition->eigenvals,
                                 partition->eigenvecs,
                                 partition->inv_eigenvecs,
                                 count,
                                 partition->attributes);
}

PLL_EXPORT void pll_set_frequencies(pll_partition_t * partition,
                                    unsigned int freqs_index,
                                    const double * frequencies)
{
  memcpy(partition->frequencies[freqs_index],
         frequencies,
         partition->states*sizeof(double));
  partition->eigen_decomp_valid[freqs_index] = 0;
}

PLL_EXPORT void pll_set_category_rates(pll_partition_t * partition,
                                       const double * rates)
{
  memcpy(partition->rates, rates, partition->rate_cats*sizeof(double));
}

PLL_EXPORT void pll_set_category_weights(pll_partition_t * partition,
                                         const double * rate_weights)
{
  memcpy(partition->rate_weights, rate_weights,
         partition->rate_cats*sizeof(double));
}

PLL_EXPORT void pll_set_subst_params(pll_partition_t * partition,
                                     unsigned int params_index,
                                     const double * params)
{
  unsigned int count = partition->states * (partition->states-1) / 2;

  memcpy(partition->subst_params[params_index],
         params, count*sizeof(double));
  partition->eigen_decomp_valid[params_index] = 0;

  /* NOTE: For protein models PLL/RAxML do a rate scaling by 10.0/max_rate */
}

PLL_EXPORT int pll_update_invariant_sites_proportion(pll_partition_t * partition,
                                                     unsigned int params_index,
                                                     double prop_invar)
{

  /* check that there is no ascertainment bias correction */
  if (prop_invar != 0.0 && (partition->attributes & PLL_ATTRIB_AB_MASK))
  {
    pll_errno = PLL_ERROR_INVAR_INCOMPAT;
    snprintf(pll_errmsg,
             200,
             "Invariant sites are not compatible with asc bias correction");
    return PLL_FAILURE;
  }

  /* validate new invariant sites proportion */
  if (prop_invar < 0 || prop_invar >= 1)
  {
    pll_errno = PLL_ERROR_INVAR_PROPORTION;
    snprintf(pll_errmsg,
             200,
             "Invalid proportion of invariant sites (%f)", prop_invar);
    return PLL_FAILURE;
  }

  if (params_index > partition->rate_matrices)
  {
    pll_errno = PLL_ERROR_INVAR_PARAMINDEX;
    snprintf(pll_errmsg,
             200,
             "Invalid params index (%d)", params_index);
    return PLL_FAILURE;
  }

  if (prop_invar > 0.0 && !partition->invariant)
  {
    if (!pll_update_invariant_sites(partition))
    {
      pll_errno = PLL_ERROR_INVAR_NONEFOUND;
      snprintf(pll_errmsg,
               200,
               "No invariant sites found");
      return PLL_FAILURE;
    }
  }

  partition->prop_invar[params_index] = prop_invar;

  return PLL_SUCCESS;
}

PLL_EXPORT unsigned int pll_count_invariant_sites(pll_partition_t * partition,
                                                  unsigned int * state_inv_count)
{
  unsigned int i,j,k;
  unsigned int invariant_count = 0;
  unsigned int tips = partition->tips;
  unsigned int sites = partition->sites;
  unsigned int states = partition->states;
  unsigned int gap_state = 0;
  unsigned int cur_state;
  int * invariant = partition->invariant;
  double * tipclv;

  /* gap state has always all bits set to one */
  for (i = 0; i < states; ++i)
  {
    gap_state <<= 1;
    gap_state |= 1;
  }

  if (state_inv_count)
    memset(state_inv_count, 0, states*sizeof(unsigned int));

  if (invariant)
  {
    /* count the invariant sites for each state */
    for (i=0; i<sites; ++i)
    {
      if (invariant[i] > -1)
      {
        cur_state = invariant[i];
        /* since the invariant sites array is generated in the library,
           it should not contain invalid values */
        assert (cur_state < states);

        /* increase the counter and per-state count */
        invariant_count += partition->pattern_weights[i];
        if (state_inv_count)
          state_inv_count[cur_state]++;
      }
    }
  }
  else
  {
    if (partition->attributes & PLL_ATTRIB_PATTERN_TIP)
    {
      for (j = 0; j < sites; ++j)
      {
        cur_state = gap_state;
        for (i = 0; i < tips; ++i)
        {
          cur_state &= ((unsigned int)(partition->tipchars[i][j]));
          if  (!cur_state)
          {
            break;
          }
        }
        if (__builtin_popcount(cur_state) == 1)
        {
          invariant_count += partition->pattern_weights[j];
          if (state_inv_count)
            state_inv_count[__builtin_ctz(cur_state)]++;
        }
      }
    }
    else
    {
      /* warning: note that this operation traverses the clvs by columns, and
         hence it may be slow. If PLL_ATTRIB_PATTERN_TIP is not set, I suggest
         to call pll_update_invariant_sites() before calling this function in
         order to populate partition->invariant beforehand. It can be freed
         afterwards. */
      unsigned int span_padded = partition->rate_cats * partition->states_padded;

      for (j = 0; j < sites; ++j)
      {
        unsigned int clv_shift = j*span_padded;
        tipclv = partition->clv[0] + clv_shift;
        unsigned int state = gap_state;
        for (i = 0; i < tips; ++i)
        {
          tipclv = partition->clv[i] + clv_shift;
          cur_state = 0;
          for (k = 0; k < states; ++k)
          {
            cur_state |= ((unsigned int)tipclv[k] << k);
          }
          state &= cur_state;
          if (!state)
          {
            break;
          }
        }
        if (__builtin_popcount(state) == 1)
        {
          invariant_count += partition->pattern_weights[j];
          if (state_inv_count)
            state_inv_count[__builtin_ctz(state)]++;
        }
      }
    }
  }
  return invariant_count;
}

PLL_EXPORT int pll_update_invariant_sites(pll_partition_t * partition)
{
  unsigned int i,j,k;
  unsigned int state;
  unsigned int states = partition->states;
  unsigned int states_padded = partition->states_padded;
  unsigned int sites = partition->sites;
  unsigned int tips = partition->tips;
  unsigned int rate_cats = partition->rate_cats;
  unsigned int gap_state = 0;
  int * invariant;
  double * tipclv;

  /* gap state has always all bits set to one */
  for (i = 0; i < states; ++i)
  {
    gap_state <<= 1;
    gap_state |= 1;
  }

  /* allocate array (on first call) denoting the frequency index for invariant
     sites, or -1 for variant sites */
  if (!partition->invariant)
  {
    if (!(partition->invariant = (int *)malloc(sites * sizeof(int))))
    {
      pll_errno = PLL_ERROR_MEM_ALLOC;
      snprintf(pll_errmsg, 200,
          "Cannot allocate charmap for invariant sites array.");
      return PLL_FAILURE;
    }
  }
  invariant = partition->invariant;

  /* initialize all elements to zero */
  memset(partition->invariant, gap_state, sites*sizeof(int));

  /* depending on the attribute flag, fill each element of the invariant array
     with the bitwise AND of gap and all states in the corresponding site */
  if (partition->attributes & PLL_ATTRIB_PATTERN_TIP)
  {
    if (states == 4)
    {
      for (i = 0; i < tips; ++i)
        for (j = 0; j < sites; ++j)
        {
          state = (unsigned int)(partition->tipchars[i][j]);
            invariant[j] &= state;
        }
    }
    else
    {
      for (i = 0; i < tips; ++i)
        for (j = 0; j < sites; ++j)
        {
          state = partition->tipmap[(int)(partition->tipchars[i][j])];
            invariant[j] &= state;
        }
    }
  }
  else
  {
    unsigned int span_padded = rate_cats * states_padded;
    for (i = 0; i < tips; ++i)
    {
      tipclv = partition->clv[i];
      for (j = 0; j < sites; ++j)
      {
        state = 0;
        for (k = 0; k < states; ++k)
        {
          state |= ((unsigned int)tipclv[k] << k);
        }
          invariant[j] &= state;
        tipclv += span_padded;
      }
    }
  }

  /* if all basecalls at current site are the same and not degenerate set the
     index in invariant to the frequency index of the basecall, otherwise -1 */
  for (i = 0; i < partition->sites; ++i)
  {
    if (partition->invariant[i] == 0 || __builtin_popcount((unsigned int)(partition->invariant[i])) > 1)
      partition->invariant[i] = -1;
    else
      partition->invariant[i] = __builtin_ctz((unsigned int)(partition->invariant[i]));
  }
  return PLL_SUCCESS;
}
