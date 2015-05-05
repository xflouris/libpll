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

static int mytqli(double *d, double *e, const int n, double **z)
{
  int     m, l, iter, i, k;
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


static void mytred2(double **a, const int n, double *d, double *e)
{
  int     l, k, j, i;
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
                                   int states)
{
  int i,j,k,success;

  double ** qmatrix;

  /* normalize substitution parameters */
  int params_count = (states*states - states) / 2;
  double * params_normalized = (double *)malloc(sizeof(double) * params_count);
  if (!params_normalized) return NULL;
  memcpy(params_normalized,params,params_count*sizeof(double));
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
  printf("Here!\n");

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

void pll_update_prob_matrices(pll_partition_t * partition, 
                              int params_index, 
                              int * matrix_indices, 
                              double * branch_lengths, 
                              int count)
{
  int i,j,k,m,n;
  double *e, *d;
  double ** a;
  double * expd;
  double * temp;

  double * eigenvecs = partition->eigenvecs[params_index];
  double * inv_eigenvecs = partition->inv_eigenvecs[params_index];
  double * eigenvals = partition->eigenvals[params_index];
  double * rates = partition->rates;
  double * freqs = partition->frequencies[params_index];

  double * pmatrix;

  int states = partition->states;

  /* check whether we have cached an eigen decomposition. If not, compute it */
  if (!partition->eigen_decomp_valid[params_index]) 
  {
    a = create_ratematrix(partition->subst_params[params_index],
                          partition->frequencies[params_index],
                          partition->states);
    
    d = (double *)malloc(states*sizeof(double));
    e = (double *)malloc(states*sizeof(double));

    mytred2(a, states, d, e);
    mytqli(d, e, states, a);

    /* store eigen vectors */
    for (i = 0; i < states; ++i)
      memcpy(eigenvecs + i*states, a[i], states*sizeof(double));

    /* store eigen values */
    memcpy(eigenvals, d, states*sizeof(double));

    /* store inverse eigen vectors */
    for (k = 0, i = 0; i < states; ++i)
      for (j = i; j < states*states; j += states)
        inv_eigenvecs[k++] = eigenvecs[j];

    /* multiply the inverse eigen vectors from the left with sqrt(pi)^-1 */
    for (i = 0; i < states; ++i)
      for (j = 0; j < states; ++j)
        inv_eigenvecs[i*states+j] /= sqrt(freqs[i]);

    /* multiply the eigen vectors from the right with sqrt(pi) */
    for (i = 0; i < states; ++i)
      for (j = 0; j < states; ++j)
        eigenvecs[i*states+j] *= sqrt(freqs[j]);
  
    partition->eigen_decomp_valid[params_index] = 1;
    free(d);
    free(e);
    for (i = 0; i < states; ++i)
      free(a[i]);
    free(a);
  }

  expd = (double *)malloc(states * sizeof(double));
  temp = (double *)malloc(states*states * sizeof(double));
 
  for (i = 0; i < count; ++i)
  {
    for (n = 0; n < partition->rate_cats; ++n)
    {
      pmatrix = partition->pmatrix[matrix_indices[i]] + n*states*states;

      /* exponentiate eigenvalues */
      for (j = 0; j < states; ++j)
        expd[j] = exp(eigenvals[j] * rates[n] * branch_lengths[i] 
                                   / (1.0 - partition->prop_invar));

      for (j = 0; j < states; ++j)
        for (k = 0; k < states; ++k)
          temp[j*states+k] = inv_eigenvecs[j*states+k] * expd[k];
      
      for (j = 0; j < states; ++j)
        for (k = 0; k < states; ++k)
        {
          pmatrix[j*states+k] = 0;
          for (m = 0; m < states; ++m)
          {
            pmatrix[j*states+k] += temp[j*states+m] * eigenvecs[m*states+k];
          }
        }
    }
  }
  free(expd);
  free(temp);
}


void pll_set_frequencies(pll_partition_t * partition, 
                         int params_index, 
                         double * frequencies)
{
  memcpy(partition->frequencies[params_index], 
         frequencies, 
         partition->states*sizeof(double));
}

void pll_set_category_rates(pll_partition_t * partition, double * rates)
{
  memcpy(partition->rates, rates, partition->rate_cats*sizeof(double));
}

void pll_set_subst_params(pll_partition_t * partition, 
                          int params_index, 
                          double * params, 
                          int count)
{
  memcpy(partition->subst_params[params_index], params, count*sizeof(double));
}

PLL_EXPORT int pll_update_invariant_sites(pll_partition_t * partition, 
                                          double prop_invar)
{
  int i,j,k;
  double * tipclv;
  int state;

  if (prop_invar < 0 || prop_invar >= 1) 
    return PLL_FAILURE;

  partition->prop_invar = prop_invar;

  if (!partition->invariant)
  {
    partition->invariant = (int *)calloc(partition->sites, sizeof(int));
    if (!partition->invariant) return PLL_FAILURE;
  }
  for (i = 0; i < partition->tips; ++i)
  {
    tipclv = partition->clv[i];
    for (j = 0; j < partition->sites; ++j)
    {
      state = 0;
      for (k = 0; k < partition->states; ++k)
        state |= ((int)tipclv[k] << k);
      partition->invariant[j] |= state;
      tipclv += (partition->rate_cats * partition->sites);
    }
  }

  /* if all basecalls at current site are the same and not degenerate set the 
     index in invariant to the frequency index of the basecall, otherwise -1 */
  for (i = 0; i < partition->sites; ++i)
    if (__builtin_popcount(partition->invariant[i]) > 1)
      partition->invariant[i] = -1;
    else
      partition->invariant[i] = __builtin_ctz(partition->invariant[i]);
  
  return PLL_SUCCESS;
}


