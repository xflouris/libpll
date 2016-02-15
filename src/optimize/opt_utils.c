/*
 Copyright (C) 2015 Diego Darriba

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

 Contact: Diego Darriba <Diego.Darriba@h-its.org>,
 Exelixis Lab, Heidelberg Instutute for Theoretical Studies
 Schloss-Wolfsbrunnenweg 35, D-69118 Heidelberg, Germany
 */
#include "pll_optimize.h"

PLL_EXPORT double * pll_compute_empirical_frequencies(pll_partition_t * partition)
{
  unsigned int i, j, k, n;
  double sum_test = 0.0;
  unsigned int states         = partition->states;
  unsigned int sites          = partition->sites;
  unsigned int cats           = partition->rate_cats;
  unsigned int tips           = partition->tips;
  const unsigned int * revmap = partition->revmap;
  const unsigned int * w      = partition->pattern_weights;
  double * frequencies = (double *) calloc (states, sizeof(double));

  if (partition->attributes & PLL_ATTRIB_PATTERN_TIP)
  {
    for (i = 0; i < tips; ++i)
    {
      const char *tipchars = partition->tipchars[i];
      for (n = 0; n < sites; ++n)
      {
        unsigned int state = revmap[(int) tipchars[n]];
        double sum_site = (double)__builtin_popcount(state);
        for (k = 0; k < states; ++k)
        {
          if (state & 1)
            frequencies[k] += w[n] / sum_site;
          state >>= 1;
        }
      }
    }
  }
  else
  {
    for (i = 0; i < tips; ++i)
    {
      for (n = 0, j = 0; j < sites * states * cats; j += (states * cats), ++n)
      {
        double sum_site = 0.0;
        for (k = 0; k < states; ++k)
          sum_site += partition->clv[i][j + k];
        for (k = 0; k < states; ++k)
          frequencies[k] += w[n] * partition->clv[i][j + k] / sum_site;
      }
    }
  }

#ifndef NDEBUG
  for (k = 0; k < states; ++k)
  {
    frequencies[k] /= sites * tips;
    printf("%.4f ", frequencies[k]);
    sum_test += frequencies[k];
  }
  printf("\n");
  assert(fabs (sum_test - 1) < 1e-6);
#endif

  return frequencies;
}

PLL_EXPORT double * pll_compute_empirical_subst_rates(pll_partition_t * partition)
{
  unsigned int i, j, k, n;
  unsigned int states         = partition->states;
  unsigned int sites          = partition->sites;
  unsigned int tips           = partition->tips;
  unsigned int cats           = partition->rate_cats;
  const unsigned int * revmap = partition->revmap;
  const unsigned int * w      = partition->pattern_weights;
  char * const * tipchars     = partition->tipchars;

  unsigned int n_subst_rates  = (states * (states - 1) / 2);
  double * subst_rates        = (double *) calloc (n_subst_rates, sizeof(double));

  unsigned *pair_rates = (unsigned *) alloca(
      states * states * sizeof(unsigned));
  memset (pair_rates, 0, sizeof(unsigned) * states * states);
  unsigned *state_freq = (unsigned *) alloca(states * sizeof(unsigned));

  unsigned int undef_state = pow (2, states) - 1;
  if (partition->attributes & PLL_ATTRIB_PATTERN_TIP)
  {
    for (n = 0; n < sites; ++n)
    {
      memset (state_freq, 0, sizeof(unsigned) * (states));
      for (i = 0; i < tips; ++i)
      {
        unsigned int state = revmap[(int) tipchars[i][n]];
        if (state == undef_state)
          continue;
        for (k = 0; k < states; ++k)
        {
          if (state & 1)
            state_freq[k]++;
          state >>= 1;
        }
      }

      for (i = 0; i < states; i++)
      {
        if (state_freq[i] == 0)
          continue;
        for (j = i + 1; j < states; j++)
        {
          pair_rates[i * states + j] += state_freq[i] * state_freq[j]
              * w[n];
        }
      }
    }
  }
  else
  {
    unsigned int cur_site = 0;
    for (n = 0; n < sites * states * cats; n += (states * cats))
    {
      memset (state_freq, 0, sizeof(unsigned) * (states));
      for (i = 0; i < tips; ++i)
      {
        int unstate = 1;
        for (k = 0; k < states; ++k)
          if (partition->clv[i][n + k] == 0)
          {
            unstate = 0;
            break;
          }
        if (unstate)
          continue;
        for (k = 0; k < states; ++k)
        {
          if (partition->clv[i][n + k])
          {
            state_freq[k]++;
          }
        }
      }

      for (i = 0; i < states; i++)
      {
        if (state_freq[i] == 0)
          continue;
        for (j = i + 1; j < states; j++)
        {
          pair_rates[i * states + j] += state_freq[i] * state_freq[j]
              * w[cur_site];
        }
      }
      cur_site++;
    }
  }

  k = 0;
  double last_rate = pair_rates[(states - 2) * states + states - 1];
  if (last_rate == 0)
    last_rate = 1;
  for (i = 0; i < states - 1; i++)
  {
    for (j = i + 1; j < states; j++)
    {
      subst_rates[k++] = pair_rates[i * states + j] / last_rate;
      if (subst_rates[k - 1] < 0.01)
        subst_rates[k - 1] = 0.01;
      if (subst_rates[k - 1] > 50.0)
        subst_rates[k - 1] = 50.0;
    }
  }

  subst_rates[k - 1] = 1.0;
  return subst_rates;
}

PLL_EXPORT double pll_compute_empirical_invariant_sites(pll_partition_t *partition)
{
  unsigned int n;
  unsigned int n_inv   = 0;
  unsigned int sites   = partition->sites;

  /* reset errno */
  pll_errno = 0;

  if (!partition->invariant)
    if (!pll_update_invariant_sites(partition))
      return NAN;

  const int * invariant = partition->invariant;

  for (n=0; n<sites; ++n)
    if (invariant[n] > -1) n_inv++;

  double empirical_pinv = (double)1.0*n_inv/sites;
  return empirical_pinv;
}
