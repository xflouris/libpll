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

#include <limits.h>
#include "pll.h"

/*
    TODO: Vectorize these functions which were created to make derivatives
          work with SSE (states_padded is set to the corresponding value)
*/

static int core_update_sumtable_ti_4x4_sse(unsigned int sites,
                                           unsigned int rate_cats,
                                           const double * parent_clv,
                                           const unsigned char * left_tipchars,
                                           const unsigned int * parent_scaler,
                                           double * const * eigenvecs,
                                           double * const * inv_eigenvecs,
                                           double * const * freqs_indices,
                                           const unsigned int * tipmap,
                                           double * sumtable,
                                           unsigned int attrib)
{
  unsigned int i,j,k,n;
  unsigned int tipstate;
  unsigned int states = 4;
  double lterm = 0;
  double rterm = 0;

  const double * clvc = parent_clv;
  const double * ev;
  const double * invei;
  const double * freqs;

  double * sum = sumtable;

  unsigned int min_scaler;
  unsigned int * rate_scalings = NULL;
  int per_rate_scaling = (attrib & PLL_ATTRIB_RATE_SCALERS) ? 1 : 0;

  /* powers of scale threshold for undoing the scaling */
  double scale_minlh[PLL_SCALE_RATE_MAXDIFF];
  if (per_rate_scaling)
  {
    rate_scalings = (unsigned int*) calloc(rate_cats, sizeof(unsigned int));

    double scale_factor = 1.0;
    for (i = 0; i < PLL_SCALE_RATE_MAXDIFF; ++i)
    {
      scale_factor *= PLL_SCALE_THRESHOLD;
      scale_minlh[i] = scale_factor;
    }
  }

  /* build sumtable */
  for (n = 0; n < sites; n++)
  {
    if (per_rate_scaling)
    {
      /* compute minimum per-rate scaler -> common per-site scaler */
      min_scaler = UINT_MAX;
      for (i = 0; i < rate_cats; ++i)
      {
        rate_scalings[i] = (parent_scaler) ? parent_scaler[n*rate_cats+i] : 0;
        if (rate_scalings[i] < min_scaler)
          min_scaler = rate_scalings[i];
      }

      /* compute relative capped per-rate scalers */
      for (i = 0; i < rate_cats; ++i)
      {
        rate_scalings[i] = PLL_MIN(rate_scalings[i] - min_scaler,
                                   PLL_SCALE_RATE_MAXDIFF);
      }
    }

    for (i = 0; i < rate_cats; ++i)
    {
      ev    = eigenvecs[i];
      invei = inv_eigenvecs[i];
      freqs = freqs_indices[i];

      for (j = 0; j < states; ++j)
      {
        tipstate = (unsigned int)left_tipchars[n];
        lterm = 0;
        rterm = 0;

        for (k = 0; k < states; ++k)
        {
          lterm += (tipstate & 1) * freqs[k] * invei[k*states+j];
          rterm += ev[j*states+k] * clvc[k];
          tipstate >>= 1;
        }
        sum[j] = lterm*rterm;

        if (rate_scalings && rate_scalings[i] > 0)
          sum[j] *= scale_minlh[rate_scalings[i]-1];
      }

      clvc += states;
      sum += states;
    }
  }

  if (rate_scalings)
    free(rate_scalings);

  return PLL_SUCCESS;
}

PLL_EXPORT int pll_core_update_sumtable_ii_sse(unsigned int states,
                                               unsigned int sites,
                                               unsigned int rate_cats,
                                               const double * parent_clv,
                                               const double * child_clv,
                                               const unsigned int * parent_scaler,
                                               const unsigned int * child_scaler,
                                               double * const * eigenvecs,
                                               double * const * inv_eigenvecs,
                                               double * const * freqs_indices,
                                               double * sumtable,
                                               unsigned int attrib)
{
  unsigned int i, j, k, n;
  double lterm = 0;
  double rterm = 0;

  const double * clvp = parent_clv;
  const double * clvc = child_clv;
  const double * ev;
  const double * invev;
  const double * freqs;

  double * sum = sumtable;

  unsigned int states_padded = (states+1) & 0xFFFFFFFE;

  unsigned int min_scaler;
  unsigned int * rate_scalings = NULL;
  int per_rate_scaling = (attrib & PLL_ATTRIB_RATE_SCALERS) ? 1 : 0;

  /* powers of scale threshold for undoing the scaling */
  double scale_minlh[PLL_SCALE_RATE_MAXDIFF];
  if (per_rate_scaling)
  {
    rate_scalings = (unsigned int*) calloc(rate_cats, sizeof(unsigned int));

    double scale_factor = 1.0;
    for (i = 0; i < PLL_SCALE_RATE_MAXDIFF; ++i)
    {
      scale_factor *= PLL_SCALE_THRESHOLD;
      scale_minlh[i] = scale_factor;
    }
  }

  /* build sumtable */
  for (n = 0; n < sites; n++)
  {
    if (per_rate_scaling)
    {
      /* compute minimum per-rate scaler -> common per-site scaler */
      min_scaler = UINT_MAX;
      for (i = 0; i < rate_cats; ++i)
      {
        rate_scalings[i] = (parent_scaler) ? parent_scaler[n*rate_cats+i] : 0;
        rate_scalings[i] += (child_scaler) ? child_scaler[n*rate_cats+i] : 0;
        if (rate_scalings[i] < min_scaler)
          min_scaler = rate_scalings[i];
      }

      /* compute relative capped per-rate scalers */
      for (i = 0; i < rate_cats; ++i)
      {
        rate_scalings[i] = PLL_MIN(rate_scalings[i] - min_scaler,
                                   PLL_SCALE_RATE_MAXDIFF);
      }
    }

    for (i = 0; i < rate_cats; ++i)
    {
      ev = eigenvecs[i];
      invev = inv_eigenvecs[i];
      freqs = freqs_indices[i];

      for (j = 0; j < states; ++j)
      {
        lterm = 0;
        rterm = 0;

        for (k = 0; k < states; ++k)
        {
          lterm += clvp[k] * freqs[k] * invev[k*states+j];
          rterm += ev[j*states+k] * clvc[k];
        }

        sum[j] = lterm*rterm;

        if (rate_scalings && rate_scalings[i] > 0)
          sum[j] *= scale_minlh[rate_scalings[i]-1];
      }

      clvc += states_padded;
      clvp += states_padded;
      sum += states_padded;
    }
  }

  if (rate_scalings)
    free(rate_scalings);

  return PLL_SUCCESS;
}

PLL_EXPORT int pll_core_update_sumtable_ti_sse(unsigned int states,
                                               unsigned int sites,
                                               unsigned int rate_cats,
                                               const double * parent_clv,
                                               const unsigned char * left_tipchars,
                                               const unsigned int * parent_scaler,
                                               double * const * eigenvecs,
                                               double * const * inv_eigenvecs,
                                               double * const * freqs_indices,
                                               const unsigned int * tipmap,
                                               double * sumtable,
                                               unsigned int attrib)
{
  unsigned int i,j,k,n;
  unsigned int tipstate;
  double lterm = 0;
  double rterm = 0;

  double * sum = sumtable;
  const double * clvc = parent_clv;
  const double * ev;
  const double * invev;
  const double * freqs;

  if (states == 4)
  {
    return core_update_sumtable_ti_4x4_sse(sites,
                                           rate_cats,
                                           parent_clv,
                                           left_tipchars,
                                           parent_scaler,
                                           eigenvecs,
                                           inv_eigenvecs,
                                           freqs_indices,
                                           tipmap,
                                           sumtable,
                                           attrib);
  }

  unsigned int states_padded = (states+1) & 0xFFFFFFFE;

  /* build sumtable: non-vectorized version, general case */
  for (n = 0; n < sites; n++)
  {
    for (i = 0; i < rate_cats; ++i)
    {
      ev    = eigenvecs[i];
      invev = inv_eigenvecs[i];
      freqs = freqs_indices[i];

      for (j = 0; j < states; ++j)
      {
        tipstate = tipmap[(unsigned int)left_tipchars[n]];
        lterm = 0;
        rterm = 0;

        for (k = 0; k < states; ++k)
        {
          lterm += (tipstate & 1) * freqs[k] * invev[k*states+j];
          rterm += ev[j*states+k] * clvc[k];
          tipstate >>= 1;
        }
        sum[j] = lterm*rterm;
      }

      clvc += states_padded;
      sum += states_padded;
    }
  }
  return PLL_SUCCESS;
}
