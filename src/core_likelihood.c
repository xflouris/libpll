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

static void fill_parent_scaler(unsigned int sites,
                               unsigned int * parent_scaler,
                               const unsigned int * left_scaler,
                               const unsigned int * right_scaler)
{
  unsigned int i;

  if (!left_scaler && !right_scaler)
    memset(parent_scaler, 0, sizeof(unsigned int) * sites);
  else if (left_scaler && right_scaler)
  {
    memcpy(parent_scaler, left_scaler, sizeof(unsigned int) * sites);
    for (i = 0; i < sites; ++i)
      parent_scaler[i] += right_scaler[i];
  }
  else
  {
    if (left_scaler)
      memcpy(parent_scaler, left_scaler, sizeof(unsigned int) * sites);
    else
      memcpy(parent_scaler, right_scaler, sizeof(unsigned int) * sites);
  }
}

PLL_EXPORT void pll_core_update_partial_tt_4x4(unsigned int sites,
                                               unsigned int rate_cats,
                                               double * parent_clv,
                                               unsigned int * parent_scaler,
                                               const unsigned char * left_tipchars,
                                               const unsigned char * right_tipchars,
                                               const double * lookup)
{
  unsigned int j,k,n;
  unsigned int states = 4;
  unsigned int span = states * rate_cats;
  const double * offset;

  if (parent_scaler)
    memset(parent_scaler, 0, sizeof(unsigned int) * sites);

  for (n = 0; n < sites; ++n)
  {
    j = (unsigned int)(left_tipchars[n]);
    k = (unsigned int)(right_tipchars[n]);

    offset = lookup;
    offset += ((j << 4) + k)*span;

    memcpy(parent_clv, offset, span*sizeof(double));

    parent_clv += span;
  }
}

PLL_EXPORT void pll_core_update_partial_tt(unsigned int states,
                                           unsigned int sites,
                                           unsigned int rate_cats,
                                           double * parent_clv,
                                           unsigned int * parent_scaler,
                                           const unsigned char * left_tipchars,
                                           const unsigned char * right_tipchars,
                                           const unsigned int * tipmap,
                                           unsigned int tipmap_size,
                                           const double * lookup,
                                           unsigned int attrib)
{
  unsigned int j,k,n;
  const double * offset;

  #ifdef HAVE_AVX
  if (attrib & PLL_ATTRIB_ARCH_AVX)
  {
    if (states == 4)
      pll_core_update_partial_tt_4x4_avx(sites,
                                         rate_cats,
                                         parent_clv,
                                         parent_scaler,
                                         left_tipchars,
                                         right_tipchars,
                                         lookup);
    else
      pll_core_update_partial_tt_avx(states,
                                     sites,
                                     rate_cats,
                                     parent_clv,
                                     parent_scaler,
                                     left_tipchars,
                                     right_tipchars,
                                     lookup,
                                     tipmap_size);

    return;
  }
  #endif

  unsigned int span = states * rate_cats;
  unsigned int log2_maxstates = (unsigned int)ceil(log2(tipmap_size));

  if (parent_scaler)
    memset(parent_scaler, 0, sizeof(unsigned int) * sites);

  for (n = 0; n < sites; ++n)
  {
    j = (unsigned int)(left_tipchars[n]);
    k = (unsigned int)(right_tipchars[n]);

    offset = lookup;
    offset += ((j << log2_maxstates) + k)*span;

    memcpy(parent_clv, offset, span*sizeof(double));

    parent_clv += span;
  }
}

PLL_EXPORT void pll_core_update_partial_ti_4x4(unsigned int sites,
                                               unsigned int rate_cats,
                                               double * parent_clv,
                                               unsigned int * parent_scaler,
                                               const unsigned char * left_tipchars,
                                               const double * right_clv,
                                               const double * left_matrix,
                                               const double * right_matrix,
                                               const unsigned int * right_scaler,
                                               unsigned int attrib)
{
  unsigned int states = 4;
  unsigned int scaling;
  unsigned int i,j,k,n;
  unsigned int span = states * rate_cats;

  const double * lmat;
  const double * rmat;

  #ifdef HAVE_AVX
  if (attrib & PLL_ATTRIB_ARCH_AVX)
  {
    pll_core_update_partial_ti_4x4_avx(sites,
                                       rate_cats,
                                       parent_clv,
                                       parent_scaler,
                                       left_tipchars,
                                       right_clv,
                                       left_matrix,
                                       right_matrix,
                                       right_scaler);
    return;
  }
  #endif

  if (parent_scaler)
    fill_parent_scaler(sites, parent_scaler, NULL, right_scaler);

  for (n = 0; n < sites; ++n)
  {
    lmat = left_matrix;
    rmat = right_matrix;

    scaling = (parent_scaler) ? 1 : 0;

    for (k = 0; k < rate_cats; ++k)
    {
      for (i = 0; i < states; ++i)
      {
        double terma = 0;
        double termb = 0;
        unsigned int lstate = left_tipchars[n];
        for (j = 0; j < states; ++j)
        {
          if (lstate & 1)
            terma += lmat[j];

          termb += rmat[j] * right_clv[j];

          lstate >>= 1;
        }
        parent_clv[i] = terma*termb;
        lmat += states;
        rmat += states;

        scaling = scaling && (parent_clv[i] < PLL_SCALE_THRESHOLD);
      }
      parent_clv += states;
      right_clv  += states;
    }
    /* if *all* entries of the site CLV were below the threshold then scale
       (all) entries by PLL_SCALE_FACTOR */
    if (scaling)
    {
      parent_clv -= span;
      for (i = 0; i < span; ++i)
        parent_clv[i] *= PLL_SCALE_FACTOR;
      parent_clv += span;
      parent_scaler[n] += 1;
    }
  }
}

PLL_EXPORT void pll_core_update_partial_ti(unsigned int states,
                                           unsigned int sites,
                                           unsigned int rate_cats,
                                           double * parent_clv,
                                           unsigned int * parent_scaler,
                                           const unsigned char * left_tipchars,
                                           const double * right_clv,
                                           const double * left_matrix,
                                           const double * right_matrix,
                                           const unsigned int * right_scaler,
                                           const unsigned int * tipmap,
                                           unsigned int attrib)
{
  int scaling;
  unsigned int i,j,k,n;
  unsigned int span = states * rate_cats;

  const double * lmat;
  const double * rmat;

#ifdef HAVE_AVX
    if ((attrib & PLL_ATTRIB_ARCH_AVX))
    {
      if (states == 4)
        pll_core_update_partial_ti_4x4_avx(sites,
                                           rate_cats,
                                           parent_clv,
                                           parent_scaler,
                                           left_tipchars,
                                           right_clv,
                                           left_matrix,
                                           right_matrix,
                                           right_scaler);
      else
        pll_core_update_partial_ti_avx(states,
                                       sites,
                                       rate_cats,
                                       parent_clv,
                                       parent_scaler,
                                       left_tipchars,
                                       right_clv,
                                       left_matrix,
                                       right_matrix,
                                       right_scaler,
                                       tipmap);
      return;
    }
#endif
    if (states == 4)
    {
      pll_core_update_partial_ti_4x4(sites,
                                     rate_cats,
                                     parent_clv,
                                     parent_scaler,
                                     left_tipchars,
                                     right_clv,
                                     left_matrix,
                                     right_matrix,
                                     right_scaler,
                                     attrib);
      return;
    }

  if (parent_scaler)
    fill_parent_scaler(sites, parent_scaler, NULL, right_scaler);

  for (n = 0; n < sites; ++n)
  {
    lmat = left_matrix;
    rmat = right_matrix;

    scaling = (parent_scaler) ? 1 : 0;

    for (k = 0; k < rate_cats; ++k)
    {
      for (i = 0; i < states; ++i)
      {
        double terma = 0;
        double termb = 0;
        unsigned int lstate = tipmap[(unsigned int)left_tipchars[n]];
        for (j = 0; j < states; ++j)
        {
          if (lstate & 1)
            terma += lmat[j];

          termb += rmat[j] * right_clv[j];

          lstate >>= 1;
        }
        parent_clv[i] = terma*termb;
        lmat += states;
        rmat += states;

        scaling = scaling && (parent_clv[i] < PLL_SCALE_THRESHOLD);
      }
      parent_clv += states;
      right_clv  += states;
    }
    /* if *all* entries of the site CLV were below the threshold then scale
       (all) entries by PLL_SCALE_FACTOR */
    if (scaling)
    {
      parent_clv -= span;
      for (i = 0; i < span; ++i)
        parent_clv[i] *= PLL_SCALE_FACTOR;
      parent_clv += span;
      parent_scaler[n] += 1;
    }
  }
}

PLL_EXPORT void pll_core_update_partial_ii(unsigned int states,
                                           unsigned int sites,
                                           unsigned int rate_cats,
                                           double * parent_clv,
                                           unsigned int * parent_scaler,
                                           const double * left_clv,
                                           const double * right_clv,
                                           const double * left_matrix,
                                           const double * right_matrix,
                                           const unsigned int * left_scaler,
                                           const unsigned int * right_scaler,
                                           unsigned int attrib)
{
  unsigned int i,j,k,n;
  unsigned int scaling;

  const double * lmat;
  const double * rmat;

  unsigned int span = states * rate_cats;

#ifdef HAVE_AVX
  if (attrib & PLL_ATTRIB_ARCH_AVX)
  {
    pll_core_update_partial_ii_avx(states,
                                   sites,
                                   rate_cats,
                                   parent_clv,
                                   parent_scaler,
                                   left_clv,
                                   right_clv,
                                   left_matrix,
                                   right_matrix,
                                   left_scaler,
                                   right_scaler);
    return;
  }
#endif
#ifdef HAVE_SSE
#endif

  /* add up the scale vectors of the two children if available */
  if (parent_scaler)
    fill_parent_scaler(sites, parent_scaler, left_scaler, right_scaler);

  /* compute CLV */
  for (n = 0; n < sites; ++n)
  {
    lmat = left_matrix;
    rmat = right_matrix;
    scaling = (parent_scaler) ? 1 : 0;

    for (k = 0; k < rate_cats; ++k)
    {
      for (i = 0; i < states; ++i)
      {
        double terma = 0;
        double termb = 0;
        for (j = 0; j < states; ++j)
        {
          terma += lmat[j] * left_clv[j];
          termb += rmat[j] * right_clv[j];
        }
        parent_clv[i] = terma*termb;
        lmat += states;
        rmat += states;

        scaling = scaling && (parent_clv[i] < PLL_SCALE_THRESHOLD);
      }
      parent_clv += states;
      left_clv   += states;
      right_clv  += states;
    }
    /* if *all* entries of the site CLV were below the threshold then scale
       (all) entries by PLL_SCALE_FACTOR */
    if (scaling)
    {
      parent_clv -= span;
      for (i = 0; i < span; ++i)
        parent_clv[i] *= PLL_SCALE_FACTOR;
      parent_clv += span;
      parent_scaler[n] += 1;
    }
  }
}

PLL_EXPORT void pll_core_create_lookup_4x4(unsigned int rate_cats,
                                           double * lookup,
                                           const double * left_matrix,
                                           const double * right_matrix)
{
  unsigned int i,j,k,n,m;
  unsigned int maxstates = 16;
  unsigned int states = 4;
  unsigned int index = 0;

  /* precompute first the entries that contain only one 1 */
  double termj = 0;
  double termk = 0;

  const double * jmat;
  const double * kmat;

  /* go through all pairs j,k of states for the two tips; i is the inner
     node state */
  for (j = 0; j < maxstates; ++j)
  {
    for (k = 0; k < maxstates; ++k)
    {
      jmat = left_matrix;
      kmat = right_matrix;

      /* precompute the likelihood for each state and each rate */
      for (n = 0; n < rate_cats; ++n)
      {
        for (i = 0; i < states; ++i)
        {
          termj = 0;
          termk = 0;

          unsigned int jstate = j;
          unsigned int kstate = k;

          /* decompose basecall into the encoded residues and set the appropriate
             positions in the tip vector */
          for (m = 0; m < states; ++m)
          {
            if (jstate & 1)
              termj += jmat[m];

            if (kstate & 1)
              termk += kmat[m];

            jstate >>= 1;
            kstate >>= 1;
          }

          jmat += states;
          kmat += states;
          lookup[index++] = termj*termk;
        }
      }
    }
  }
}

PLL_EXPORT void pll_core_create_lookup(unsigned int states,
                                       unsigned int rate_cats,
                                       double * lookup,
                                       const double * left_matrix,
                                       const double * right_matrix,
                                       unsigned int * tipmap,
                                       unsigned int tipmap_size,
                                       unsigned int attrib)
{

  #ifdef HAVE_AVX
  if (attrib & PLL_ATTRIB_ARCH_AVX)
  {
    if (states == 4)
      pll_core_create_lookup_4x4_avx(rate_cats,
                                     lookup,
                                     left_matrix,
                                     right_matrix);
    else
      pll_core_create_lookup_avx(states,
                                 rate_cats,
                                 lookup,
                                 left_matrix,
                                 right_matrix,
                                 tipmap,
                                 tipmap_size);
    return;
  }
  #endif
  if (states == 4)
  {
    pll_core_create_lookup_4x4(rate_cats,
                               lookup,
                               left_matrix,
                               right_matrix);
    return;
  }

  unsigned int i,j,k,n,m;
  unsigned int index = 0;
  unsigned int maxstates = tipmap_size;

  unsigned int log2_maxstates = (unsigned int)ceil(log2(maxstates));
  unsigned int span = states*rate_cats;

  /* precompute first the entries that contain only one 1 */
  double termj = 0;
  double termk = 0;

  const double * jmat;
  const double * kmat;
  double * lh_statepair;

  /* go through all pairs j,k of states for the two tips; i is the inner
     node state */
  for (j = 0; j < maxstates; ++j)
  {
    for (k = 0; k < maxstates; ++k)
    {
      jmat = left_matrix;
      kmat = right_matrix;
      index = 0;

      /* find offset of state-pair in the precomputation table */
      lh_statepair = lookup;
      lh_statepair += ((j << log2_maxstates) + k)*span;

      /* precompute the likelihood for each state and each rate */
      for (n = 0; n < rate_cats; ++n)
      {
        for (i = 0; i < states; ++i)
        {
          termj = 0;
          termk = 0;

          unsigned int jstate = tipmap[j];
          unsigned int kstate = tipmap[k];

          /* decompose basecall into the encoded residues and set the appropriate
             positions in the tip vector */
          for (m = 0; m < states; ++m)
          {
            if (jstate & 1)
              termj += jmat[m];

            if (kstate & 1)
              termk += kmat[m];

            jstate >>= 1;
            kstate >>= 1;
          }

          jmat += states;
          kmat += states;
          lh_statepair[index++] = termj*termk;
        }
      }
    }
  }
}

PLL_EXPORT int pll_core_update_sumtable_ti_4x4(unsigned int sites,
                                               unsigned int rate_cats,
                                               const double * parent_clv,
                                               const unsigned char * left_tipchars,
                                               double ** eigenvecs,
                                               double ** inv_eigenvecs,
                                               double ** freqs,
                                               unsigned int * tipmap,
                                               double *sumtable,
                                               unsigned int attrib)
{
  unsigned int i, j, k, n;
  unsigned int tipstate;
  double lefterm = 0;
  double righterm = 0;

  double * sum = sumtable;
  const double * t_clvc = parent_clv;
  const double * t_eigenvecs;
  const double * t_inv_eigenvecs;
  const double * t_freqs;

  unsigned int states = 4;

  /* build sumtable */
  for (n = 0; n < sites; n++)
  {
    for (i = 0; i < rate_cats; ++i)
    {
      t_eigenvecs = eigenvecs[i];
      t_inv_eigenvecs = inv_eigenvecs[i];
      t_freqs = freqs[i];

      for (j = 0; j < states; ++j)
      {
        tipstate = (unsigned int) left_tipchars[n];
        lefterm = 0;
        righterm = 0;
        for (k = 0; k < states; ++k)
        {
          lefterm += (tipstate & 1) * t_freqs[k]
              * t_inv_eigenvecs[k * states + j];
          righterm += t_eigenvecs[j * states + k] * t_clvc[k];
          tipstate >>= 1;
        }
        sum[j] = lefterm * righterm;
      }

      t_clvc += states;
      sum += states;
    }
  }

  return PLL_SUCCESS;
}

PLL_EXPORT int pll_core_update_sumtable_ii(unsigned int states,
                                           unsigned int sites,
                                           unsigned int rate_cats,
                                           const double * parent_clv,
                                           const double * child_clv,
                                           double ** eigenvecs,
                                           double ** inv_eigenvecs,
                                           double ** freqs,
                                           double *sumtable,
                                           unsigned int attrib)
{
  unsigned int i, j, k, n;
  double lefterm  = 0;
  double righterm = 0;

  double * sum                   = sumtable;
  const double * t_clvp          = parent_clv;
  const double * t_clvc          = child_clv;
  const double * t_eigenvecs;
  const double * t_inv_eigenvecs;
  const double * t_freqs;

#ifdef HAVE_AVX
  if (attrib & PLL_ATTRIB_ARCH_AVX)
  {
    return pll_core_update_sumtable_ii_avx(states,
                                           sites,
                                           rate_cats,
                                           parent_clv,
                                           child_clv,
                                           eigenvecs,
                                           inv_eigenvecs,
                                           freqs,
                                           sumtable);
  }
#endif
#ifdef HAVE_SSE
#endif

  /* build sumtable */
  for (n = 0; n < sites; n++)
  {
    for (i = 0; i < rate_cats; ++i)
    {
      t_eigenvecs     = eigenvecs[i];
      t_inv_eigenvecs = inv_eigenvecs[i];
      t_freqs         = freqs[i];

      for (j = 0; j < states; ++j)
      {
        lefterm = 0;
        righterm = 0;
        for (k = 0; k < states; ++k)
        {
          lefterm  += t_clvp[k] * t_freqs[k] * t_inv_eigenvecs[k * states + j];
          righterm += t_eigenvecs[j * states + k] * t_clvc[k];
        }
        sum[j] = lefterm * righterm;
      }
      t_clvc += states;
      t_clvp += states;
      sum += states;
    }
  }

  return PLL_SUCCESS;
}

PLL_EXPORT int pll_core_update_sumtable_ti(unsigned int states,
                                           unsigned int sites,
                                           unsigned int rate_cats,
                                           const double * parent_clv,
                                           const unsigned char * left_tipchars,
                                           double ** eigenvecs,
                                           double ** inv_eigenvecs,
                                           double ** freqs,
                                           unsigned int * tipmap,
                                           double *sumtable,
                                           unsigned int attrib)
{
  unsigned int i, j, k, n;
  unsigned int tipstate;
  double lefterm = 0;
  double righterm = 0;

  double * sum = sumtable;
  const double * t_clvc = parent_clv;
  const double * t_eigenvecs;
  const double * t_inv_eigenvecs;
  const double * t_freqs;

  unsigned int states_padded = states;

  if (states == 4)
  {
    return pll_core_update_sumtable_ti_4x4(sites,
                                    rate_cats,
                                    parent_clv,
                                    left_tipchars,
                                    eigenvecs,
                                    inv_eigenvecs,
                                    freqs,
                                    tipmap,
                                    sumtable,
                                    attrib);
  }

#ifdef HAVE_AVX
  if (attrib & PLL_ATTRIB_ARCH_AVX)
  {
    states_padded = (states+3) & 0xFFFFFFFC;
  }
#endif
#ifdef HAVE_SSE
#endif

  /* build sumtable */
  for (n = 0; n < sites; n++)
  {
    for (i = 0; i < rate_cats; ++i)
    {
      t_eigenvecs = eigenvecs[i];
      t_inv_eigenvecs = inv_eigenvecs[i];
      t_freqs = freqs[i];

      for (j = 0; j < states; ++j)
      {
        tipstate = tipmap[(unsigned int)left_tipchars[n]];
        lefterm = 0;
        righterm = 0;
        for (k = 0; k < states; ++k)
        {
          lefterm += (tipstate & 1) * t_freqs[k]
              * t_inv_eigenvecs[k * states + j];
          righterm += t_eigenvecs[j * states + k] * t_clvc[k];
          tipstate >>= 1;
        }
        sum[j] = lefterm * righterm;
      }
      t_clvc += states_padded;
      sum += states_padded;
    }
  }

  return PLL_SUCCESS;
}

static void core_site_likelihood_derivatives(unsigned int states,
                                             unsigned int states_padded,
                                             unsigned int rate_cats,
                                             const double * rate_weights,
                                             const int * invariant,
                                             const double * prop_invar,
                                             double ** freqs,
                                             const double * sumtable,
                                             const double * diagptable,
                                             double * site_lk)
{
  unsigned int i,j;
  double inv_site_lk = 0.0;
  double cat_sitelk[3];
  const double *sum = sumtable;
  const double * diagp = diagptable;
  const double * t_freqs;
  double t_prop_invar;

  site_lk[0] = site_lk[1] = site_lk[2] = 0;
  for (i = 0; i < rate_cats; ++i)
  {
    t_freqs = freqs[i];

    cat_sitelk[0] = cat_sitelk[1] = cat_sitelk[2] = 0;
    for (j = 0; j < states; ++j)
    {
      cat_sitelk[0] += sum[j] * diagp[0];
      cat_sitelk[1] += sum[j] * diagp[1];
      cat_sitelk[2] += sum[j] * diagp[2];
      diagp += 3;
    }

    /* account for invariant sites */
    t_prop_invar = prop_invar[i];
    if (t_prop_invar > 0)
    {
      inv_site_lk =
          (invariant[0] == -1) ? 0 : t_freqs[invariant[0]] * t_prop_invar;
      cat_sitelk[0] = cat_sitelk[0] * (1. - t_prop_invar) + inv_site_lk;
      cat_sitelk[1] = cat_sitelk[1] * (1. - t_prop_invar);
      cat_sitelk[2] = cat_sitelk[2] * (1. - t_prop_invar);
    }

    site_lk[0] += cat_sitelk[0] * rate_weights[i];
    site_lk[1] += cat_sitelk[1] * rate_weights[i];
    site_lk[2] += cat_sitelk[2] * rate_weights[i];

    sum += states_padded;
  }
}

PLL_EXPORT double pll_core_likelihood_derivatives(unsigned int states,
                                                  unsigned int sites,
                                                  unsigned int rate_cats,
                                                  const double * rate_weights,
                                                  const unsigned int * parent_scaler,
                                                  const unsigned int * child_scaler,
                                                  const int * invariant,
                                                  const unsigned int * pattern_weights,
                                                  double branch_length,
                                                  const double * prop_invar,
                                                  double ** freqs,
                                                  const const double * rates,
                                                  double ** eigenvals,
                                                  const double * sumtable,
                                                  double * d_f,
                                                  double * dd_f,
                                                  unsigned int attrib)
{
  unsigned int n, i, j;
  unsigned int ef_sites;
  double site_lk[3];
  double logLK = 0.0;

  const double * sum;
  double deriv1, deriv2;

  const double * t_eigenvals;
  double t_branch_length;
  unsigned int scale_factors;

  double *diagptable, *diagp;
  const int * invariant_ptr;
  double ki;

  unsigned int states_padded = states;
  unsigned int pattern_weight_sum = 0;

#ifdef HAVE_AVX
  if (attrib & PLL_ATTRIB_ARCH_AVX)
  {
    states_padded = (states+3) & 0xFFFFFFFC;
  }
#endif
#ifdef HAVE_SSE
#endif

  /* For Stamatakis correction, the likelihood derivatives are computed in
     the usual way for the additional per-state sites. */
  if ((attrib & PLL_ATTRIB_ASC_BIAS_MASK) == PLL_ATTRIB_ASC_BIAS_STAMATAKIS)
  {
    ef_sites = sites + states;
  }
  else
  {
    ef_sites = sites;
  }

  *d_f = 0.0;
  *dd_f = 0.0;

  diagptable = (double *) calloc (rate_cats * states * 3, sizeof(double));
  if (!diagptable)
  {
    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf (pll_errmsg, 200, "Cannot allocate memory for diagptable");
    return -INFINITY;
  }

  /* pre-compute the derivatives of the P matrix for all discrete GAMMA rates */
  diagp = diagptable;
  for(i = 0; i < rate_cats; ++i)
  {
    t_eigenvals = eigenvals[i];
    ki = rates[i];
    t_branch_length = branch_length/(1.0 - prop_invar[i]);
    for(j = 0; j < states; ++j)
    {
      diagp[0] = exp(t_eigenvals[j] * ki * t_branch_length);
      diagp[1] = t_eigenvals[j] * ki * diagp[0];
      diagp[2] = t_eigenvals[j] * ki * t_eigenvals[j] * ki * diagp[0];
      diagp += 3;
    }
  }

  sum = sumtable;
  invariant_ptr = invariant;
  for (n = 0; n < ef_sites; ++n)
  {
    core_site_likelihood_derivatives(states,
                                     states_padded,
                                     rate_cats,
                                     rate_weights,
                                     invariant_ptr,
                                     prop_invar,
                                     freqs,
                                     sum,
                                     diagptable,
                                     site_lk);
    invariant_ptr++;
    sum += rate_cats * states_padded;

    scale_factors = (parent_scaler) ? parent_scaler[n] : 0;
    scale_factors += (child_scaler) ? child_scaler[n] : 0;

    logLK += log (site_lk[0]) * pattern_weights[n];
    if (scale_factors)
    {
      logLK += scale_factors * log (PLL_SCALE_THRESHOLD);
    }

    /* build derivatives */
    deriv1 = (-site_lk[1] / site_lk[0]);
    deriv2 = (deriv1 * deriv1 - (site_lk[2] / site_lk[0]));
    *d_f += pattern_weights[n] * deriv1;
    *dd_f += pattern_weights[n] * deriv2;
    pattern_weight_sum += pattern_weights[n];
  }

  /* account for ascertainment bias correction */
  if (attrib & PLL_ATTRIB_ASC_BIAS_MASK)
  {
    double asc_Lk[3] = {0.0, 0.0, 0.0};
    unsigned int sum_w_inv = 0;
    double asc_scaling;
    int asc_bias_type = attrib & PLL_ATTRIB_ASC_BIAS_MASK;

    if (asc_bias_type != PLL_ATTRIB_ASC_BIAS_STAMATAKIS)
    {
      for (n=0; n<states; ++n)
      {
        /* compute the site LK derivatives for the additional per-state sites */
        core_site_likelihood_derivatives(states,
                                         states_padded,
                                         rate_cats,
                                         rate_weights,
                                         0, /* prop invar disallowed */
                                         prop_invar,
                                         freqs,
                                         sum,
                                         diagptable,
                                         site_lk);
        sum += rate_cats * states_padded;

        scale_factors = (parent_scaler) ? parent_scaler[n] : 0;
        scale_factors += (child_scaler) ? child_scaler[n] : 0;
        asc_scaling = pow(PLL_SCALE_THRESHOLD, (double)scale_factors);

        /* sum over likelihood and 1st and 2nd derivative / apply scaling */
        asc_Lk[0] += site_lk[0] * asc_scaling;
        asc_Lk[1] += site_lk[1] * asc_scaling;
        asc_Lk[2] += site_lk[2] * asc_scaling;

        sum_w_inv += pattern_weights[sites + n];
      }

      switch(asc_bias_type)
      {
        case PLL_ATTRIB_ASC_BIAS_LEWIS:
          /* correct log-likelihood */
          logLK -= pattern_weight_sum * log(1.0 - asc_Lk[0]);

          /* derivatives of log(1.0 - (sum Li(s) over states 's')) */
      		*d_f  -= pattern_weight_sum * (asc_Lk[1] / (asc_Lk[0] - 1.0));
      		*dd_f -= pattern_weight_sum *
                     (((asc_Lk[0] - 1.0) * asc_Lk[2] - asc_Lk[1] * asc_Lk[1]) /
                     ((asc_Lk[0] - 1.0) * (asc_Lk[0] - 1.0)));
          break;
        case PLL_ATTRIB_ASC_BIAS_FELSENSTEIN:
          /* correct log-likelihood */
          logLK += sum_w_inv * log(asc_Lk[0]);

          /* derivatives of log(sum Li(s) over states 's') */
      		*d_f  += sum_w_inv * (asc_Lk[1] / asc_Lk[0]);
      		*dd_f += sum_w_inv *
                     (((asc_Lk[2] * asc_Lk[0]) - asc_Lk[1] * asc_Lk[1]) /
                     (asc_Lk[0] * asc_Lk[0]));
        break;
        default:
          pll_errno = PLL_ERROR_ASC_BIAS;
          snprintf(pll_errmsg, 200, "Illegal ascertainment bias algorithm");
          return -INFINITY;
      }
    }
  }

  free (diagptable);
  return logLK;
}
