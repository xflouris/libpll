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

  #ifdef HAVE_SSE3
  if (attrib & PLL_ATTRIB_ARCH_SSE)
  {
    if (states == 4)
      pll_core_update_partial_tt_4x4_sse(sites,
                                         rate_cats,
                                         parent_clv,
                                         parent_scaler,
                                         left_tipchars,
                                         right_tipchars,
                                         lookup);
    else
      pll_core_update_partial_tt_sse(states,
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
                                         lookup,
                                         attrib);
    else
      pll_core_update_partial_tt_avx(states,
                                     sites,
                                     rate_cats,
                                     parent_clv,
                                     parent_scaler,
                                     left_tipchars,
                                     right_tipchars,
                                     lookup,
                                     tipmap_size,
                                     attrib);

    return;
  }
  #endif
  #ifdef HAVE_AVX2
  if (attrib & PLL_ATTRIB_ARCH_AVX2)
  {
    if (states == 4)
      pll_core_update_partial_tt_4x4_avx(sites,
                                         rate_cats,
                                         parent_clv,
                                         parent_scaler,
                                         left_tipchars,
                                         right_tipchars,
                                         lookup,
                                         attrib);
    else
      pll_core_update_partial_tt_avx(states,
                                     sites,
                                     rate_cats,
                                     parent_clv,
                                     parent_scaler,
                                     left_tipchars,
                                     right_tipchars,
                                     lookup,
                                     tipmap_size,
                                     attrib);

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

  #ifdef HAVE_SSE3
  if (attrib & PLL_ATTRIB_ARCH_SSE)
  {
    pll_core_update_partial_ti_4x4_sse(sites,
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
                                       right_scaler,
                                       attrib);
    return;
  }
  #endif
  #ifdef HAVE_AVX2
  if (attrib & PLL_ATTRIB_ARCH_AVX2)
  {
    pll_core_update_partial_ti_4x4_avx(sites,
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
                                           unsigned int tipmap_size,
                                           unsigned int attrib)
{
  int scaling;
  unsigned int i,j,k,n;
  unsigned int span = states * rate_cats;

  const double * lmat;
  const double * rmat;

#ifdef HAVE_SSE3
  if ((attrib & PLL_ATTRIB_ARCH_SSE))
  {
    if (states == 4)
      pll_core_update_partial_ti_4x4_sse(sites,
                                         rate_cats,
                                         parent_clv,
                                         parent_scaler,
                                         left_tipchars,
                                         right_clv,
                                         left_matrix,
                                         right_matrix,
                                         right_scaler);
    else
      pll_core_update_partial_ti_sse(states,
                                     sites,
                                     rate_cats,
                                     parent_clv,
                                     parent_scaler,
                                     left_tipchars,
                                     right_clv,
                                     left_matrix,
                                     right_matrix,
                                     right_scaler,
                                     tipmap,
                                     tipmap_size);
    return;
  }
#endif
#ifdef HAVE_AVX
  if ((attrib & PLL_ATTRIB_ARCH_AVX))
  {
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
                                   tipmap,
                                   tipmap_size,
                                   attrib);
    return;
  }
#endif
#ifdef HAVE_AVX2
  if ((attrib & PLL_ATTRIB_ARCH_AVX2))
  {
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
                                   tipmap,
                                   tipmap_size,
                                   attrib);
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

#ifdef HAVE_SSE3
  if (attrib & PLL_ATTRIB_ARCH_SSE)
  {
    pll_core_update_partial_ii_sse(states,
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
                                   right_scaler,
                                   attrib);
    return;
  }
#endif
#ifdef HAVE_AVX2
  if (attrib & PLL_ATTRIB_ARCH_AVX2)
  {
    pll_core_update_partial_ii_avx2(states,
                                   sites,
                                   rate_cats,
                                   parent_clv,
                                   parent_scaler,
                                   left_clv,
                                   right_clv,
                                   left_matrix,
                                   right_matrix,
                                   left_scaler,
                                   right_scaler,
                                   attrib);
    return;
  }
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

  #ifdef HAVE_SSE3
  if (attrib & PLL_ATTRIB_ARCH_SSE)
  {
    if (states == 4)
      pll_core_create_lookup_4x4_sse(rate_cats,
                                     lookup,
                                     left_matrix,
                                     right_matrix);
    else
      pll_core_create_lookup_sse(states,
                                 rate_cats,
                                 lookup,
                                 left_matrix,
                                 right_matrix,
                                 tipmap,
                                 tipmap_size);
    return;
  }
  #endif
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
  #ifdef HAVE_AVX2
  if (attrib & PLL_ATTRIB_ARCH_AVX2)
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

