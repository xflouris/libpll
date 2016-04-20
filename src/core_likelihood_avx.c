/*
    Copyright (C) 2016 Tomas Flouri, Kassian Kobert

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

PLL_EXPORT void pll_core_create_lookup_avx(unsigned int states,
                                           unsigned int rate_cats,
                                           double * ttlookup,
                                           const double * left_matrix,
                                           const double * right_matrix,
                                           unsigned int * tipmap,
                                           unsigned int tipmap_size)
{
  if (states == 4)
  {
    pll_core_create_lookup_4x4_avx(rate_cats,
                                   ttlookup,
                                   left_matrix,
                                   right_matrix);
    return;
  }

  unsigned int i,j,k,n,m;
  unsigned int states_padded = (states+3) & 0xFFFFFFFC;
  unsigned int maxstates = tipmap_size;
  unsigned int index = 0;

  unsigned int log2_maxstates = (unsigned int)ceil(log2(maxstates));
  unsigned int log2_states_padded = (unsigned int)ceil(log2(states_padded));
  unsigned int log2_rates = (unsigned int)ceil(log2(rate_cats));

  /* precompute first the entries that contain only one 1 */
  double termj = 0;
  double termk = 0;

  const double * jmat;
  const double * kmat;
  double * lookup;

  /* go through all pairs j,k of states for the two tips; i is the inner
     node state */
  for (j = 0; j < maxstates; ++j)
  {
    for (k = 0; k < maxstates; ++k)
    {
      jmat = left_matrix;
      kmat = right_matrix;

      /* find offset of state-pair in the precomputation table */
      lookup = ttlookup;
      lookup += ((j << log2_maxstates) + k) << (log2_states_padded+log2_rates);

      /* precompute the likelihood for each state and each rate */
      for (n = 0; n < rate_cats; ++n)
      {
        index = 0;
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

          jmat += states_padded;
          kmat += states_padded;
          lookup[index++] = termj*termk;
        }
        lookup += states_padded;
      }
    }
  }
}

PLL_EXPORT void pll_core_create_lookup_4x4_avx(unsigned int rate_cats,
                                               double * lookup,
                                               const double * left_matrix,
                                               const double * right_matrix)
{
  unsigned int j,k,n;
  unsigned int maxstates = 16;
  unsigned int states = 4;

  __m256d ymm0,ymm1,ymm2,ymm3,ymm4,ymm5,ymm6,ymm7;
  __m256d xmm0,xmm1,xmm2,xmm3,xmm4,xmm5,xmm6,xmm7;
  __m256i jmask,kmask;

  const double * jmat;
  const double * kmat;

  for (j = 0; j < maxstates; ++j)
  {
    for (k = 0; k < maxstates; ++k)
    {
      jmat = left_matrix;
      kmat = right_matrix;
    
      jmask = _mm256_set_epi64x(
                 ((j >> 3) & 1) ? ~0 : 0,
                 ((j >> 2) & 1) ? ~0 : 0,
                 ((j >> 1) & 1) ? ~0 : 0,
                 (j & 1) ? ~0 : 0);

      kmask = _mm256_set_epi64x(
                 ((k >> 3) & 1) ? ~0 : 0,
                 ((k >> 2) & 1) ? ~0 : 0,
                 ((k >> 1) & 1) ? ~0 : 0,
                 (k & 1) ? ~0 : 0);

      for (n = 0; n < rate_cats; ++n)
      {
        xmm0 = _mm256_maskload_pd(jmat,jmask);
        ymm0 = _mm256_maskload_pd(kmat,kmask);

        jmat += states;
        kmat += states;

        xmm1 = _mm256_maskload_pd(jmat,jmask);
        ymm1 = _mm256_maskload_pd(kmat,kmask);

        jmat += states;
        kmat += states;

        xmm2 = _mm256_maskload_pd(jmat,jmask);
        ymm2 = _mm256_maskload_pd(kmat,kmask);

        jmat += states;
        kmat += states;

        xmm3 = _mm256_maskload_pd(jmat,jmask);
        ymm3 = _mm256_maskload_pd(kmat,kmask);

        jmat += states;
        kmat += states;

        /* compute x */
        xmm4 = _mm256_unpackhi_pd(xmm0,xmm1);
        xmm5 = _mm256_unpacklo_pd(xmm0,xmm1);

        xmm6 = _mm256_unpackhi_pd(xmm2,xmm3);
        xmm7 = _mm256_unpacklo_pd(xmm2,xmm3);

        xmm0 = _mm256_add_pd(xmm4,xmm5);
        xmm1 = _mm256_add_pd(xmm6,xmm7);

        xmm2 = _mm256_permute2f128_pd(xmm0,xmm1, _MM_SHUFFLE(0,2,0,1));
        xmm3 = _mm256_blend_pd(xmm0,xmm1,12);
        xmm4 = _mm256_add_pd(xmm2,xmm3);

        /* compute y */
        ymm4 = _mm256_unpackhi_pd(ymm0,ymm1);
        ymm5 = _mm256_unpacklo_pd(ymm0,ymm1);

        ymm6 = _mm256_unpackhi_pd(ymm2,ymm3);
        ymm7 = _mm256_unpacklo_pd(ymm2,ymm3);

        ymm0 = _mm256_add_pd(ymm4,ymm5);
        ymm1 = _mm256_add_pd(ymm6,ymm7);

        ymm2 = _mm256_permute2f128_pd(ymm0,ymm1, _MM_SHUFFLE(0,2,0,1));
        ymm3 = _mm256_blend_pd(ymm0,ymm1,12);
        ymm4 = _mm256_add_pd(ymm2,ymm3);

        /* compute x*y */
        xmm0 = _mm256_mul_pd(xmm4,ymm4);

        _mm256_store_pd(lookup, xmm0);

        lookup += states;
      }
    }
  }
}

PLL_EXPORT void pll_core_update_partial_ii_4x4_avx(unsigned int sites,
                                                   unsigned int rate_cats,
                                                   double * parent_clv,
                                                   unsigned int * parent_scaler,
                                                   const double * left_clv,
                                                   const double * right_clv,
                                                   const double * left_matrix,
                                                   const double * right_matrix,
                                                   const unsigned int * left_scaler,
                                                   const unsigned int * right_scaler)
{
  unsigned int states = 4;
  unsigned int n,k,i;
  unsigned int scaling;

  const double * lmat;
  const double * rmat;

  __m256d ymm0,ymm1,ymm2,ymm3,ymm4,ymm5,ymm6,ymm7;
  __m256d xmm0,xmm1,xmm2,xmm3,xmm4,xmm5,xmm6,xmm7;

  unsigned int span = states * rate_cats;

  /* add up the scale vector of the two children if available */
  if (parent_scaler)
    fill_parent_scaler(sites, parent_scaler, left_scaler, right_scaler);

  for (n = 0; n < sites; ++n)
  {
    lmat = left_matrix;
    rmat = right_matrix;
    scaling = (parent_scaler) ? 1 : 0;

    for (k = 0; k < rate_cats; ++k)
    {
      /* compute vector of x */
      xmm4 = _mm256_load_pd(lmat);
      xmm5 = _mm256_load_pd(left_clv);
      xmm0 = _mm256_mul_pd(xmm4,xmm5);

      ymm4 = _mm256_load_pd(rmat);
      ymm5 = _mm256_load_pd(right_clv);
      ymm0 = _mm256_mul_pd(ymm4,ymm5);

      lmat += states;
      rmat += states;

      xmm4 = _mm256_load_pd(lmat);
      xmm1 = _mm256_mul_pd(xmm4,xmm5);

      ymm4 = _mm256_load_pd(rmat);
      ymm1 = _mm256_mul_pd(ymm4,ymm5);

      lmat += states;
      rmat += states;

      xmm4 = _mm256_load_pd(lmat);
      xmm2 = _mm256_mul_pd(xmm4,xmm5);

      ymm4 = _mm256_load_pd(rmat);
      ymm2 = _mm256_mul_pd(ymm4,ymm5);

      lmat += states;
      rmat += states;

      xmm4 = _mm256_load_pd(lmat);
      xmm3 = _mm256_mul_pd(xmm4,xmm5);

      ymm4 = _mm256_load_pd(rmat);
      ymm3 = _mm256_mul_pd(ymm4,ymm5);

      lmat += states;
      rmat += states;

      /* compute x */
      xmm4 = _mm256_unpackhi_pd(xmm0,xmm1);
      xmm5 = _mm256_unpacklo_pd(xmm0,xmm1);

      xmm6 = _mm256_unpackhi_pd(xmm2,xmm3);
      xmm7 = _mm256_unpacklo_pd(xmm2,xmm3);

      xmm0 = _mm256_add_pd(xmm4,xmm5);
      xmm1 = _mm256_add_pd(xmm6,xmm7);

      xmm2 = _mm256_permute2f128_pd(xmm0,xmm1, _MM_SHUFFLE(0,2,0,1));
      xmm3 = _mm256_blend_pd(xmm0,xmm1,12);
      xmm4 = _mm256_add_pd(xmm2,xmm3);

      /* compute y */
      ymm4 = _mm256_unpackhi_pd(ymm0,ymm1);
      ymm5 = _mm256_unpacklo_pd(ymm0,ymm1);

      ymm6 = _mm256_unpackhi_pd(ymm2,ymm3);
      ymm7 = _mm256_unpacklo_pd(ymm2,ymm3);

      ymm0 = _mm256_add_pd(ymm4,ymm5);
      ymm1 = _mm256_add_pd(ymm6,ymm7);

      ymm2 = _mm256_permute2f128_pd(ymm0,ymm1, _MM_SHUFFLE(0,2,0,1));
      ymm3 = _mm256_blend_pd(ymm0,ymm1,12);
      ymm4 = _mm256_add_pd(ymm2,ymm3);

      /* compute x*y */
      xmm0 = _mm256_mul_pd(xmm4,ymm4);

      _mm256_store_pd(parent_clv, xmm0);

      for (i = 0; i < states; ++i)
        scaling = scaling && (parent_clv[i] < PLL_SCALE_THRESHOLD);

      parent_clv += states;
      left_clv   += states;
      right_clv  += states;
    }
    /* if *all* entries of the site CLV were below the threshold then scale
       (all) entries by PLL_SCALE_FACTOR */
    if (scaling)
    {
      __m256d v_scale_factor = _mm256_set_pd(PLL_SCALE_FACTOR,
                                             PLL_SCALE_FACTOR,
                                             PLL_SCALE_FACTOR,
                                             PLL_SCALE_FACTOR);

      parent_clv -= span;
      for (i = 0; i < span; i += 4)
      {
        __m256d v_prod = _mm256_load_pd(parent_clv + i);
        v_prod = _mm256_mul_pd(v_prod,v_scale_factor);
        _mm256_store_pd(parent_clv + i, v_prod);
      }
      parent_clv += span;
      parent_scaler[n] += 1;
    }
  }
}

PLL_EXPORT void pll_core_update_partial_tt_avx(unsigned int states,
                                               unsigned int sites,
                                               unsigned int rate_cats,
                                               double * parent_clv,
                                               unsigned int * parent_scaler,
                                               const char * left_tipchars,
                                               const char * right_tipchars,
                                               const double * lookup,
                                               unsigned int tipstates_count)
{
  unsigned int j,k,n;
  unsigned int log2_rates = (unsigned int)ceil(log2(rate_cats));
  unsigned int log2_maxstates = (unsigned int)ceil(log2(tipstates_count));
  unsigned int states_padded = (states+3) & 0xFFFFFFFC;
  unsigned int log2_states_padded = (unsigned int)ceil(log2(states_padded));
  unsigned int span_padded = states_padded * rate_cats;
  const double * offset;

  if (states == 4)
  {
    pll_core_update_partial_tt_4x4_avx(sites,
                                       rate_cats,
                                       parent_clv,
                                       parent_scaler,
                                       left_tipchars,
                                       right_tipchars,
                                       lookup);
    return;
  }

  if (parent_scaler)
    memset(parent_scaler, 0, sizeof(unsigned int) * sites);

  for (n = 0; n < sites; ++n)
  {
    j = (unsigned int)(left_tipchars[n]);
    k = (unsigned int)(right_tipchars[n]);

    offset = lookup;
    offset += ((j << log2_maxstates) + k) << (log2_states_padded+log2_rates);

    memcpy(parent_clv, offset, span_padded*sizeof(double));

    parent_clv += span_padded;
  }
}

PLL_EXPORT void pll_core_update_partial_tt_4x4_avx(unsigned int sites,
                                                   unsigned int rate_cats,
                                                   double * parent_clv,
                                                   unsigned int * parent_scaler,
                                                   const char * left_tipchars,
                                                   const char * right_tipchars,
                                                   const double * lookup)
{
  unsigned int j,k,n;
  unsigned int log2_rates = (unsigned int)ceil(log2(rate_cats));
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
    offset += (( j << 4) + k) << (2+log2_rates);

    memcpy(parent_clv, offset, span*sizeof(double));

    parent_clv += span;
  }
}

/* not vectorized for a general number of states */
PLL_EXPORT void pll_core_update_partial_ti_avx(unsigned int states,
                                               unsigned int sites,                                   
                                               unsigned int rate_cats,
                                               double * parent_clv,
                                               unsigned int * parent_scaler,
                                               const char * left_tipchars,
                                               const double * right_clv,
                                               const double * left_matrix,
                                               const double * right_matrix,
                                               const unsigned int * right_scaler,
                                               const unsigned int * tipmap)
{
  int scaling;
  unsigned int i,j,k,n;

  const double * lmat;
  const double * rmat;

  if (states == 4)
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

  unsigned int states_padded = (states+3) & 0xFFFFFFFC;
  unsigned int span = states * rate_cats;
  unsigned int span_padded = states_padded * rate_cats;

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
        unsigned int lstate = tipmap[(int)left_tipchars[n]];
        for (j = 0; j < states; ++j)
        {
          if (lstate & 1)
            terma += lmat[j];

          termb += rmat[j] * right_clv[j];

          lstate >>= 1;
        }
        parent_clv[i] = terma*termb;
        lmat += states_padded;
        rmat += states_padded;

        scaling = scaling && (parent_clv[i] < PLL_SCALE_THRESHOLD);
      }
      parent_clv += states_padded;
      right_clv  += states_padded;
    }
    /* if *all* entries of the site CLV were below the threshold then scale
       (all) entries by PLL_SCALE_FACTOR */
    if (scaling)
    {
      parent_clv -= span_padded;
      for (i = 0; i < span; ++i)
        parent_clv[i] *= PLL_SCALE_FACTOR;
      parent_clv += span_padded;
      parent_scaler[n] += 1;
    }
  }
}

PLL_EXPORT void pll_core_update_partial_ti_4x4_avx(unsigned int sites,
                                                   unsigned int rate_cats,
                                                   double * parent_clv,
                                                   unsigned int * parent_scaler,
                                                   const char * left_tipchar,
                                                   const double * right_clv,
                                                   const double * left_matrix,
                                                   const double * right_matrix,
                                                   const unsigned int * right_scaler)
{
  unsigned int states = 4;
  unsigned int scaling;
  unsigned int i,k,n;
  
  const double * lmat;
  const double * rmat;

  unsigned int span = states * rate_cats;
  unsigned int lstate;

  __m256d ymm0,ymm1,ymm2,ymm3,ymm4,ymm5,ymm6,ymm7;
  __m256d xmm0,xmm1,xmm2,xmm3,xmm4,xmm5,xmm6,xmm7;
  __m256i mask;

  if (parent_scaler)
    fill_parent_scaler(sites, parent_scaler, NULL, right_scaler);

  for (n = 0; n < sites; ++n)
  {
    lmat = left_matrix;
    rmat = right_matrix;

    scaling = (parent_scaler) ? 1 : 0;
    
    lstate = left_tipchar[n];

    mask = _mm256_set_epi64x(
              ((lstate >> 3) & 1) ? ~0 : 0,
              ((lstate >> 2) & 1) ? ~0 : 0,
              ((lstate >> 1) & 1) ? ~0 : 0,
              ((lstate >> 0) & 1) ? ~0 : 0);

    for (k = 0; k < rate_cats; ++k)
    {
      xmm0 = _mm256_maskload_pd(lmat,mask);

      ymm4 = _mm256_load_pd(rmat);
      ymm5 = _mm256_load_pd(right_clv);
      ymm0 = _mm256_mul_pd(ymm4,ymm5);

      lmat += states;
      rmat += states;

      xmm1 = _mm256_maskload_pd(lmat,mask);

      ymm4 = _mm256_load_pd(rmat);
      ymm1 = _mm256_mul_pd(ymm4,ymm5);

      lmat += states;
      rmat += states;

      xmm2 = _mm256_maskload_pd(lmat,mask);

      ymm4 = _mm256_load_pd(rmat);
      ymm2 = _mm256_mul_pd(ymm4,ymm5);

      lmat += states;
      rmat += states;

      xmm3 = _mm256_maskload_pd(lmat,mask);

      ymm4 = _mm256_load_pd(rmat);
      ymm3 = _mm256_mul_pd(ymm4,ymm5);

      lmat += states;
      rmat += states;

      /* compute x */
      xmm4 = _mm256_unpackhi_pd(xmm0,xmm1);
      xmm5 = _mm256_unpacklo_pd(xmm0,xmm1);

      xmm6 = _mm256_unpackhi_pd(xmm2,xmm3);
      xmm7 = _mm256_unpacklo_pd(xmm2,xmm3);

      xmm0 = _mm256_add_pd(xmm4,xmm5);
      xmm1 = _mm256_add_pd(xmm6,xmm7);

      xmm2 = _mm256_permute2f128_pd(xmm0,xmm1, _MM_SHUFFLE(0,2,0,1));
      xmm3 = _mm256_blend_pd(xmm0,xmm1,12);
      xmm4 = _mm256_add_pd(xmm2,xmm3);

      /* compute y */
      ymm4 = _mm256_unpackhi_pd(ymm0,ymm1);
      ymm5 = _mm256_unpacklo_pd(ymm0,ymm1);

      ymm6 = _mm256_unpackhi_pd(ymm2,ymm3);
      ymm7 = _mm256_unpacklo_pd(ymm2,ymm3);

      ymm0 = _mm256_add_pd(ymm4,ymm5);
      ymm1 = _mm256_add_pd(ymm6,ymm7);

      ymm2 = _mm256_permute2f128_pd(ymm0,ymm1, _MM_SHUFFLE(0,2,0,1));
      ymm3 = _mm256_blend_pd(ymm0,ymm1,12);
      ymm4 = _mm256_add_pd(ymm2,ymm3);

      /* compute x*y */
      xmm0 = _mm256_mul_pd(xmm4,ymm4);

      _mm256_store_pd(parent_clv, xmm0);

      for (i = 0; i < states; ++i)
        scaling = scaling && (parent_clv[i] < PLL_SCALE_THRESHOLD);

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


PLL_EXPORT void pll_core_update_partial_ii_avx(unsigned int states,
                                               unsigned int sites,
                                               unsigned int rate_cats,
                                               double * parent_clv,
                                               unsigned int * parent_scaler,
                                               const double * left_clv,
                                               const double * right_clv,
                                               const double * left_matrix,
                                               const double * right_matrix,
                                               const unsigned int * left_scaler,
                                               const unsigned int * right_scaler)
{
  unsigned int i,j,k,n;
  unsigned int scaling;

  const double * lmat;
  const double * rmat;

  unsigned int states_padded = (states+3) & 0xFFFFFFFC;
  unsigned int span_padded = states_padded * rate_cats;

  /* dedicated functions for 4x4 matrices */
  if (states == 4)
  {
    pll_core_update_partial_ii_4x4_avx(sites,
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

  /* add up the scale vector of the two children if available */
  if (parent_scaler)
    fill_parent_scaler(sites, parent_scaler, left_scaler, right_scaler);

  size_t displacement = (states_padded - states) * (states_padded);

  /* compute CLV */
  for (n = 0; n < sites; ++n)
  {
    lmat = left_matrix;
    rmat = right_matrix;
    scaling = (parent_scaler) ? 1 : 0;
    for (k = 0; k < rate_cats; ++k)
    {
      for (i = 0; i < states_padded; i += 4)
      {

        __m256d v_terma0 = _mm256_setzero_pd();
        __m256d v_termb0 = _mm256_setzero_pd();
        __m256d v_terma1 = _mm256_setzero_pd();
        __m256d v_termb1 = _mm256_setzero_pd();
        __m256d v_terma2 = _mm256_setzero_pd();
        __m256d v_termb2 = _mm256_setzero_pd();
        __m256d v_terma3 = _mm256_setzero_pd();
        __m256d v_termb3 = _mm256_setzero_pd();

        __m256d v_mat;
        __m256d v_clv;

        for (j = 0; j < states_padded; j += 4)
        {
          v_mat    = _mm256_load_pd(lmat);
          v_clv    = _mm256_load_pd(left_clv+j);
          v_terma0 = _mm256_add_pd(v_terma0,
                                   _mm256_mul_pd(v_mat,v_clv));

          v_mat    = _mm256_load_pd(rmat);
          v_clv    = _mm256_load_pd(right_clv+j);
          v_termb0 = _mm256_add_pd(v_termb0,
                                   _mm256_mul_pd(v_mat,v_clv));

          lmat += 4;
          rmat += 4;
        }

        for (j = 0; j < states_padded; j += 4)
        {
          v_mat    = _mm256_load_pd(lmat);
          v_clv    = _mm256_load_pd(left_clv+j);
          v_terma1 = _mm256_add_pd(v_terma1,
                                   _mm256_mul_pd(v_mat,v_clv));

          v_mat    = _mm256_load_pd(rmat);
          v_clv    = _mm256_load_pd(right_clv+j);
          v_termb1 = _mm256_add_pd(v_termb1,
                                   _mm256_mul_pd(v_mat,v_clv));

          lmat += 4;
          rmat += 4;
        }

        for (j = 0; j < states_padded; j += 4)
        {
          v_mat    = _mm256_load_pd(lmat);
          v_clv    = _mm256_load_pd(left_clv+j);
          v_terma2 = _mm256_add_pd(v_terma2,
                                   _mm256_mul_pd(v_mat,v_clv));

          v_mat    = _mm256_load_pd(rmat);
          v_clv    = _mm256_load_pd(right_clv+j);
          v_termb2 = _mm256_add_pd(v_termb2,
                                   _mm256_mul_pd(v_mat,v_clv));

          lmat += 4;
          rmat += 4;
        }

        for (j = 0; j < states_padded; j += 4)
        {
          v_mat    = _mm256_load_pd(lmat);
          v_clv    = _mm256_load_pd(left_clv+j);
          v_terma3 = _mm256_add_pd(v_terma3,
                                   _mm256_mul_pd(v_mat,v_clv));

          v_mat    = _mm256_load_pd(rmat);
          v_clv    = _mm256_load_pd(right_clv+j);
          v_termb3 = _mm256_add_pd(v_termb3,
                                   _mm256_mul_pd(v_mat,v_clv));

          lmat += 4;
          rmat += 4;
        }

        __m256d xmm0 = _mm256_unpackhi_pd(v_terma0,v_terma1);
        __m256d xmm1 = _mm256_unpacklo_pd(v_terma0,v_terma1);

        __m256d xmm2 = _mm256_unpackhi_pd(v_terma2,v_terma3);
        __m256d xmm3 = _mm256_unpacklo_pd(v_terma2,v_terma3);

        xmm0 = _mm256_add_pd(xmm0,xmm1);
        xmm1 = _mm256_add_pd(xmm2,xmm3);

        xmm2 = _mm256_permute2f128_pd(xmm0,xmm1, _MM_SHUFFLE(0,2,0,1));
        
        xmm3 = _mm256_blend_pd(xmm0,xmm1,12);

        __m256d v_terma_sum = _mm256_add_pd(xmm2,xmm3);

        /* compute termb */

        xmm0 = _mm256_unpackhi_pd(v_termb0,v_termb1);
        xmm1 = _mm256_unpacklo_pd(v_termb0,v_termb1);

        xmm2 = _mm256_unpackhi_pd(v_termb2,v_termb3);
        xmm3 = _mm256_unpacklo_pd(v_termb2,v_termb3);

        xmm0 = _mm256_add_pd(xmm0,xmm1);
        xmm1 = _mm256_add_pd(xmm2,xmm3);

        xmm2 = _mm256_permute2f128_pd(xmm0,xmm1, _MM_SHUFFLE(0,2,0,1));
        
        xmm3 = _mm256_blend_pd(xmm0,xmm1,12);

        __m256d v_termb_sum = _mm256_add_pd(xmm2,xmm3);

        __m256d v_prod = _mm256_mul_pd(v_terma_sum,v_termb_sum);

        _mm256_store_pd(parent_clv+i, v_prod);

      }

      lmat -= displacement;
      rmat -= displacement;

      for (j = 0; j < states; ++j)
        scaling = scaling && (parent_clv[j] < PLL_SCALE_THRESHOLD);

      parent_clv += states_padded;
      left_clv   += states_padded;
      right_clv  += states_padded;
    }
    /* if *all* entries of the site CLV were below the threshold then scale
       (all) entries by PLL_SCALE_FACTOR */
    if (scaling)
    {
      __m256d v_scale_factor = _mm256_set_pd(PLL_SCALE_FACTOR,
                                             PLL_SCALE_FACTOR,
                                             PLL_SCALE_FACTOR,
                                             PLL_SCALE_FACTOR);

      parent_clv -= span_padded;
      for (i = 0; i < span_padded; i += 4)
      {
        __m256d v_prod = _mm256_load_pd(parent_clv + i);
        v_prod = _mm256_mul_pd(v_prod,v_scale_factor);
        _mm256_store_pd(parent_clv + i, v_prod);
      }
      parent_clv += span_padded;
      parent_scaler[n] += 1;
    }
  }
}

PLL_EXPORT int pll_core_update_sumtable_ii_4x4_avx(unsigned int sites,
                                                   unsigned int rate_cats,
                                                   const double * clvp,
                                                   const double * clvc,
                                                   double ** eigenvecs,
                                                   double ** inv_eigenvecs,
                                                   double ** freqs,
                                                   double *sumtable)
{
  unsigned int i, j, k, n;

  /* build sumtable */
  double * sum = sumtable;

  const double * t_clvp = clvp;
  const double * t_clvc = clvc;
  double * t_eigenvecs;
  double * t_freqs;

  unsigned int states = 4;

  /* transposed inv_eigenvecs */
  double * tt_inv_eigenvecs = (double *) pll_aligned_alloc (
      (states * states * rate_cats) * sizeof(double),
      PLL_ALIGNMENT_AVX);

  if (!tt_inv_eigenvecs)
  {
    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf (pll_errmsg, 200, "Cannot allocate memory for tt_inv_eigenvecs");
    return PLL_FAILURE;
  }

  for (i = 0; i < rate_cats; ++i)
  {
    for (j = 0; j < states; ++j)
      for (k = 0; k < states; ++k)
      {
        tt_inv_eigenvecs[i * states * states + j * states + k] =
            inv_eigenvecs[i][states * states + k * states + j];
      }
  }

  /* vectorized loop from update_sumtable() */
  for (n = 0; n < sites; n++)
  {
    for (i = 0; i < rate_cats; ++i)
    {
      t_eigenvecs = eigenvecs[i];
      t_freqs = freqs[i];

      const double * c_eigenvecs = t_eigenvecs;
      const double * ct_inv_eigenvecs = tt_inv_eigenvecs;

      __m256d v_lefterm[4], v_righterm[4];
      v_lefterm[0] = v_lefterm[1] = v_lefterm[2] = v_lefterm[3] = _mm256_setzero_pd ();
      v_righterm[0] = v_righterm[1] = v_righterm[2] = v_righterm[3] = _mm256_setzero_pd ();

      __m256d v_eigen;
      __m256d v_freqs;
      __m256d v_clvp, v_clvc;

      v_clvp = _mm256_load_pd (t_clvp);
      v_clvc = _mm256_load_pd (t_clvc);
      v_freqs = _mm256_load_pd (t_freqs);

      v_eigen = _mm256_load_pd (ct_inv_eigenvecs);
      v_lefterm[0] = _mm256_add_pd (
          v_lefterm[0],
          _mm256_mul_pd (v_freqs, _mm256_mul_pd (v_eigen, v_clvp)));
      v_eigen = _mm256_load_pd (c_eigenvecs);
      v_righterm[0] = _mm256_add_pd (v_righterm[0],
                                     _mm256_mul_pd (v_eigen, v_clvc));
      c_eigenvecs += 4;
      ct_inv_eigenvecs += 4;

      v_eigen = _mm256_load_pd (ct_inv_eigenvecs);
      v_lefterm[1] = _mm256_add_pd (
          v_lefterm[1],
          _mm256_mul_pd (v_freqs, _mm256_mul_pd (v_eigen, v_clvp)));
      v_eigen = _mm256_load_pd (c_eigenvecs);
      v_righterm[1] = _mm256_add_pd (v_righterm[1],
                                     _mm256_mul_pd (v_eigen, v_clvc));
      c_eigenvecs += 4;
      ct_inv_eigenvecs += 4;

      v_eigen = _mm256_load_pd (ct_inv_eigenvecs);
      v_lefterm[2] = _mm256_add_pd (
          v_lefterm[2],
          _mm256_mul_pd (v_freqs, _mm256_mul_pd (v_eigen, v_clvp)));
      v_eigen = _mm256_load_pd (c_eigenvecs);
      v_righterm[2] = _mm256_add_pd (v_righterm[2],
                                     _mm256_mul_pd (v_eigen, v_clvc));
      c_eigenvecs += 4;
      ct_inv_eigenvecs += 4;

      v_eigen = _mm256_load_pd (ct_inv_eigenvecs);
      v_lefterm[3] = _mm256_add_pd (
          v_lefterm[3],
          _mm256_mul_pd (v_freqs, _mm256_mul_pd (v_eigen, v_clvp)));
      v_eigen = _mm256_load_pd (c_eigenvecs);
      v_righterm[3] = _mm256_add_pd (v_righterm[3],
                                     _mm256_mul_pd (v_eigen, v_clvc));
      c_eigenvecs += 4;
      ct_inv_eigenvecs += 4;

      /* compute lefterm */
      __m256d xmm0 = _mm256_unpackhi_pd (v_lefterm[0], v_lefterm[1]);
      __m256d xmm1 = _mm256_unpacklo_pd (v_lefterm[0], v_lefterm[1]);
      __m256d xmm2 = _mm256_unpackhi_pd (v_lefterm[2], v_lefterm[3]);
      __m256d xmm3 = _mm256_unpacklo_pd (v_lefterm[2], v_lefterm[3]);
      xmm0 = _mm256_add_pd (xmm0, xmm1);
      xmm1 = _mm256_add_pd (xmm2, xmm3);
      xmm2 = _mm256_permute2f128_pd (xmm0, xmm1, _MM_SHUFFLE(0, 2, 0, 1));
      xmm3 = _mm256_blend_pd (xmm0, xmm1, 12);
      __m256d v_lefterm_sum = _mm256_add_pd (xmm2, xmm3);

      /* compute righterm */
      xmm0 = _mm256_unpackhi_pd (v_righterm[0], v_righterm[1]);
      xmm1 = _mm256_unpacklo_pd (v_righterm[0], v_righterm[1]);
      xmm2 = _mm256_unpackhi_pd (v_righterm[2], v_righterm[3]);
      xmm3 = _mm256_unpacklo_pd (v_righterm[2], v_righterm[3]);
      xmm0 = _mm256_add_pd (xmm0, xmm1);
      xmm1 = _mm256_add_pd (xmm2, xmm3);
      xmm2 = _mm256_permute2f128_pd (xmm0, xmm1, _MM_SHUFFLE(0, 2, 0, 1));
      xmm3 = _mm256_blend_pd (xmm0, xmm1, 12);
      __m256d v_righterm_sum = _mm256_add_pd (xmm2, xmm3);

      /* update sum */
      __m256d v_prod = _mm256_mul_pd (v_lefterm_sum, v_righterm_sum);
      _mm256_store_pd (sum, v_prod);

      t_clvc += states;
      t_clvp += states;
      sum += states;
    }
  }

  pll_aligned_free (tt_inv_eigenvecs);

  return PLL_SUCCESS;
}

PLL_EXPORT int pll_core_update_sumtable_ii_avx(unsigned int states,
                                               unsigned int sites,
                                               unsigned int rate_cats,
                                               const double * clvp,
                                               const double * clvc,
                                               double ** eigenvecs,
                                               double ** inv_eigenvecs,
                                               double ** freqs,
                                               double *sumtable)
{
  unsigned int i, j, k, n;

  /* build sumtable */
  double * sum = sumtable;

  const double * t_clvp = clvp;
  const double * t_clvc = clvc;
  double * t_eigenvecs;
  double * t_freqs;

  /* dedicated functions for 4x4 matrices */
  if (states == 4)
  {
    return pll_core_update_sumtable_ii_4x4_avx(sites,
                                       rate_cats,
                                       clvp,
                                       clvc,
                                       eigenvecs,
                                       inv_eigenvecs,
                                       freqs,
                                       sumtable);
  }

  unsigned int states_padded = (states+3) & 0xFFFFFFFC;

  /* padded eigenvecs */
  double * tt_eigenvecs = (double *) pll_aligned_alloc (
        (states_padded * states_padded * rate_cats) * sizeof(double),
        PLL_ALIGNMENT_AVX);

  if (!tt_eigenvecs)
  {
    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf (pll_errmsg, 200, "Cannot allocate memory for tt_eigenvecs");
    return PLL_FAILURE;
  }

  /* transposed padded inv_eigenvecs */
  double * tt_inv_eigenvecs = (double *) pll_aligned_alloc (
      (states_padded * states_padded * rate_cats) * sizeof(double),
      PLL_ALIGNMENT_AVX);

  if (!tt_inv_eigenvecs)
  {
    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf (pll_errmsg, 200, "Cannot allocate memory for tt_inv_eigenvecs");
    return PLL_FAILURE;
  }

  for (i = 0; i < rate_cats; ++i)
  {
    for (j = 0; j < states; ++j)
      for (k = 0; k < states; ++k)
      {
        tt_inv_eigenvecs[i * states_padded * states_padded + j * states_padded
            + k] = inv_eigenvecs[i][k * states + j];
        tt_eigenvecs[i * states_padded * states_padded + j * states_padded
            + k] = eigenvecs[i][j * states + k];
      }
  }

  /* vectorized loop from update_sumtable() */
  for (n = 0; n < sites; n++)
  {
    t_eigenvecs = tt_eigenvecs;
    for (i = 0; i < rate_cats; ++i)
    {
      t_freqs = freqs[i];

      const double * c_eigenvecs = t_eigenvecs;
      const double * ct_inv_eigenvecs = tt_inv_eigenvecs;

      for (j = 0; j < states_padded; j += 4)
      {
        __m256d v_lefterm0 = _mm256_setzero_pd ();
        __m256d v_righterm0 = _mm256_setzero_pd ();
        __m256d v_lefterm1 = _mm256_setzero_pd ();
        __m256d v_righterm1 = _mm256_setzero_pd ();
        __m256d v_lefterm2 = _mm256_setzero_pd ();
        __m256d v_righterm2 = _mm256_setzero_pd ();
        __m256d v_lefterm3 = _mm256_setzero_pd ();
        __m256d v_righterm3 = _mm256_setzero_pd ();

        __m256d v_eigen;
        __m256d v_freqs;
        __m256d v_clv;

        for (k = 0; k < states_padded; k += 4)
        {
          v_freqs = _mm256_load_pd (t_freqs + k);
          v_eigen = _mm256_load_pd (ct_inv_eigenvecs);
          v_clv = _mm256_load_pd (t_clvp + k);
          v_lefterm0 = _mm256_add_pd (
              v_lefterm0,
              _mm256_mul_pd (v_freqs, _mm256_mul_pd (v_eigen, v_clv)));

          v_clv = _mm256_load_pd (t_clvc + k);
          v_eigen = _mm256_load_pd (c_eigenvecs);
          v_righterm0 = _mm256_add_pd (v_righterm0,
                                       _mm256_mul_pd (v_eigen, v_clv));

          c_eigenvecs += 4;
          ct_inv_eigenvecs += 4;
        }

        for (k = 0; k < states_padded; k += 4)
        {
          v_freqs = _mm256_load_pd (t_freqs + k);
          v_eigen = _mm256_load_pd (ct_inv_eigenvecs);
          v_clv = _mm256_load_pd (t_clvp + k);
          v_lefterm1 = _mm256_add_pd (
              v_lefterm1,
              _mm256_mul_pd (v_freqs, _mm256_mul_pd (v_eigen, v_clv)));

          v_clv = _mm256_load_pd (t_clvc + k);
          v_eigen = _mm256_load_pd (c_eigenvecs);
          v_righterm1 = _mm256_add_pd (v_righterm1,
                                       _mm256_mul_pd (v_eigen, v_clv));

          c_eigenvecs += 4;
          ct_inv_eigenvecs += 4;
        }

        for (k = 0; k < states_padded; k += 4)
        {
          v_freqs = _mm256_load_pd (t_freqs + k);
          v_eigen = _mm256_load_pd (ct_inv_eigenvecs);
          v_clv = _mm256_load_pd (t_clvp + k);
          v_lefterm2 = _mm256_add_pd (
              v_lefterm2,
              _mm256_mul_pd (v_freqs, _mm256_mul_pd (v_eigen, v_clv)));

          v_clv = _mm256_load_pd (t_clvc + k);
          v_eigen = _mm256_load_pd (c_eigenvecs);
          v_righterm2 = _mm256_add_pd (v_righterm2,
                                       _mm256_mul_pd (v_eigen, v_clv));

          c_eigenvecs += 4;
          ct_inv_eigenvecs += 4;
        }

        for (k = 0; k < states_padded; k += 4)
        {
          v_freqs = _mm256_load_pd (t_freqs + k);
          v_eigen = _mm256_load_pd (ct_inv_eigenvecs);
          v_clv = _mm256_load_pd (t_clvp + k);
          v_lefterm3 = _mm256_add_pd (
              v_lefterm3,
              _mm256_mul_pd (v_freqs, _mm256_mul_pd (v_eigen, v_clv)));

          v_clv = _mm256_load_pd (t_clvc + k);
          v_eigen = _mm256_load_pd (c_eigenvecs);
          v_righterm3 = _mm256_add_pd (v_righterm3,
                                       _mm256_mul_pd (v_eigen, v_clv));

          c_eigenvecs += 4;
          ct_inv_eigenvecs += 4;
        }

        /* compute lefterm */
        __m256d xmm0 = _mm256_unpackhi_pd (v_lefterm0, v_lefterm1);
        __m256d xmm1 = _mm256_unpacklo_pd (v_lefterm0, v_lefterm1);
        __m256d xmm2 = _mm256_unpackhi_pd (v_lefterm2, v_lefterm3);
        __m256d xmm3 = _mm256_unpacklo_pd (v_lefterm2, v_lefterm3);
        xmm0 = _mm256_add_pd (xmm0, xmm1);
        xmm1 = _mm256_add_pd (xmm2, xmm3);
        xmm2 = _mm256_permute2f128_pd (xmm0, xmm1, _MM_SHUFFLE(0, 2, 0, 1));
        xmm3 = _mm256_blend_pd (xmm0, xmm1, 12);
        __m256d v_lefterm_sum = _mm256_add_pd (xmm2, xmm3);

        /* compute righterm */
        xmm0 = _mm256_unpackhi_pd (v_righterm0, v_righterm1);
        xmm1 = _mm256_unpacklo_pd (v_righterm0, v_righterm1);
        xmm2 = _mm256_unpackhi_pd (v_righterm2, v_righterm3);
        xmm3 = _mm256_unpacklo_pd (v_righterm2, v_righterm3);
        xmm0 = _mm256_add_pd (xmm0, xmm1);
        xmm1 = _mm256_add_pd (xmm2, xmm3);
        xmm2 = _mm256_permute2f128_pd (xmm0, xmm1, _MM_SHUFFLE(0, 2, 0, 1));
        xmm3 = _mm256_blend_pd (xmm0, xmm1, 12);
        __m256d v_righterm_sum = _mm256_add_pd (xmm2, xmm3);

        /* update sum */
        __m256d v_prod = _mm256_mul_pd (v_lefterm_sum, v_righterm_sum);
        _mm256_store_pd (sum + j, v_prod);
      }

      t_clvc += states_padded;
      t_clvp += states_padded;
      sum    += states_padded;

      t_eigenvecs += states_padded * states_padded;
    }
  }

  pll_aligned_free (tt_inv_eigenvecs);
  pll_aligned_free (tt_eigenvecs);

  return PLL_SUCCESS;
}

PLL_EXPORT double pll_core_likelihood_derivatives_avx(unsigned int states,
                                                  unsigned int sites,
                                                  unsigned int rate_cats,
                                                  double * rate_weights,
                                                  unsigned int * parent_scaler,
                                                  unsigned int * child_scaler,
                                                  int * invariant,
                                                  unsigned int * pattern_weights,
                                                  double branch_length,
                                                  double * prop_invar,
                                                  double ** freqs,
                                                  const double * rates,
                                                  double ** eigenvals,
                                                  double * sumtable,
                                                  double * d_f,
                                                  double * dd_f)
{
  unsigned int n, i, j;
  double site_lk[3];
  const double * sum;
  double logLK = 0.0;

  double deriv1, deriv2;

  const double * t_eigenvals;
  const double * t_freqs;
  double t_prop_invar;
  double t_branch_length;

  unsigned int states_padded = (states + 3) & 0xFFFFFFFC;

  *d_f = 0.0;
  *dd_f = 0.0;

  double *diagptable = (double *) calloc (rate_cats * states * 3,
                                          sizeof(double)), *diagp, ki;

  if (!diagptable)
  {
    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf (pll_errmsg, 200, "Cannot allocate memory for diagptable");
    return -INFINITY;
  }

  /* pre-compute the derivatives of the P matrix for all discrete GAMMA rates */
  diagp = diagptable;
  for (i = 0; i < rate_cats; ++i)
  {
    t_eigenvals = eigenvals[i];
    ki = rates[i];
    t_branch_length = branch_length / (1.0 - prop_invar[i]);
    for (j = 0; j < states; ++j)
    {
      diagp[0] = exp (t_eigenvals[j] * ki * t_branch_length);
      diagp[1] = t_eigenvals[j] * ki * diagp[0];
      diagp[2] = t_eigenvals[j] * ki * t_eigenvals[j] * ki * diagp[0];
      diagp += 3;
    }
  }

  sum = sumtable;
  for (n = 0; n < sites; ++n)
  {
    double inv_site_lk = 0.0;
    site_lk[0] = site_lk[1] = site_lk[2] = 0;
    diagp = diagptable;
    for (i = 0; i < rate_cats; ++i)
    {
      t_freqs = freqs[i];

      ki = rates[i];
      double cat_sitelk0 = 0, cat_sitelk1 = 0, cat_sitelk2 = 0;
      for (j = 0; j < states; ++j)
      {

        cat_sitelk0 += sum[j] * diagp[0];
        cat_sitelk1 += sum[j] * diagp[1];
        cat_sitelk2 += sum[j] * diagp[2];
        diagp += 3;
      }

      /* account for invariant sites */
      t_prop_invar = prop_invar[i];
      if (t_prop_invar > 0)
      {
        inv_site_lk =
            (invariant[n] == -1) ? 0 : t_freqs[invariant[n]] * t_prop_invar;
        cat_sitelk0 = cat_sitelk0 * (1. - t_prop_invar) + inv_site_lk;
        cat_sitelk1 = cat_sitelk1 * (1. - t_prop_invar);
        cat_sitelk2 = cat_sitelk2 * (1. - t_prop_invar);
      }

      site_lk[0] += cat_sitelk0 * rate_weights[i];
      site_lk[1] += cat_sitelk1 * rate_weights[i];
      site_lk[2] += cat_sitelk2 * rate_weights[i];

      sum += states_padded;
    }

    unsigned int scale_factors;
    scale_factors = (parent_scaler) ? parent_scaler[n] : 0;
    scale_factors += (child_scaler) ? child_scaler[n] : 0;

    //    if (site_lk[0] < PLL_SCALE_THRESHOLD_SQRT)
    //    {
    //      /* correct for underflow */
    //      scale_factors += 1;
    //      double lk_div = PLL_SCALE_THRESHOLD;
    //      site_lk[0] /= lk_div;
    //      site_lk[1] /= lk_div;
    //      site_lk[2] /= lk_div;
    //    }

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
  }

  free (diagptable);
  return logLK;
}
