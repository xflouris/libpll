/*
    Copyright (C) 2015 Tomas Flouri

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

PLL_EXPORT double pll_core_root_loglikelihood_avx(unsigned int states,
                                                  unsigned int sites,
                                                  unsigned int rate_cats,
                                                  const double * clv,
                                                  const unsigned int * scaler,
                                                  double * const * frequencies,
                                                  const double * rate_weights,
                                                  const unsigned int * pattern_weights,
                                                  const double * invar_proportion,
                                                  const int * invar_indices,
                                                  const unsigned int * freqs_indices,
                                                  double * persite_lnl)
{
  unsigned int i,j,k;
  double logl = 0;
  double prop_invar = 0;

  const double * freqs = NULL;

  double term, term_r;
  double inv_site_lk;

  unsigned int states_padded = (states+3) & 0xFFFFFFFC;

  __m256d xmm0, xmm1, xmm2, xmm3;

  for (i = 0; i < sites; ++i)
  {
    term = 0;
    for (j = 0; j < rate_cats; ++j)
    {
      freqs = frequencies[freqs_indices[j]];
      xmm3 = _mm256_setzero_pd();

      for (k = 0; k < states_padded; k += 4)
      {
        /* load frequencies for current rate matrix */
        xmm0 = _mm256_load_pd(freqs);

        /* load clv */
        xmm1 = _mm256_load_pd(clv);

        /* multiply with frequencies */
        xmm2 = _mm256_mul_pd(xmm0,xmm1);

        xmm3 = _mm256_add_pd(xmm3,xmm2);

        freqs += 4;
        clv += 4;
      }

      /* add up the elements of xmm2 */
      xmm1 = _mm256_hadd_pd(xmm3,xmm3);

      term_r = ((double *)&xmm1)[0] + ((double *)&xmm1)[2];

      /* account for invariant sites */
      prop_invar = invar_proportion ? invar_proportion[freqs_indices[j]] : 0;
      if (prop_invar > 0)
      {
        freqs = frequencies[freqs_indices[j]];
        inv_site_lk = (invar_indices[i] == -1) ?
                           0 : freqs[invar_indices[i]];
        term += rate_weights[j] * (term_r * (1 - prop_invar) +
                                   inv_site_lk*prop_invar);
      }
      else
      {
        term += term_r * rate_weights[j];
      }
    }

    /* compute site log-likelihood and scale if necessary */
    term = log(term);
    if (scaler && scaler[i])
      term += scaler[i] * log(PLL_SCALE_THRESHOLD);

    term *= pattern_weights[i];

    /* store per-site log-likelihood */
    if (persite_lnl)
      persite_lnl[i] = term;

    logl += term;
  }
  return logl;
}

PLL_EXPORT
double pll_core_root_loglikelihood_4x4_avx(unsigned int sites,
                                           unsigned int rate_cats,
                                           const double * clv,
                                           const unsigned int * scaler,
                                           double * const * frequencies,
                                           const double * rate_weights,
                                           const unsigned int * pattern_weights,
                                           const double * invar_proportion,
                                           const int * invar_indices,
                                           const unsigned int * freqs_indices,
                                           double * persite_lnl)
{
  unsigned int i,j;
  double logl = 0;
  double prop_invar = 0;

  const double * freqs = NULL;

  double term, term_r;
  double inv_site_lk;

  __m256d xmm0, xmm1, xmm2;

  for (i = 0; i < sites; ++i)
  {
    term = 0;
    for (j = 0; j < rate_cats; ++j)
    {
      freqs = frequencies[freqs_indices[j]];

      /* load frequencies for current rate matrix */
      xmm0 = _mm256_load_pd(freqs);

      /* load clv */
      xmm1 = _mm256_load_pd(clv);

      /* multiply with frequencies */
      xmm2 = _mm256_mul_pd(xmm0,xmm1);

      /* add up the elements of xmm2 */
      xmm1 = _mm256_hadd_pd(xmm2,xmm2);

      term_r = ((double *)&xmm1)[0] + ((double *)&xmm1)[2];

      /* account for invariant sites */
      prop_invar = invar_proportion ? invar_proportion[freqs_indices[j]] : 0;
      if (prop_invar > 0)
      {
        inv_site_lk = (invar_indices[i] == -1) ?
                           0 : freqs[invar_indices[i]];
        term += rate_weights[j] * (term_r * (1 - prop_invar) +
                                   inv_site_lk*prop_invar);
      }
      else
      {
        term += term_r * rate_weights[j];
      }

      clv += 4;
    }

    /* compute site log-likelihood and scale if necessary */
    term = log(term);
    if (scaler && scaler[i])
      term += scaler[i] * log(PLL_SCALE_THRESHOLD);

    term *= pattern_weights[i];

    /* store per-site log-likelihood */
    if (persite_lnl)
      persite_lnl[i] = term;

    logl += term;
  }
  return logl;
}

PLL_EXPORT
double pll_core_edge_loglikelihood_ti_4x4_avx(unsigned int sites,
                                              unsigned int rate_cats,
                                              const double * parent_clv,
                                              const unsigned int * parent_scaler,
                                              const unsigned char * tipchars,
                                              const double * pmatrix,
                                              double * const * frequencies,
                                              const double * rate_weights,
                                              const unsigned int * pattern_weights,
                                              const double * invar_proportion,
                                              const int * invar_indices,
                                              const unsigned int * freqs_indices,
                                              double * persite_lnl,
                                              unsigned int attrib)
{
  unsigned int n,i;
  double logl = 0;
  double prop_invar = 0;

  const double * clvp = parent_clv;
  const double * pmat;
  const double * freqs = NULL;

  double terma, terma_r;
  double site_lk, inv_site_lk;

  unsigned int cstate;
  unsigned int states_padded = 4;
  unsigned int span = rate_cats*states_padded;

  __m256d xmm0, xmm1, xmm2, xmm3, xmm4, xmm5, xmm6, xmm7;
  __m256i mask;

  unsigned int site_scalings;
  unsigned int * rate_scalings = NULL;
  int per_rate_scaling = parent_scaler && (attrib & PLL_ATTRIB_RATE_SCALERS) ? 1 : 0;

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

  /* precompute a lookup table of four values per entry (one for each state),
     for all 16 states (including ambiguities) and for each rate category. */
  double * lookup = pll_aligned_alloc(64*rate_cats*sizeof(double),
                                      PLL_ALIGNMENT_AVX);
  if (!lookup)
  {
    /* TODO: in the highly unlikely event that allocation fails, we should
       resort to a non-lookup-precomputation version of this function,
       available at commit e.g.  a4fc873fdc65741e402cdc1c59919375143d97d1 */
    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf(pll_errmsg, 200, "Cannot allocate space for precomputation.");
    return 1;
  }

  /* skip first entry of lookup table as it is never used */
  double * ptr = lookup + span;

  /* iterate all ambiguities skipping 0 */
  for (i = 1; i < 16; ++i)
  {
    pmat = pmatrix;

    /* mask the entries of pmatrix row to be loaded */
    mask = _mm256_set_epi64x(
              ((i >> 3) & 1) ? ~0 : 0,
              ((i >> 2) & 1) ? ~0 : 0,
              ((i >> 1) & 1) ? ~0 : 0,
              ((i >> 0) & 1) ? ~0 : 0);

    for (n = 0; n < rate_cats; ++n)
    {
      freqs = frequencies[freqs_indices[n]];

      /* masked load of the four matrix rows */
      xmm0 = _mm256_maskload_pd(pmat+0, mask);
      xmm1 = _mm256_maskload_pd(pmat+4, mask);
      xmm2 = _mm256_maskload_pd(pmat+8, mask);
      xmm3 = _mm256_maskload_pd(pmat+12,mask);

      /* create a vector containing the sums of xmm0, xmm1, xmm2, xmm3 */
      xmm4 = _mm256_unpackhi_pd(xmm0,xmm1);
      xmm5 = _mm256_unpacklo_pd(xmm0,xmm1);

      xmm6 = _mm256_unpackhi_pd(xmm2,xmm3);
      xmm7 = _mm256_unpacklo_pd(xmm2,xmm3);

      xmm0 = _mm256_add_pd(xmm4,xmm5);
      xmm1 = _mm256_add_pd(xmm6,xmm7);

      xmm2 = _mm256_permute2f128_pd(xmm0,xmm1, _MM_SHUFFLE(0,2,0,1));
      xmm3 = _mm256_blend_pd(xmm0,xmm1,12);
      xmm4 = _mm256_add_pd(xmm2,xmm3);

      /* load frequencies */
      xmm0 = _mm256_load_pd(freqs);

      /* multiply sums with frequencies */
      xmm1 = _mm256_mul_pd(xmm0,xmm4);

      /* store the result in the lookup table */
      _mm256_store_pd(ptr,xmm1);

      /* move pointers */
      ptr  += 4;
      pmat += 16;
    }
  }

  for (n = 0; n < sites; ++n)
  {
    pmat = pmatrix;
    terma = 0;

    cstate = tipchars[n];

    unsigned int coffset = cstate*span;

    if (per_rate_scaling)
    {
      const unsigned int * pscaler = parent_scaler + n*rate_cats;

      /* compute minimum per-rate scaler -> common per-site scaler */
      site_scalings = UINT_MAX;
      for (i = 0; i < rate_cats; ++i)
      {
        if (pscaler[i] < site_scalings)
          site_scalings = pscaler[i];
      }

      /* compute relative capped per-rate scalers */
      for (i = 0; i < rate_cats; ++i)
      {
        rate_scalings[i] = PLL_MIN(pscaler[i] - site_scalings,
                                   PLL_SCALE_RATE_MAXDIFF);
      }
    }
    else
    {
      /* count number of scaling factors to account for */
      site_scalings = (parent_scaler) ? parent_scaler[n] : 0;
    }

    for (i = 0; i < rate_cats; ++i)
    {
      /* load precomputed lookup table into xmm3 */
      xmm1 = _mm256_load_pd(lookup+coffset);

      /* multiply with clvp */
      xmm2 = _mm256_load_pd(clvp);
      xmm0 = _mm256_mul_pd(xmm1,xmm2);

      /* add up the elements of xmm0 */
      xmm1 = _mm256_hadd_pd(xmm0,xmm0);
      terma_r = ((double *)&xmm1)[0] + ((double *)&xmm1)[2];

      /* apply per-rate scalers, if necessary */
      if (rate_scalings && rate_scalings[i] > 0)
      {
        terma_r *= scale_minlh[rate_scalings[i]-1];
      }

      if (terma_r > 0.)
      {
        /* account for invariant sites */
        prop_invar = invar_proportion ? invar_proportion[freqs_indices[i]] : 0;
        if (prop_invar > 0)
        {
          inv_site_lk = (invar_indices[n] == -1) ?
                            0 : freqs[invar_indices[n]];
          terma += rate_weights[i] * (terma_r * (1 - prop_invar) +
                   inv_site_lk * prop_invar);
        }
        else
        {
          terma += terma_r * rate_weights[i];
        }
      }

      clvp += states_padded;
      coffset += 4;
    }

    /* compute site log-likelihood and scale if necessary */
    site_lk = log(terma);

    /* apply per-site scaler, if necessary */
    if (site_scalings)
      site_lk += site_scalings * log(PLL_SCALE_THRESHOLD);

    site_lk *= pattern_weights[n];

    /* store per-site log-likelihood */
    if (persite_lnl)
      persite_lnl[n] = site_lk;

    logl += site_lk;
  }

  pll_aligned_free(lookup);
  if (rate_scalings)
    free(rate_scalings);

  return logl;
}

PLL_EXPORT
double pll_core_edge_loglikelihood_ti_20x20_avx(unsigned int sites,
                                                unsigned int rate_cats,
                                                const double * parent_clv,
                                                const unsigned int * parent_scaler,
                                                const unsigned char * tipchars,
                                                const unsigned int * tipmap,
                                                unsigned int tipmap_size,
                                                const double * pmatrix,
                                                double * const * frequencies,
                                                const double * rate_weights,
                                                const unsigned int * pattern_weights,
                                                const double * invar_proportion,
                                                const int * invar_indices,
                                                const unsigned int * freqs_indices,
                                                double * persite_lnl,
                                                unsigned int attrib)
{
  unsigned int n,i,j,m;
  double logl = 0;
  double prop_invar = 0;

  const double * clvp = parent_clv;
  const double * pmat;
  const double * freqs = NULL;

  double terma, terma_r;
  double site_lk, inv_site_lk;

  unsigned int cstate;
  unsigned int states = 20;
  unsigned int states_padded = states;

  __m256d xmm0, xmm1, xmm2;

  size_t displacement = (states_padded - states) * (states_padded);

  unsigned int span = states_padded * rate_cats;
  unsigned int maxstates = tipmap_size;

  /* scaling stuff */
  unsigned int site_scalings;
  unsigned int * rate_scalings = NULL;
  int per_rate_scaling = (attrib & PLL_ATTRIB_RATE_SCALERS) ? 1 : 0;

  /* powers of scale threshold for undoing the scaling */
  double scale_minlh[PLL_SCALE_RATE_MAXDIFF];
  if (per_rate_scaling)
  {
    rate_scalings = (unsigned int*) calloc(rate_cats, sizeof(unsigned int));

    if (!rate_scalings)
    {
      pll_errno = PLL_ERROR_MEM_ALLOC;
      snprintf(pll_errmsg, 200, "Cannot allocate space for rate scalers.");
      return -INFINITY;
    }

    double scale_factor = 1.0;
    for (i = 0; i < PLL_SCALE_RATE_MAXDIFF; ++i)
    {
      scale_factor *= PLL_SCALE_THRESHOLD;
      scale_minlh[i] = scale_factor;
    }
  }

  /* precompute a lookup table of four values per entry (one for each state),
     for all 16 states (including ambiguities) and for each rate category. */
  double * lookup = pll_aligned_alloc(maxstates*span*sizeof(double),
                                      PLL_ALIGNMENT_AVX);
  if (!lookup)
  {
    /* TODO: in the highly unlikely event that allocation fails, we should
       resort to a non-lookup-precomputation version of this function,
       available at commit e.g.  a4fc873fdc65741e402cdc1c59919375143d97d1 */
    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf(pll_errmsg, 200, "Cannot allocate space for precomputation.");
    return 0.;
  }

  double * ptr = lookup;

  /* precompute left-side values and store them in lookup table */
  for (j = 0; j < maxstates; ++j)
  {
    pmat = pmatrix;

    unsigned int state = tipmap[j];

    int ss = __builtin_popcount(state) == 1 ? __builtin_ctz(state) : -1;

    for (n = 0; n < rate_cats; ++n)
    {
      freqs = frequencies[freqs_indices[n]];

      for (i = 0; i < states; ++i)
      {
        double terml;
        if (ss != -1)
        {
          /* special case for non-ambiguous states */
          terml = pmat[ss];
        }
        else
        {
          terml = 0;
          for (m = 0; m < states; ++m)
          {
            if ((state>>m) & 1)
            {
              terml += pmat[m];
            }
          }
        }

        pmat += states;

        ptr[i] = terml * freqs[i];
      }

      ptr += states;
    }
  }

  for (n = 0; n < sites; ++n)
  {
    terma = 0;

    cstate = (unsigned int) tipchars[n];
    unsigned int loffset = cstate*span;

    if (per_rate_scaling)
    {
      /* compute minimum per-rate scaler -> common per-site scaler */
      site_scalings = UINT_MAX;
      for (i = 0; i < rate_cats; ++i)
      {
        rate_scalings[i] = (parent_scaler) ? parent_scaler[n*rate_cats+i] : 0;
        if (rate_scalings[i] < site_scalings)
          site_scalings = rate_scalings[i];
      }

      /* compute relative capped per-rate scalers */
      for (i = 0; i < rate_cats; ++i)
      {
        rate_scalings[i] = PLL_MIN(rate_scalings[i] - site_scalings,
                                   PLL_SCALE_RATE_MAXDIFF);
      }
    }
    else
    {
      /* count number of scaling factors to account for */
      site_scalings =  (parent_scaler) ? parent_scaler[n] : 0;
    }

    for (i = 0; i < rate_cats; ++i)
    {
      xmm1 = _mm256_setzero_pd();

      /* iterate over quadruples of rows */
      for (j = 0; j < states_padded; j += 4)
      {
        /* load value from lookup table */
        xmm2 = _mm256_load_pd(lookup+loffset);

        /* multiply with clvp */
        xmm0 = _mm256_load_pd(clvp);
        xmm0 = _mm256_mul_pd(xmm2,xmm0);

        /* add up the elements of xmm1 to the accumulator */
        xmm1 = _mm256_add_pd(xmm1,xmm0);

        clvp += 4;
        loffset += 4;
      }

      xmm0 = _mm256_hadd_pd(xmm1,xmm1);
      terma_r = ((double *)&xmm0)[0] + ((double *)&xmm0)[2];

      /* apply per-rate scalers, if necessary */
      if (rate_scalings && rate_scalings[i] > 0)
      {
        terma_r *= scale_minlh[rate_scalings[i]-1];
      }

      /* account for invariant sites */
      prop_invar = invar_proportion ? invar_proportion[freqs_indices[i]] : 0;
      if (prop_invar > 0)
      {
        freqs = frequencies[freqs_indices[i]];
        inv_site_lk = (invar_indices[n] == -1) ?
                          0 : freqs[invar_indices[n]];
        terma += rate_weights[i] * (terma_r * (1 - prop_invar) +
                 inv_site_lk * prop_invar);
      }
      else
      {
        terma += terma_r * rate_weights[i];
      }

      pmat -= displacement;
    }

    /* compute site log-likelihood and scale if necessary */
    site_lk = log(terma);
    if (site_scalings)
      site_lk += site_scalings * log(PLL_SCALE_THRESHOLD);

    site_lk *= pattern_weights[n];

    /* store per-site log-likelihood */
    if (persite_lnl)
      persite_lnl[n] = site_lk;

    logl += site_lk;
  }

  pll_aligned_free(lookup);
  if (rate_scalings)
    free(rate_scalings);

  return logl;
}

PLL_EXPORT
double pll_core_edge_loglikelihood_ti_avx(unsigned int states,
                                          unsigned int sites,
                                          unsigned int rate_cats,
                                          const double * parent_clv,
                                          const unsigned int * parent_scaler,
                                          const unsigned char * tipchars,
                                          const unsigned int * tipmap,
                                          const double * pmatrix,
                                          double * const * frequencies,
                                          const double * rate_weights,
                                          const unsigned int * pattern_weights,
                                          const double * invar_proportion,
                                          const int * invar_indices,
                                          const unsigned int * freqs_indices,
                                          double * persite_lnl,
                                          unsigned int attrib)
{
  unsigned int n,i,j,k;
  double logl = 0;
  double prop_invar = 0;

  const double * clvp = parent_clv;
  const double * pmat;
  const double * freqs = NULL;

  double terma, terma_r;
  double site_lk, inv_site_lk;

  unsigned int cstate;
  unsigned int states_padded = (states+3) & 0xFFFFFFFC;

  __m256d xmm0, xmm1, xmm2, xmm3, xmm4, xmm5, xmm6, xmm7;
  __m256i mask;

  size_t displacement = (states_padded - states) * (states_padded);

  /* scaling stuff */
  unsigned int site_scalings;
  unsigned int * rate_scalings = NULL;
  int per_rate_scaling = (attrib & PLL_ATTRIB_RATE_SCALERS) ? 1 : 0;

  /* powers of scale threshold for undoing the scaling */
  double scale_minlh[PLL_SCALE_RATE_MAXDIFF];
  if (per_rate_scaling)
  {
    rate_scalings = (unsigned int*) calloc(rate_cats, sizeof(unsigned int));

    if (!rate_scalings)
    {
      pll_errno = PLL_ERROR_MEM_ALLOC;
      snprintf(pll_errmsg, 200, "Cannot allocate space for rate scalers.");
      return -INFINITY;
    }

    double scale_factor = 1.0;
    for (i = 0; i < PLL_SCALE_RATE_MAXDIFF; ++i)
    {
      scale_factor *= PLL_SCALE_THRESHOLD;
      scale_minlh[i] = scale_factor;
    }
  }

  for (n = 0; n < sites; ++n)
  {
    pmat = pmatrix;
    terma = 0;

    cstate = tipmap[tipchars[n]];

    if (per_rate_scaling)
    {
      /* compute minimum per-rate scaler -> common per-site scaler */
      site_scalings = UINT_MAX;
      for (i = 0; i < rate_cats; ++i)
      {
        rate_scalings[i] = (parent_scaler) ? parent_scaler[n*rate_cats+i] : 0;
        if (rate_scalings[i] < site_scalings)
          site_scalings = rate_scalings[i];
      }

      /* compute relative capped per-rate scalers */
      for (i = 0; i < rate_cats; ++i)
      {
        rate_scalings[i] = PLL_MIN(rate_scalings[i] - site_scalings,
                                   PLL_SCALE_RATE_MAXDIFF);
      }
    }
    else
    {
      /* count number of scaling factors to account for */
      site_scalings =  (parent_scaler) ? parent_scaler[n] : 0;
    }

    for (i = 0; i < rate_cats; ++i)
    {
      freqs = frequencies[freqs_indices[i]];
      terma_r = 0;

      /* iterate over quadruples of rows */
      for (j = 0; j < states_padded; j += 4)
      {
        xmm0 = _mm256_setzero_pd();
        xmm1 = _mm256_setzero_pd();
        xmm2 = _mm256_setzero_pd();
        xmm3 = _mm256_setzero_pd();

        /* point to the four rows */
        const double * row0 = pmat;
        const double * row1 = row0 + states_padded;
        const double * row2 = row1 + states_padded;
        const double * row3 = row2 + states_padded;

        /* set position of least significant bit in character state */
        register int lsb = 0;

        /* iterate quadruples of columns */
        for (k = 0; k < states_padded; k += 4)
        {

          /* set mask */
          mask = _mm256_set_epi64x(
                    ((cstate >> (lsb+3)) & 1) ? ~0 : 0,
                    ((cstate >> (lsb+2)) & 1) ? ~0 : 0,
                    ((cstate >> (lsb+1)) & 1) ? ~0 : 0,
                    ((cstate >> (lsb+0)) & 1) ? ~0 : 0);

          lsb += 4;

          /* row 0 */
          xmm4 = _mm256_maskload_pd(row0,mask);
          xmm0 = _mm256_add_pd(xmm0,xmm4);
          row0 += 4;

          /* row 1 */
          xmm4 = _mm256_maskload_pd(row1,mask);
          xmm1 = _mm256_add_pd(xmm1,xmm4);
          row1 += 4;

          /* row 2 */
          xmm4 = _mm256_maskload_pd(row2,mask);
          xmm2 = _mm256_add_pd(xmm2,xmm4);
          row2 += 4;

          /* row 3 */
          xmm4 = _mm256_maskload_pd(row3,mask);
          xmm3 = _mm256_add_pd(xmm3,xmm4);
          row3 += 4;
        }

        /* point pmatrix to the next four rows */
        pmat = row3;

        /* create a vector containing the sums of xmm0, xmm1, xmm2, xmm3 */
        xmm4 = _mm256_unpackhi_pd(xmm0,xmm1);
        xmm5 = _mm256_unpacklo_pd(xmm0,xmm1);

        xmm6 = _mm256_unpackhi_pd(xmm2,xmm3);
        xmm7 = _mm256_unpacklo_pd(xmm2,xmm3);

        xmm0 = _mm256_add_pd(xmm4,xmm5);
        xmm1 = _mm256_add_pd(xmm6,xmm7);

        xmm2 = _mm256_permute2f128_pd(xmm0,xmm1, _MM_SHUFFLE(0,2,0,1));
        xmm3 = _mm256_blend_pd(xmm0,xmm1,12);
        xmm0 = _mm256_add_pd(xmm2,xmm3);

        /* multiply with frequencies */
        xmm1 = _mm256_load_pd(freqs);
        xmm2 = _mm256_mul_pd(xmm0,xmm1);

        /* multiply with clvp */
        xmm0 = _mm256_load_pd(clvp);
        xmm1 = _mm256_mul_pd(xmm2,xmm0);

        /* add up the elements of xmm1 */
        xmm0 = _mm256_hadd_pd(xmm1,xmm1);
        terma_r += ((double *)&xmm0)[0] + ((double *)&xmm0)[2];

        freqs += 4;
        clvp += 4;
      }

      /* apply per-rate scalers, if necessary */
      if (rate_scalings && rate_scalings[i] > 0)
      {
        terma_r *= scale_minlh[rate_scalings[i]-1];
      }

      /* account for invariant sites */
      prop_invar = invar_proportion ? invar_proportion[freqs_indices[i]] : 0;
      if (prop_invar > 0)
      {
        freqs = frequencies[freqs_indices[i]];
        inv_site_lk = (invar_indices[n] == -1) ?
                          0 : freqs[invar_indices[n]];
        terma += rate_weights[i] * (terma_r * (1 - prop_invar) +
                 inv_site_lk * prop_invar);
      }
      else
      {
        terma += terma_r * rate_weights[i];
      }

      pmat -= displacement;
    }

    /* compute site log-likelihood and scale if necessary */
    site_lk = log(terma);
    if (site_scalings)
      site_lk += site_scalings * log(PLL_SCALE_THRESHOLD);

    site_lk *= pattern_weights[n];

    /* store per-site log-likelihood */
    if (persite_lnl)
      persite_lnl[n] = site_lk;

    logl += site_lk;
  }

  if (rate_scalings)
    free(rate_scalings);

  return logl;
}

PLL_EXPORT
double pll_core_edge_loglikelihood_ii_avx(unsigned int states,
                                          unsigned int sites,
                                          unsigned int rate_cats,
                                          const double * parent_clv,
                                          const unsigned int * parent_scaler,
                                          const double * child_clv,
                                          const unsigned int * child_scaler,
                                          const double * pmatrix,
                                          double * const * frequencies,
                                          const double * rate_weights,
                                          const unsigned int * pattern_weights,
                                          const double * invar_proportion,
                                          const int * invar_indices,
                                          const unsigned int * freqs_indices,
                                          double * persite_lnl,
                                          unsigned int attrib)
{
  unsigned int n,i,j,k;
  double logl = 0;
  double prop_invar = 0;

  const double * clvp = parent_clv;
  const double * clvc = child_clv;
  const double * pmat;
  const double * freqs = NULL;

  double terma, terma_r;
  double site_lk, inv_site_lk;

  unsigned int states_padded = (states+3) & 0xFFFFFFFC;

  __m256d xmm0, xmm1, xmm2, xmm3, xmm4, xmm5, xmm6, xmm7;

  size_t displacement = (states_padded - states) * (states_padded);

  /* scaling stuff */
  unsigned int site_scalings;
  unsigned int * rate_scalings = NULL;
  int per_rate_scaling = (attrib & PLL_ATTRIB_RATE_SCALERS) ? 1 : 0;

  /* powers of scale threshold for undoing the scaling */
  double scale_minlh[PLL_SCALE_RATE_MAXDIFF];
  if (per_rate_scaling)
  {
    rate_scalings = (unsigned int*) calloc(rate_cats, sizeof(unsigned int));

    if (!rate_scalings)
    {
      pll_errno = PLL_ERROR_MEM_ALLOC;
      snprintf(pll_errmsg, 200, "Cannot allocate space for precomputation.");
      return -INFINITY;
    }

    double scale_factor = 1.0;
    for (i = 0; i < PLL_SCALE_RATE_MAXDIFF; ++i)
    {
      scale_factor *= PLL_SCALE_THRESHOLD;
      scale_minlh[i] = scale_factor;
    }
  }

  for (n = 0; n < sites; ++n)
  {
    pmat = pmatrix;
    terma = 0;

    if (per_rate_scaling)
    {
      /* compute minimum per-rate scaler -> common per-site scaler */
      site_scalings = UINT_MAX;
      for (i = 0; i < rate_cats; ++i)
      {
        rate_scalings[i] = (parent_scaler) ? parent_scaler[n*rate_cats+i] : 0;
        rate_scalings[i] += (child_scaler) ? child_scaler[n*rate_cats+i] : 0;
        if (rate_scalings[i] < site_scalings)
          site_scalings = rate_scalings[i];
      }

      /* compute relative capped per-rate scalers */
      for (i = 0; i < rate_cats; ++i)
      {
        rate_scalings[i] = PLL_MIN(rate_scalings[i] - site_scalings,
                                   PLL_SCALE_RATE_MAXDIFF);
      }
    }
    else
    {
      /* count number of scaling factors to account for */
      site_scalings =  (parent_scaler) ? parent_scaler[n] : 0;
      site_scalings += (child_scaler) ? child_scaler[n] : 0;
    }

    for (i = 0; i < rate_cats; ++i)
    {
      freqs = frequencies[freqs_indices[i]];
      terma_r = 0;

      /* iterate over quadruples of rows */
      for (j = 0; j < states_padded; j += 4)
      {
        xmm0 = _mm256_setzero_pd();
        xmm1 = _mm256_setzero_pd();
        xmm2 = _mm256_setzero_pd();
        xmm3 = _mm256_setzero_pd();

        /* point to the four rows */
        const double * row0 = pmat;
        const double * row1 = row0 + states_padded;
        const double * row2 = row1 + states_padded;
        const double * row3 = row2 + states_padded;

        /* iterate quadruples of columns */
        for (k = 0; k < states_padded; k += 4)
        {
          xmm5 = _mm256_load_pd(clvc+k);

          /* row 0 */
          xmm4 = _mm256_load_pd(row0);
          xmm6 = _mm256_mul_pd(xmm4,xmm5);
          xmm0 = _mm256_add_pd(xmm0,xmm6);
          row0 += 4;

          /* row 1 */
          xmm4 = _mm256_load_pd(row1);
          xmm6 = _mm256_mul_pd(xmm4,xmm5);
          xmm1 = _mm256_add_pd(xmm1,xmm6);
          row1 += 4;

          /* row 2 */
          xmm4 = _mm256_load_pd(row2);
          xmm6 = _mm256_mul_pd(xmm4,xmm5);
          xmm2 = _mm256_add_pd(xmm2,xmm6);
          row2 += 4;

          /* row 3 */
          xmm4 = _mm256_load_pd(row3);
          xmm6 = _mm256_mul_pd(xmm4,xmm5);
          xmm3 = _mm256_add_pd(xmm3,xmm6);
          row3 += 4;
        }

        /* point pmatrix to the next four rows */
        pmat = row3;

        /* create a vector containing the sums of xmm0, xmm1, xmm2, xmm3 */
        xmm4 = _mm256_unpackhi_pd(xmm0,xmm1);
        xmm5 = _mm256_unpacklo_pd(xmm0,xmm1);

        xmm6 = _mm256_unpackhi_pd(xmm2,xmm3);
        xmm7 = _mm256_unpacklo_pd(xmm2,xmm3);

        xmm0 = _mm256_add_pd(xmm4,xmm5);
        xmm1 = _mm256_add_pd(xmm6,xmm7);

        xmm2 = _mm256_permute2f128_pd(xmm0,xmm1, _MM_SHUFFLE(0,2,0,1));
        xmm3 = _mm256_blend_pd(xmm0,xmm1,12);
        xmm0 = _mm256_add_pd(xmm2,xmm3);

        /* multiply with frequencies */
        xmm1 = _mm256_load_pd(freqs);
        xmm2 = _mm256_mul_pd(xmm0,xmm1);

        /* multiply with clvp */
        xmm0 = _mm256_load_pd(clvp);
        xmm1 = _mm256_mul_pd(xmm2,xmm0);

        /* add up the elements of xmm1 */
        xmm0 = _mm256_hadd_pd(xmm1,xmm1);
        terma_r += ((double *)&xmm0)[0] + ((double *)&xmm0)[2];

        freqs += 4;
        clvp += 4;
      }

      /* apply per-rate scalers, if necessary */
      if (rate_scalings && rate_scalings[i] > 0)
      {
        terma_r *= scale_minlh[rate_scalings[i]-1];
      }

      /* account for invariant sites */
      prop_invar = invar_proportion ? invar_proportion[freqs_indices[i]] : 0;
      if (prop_invar > 0)
      {
        freqs = frequencies[freqs_indices[i]];
        inv_site_lk = (invar_indices[n] == -1) ?
                          0 : freqs[invar_indices[n]];
        terma += rate_weights[i] * (terma_r * (1 - prop_invar) +
                 inv_site_lk * prop_invar);
      }
      else
      {
        terma += terma_r * rate_weights[i];
      }

      clvc += states_padded;
      pmat -= displacement;
    }

    /* compute site log-likelihood and scale if necessary */
    site_lk = log(terma);
    if (site_scalings)
      site_lk += site_scalings * log(PLL_SCALE_THRESHOLD);

    site_lk *= pattern_weights[n];

    /* store per-site log-likelihood */
    if (persite_lnl)
      persite_lnl[n] = site_lk;

    logl += site_lk;
  }

  if (rate_scalings)
    free(rate_scalings);

  return logl;
}

PLL_EXPORT
double pll_core_edge_loglikelihood_ii_4x4_avx(unsigned int sites,
                                              unsigned int rate_cats,
                                              const double * parent_clv,
                                              const unsigned int * parent_scaler,
                                              const double * child_clv,
                                              const unsigned int * child_scaler,
                                              const double * pmatrix,
                                              double * const * frequencies,
                                              const double * rate_weights,
                                              const unsigned int * pattern_weights,
                                              const double * invar_proportion,
                                              const int * invar_indices,
                                              const unsigned int * freqs_indices,
                                              double * persite_lnl,
                                              unsigned int attrib)
{
  unsigned int n,i;
  double logl = 0;
  double prop_invar = 0;

  const double * clvp = parent_clv;
  const double * clvc = child_clv;
  const double * pmat;
  const double * freqs = NULL;

  double terma, terma_r;
  double site_lk, inv_site_lk;

  unsigned int states = 4;
  unsigned int states_padded = 4;

  __m256d xmm0, xmm1, xmm2, xmm3, xmm4, xmm5, xmm6;

  unsigned int site_scalings;
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

  for (n = 0; n < sites; ++n)
  {
    pmat = pmatrix;
    terma = 0;

    if (per_rate_scaling)
    {
      /* compute minimum per-rate scaler -> common per-site scaler */
      site_scalings = UINT_MAX;
      for (i = 0; i < rate_cats; ++i)
      {
        rate_scalings[i] = (parent_scaler) ? parent_scaler[n*rate_cats+i] : 0;
        rate_scalings[i] += (child_scaler) ? child_scaler[n*rate_cats+i] : 0;
        if (rate_scalings[i] < site_scalings)
          site_scalings = rate_scalings[i];
      }

      /* compute relative capped per-rate scalers */
      for (i = 0; i < rate_cats; ++i)
      {
        rate_scalings[i] = PLL_MIN(rate_scalings[i] - site_scalings,
                                   PLL_SCALE_RATE_MAXDIFF);
      }
    }
    else
    {
      /* count number of scaling factors to account for */
      site_scalings =  (parent_scaler) ? parent_scaler[n] : 0;
      site_scalings += (child_scaler) ? child_scaler[n] : 0;
    }

    for (i = 0; i < rate_cats; ++i)
    {
      freqs = frequencies[freqs_indices[i]];

      /* load frequencies for current rate matrix */
      xmm0 = _mm256_load_pd(freqs);

      /* load clvc */
      xmm1 = _mm256_load_pd(clvc);

      /* load pmatrix row 1 and multiply with clvc */
      xmm2 = _mm256_load_pd(pmat);
      xmm3 = _mm256_mul_pd(xmm1,xmm2);

      /* load pmatrix row 2 and multiply with clvc */
      pmat += states;
      xmm2 = _mm256_load_pd(pmat);
      xmm4 = _mm256_mul_pd(xmm1,xmm2);

      /* load pmatrix row 3 and multiply with clvc */
      pmat += states;
      xmm2 = _mm256_load_pd(pmat);
      xmm5 = _mm256_mul_pd(xmm1,xmm2);

      /* load pmatrix row 4 and multiply with clvc */
      pmat += states;
      xmm2 = _mm256_load_pd(pmat);
      xmm6 = _mm256_mul_pd(xmm1,xmm2);

      /* point to the pmatrix for the next rate category */
      pmat += states;

      /* create a vector containing the sums of xmm3, xmm4, xmm5, xmm6 */
      xmm1 = _mm256_unpackhi_pd(xmm3,xmm4);
      xmm2 = _mm256_unpacklo_pd(xmm3,xmm4);

      xmm3 = _mm256_unpackhi_pd(xmm5,xmm6);
      xmm4 = _mm256_unpacklo_pd(xmm5,xmm6);

      xmm5 = _mm256_add_pd(xmm1,xmm2);
      xmm6 = _mm256_add_pd(xmm3,xmm4);

      xmm1 = _mm256_permute2f128_pd(xmm5,xmm6, _MM_SHUFFLE(0,2,0,1));
      xmm2 = _mm256_blend_pd(xmm5,xmm6,12);
      xmm3 = _mm256_add_pd(xmm1,xmm2);

      /* multiply with frequencies */
      xmm1 = _mm256_mul_pd(xmm0,xmm3);

      /* multiply with clvp */
      xmm2 = _mm256_load_pd(clvp);
      xmm0 = _mm256_mul_pd(xmm1,xmm2);

      /* add up the elements of xmm0 */
      xmm1 = _mm256_hadd_pd(xmm0,xmm0);
      terma_r = ((double *)&xmm1)[0] + ((double *)&xmm1)[2];

      /* apply per-rate scalers, if necessary */
      if (rate_scalings && rate_scalings[i] > 0)
      {
        terma_r *= scale_minlh[rate_scalings[i]-1];
      }

      if (terma_r > 0.)
      {
        /* account for invariant sites */
        prop_invar = invar_proportion ? invar_proportion[freqs_indices[i]] : 0;
        if (prop_invar > 0)
        {
          inv_site_lk = (invar_indices[n] == -1) ?
                            0 : freqs[invar_indices[n]];
          terma += rate_weights[i] * (terma_r * (1 - prop_invar) +
                   inv_site_lk * prop_invar);
        }
        else
        {
          terma += terma_r * rate_weights[i];
        }
      }

      clvp += states_padded;
      clvc += states_padded;
    }

    /* compute site log-likelihood */
    site_lk = log(terma);

    /* apply per-site scaler, if necessary */
    if (site_scalings)
      site_lk += site_scalings * log(PLL_SCALE_THRESHOLD);

    site_lk *= pattern_weights[n];

    /* store per-site log-likelihood */
    if (persite_lnl)
      persite_lnl[n] = site_lk;

    logl += site_lk;
  }

  if (rate_scalings)
    free(rate_scalings);

  return logl;
}
