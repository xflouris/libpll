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

#include "pll.h"

PLL_EXPORT double pll_core_root_loglikelihood_sse(unsigned int states,
                                                  unsigned int sites,
                                                  unsigned int rate_cats,
                                                  const double * clv,
                                                  const unsigned int * scaler,
                                                  double ** frequencies,
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

  __m128d xmm0, xmm1, xmm2, xmm3;

  for (i = 0; i < sites; ++i)
  {
    term = 0;
    for (j = 0; j < rate_cats; ++j)
    {
      freqs = frequencies[freqs_indices[j]];
      xmm3 = _mm_setzero_pd();

      for (k = 0; k < states_padded; k += 2)
      {
        /* load frequencies for current rate matrix */
        xmm0 = _mm_load_pd(freqs);

        /* load clv */
        xmm1 = _mm_load_pd(clv);

        /* multiply with frequencies */
        xmm2 = _mm_mul_pd(xmm0,xmm1);

        xmm3 = _mm_add_pd(xmm3,xmm2);

        freqs += 2;
        clv += 2;
      }

      term_r = ((double *)&xmm3)[0] + ((double *)&xmm3)[1];

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

PLL_EXPORT double pll_core_root_loglikelihood_4x4_sse(unsigned int sites,
                                                      unsigned int rate_cats,
                                                      const double * clv,
                                                      const unsigned int * scaler,
                                                      double ** frequencies,
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

  __m128d xmm0, xmm1, xmm2, xmm3, xmm4, xmm5;

  for (i = 0; i < sites; ++i)
  {
    term = 0;
    for (j = 0; j < rate_cats; ++j)
    {
      freqs = frequencies[freqs_indices[j]];

      /* load frequencies for current rate matrix */
      xmm0 = _mm_load_pd(freqs+0);
      xmm1 = _mm_load_pd(freqs+2);

      /* load clv */
      xmm2 = _mm_load_pd(clv+0);
      xmm3 = _mm_load_pd(clv+2);

      /* multiply with frequencies */
      xmm4 = _mm_mul_pd(xmm0,xmm2);
      xmm5 = _mm_mul_pd(xmm1,xmm3);

      /* add up the elements of xmm2 */
      xmm1 = _mm_add_pd(xmm4,xmm5);

      term_r = ((double *)&xmm1)[0] + ((double *)&xmm1)[1];

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
double pll_core_edge_loglikelihood_ti_sse(unsigned int states,
                                          unsigned int sites,
                                          unsigned int rate_cats,
                                          const double * parent_clv,
                                          const unsigned int * parent_scaler,
                                          const unsigned char * tipchars,
                                          const unsigned int * tipmap,
                                          const double * pmatrix,
                                          double ** frequencies,
                                          const double * rate_weights,
                                          const unsigned int * pattern_weights,
                                          const double * invar_proportion,
                                          const int * invar_indices,
                                          const unsigned int * freqs_indices,
                                          double * persite_lnl)
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
  unsigned int scale_factors;
  unsigned int states_padded = (states+1) & 0xFFFFFFFE;

  __m128d xmm0, xmm1, xmm2, xmm3, xmm4, xmm5;

  size_t displacement = (states_padded - states) * (states_padded);

  xmm5 = _mm_setzero_pd();

  for (n = 0; n < sites; ++n)
  {
    pmat = pmatrix;
    terma = 0;

    cstate = tipmap[tipchars[n]];

    for (i = 0; i < rate_cats; ++i)
    {
      freqs = frequencies[freqs_indices[i]];
      terma_r = 0;

      /* iterate over pairs of rows */
      for (j = 0; j < states_padded; j += 2)
      {
        xmm0 = _mm_setzero_pd();
        xmm1 = _mm_setzero_pd();

        /* point to the two rows */
        const double * row0 = pmat;
        const double * row1 = row0 + states_padded;

        /* set position of least significant bit in character state */
        register int lsb = 0;

        /* iterate pairs of columns */
        for (k = 0; k < states_padded; k += 2)
        {
          /* set mask */
          xmm2 = _mm_set_pd(((cstate >> (lsb+1)) & 1) ? 1 : 0,
                            ((cstate >> (lsb+0)) & 1) ? 1 : 0);
          xmm3 = _mm_cmpgt_pd(xmm2,xmm5);

          lsb += 2;

          /* row 0 */
          xmm2 = _mm_load_pd(row0);
          xmm4 = _mm_and_pd(xmm2,xmm3);
          xmm0 = _mm_add_pd(xmm0,xmm4);
          row0 += 2;

          /* row 1 */
          xmm2 = _mm_load_pd(row1);
          xmm4 = _mm_and_pd(xmm2,xmm3);
          xmm1 = _mm_add_pd(xmm1,xmm4);
          row1 += 2;
        }

        /* point pmatrix to the next four rows */
        pmat = row1;

        /* create a vector containing the sums of xmm0 and xmm1 */
        xmm3 = _mm_hadd_pd(xmm0,xmm1);

        /* multiply with frequencies */
        xmm1 = _mm_load_pd(freqs);
        xmm2 = _mm_mul_pd(xmm3,xmm1);

        /* multiply with clvp */
        xmm0 = _mm_load_pd(clvp);
        xmm1 = _mm_mul_pd(xmm2,xmm0);

        /* add up the elements of xmm1 */
        terma_r += ((double *)&xmm1)[0] + ((double *)&xmm1)[1];

        freqs += 2;
        clvp += 2;
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
    /* count number of scaling factors to acount for */
    scale_factors = (parent_scaler) ? parent_scaler[n] : 0;

    /* compute site log-likelihood and scale if necessary */
    site_lk = log(terma);
    if (scale_factors)
      site_lk += scale_factors * log(PLL_SCALE_THRESHOLD);

    site_lk *= pattern_weights[n];

    /* store per-site log-likelihood */
    if (persite_lnl)
      persite_lnl[n] = site_lk;

    logl += site_lk;
  }
  return logl;
}

PLL_EXPORT
double pll_core_edge_loglikelihood_ii_sse(unsigned int states,
                                          unsigned int sites,
                                          unsigned int rate_cats,
                                          const double * parent_clv,
                                          const unsigned int * parent_scaler,
                                          const double * child_clv,
                                          const unsigned int * child_scaler,
                                          const double * pmatrix,
                                          double ** frequencies,
                                          const double * rate_weights,
                                          const unsigned int * pattern_weights,
                                          const double * invar_proportion,
                                          const int * invar_indices,
                                          const unsigned int * freqs_indices,
                                          double * persite_lnl)
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

  unsigned int scale_factors;
  unsigned int states_padded = (states+1) & 0xFFFFFFFE;

  __m128d xmm0, xmm1, xmm2, xmm3, xmm4, xmm5, xmm6;

  size_t displacement = (states_padded - states) * (states_padded);

  for (n = 0; n < sites; ++n)
  {
    pmat = pmatrix;
    terma = 0;
    for (i = 0; i < rate_cats; ++i)
    {
      freqs = frequencies[freqs_indices[i]];
      terma_r = 0;

      /* iterate over pairs of rows */
      for (j = 0; j < states_padded; j += 2)
      {
        xmm0 = _mm_setzero_pd();
        xmm1 = _mm_setzero_pd();

        /* point to the two rows */
        const double * row0 = pmat;
        const double * row1 = row0 + states_padded;

        /* iterate pairs of columns */
        for (k = 0; k < states_padded; k += 2)
        {
          xmm5 = _mm_load_pd(clvc+k);

          /* row 0 */
          xmm4 = _mm_load_pd(row0);
          xmm6 = _mm_mul_pd(xmm4,xmm5);
          xmm0 = _mm_add_pd(xmm0,xmm6);
          row0 += 2;

          /* row 1 */
          xmm4 = _mm_load_pd(row1);
          xmm6 = _mm_mul_pd(xmm4,xmm5);
          xmm1 = _mm_add_pd(xmm1,xmm6);
          row1 += 2;
        }

        /* point pmatrix to the next two rows */
        pmat = row1;

        /* create a vector containing the sums of xmm0 and xmm1 */
        xmm3 = _mm_hadd_pd(xmm0,xmm1);

        /* multiply with frequencies */
        xmm1 = _mm_load_pd(freqs);
        xmm2 = _mm_mul_pd(xmm3,xmm1);

        /* multiply with clvp */
        xmm0 = _mm_load_pd(clvp);
        xmm1 = _mm_mul_pd(xmm2,xmm0);

        /* add up the elements of xmm1 */
        terma_r += ((double *)&xmm1)[0] + ((double *)&xmm1)[1];

        freqs += 2;
        clvp += 2;
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
    /* count number of scaling factors to acount for */
    scale_factors = (parent_scaler) ? parent_scaler[n] : 0;
    scale_factors += (child_scaler) ? child_scaler[n] : 0;

    /* compute site log-likelihood and scale if necessary */
    site_lk = log(terma);
    if (scale_factors)
      site_lk += scale_factors * log(PLL_SCALE_THRESHOLD);

    site_lk *= pattern_weights[n];

    /* store per-site log-likelihood */
    if (persite_lnl)
      persite_lnl[n] = site_lk;

    logl += site_lk;
  }
  return logl;
}

PLL_EXPORT
double pll_core_edge_loglikelihood_ii_4x4_sse(unsigned int sites,
                                              unsigned int rate_cats,
                                              const double * parent_clv,
                                              const unsigned int * parent_scaler,
                                              const double * child_clv,
                                              const unsigned int * child_scaler,
                                              const double * pmatrix,
                                              double ** frequencies,
                                              const double * rate_weights,
                                              const unsigned int * pattern_weights,
                                              const double * invar_proportion,
                                              const int * invar_indices,
                                              const unsigned int * freqs_indices,
                                              double * persite_lnl)
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

  unsigned int scale_factors;
  unsigned int states = 4;
  unsigned int states_padded = 4;

  __m128d xmm0, xmm1, xmm2, xmm3, xmm4, xmm5, xmm6, xmm7;

  for (n = 0; n < sites; ++n)
  {
    pmat = pmatrix;
    terma = 0;
    for (i = 0; i < rate_cats; ++i)
    {
      freqs = frequencies[freqs_indices[i]];

      /* load first two frequencies for current rate matrix */
      xmm0 = _mm_load_pd(freqs);

      /* load clvc */
      xmm1 = _mm_load_pd(clvc);
      xmm2 = _mm_load_pd(clvc+2);

      /* load pmatrix row 1 */
      xmm3 = _mm_load_pd(pmat);
      xmm4 = _mm_load_pd(pmat+2);

      /* load pmatrix row 2 */
      pmat += states;
      xmm5 = _mm_load_pd(pmat);
      xmm6 = _mm_load_pd(pmat+2);

      /* multiply row 1 with clvc */
      xmm7 = _mm_mul_pd(xmm3,xmm1);
      xmm3 = _mm_mul_pd(xmm4,xmm2);

      /* multiply row 2 with clvc */
      xmm4 = _mm_mul_pd(xmm5,xmm1);
      xmm5 = _mm_mul_pd(xmm6,xmm2);

      /* add the four values of row 1 into a pair */
      xmm6 = _mm_add_pd(xmm7,xmm3);

      /* add the four values  of row 2 into a pair */
      xmm7 = _mm_add_pd(xmm4,xmm5);

      /* horizontally add the two pairs */
      xmm4 = _mm_hadd_pd(xmm6,xmm7);

      /* multiply with frequencies */
      xmm5 = _mm_mul_pd(xmm4,xmm0);

      /* multiply with clvp */
      xmm4 = _mm_load_pd(clvp);
      xmm6 = _mm_mul_pd(xmm5,xmm4);     /* needed */

      /* load pmatrix row 3 */
      pmat += states;
      xmm3 = _mm_load_pd(pmat);
      xmm4 = _mm_load_pd(pmat+2);

      /* load pmatrix row 4 */
      pmat += states;
      xmm5 = _mm_load_pd(pmat);
      xmm0 = _mm_load_pd(pmat+2);

      /* point pmatrix to next rate category */
      pmat += states;

      /* multiply row 3 with clvc */
      xmm7 = _mm_mul_pd(xmm3,xmm1);
      xmm3 = _mm_mul_pd(xmm4,xmm2);

      /* multiply row 4 with clvc */
      xmm4 = _mm_mul_pd(xmm5,xmm1);
      xmm5 = _mm_mul_pd(xmm0,xmm2);

      /* add the four values of row 3 into a pair */
      xmm0 = _mm_add_pd(xmm7,xmm3);

      /* add the four values  of row 4 into a pair */
      xmm7 = _mm_add_pd(xmm4,xmm5);

      /* horizontally add the two pairs */
      xmm4 = _mm_hadd_pd(xmm0,xmm7);

      /* multiply with frequencies */
      xmm0 = _mm_load_pd(freqs+2);
      xmm5 = _mm_mul_pd(xmm4,xmm0);

      /* multiply with clvp */
      xmm4 = _mm_load_pd(clvp+2);
      xmm7 = _mm_mul_pd(xmm5,xmm4);     /* needed */

      /* add the two results together */
      xmm0 = _mm_hadd_pd(xmm6,xmm7);
      terma_r = ((double *)&xmm0)[0] + ((double *)&xmm0)[1];

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

      clvp += states_padded;
      clvc += states_padded;
    }

    /* count number of scaling factors to acount for */
    scale_factors = (parent_scaler) ? parent_scaler[n] : 0;
    scale_factors += (child_scaler) ? child_scaler[n] : 0;

    /* compute site log-likelihood and scale if necessary */
    site_lk = log(terma);
    if (scale_factors)
      site_lk += scale_factors * log(PLL_SCALE_THRESHOLD);

    site_lk *= pattern_weights[n];

    /* store per-site log-likelihood */
    if (persite_lnl)
      persite_lnl[n] = site_lk;

    logl += site_lk;
  }
  return logl;
}

PLL_EXPORT
double pll_core_edge_loglikelihood_ti_4x4_sse(unsigned int sites,
                                              unsigned int rate_cats,
                                              const double * parent_clv,
                                              const unsigned int * parent_scaler,
                                              const unsigned char * tipchars,
                                              const double * pmatrix,
                                              double ** frequencies,
                                              const double * rate_weights,
                                              const unsigned int * pattern_weights,
                                              const double * invar_proportion,
                                              const int * invar_indices,
                                              const unsigned int * freqs_indices,
                                              double * persite_lnl)
{
  unsigned int i,k,n;
  double logl = 0;
  double prop_invar = 0;

  const double * clvp = parent_clv;
  const double * pmat;
  const double * freqs = NULL;

  double terma, terma_r;
  double site_lk, inv_site_lk;

  unsigned int scale_factors;
  unsigned int cstate;
  unsigned int states_padded = 4;
  unsigned int span = rate_cats*states_padded;

  __m128d xmm0, xmm1, xmm2, xmm3, xmm4, xmm5, xmm6, xmm7;
  __m128d ymm0,ymm1,ymm2,ymm3;

  /* precompute a lookup table of four values per entry (one for each state),
     for all 16 states (including ambiguities) and for each rate category. */
  double * lookup = pll_aligned_alloc(64*rate_cats*sizeof(double),
                                      PLL_ALIGNMENT_SSE);
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

  ymm3 = _mm_setzero_pd();

  /* iterate all ambiguities skipping 0 */
  for (i = 1; i < 16; ++i)
  {
    pmat = pmatrix;

    /* mask the entries of pmatrix row to be loaded */
    xmm1 = _mm_set_pd(((i >> 1) & 1) ? 1 : 0,
                      ((i >> 0) & 1) ? 1 : 0);
    xmm0 = _mm_cmpgt_pd(xmm1,ymm3);

    xmm2 = _mm_set_pd(((i >> 3) & 1) ? 1 : 0,
                      ((i >> 2) & 1) ? 1 : 0);
    xmm1 = _mm_cmpgt_pd(xmm2,ymm3);

    for (k = 0; k < rate_cats; ++k)
    {
      freqs = frequencies[freqs_indices[k]];

      /* load row0 from matrix */
      ymm0 = _mm_load_pd(pmat);
      ymm1 = _mm_load_pd(pmat+2);

      /* mask row0 from matrix */
      xmm4 = _mm_and_pd(ymm0,xmm0);
      xmm5 = _mm_and_pd(ymm1,xmm1);
      xmm6 = _mm_add_pd(xmm4,xmm5);     /* a1+a3 | a2+a4 */

      pmat += 4;

      /* load row1 from left matrix */
      ymm0 = _mm_load_pd(pmat);
      ymm1 = _mm_load_pd(pmat+2);

      /* mask row */
      xmm4 = _mm_and_pd(ymm0,xmm0);
      xmm5 = _mm_and_pd(ymm1,xmm1);
      xmm7 = _mm_add_pd(xmm4,xmm5);     /* a5+a7 | a3+a8 */

      ymm2 = _mm_hadd_pd(xmm6,xmm7);

      /* laod frequencies */
      xmm6 = _mm_load_pd(freqs);
      xmm7 = _mm_mul_pd(ymm2,xmm6);
      _mm_store_pd(ptr,xmm7);

      /* point to row2 */
      pmat += 4;

      /* load row2 from left matrix */
      ymm0 = _mm_load_pd(pmat);
      ymm1 = _mm_load_pd(pmat+2);

      /* mask row */
      xmm4 = _mm_and_pd(ymm0,xmm0);
      xmm5 = _mm_and_pd(ymm1,xmm1);
      xmm6 = _mm_add_pd(xmm4,xmm5);     /* a9+a11 | a10+a12 */

      pmat += 4;

      /* load row3 from left matrix */
      ymm0 = _mm_load_pd(pmat);
      ymm1 = _mm_load_pd(pmat+2);

      /* mask row */
      xmm4 = _mm_and_pd(ymm0,xmm0);
      xmm5 = _mm_and_pd(ymm1,xmm1);
      xmm7 = _mm_add_pd(xmm4,xmm5);     /* a13+a15 | a14+a16 */

      ymm2 = _mm_hadd_pd(xmm6,xmm7);

      /* laod frequencies */
      xmm6 = _mm_load_pd(freqs+2);
      xmm7 = _mm_mul_pd(ymm2,xmm6);

      _mm_store_pd(ptr+2,xmm7);

      /* move pointers */
      ptr  += 4;
      pmat += 4;
    }
  }

  for (n = 0; n < sites; ++n)
  {
    pmat = pmatrix;
    terma = 0;

    cstate = tipchars[n];

    unsigned int coffset = cstate*span;

    for (i = 0; i < rate_cats; ++i)
    {
      freqs = frequencies[freqs_indices[i]];

      /* load CLV from clvp */
      xmm0 = _mm_load_pd(clvp);
      xmm1 = _mm_load_pd(clvp+2);

      /* load precomputed lookup table into xmm3 */
      xmm2 = _mm_load_pd(lookup+coffset);
      xmm3 = _mm_load_pd(lookup+coffset+2);

      /* multiply with clvp */
      xmm4 = _mm_mul_pd(xmm0,xmm2);
      xmm5 = _mm_mul_pd(xmm1,xmm3);

      /* add up the elements of xmm0 */
      xmm1 = _mm_hadd_pd(xmm4,xmm5);
      terma_r = ((double *)&xmm1)[0] + ((double *)&xmm1)[1];

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

      clvp += states_padded;
      coffset += 4;
    }

    /* count number of scaling factors to acount for */
    scale_factors = (parent_scaler) ? parent_scaler[n] : 0;

    /* compute site log-likelihood and scale if necessary */
    site_lk = log(terma);
    if (scale_factors)
      site_lk += scale_factors * log(PLL_SCALE_THRESHOLD);

    site_lk *= pattern_weights[n];

    /* store per-site log-likelihood */
    if (persite_lnl)
      persite_lnl[n] = site_lk;

    logl += site_lk;
  }
  pll_aligned_free(lookup);
  return logl;
}
