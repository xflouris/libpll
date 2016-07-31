/*
    Copyright (C) 2016 Tomas Flouri, Diego Darriba, Alexey Kozlov

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
            inv_eigenvecs[i][k * states + j];
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

  memset(tt_eigenvecs, 0, (states_padded * states_padded * rate_cats) * sizeof(double));
  memset(tt_inv_eigenvecs, 0, (states_padded * states_padded * rate_cats) * sizeof(double));

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
    const double * c_eigenvecs      = tt_eigenvecs;
    const double * ct_inv_eigenvecs = tt_inv_eigenvecs;
    for (i = 0; i < rate_cats; ++i)
    {
      t_freqs = freqs[i];

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

          c_eigenvecs      += 4;
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

          c_eigenvecs      += 4;
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

          c_eigenvecs      += 4;
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

          c_eigenvecs      += 4;
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
    }
  }

  pll_aligned_free (tt_inv_eigenvecs);
  pll_aligned_free (tt_eigenvecs);

  return PLL_SUCCESS;
}

PLL_EXPORT int pll_core_update_sumtable_ti_avx(unsigned int states,
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
  unsigned int states_padded = (states+3) & 0xFFFFFFFC;

  unsigned int i, j, k, n;
  unsigned int tipstate;

  double * sum = sumtable;
  const double * t_clvc = parent_clv;
  const double * t_eigenvecs_trans;

  double * eigenvecs_trans = (double *) pll_aligned_alloc (
      (states_padded * states_padded * rate_cats) * sizeof(double),
      PLL_ALIGNMENT_AVX);

  double * precomp_left = (double *) pll_aligned_alloc (
      (states_padded * states_padded * rate_cats) * sizeof(double),
      PLL_ALIGNMENT_AVX);

  if (!eigenvecs_trans || !precomp_left)
  {
    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf (pll_errmsg, 200, "Cannot allocate memory for tt_inv_eigenvecs");
    return PLL_FAILURE;
  }

  /* transpose eigenvecs matrix -> for efficient vectorization */
  for (i = 0; i < rate_cats; ++i)
  {
    for (j = 0; j < states_padded; ++j)
      for (k = 0; k < states_padded; ++k)
      {
        eigenvecs_trans[i*states_padded*states_padded + j*states_padded + k] =
            (j < states && k < states) ? eigenvecs[i][k*states + j] : 0.;
      }
  }

  /* precompute left terms since they are the same for every site */
  double * t_precomp = precomp_left;
  for (i = 0; i < rate_cats; ++i)
  {
        for (k = 0; k < states_padded; ++k)
        {
            for (j = 0; j < states_padded; j += 4)
            {
              if (k < states && j < states)
                {
                  __m256d v_freqs = _mm256_set1_pd(freqs[i][k]);
                  __m256d v_eigen = _mm256_load_pd(inv_eigenvecs[i] +
                                                               k*states + j);
                  __m256d v_lefterm =  _mm256_mul_pd(v_eigen, v_freqs);
                  _mm256_store_pd(t_precomp, v_lefterm);
                }
              else
                memset(t_precomp, 0, 4 * sizeof(double));

              t_precomp += 4;
            }
        }
  }

  /* build sumtable */
  for (n = 0; n < sites; n++)
  {
    // TODO: Fix this BS!
    if (states == 4)
      tipstate = (unsigned int) left_tipchars[n];
    else
      tipstate = tipmap[(unsigned int)left_tipchars[n]];

    for (i = 0; i < rate_cats; ++i)
    {
      t_eigenvecs_trans = eigenvecs_trans + i * states_padded * states_padded;
      t_precomp = precomp_left + i * states_padded * states_padded;

      for (j = 0; j < states_padded; j += 4)
        {
          __m256d v_lefterm = _mm256_setzero_pd();
          __m256d v_righterm = _mm256_setzero_pd();

          for (k = 0; k < states; ++k)
            {
              __m256d v_clvc = _mm256_set1_pd(t_clvc[k]);
              __m256d v_eigen = _mm256_load_pd(t_eigenvecs_trans +
                                                          k*states_padded + j);
              v_righterm =  _mm256_add_pd(v_righterm,
                                          _mm256_mul_pd(v_eigen, v_clvc));

              if ((tipstate >> k) & 1)
                {
                  v_lefterm = _mm256_add_pd(v_lefterm,
                                 _mm256_load_pd(t_precomp +
                                                         k*states_padded + j));
                }
            }

          __m256d v_sum = _mm256_mul_pd(v_lefterm, v_righterm);
          _mm256_store_pd(sum + j, v_sum);

        }

      t_clvc += states_padded;
      sum += states_padded;
    }
  }

  pll_aligned_free(eigenvecs_trans);
  pll_aligned_free(precomp_left);

  return PLL_SUCCESS;
}

PLL_EXPORT void core_site_likelihood_derivatives_avx(unsigned int states,
                                             unsigned int states_padded,
                                             unsigned int rate_cats,
                                             const double * rate_weights,
                                             const double * prop_invar,
                                             const double * lk_invar,
                                             const double * sumtable,
                                             const double * diagptable,
                                             double * site_lk)
{
  unsigned int i,j;
  const double *sum = sumtable;
  const double * diagp = diagptable;

  __m256d v_sitelk = _mm256_setzero_pd ();
  for (i = 0; i < rate_cats; ++i)
  {
    __m256d v_cat_sitelk = _mm256_setzero_pd ();
    for (j = 0; j < states_padded; j++)
    {
      __m256d v_diagp = _mm256_load_pd(diagp);
      __m256d v_sum = _mm256_set1_pd(sum[j]);
      v_cat_sitelk = _mm256_add_pd (v_cat_sitelk, _mm256_mul_pd(v_sum, v_diagp));

//      v_sitelk = _mm256_add_pd (v_sitelk, _mm256_mul_pd(v_sum, v_diagp));

      diagp += 4;
    }

    /* account for invariant sites */
    const double t_prop_invar = prop_invar[i];
    if (t_prop_invar > 0)
    {
      __m256d v_inv_prop = _mm256_set1_pd(1. - t_prop_invar);
      v_cat_sitelk = _mm256_mul_pd(v_cat_sitelk, v_inv_prop);

      if (lk_invar)
        {
          __m256d v_inv_lk = _mm256_setr_pd(lk_invar[i], 0., 0., 0.);
          v_cat_sitelk = _mm256_add_pd(v_cat_sitelk, v_inv_lk);
        }
    }

    __m256d v_weight = _mm256_set1_pd(rate_weights[i]);
    v_sitelk = _mm256_add_pd (v_sitelk, _mm256_mul_pd(v_cat_sitelk, v_weight));

    sum += states_padded;
  }

  _mm256_store_pd(site_lk, v_sitelk);
}

PLL_EXPORT void core_site_likelihood_derivatives_4x4_avx(unsigned int rate_cats,
                                             const double * rate_weights,
                                             const double * prop_invar,
                                             const double * lk_invar,
                                             const double * sumtable,
                                             const double * diagptable,
                                             double * site_lk)
{
  unsigned int i;
  const double *sum = sumtable;
  const double * diagp = diagptable;

  __m256d v_sitelk = _mm256_setzero_pd ();
  for (i = 0; i < rate_cats; ++i)
  {
      __m256d v_cat_sitelk = _mm256_setzero_pd ();

      __m256d v_diagp = _mm256_load_pd(diagp);
      __m256d v_sum = _mm256_set1_pd(sum[0]);
      v_cat_sitelk = _mm256_add_pd (v_cat_sitelk, _mm256_mul_pd(v_sum, v_diagp));

      v_diagp = _mm256_load_pd(diagp + 4);
      v_sum = _mm256_set1_pd(sum[1]);
      v_cat_sitelk = _mm256_add_pd (v_cat_sitelk, _mm256_mul_pd(v_sum, v_diagp));

      v_diagp = _mm256_load_pd(diagp + 8);
      v_sum = _mm256_set1_pd(sum[2]);
      v_cat_sitelk = _mm256_add_pd (v_cat_sitelk, _mm256_mul_pd(v_sum, v_diagp));

      v_diagp = _mm256_load_pd(diagp + 12);
      v_sum = _mm256_set1_pd(sum[3]);
      v_cat_sitelk = _mm256_add_pd (v_cat_sitelk, _mm256_mul_pd(v_sum, v_diagp));

      /* account for invariant sites */
      const double t_prop_invar = prop_invar[i];
      if (t_prop_invar > 0)
      {
        __m256d v_inv_prop = _mm256_set1_pd(1. - t_prop_invar);
        v_cat_sitelk = _mm256_mul_pd(v_cat_sitelk, v_inv_prop);

        if (lk_invar)
          {
            __m256d v_inv_lk = _mm256_setr_pd(lk_invar[i], 0., 0., 0.);
            v_cat_sitelk = _mm256_add_pd(v_cat_sitelk, v_inv_lk);
          }
      }

      __m256d v_weight = _mm256_set1_pd(rate_weights[i]);
      v_sitelk = _mm256_add_pd (v_sitelk, _mm256_mul_pd(v_cat_sitelk, v_weight));

      diagp += 16;
      sum += 4;
  }
  _mm256_store_pd(site_lk, v_sitelk);
}
