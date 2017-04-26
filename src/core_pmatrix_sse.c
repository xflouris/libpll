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

#define ONESTEP(x)                                                             \
    /* compute pmat row x/4 */                                                 \
    xmm12 = _mm_load_pd(inv_evecs+x);                                          \
    xmm13 = _mm_load_pd(inv_evecs+x+2);                                        \
    xmm12 = _mm_mul_pd(xmm12,xmm1);          /* temp row x/4 (0-1) */          \
    xmm13 = _mm_mul_pd(xmm13,xmm2);          /* temp row x/4 (2-3) */          \
                                                                               \
    /* multiply with row 0 of transposed eigenvector */                        \
    xmm14 = _mm_mul_pd(xmm12,xmm4);                                            \
    xmm15 = _mm_mul_pd(xmm13,xmm8);                                            \
    xmm14 = _mm_add_pd(xmm14,xmm15);                                           \
                                                                               \
    /* multiply with row 1 of transposed eigenvector */                        \
    xmm15 = _mm_mul_pd(xmm12,xmm5);                                            \
    xmm16 = _mm_mul_pd(xmm13,xmm9);                                            \
    xmm15 = _mm_add_pd(xmm15,xmm16);                                           \
                                                                               \
    xmm16 = _mm_hadd_pd(xmm14,xmm15);                                          \
    _mm_store_pd(pmat+x,xmm16);                                                \
                                                                               \
    /* multiply with row 2 of transposed eigenvector */                        \
    xmm14 = _mm_mul_pd(xmm12,xmm6);                                            \
    xmm15 = _mm_mul_pd(xmm13,xmm10);                                           \
    xmm14 = _mm_add_pd(xmm14,xmm15);                                           \
                                                                               \
    /* multiply with row 3 of transposed eigenvector */                        \
    xmm15 = _mm_mul_pd(xmm12,xmm7);                                            \
    xmm16 = _mm_mul_pd(xmm13,xmm11);                                           \
    xmm15 = _mm_add_pd(xmm15,xmm16);                                           \
                                                                               \
    xmm16 = _mm_hadd_pd(xmm14,xmm15);                                          \
    _mm_store_pd(pmat+x+2,xmm16);                                              \


PLL_EXPORT int pll_core_update_pmatrix_4x4_sse(double ** pmatrix,
                                               unsigned int rate_cats,
                                               const double * rates,
                                               const double * branch_lengths,
                                               const unsigned int * matrix_indices,
                                               const unsigned int * params_indices,
                                               const double * prop_invar,
                                               double * const * eigenvals,
                                               double * const * eigenvecs,
                                               double * const * inv_eigenvecs,
                                               unsigned int count)
{
  unsigned int i,n;
  double * expd;

  double pinvar;
  double * evecs;
  double * inv_evecs;
  double * evals;
  double * pmat;

  expd = (double *)pll_aligned_alloc(4*sizeof(double), PLL_ALIGNMENT_SSE);

  if (!expd)
  {
    if (expd) pll_aligned_free(expd);

    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf(pll_errmsg, 200, "Unable to allocate enough memory.");
    return PLL_FAILURE;
  }

  __m128d xmm0, xmm1, xmm2, xmm3, xmm4, xmm5, xmm6, xmm7, xmm8, xmm9;
  __m128d xmm10, xmm11, xmm12, xmm13, xmm14, xmm15, xmm16;
  __m128d v_onemin, v_onemax;

  xmm0 = _mm_setzero_pd();
  v_onemin = _mm_set1_pd(PLL_ONE_MIN);
  v_onemax = _mm_set1_pd(PLL_ONE_MAX);

  for (i = 0; i < count; ++i)
  {
    assert(branch_lengths[i] >= 0);
    
    xmm3 = _mm_set1_pd(branch_lengths[i]);
    pmat = pmatrix[matrix_indices[i]];

    /* compute effective pmatrix location */
    for (n = 0; n < rate_cats; ++n)
    {
      pinvar = prop_invar[params_indices[n]];
      evecs = eigenvecs[params_indices[n]];
      inv_evecs = inv_eigenvecs[params_indices[n]];
      evals = eigenvals[params_indices[n]];

      /* if branch length is zero then set the p-matrix to identity matrix */
      if (!branch_lengths[i])
      {
        _mm_store_pd(pmat+0, xmm0);
        _mm_store_pd(pmat+2, xmm0);
        _mm_store_pd(pmat+4, xmm0);
        _mm_store_pd(pmat+6, xmm0);
        _mm_store_pd(pmat+8, xmm0);
        _mm_store_pd(pmat+10,xmm0);
        _mm_store_pd(pmat+12,xmm0);
        _mm_store_pd(pmat+14,xmm0);

        pmat[0] = pmat[5] = pmat[10] = pmat[15] = 1;
      }
      else
      {
        /* exponentiate eigenvalues */

        /* 1) load eigenvalues and 2) load rate into all slots of register */
        xmm1 = _mm_load_pd(evals+0);
        xmm2 = _mm_load_pd(evals+2);
        xmm4 = _mm_set1_pd(rates[n]);

        /* multiply eigenvalues with rate */
        xmm5 = _mm_mul_pd(xmm1,xmm4);
        xmm6 = _mm_mul_pd(xmm2,xmm4);

        /* multiply product with  branch length */
        xmm7 = _mm_mul_pd(xmm5,xmm3);
        xmm8 = _mm_mul_pd(xmm6,xmm3);

        if (pinvar > PLL_MISC_EPSILON)
        {
          xmm1 = _mm_set1_pd(1.0 - pinvar);
          xmm7 = _mm_div_pd(xmm7,xmm1);
          xmm8 = _mm_div_pd(xmm8,xmm1);
        }
          
        /* TODO: implement a vectorized double-precision exponentiation */
        //xmm1 = _mm_exp_pd(xmm7);     /* expd */
        //xmm2 = _mm_exp_pd(xmm8);     /* expd */

        /* for now exponentiate non-vectorized */
        _mm_store_pd(expd+0,xmm7);
        _mm_store_pd(expd+2,xmm8);

        /* check if all values of expd are approximately one */
        xmm12 = _mm_set_pd(exp(expd[1]), exp(expd[0]));
        xmm13 = _mm_set_pd(exp(expd[3]), exp(expd[2]));

        /* */
        xmm1 = _mm_cmpgt_pd(xmm12,v_onemin);
        xmm2 = _mm_cmplt_pd(xmm13,v_onemax);
        xmm4 = _mm_and_pd(xmm1,xmm2);

        if (_mm_movemask_pd(xmm4) == 0x3)
        {
          _mm_store_pd(pmat+0, xmm0);
          _mm_store_pd(pmat+2, xmm0);
          _mm_store_pd(pmat+4, xmm0);
          _mm_store_pd(pmat+6, xmm0);
          _mm_store_pd(pmat+8, xmm0);
          _mm_store_pd(pmat+10,xmm0);
          _mm_store_pd(pmat+12,xmm0);
          _mm_store_pd(pmat+14,xmm0);

          pmat[0] = pmat[5] = pmat[10] = pmat[15] = 1;
        }
        else
        {
          /* transpose eigenvector */

          xmm1 = _mm_load_pd(evecs+0);
          xmm2 = _mm_load_pd(evecs+4);
          xmm4 = _mm_unpacklo_pd(xmm1,xmm2);     /* row 0 (0,1) */
          xmm5 = _mm_unpackhi_pd(xmm1,xmm2);     /* row 1 (0,1) */

          xmm1 = _mm_load_pd(evecs+2);
          xmm2 = _mm_load_pd(evecs+6);
          xmm6 = _mm_unpacklo_pd(xmm1,xmm2);     /* row 2 (0,1) */
          xmm7 = _mm_unpackhi_pd(xmm1,xmm2);     /* row 3 (0,1) */

          xmm1 = _mm_load_pd(evecs+8);
          xmm2 = _mm_load_pd(evecs+12);
          xmm8 = _mm_unpacklo_pd(xmm1,xmm2);     /* row 0 (2,3) */
          xmm9 = _mm_unpackhi_pd(xmm1,xmm2);     /* row 1 (2,3) */

          xmm1 = _mm_load_pd(evecs+10);
          xmm2 = _mm_load_pd(evecs+14);
          xmm10 = _mm_unpacklo_pd(xmm1,xmm2);    /* row 2 (2,3) */
          xmm11 = _mm_unpackhi_pd(xmm1,xmm2);    /* row 3 (2,3) */

          /* load exponentiated eigenvalues */
          /*
          xmm1 = _mm_set_pd(exp(expd[1]),
                            exp(expd[0]));
          xmm2 = _mm_set_pd(exp(expd[3]),
                            exp(expd[2]));
          */
          xmm1 = xmm12;
          xmm2 = xmm13;

          /* compute pmatrix */
          ONESTEP(0);
          ONESTEP(4);
          ONESTEP(8);
          ONESTEP(12);
        }
      }
      #ifdef DEBUG
      unsigned int j,k;
      for (j = 0; j < 4; ++j)
        for (k = 0; k < 4; ++k)
          assert(pmat[j*4+k] >= 0);
      #endif
      pmat = pmat+16;
    }
  }

  pll_aligned_free(expd);
  return PLL_SUCCESS;
}

