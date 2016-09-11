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

#define ONESTEP(x,baseptr)                                      \
            ymm0 = _mm256_load_pd(baseptr+0);                   \
            ymm1 = _mm256_load_pd(baseptr+4);                   \
            ymm2 = _mm256_load_pd(baseptr+8);                   \
            ymm3 = _mm256_load_pd(baseptr+12);                  \
            ymm4 = _mm256_load_pd(baseptr+16);                  \
                                                                \
            ymm5 = _mm256_mul_pd(xmm4,ymm0);                    \
            ymm6 = _mm256_mul_pd(xmm5,ymm1);                    \
            ymm7 = _mm256_mul_pd(xmm6,ymm2);                    \
            ymm8 = _mm256_mul_pd(xmm7,ymm3);                    \
            ymm9 = _mm256_mul_pd(xmm8,ymm4);                    \
                                                                \
            xmm1 = _mm256_add_pd(ymm5,ymm6);                    \
            ymm5 = _mm256_add_pd(xmm1,ymm7);                    \
            ymm6 = _mm256_add_pd(ymm5,ymm8);                    \
            x    = _mm256_add_pd(ymm6,ymm9);    /* row x */

PLL_EXPORT int pll_core_update_pmatrix_4x4_avx(double ** pmatrix,
                                               unsigned int rate_cats,
                                               double * rates,
                                               const double * branch_lengths,
                                               const unsigned int * matrix_indices,
                                               const unsigned int * params_indices,
                                               double * prop_invar,
                                               double ** eigenvals,
                                               double ** eigenvecs,
                                               double ** inv_eigenvecs,
                                               unsigned int count)
{
  unsigned int i,n;
  double * expd;

  double pinvar;
  double * evecs;
  double * inv_evecs;
  double * evals;
  double * pmat;

  expd = (double *)pll_aligned_alloc(4*sizeof(double), PLL_ALIGNMENT_AVX);

  if (!expd)
  {
    if (expd) pll_aligned_free(expd);

    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf(pll_errmsg, 200, "Unable to allocate enough memory.");
    return PLL_FAILURE;
  }

  __m256d xmm0,xmm1,xmm2,xmm3,xmm4,xmm5,xmm6,xmm7,xmm8,xmm9;
  __m256d ymm0,ymm1,ymm2,ymm3,ymm4,ymm5,ymm6,ymm7;
  xmm0 = _mm256_setzero_pd();

  for (i = 0; i < count; ++i)
  {
    assert(branch_lengths[i] >= 0);
    
    xmm3 = _mm256_set1_pd(branch_lengths[i]);
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
        _mm256_store_pd(pmat+0,xmm0);
        _mm256_store_pd(pmat+4,xmm0);
        _mm256_store_pd(pmat+8,xmm0);
        _mm256_store_pd(pmat+12,xmm0);

        pmat[0] = pmat[5] = pmat[10] = pmat[15] = 1;
      }
      else
      {
        /* exponentiate eigenvalues */

        /* 1) load eigenvalues and 2) load rate into all slots of register */
        xmm1 = _mm256_load_pd(evals);
        xmm2 = _mm256_set1_pd(rates[n]);

        /* multiply eigenvalues with rate */
        xmm4 = _mm256_mul_pd(xmm1,xmm2);

        /* multiply product with  branch length */
        xmm2 = _mm256_mul_pd(xmm4,xmm3);

        if (pinvar > PLL_MISC_EPSILON)
        {
          xmm1 = _mm256_set1_pd(1.0 - pinvar);
          xmm2 = _mm256_div_pd(xmm2,xmm1);
        }
          
        /* TODO: implement a vectorized double-precision exponentiation */
        //xmm1 = _mm256_exp_pd(xmm2);     /* expd */

        /* for now exponentiate non-vectorized */
        _mm256_store_pd(expd,xmm2);
        xmm1 = _mm256_set_pd(exp(expd[3]),
                             exp(expd[2]),
                             exp(expd[1]),
                             exp(expd[0]));

        /* multiply inverse eigenvectors with computed result */
        xmm2 = _mm256_load_pd(inv_evecs+0);
        xmm4 = _mm256_load_pd(inv_evecs+4);
        xmm5 = _mm256_load_pd(inv_evecs+8);
        xmm6 = _mm256_load_pd(inv_evecs+12);

        xmm7 = _mm256_mul_pd(xmm2,xmm1);       /* temp row 0 */
        xmm2 = _mm256_mul_pd(xmm4,xmm1);       /* temp row 1 */
        xmm4 = _mm256_mul_pd(xmm5,xmm1);       /* temp row 2 */
        xmm5 = _mm256_mul_pd(xmm6,xmm1);       /* temp row 3 */

        /* transpose eigenvector */
        xmm6 = _mm256_load_pd(evecs+0); 
        xmm1 = _mm256_load_pd(evecs+4);
        xmm8 = _mm256_load_pd(evecs+8);
        xmm9 = _mm256_load_pd(evecs+12);

        /* transpose eigenvectors */
        ymm0 = _mm256_unpacklo_pd(xmm6,xmm1);
        ymm1 = _mm256_unpackhi_pd(xmm6,xmm1);
        ymm2 = _mm256_unpacklo_pd(xmm8,xmm9);
        ymm3 = _mm256_unpackhi_pd(xmm8,xmm9);

        xmm6 = _mm256_permute2f128_pd(ymm0,ymm2,_MM_SHUFFLE(0,2,0,0));
        xmm1 = _mm256_permute2f128_pd(ymm1,ymm3,_MM_SHUFFLE(0,2,0,0));
        xmm8 = _mm256_permute2f128_pd(ymm0,ymm2,_MM_SHUFFLE(0,3,0,1));
        xmm9 = _mm256_permute2f128_pd(ymm1,ymm3,_MM_SHUFFLE(0,3,0,1));

        /* pmat row 0*/

        ymm0 = _mm256_mul_pd(xmm7,xmm6);
        ymm1 = _mm256_mul_pd(xmm7,xmm1);
        ymm2 = _mm256_mul_pd(xmm7,xmm8);
        ymm3 = _mm256_mul_pd(xmm7,xmm9);

        /* create a vector with the sums of ymm0, ymm1, ymm2, ymm3 */
        ymm4 = _mm256_unpackhi_pd(ymm0,ymm1);
        ymm5 = _mm256_unpacklo_pd(ymm0,ymm1);

        ymm6 = _mm256_unpackhi_pd(ymm2,ymm3);
        ymm7 = _mm256_unpacklo_pd(ymm2,ymm3);

        ymm0 = _mm256_add_pd(ymm4,ymm5);
        ymm1 = _mm256_add_pd(ymm6,ymm7);

        ymm2 = _mm256_permute2f128_pd(ymm0,ymm1, _MM_SHUFFLE(0,2,0,1));
        ymm3 = _mm256_blend_pd(ymm0,ymm1,12);
        ymm0 = _mm256_add_pd(ymm2,ymm3);

        _mm256_store_pd(pmat+0,ymm0);

        /* pmat row 1 */

        ymm0 = _mm256_mul_pd(xmm2,xmm6);
        ymm1 = _mm256_mul_pd(xmm2,xmm1);
        ymm2 = _mm256_mul_pd(xmm2,xmm8);
        ymm3 = _mm256_mul_pd(xmm2,xmm9);

        /* create a vector with the sums of ymm0, ymm1, ymm2, ymm3 */
        ymm4 = _mm256_unpackhi_pd(ymm0,ymm1);
        ymm5 = _mm256_unpacklo_pd(ymm0,ymm1);

        ymm6 = _mm256_unpackhi_pd(ymm2,ymm3);
        ymm7 = _mm256_unpacklo_pd(ymm2,ymm3);

        ymm0 = _mm256_add_pd(ymm4,ymm5);
        ymm1 = _mm256_add_pd(ymm6,ymm7);

        ymm2 = _mm256_permute2f128_pd(ymm0,ymm1, _MM_SHUFFLE(0,2,0,1));
        ymm3 = _mm256_blend_pd(ymm0,ymm1,12);
        ymm0 = _mm256_add_pd(ymm2,ymm3);

        _mm256_store_pd(pmat+4,ymm0);

        /* pmat row 2 */

        ymm0 = _mm256_mul_pd(xmm4,xmm6);
        ymm1 = _mm256_mul_pd(xmm4,xmm1);
        ymm2 = _mm256_mul_pd(xmm4,xmm8);
        ymm3 = _mm256_mul_pd(xmm4,xmm9);

        /* create a vector with the sums of ymm0, ymm1, ymm2, ymm3 */
        ymm4 = _mm256_unpackhi_pd(ymm0,ymm1);
        ymm5 = _mm256_unpacklo_pd(ymm0,ymm1);

        ymm6 = _mm256_unpackhi_pd(ymm2,ymm3);
        ymm7 = _mm256_unpacklo_pd(ymm2,ymm3);

        ymm0 = _mm256_add_pd(ymm4,ymm5);
        ymm1 = _mm256_add_pd(ymm6,ymm7);

        ymm2 = _mm256_permute2f128_pd(ymm0,ymm1, _MM_SHUFFLE(0,2,0,1));
        ymm3 = _mm256_blend_pd(ymm0,ymm1,12);
        ymm0 = _mm256_add_pd(ymm2,ymm3);

        _mm256_store_pd(pmat+8,ymm0);

        /* pmat row 3 */

        ymm0 = _mm256_mul_pd(xmm5,xmm6);
        ymm1 = _mm256_mul_pd(xmm5,xmm1);
        ymm2 = _mm256_mul_pd(xmm5,xmm8);
        ymm3 = _mm256_mul_pd(xmm5,xmm9);

        /* create a vector with the sums of ymm0, ymm1, ymm2, ymm3 */
        ymm4 = _mm256_unpackhi_pd(ymm0,ymm1);
        ymm5 = _mm256_unpacklo_pd(ymm0,ymm1);

        ymm6 = _mm256_unpackhi_pd(ymm2,ymm3);
        ymm7 = _mm256_unpacklo_pd(ymm2,ymm3);

        ymm0 = _mm256_add_pd(ymm4,ymm5);
        ymm1 = _mm256_add_pd(ymm6,ymm7);

        ymm2 = _mm256_permute2f128_pd(ymm0,ymm1, _MM_SHUFFLE(0,2,0,1));
        ymm3 = _mm256_blend_pd(ymm0,ymm1,12);
        ymm0 = _mm256_add_pd(ymm2,ymm3);

        _mm256_store_pd(pmat+12,ymm0);

      }
      pmat = pmat+16;
    }
  }

  pll_aligned_free(expd);
  return PLL_SUCCESS;
}

PLL_EXPORT int pll_core_update_pmatrix_20x20_avx(double ** pmatrix,
                                                 unsigned int rate_cats,
                                                 double * rates,
                                                 const double * branch_lengths,
                                                 const unsigned int * matrix_indices,
                                                 const unsigned int * params_indices,
                                                 double * prop_invar,
                                                 double ** eigenvals,
                                                 double ** eigenvecs,
                                                 double ** inv_eigenvecs,
                                                 unsigned int count)
{
  unsigned int i,n,j,k;
  double pinvar;

  int * transposed;
  double * evecs;
  double * inv_evecs;
  double * evals;
  double * pmat;
  double * expd;
  double * temp;
  double ** tran_evecs;

  expd = (double *)pll_aligned_alloc(20*sizeof(double), PLL_ALIGNMENT_AVX);
  temp = (double *)pll_aligned_alloc(400*sizeof(double), PLL_ALIGNMENT_AVX);

  /* transposed eigen vectors */
  transposed = (int *)calloc((size_t)rate_cats, sizeof(int));
  tran_evecs= (double **)calloc((size_t)rate_cats, sizeof(double *));

  if (!expd || !temp || !transposed || !tran_evecs)
  {
    if (expd) pll_aligned_free(expd);
    if (temp) pll_aligned_free(temp);
    if (transposed) free(transposed);
    if (tran_evecs) free(tran_evecs);

    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf(pll_errmsg, 200, "Unable to allocate enough memory.");
    return PLL_FAILURE;
  }

  /* transpose eigenvectors */
  /* TODO: The same trick can be applied for exponentiations */
  for (n = 0; n < rate_cats; ++n)
  {
    int index = params_indices[n];

    if (!transposed[index])
    {
      /* allocate space for transposed eigenvectors and check that
         allocation succeeds */
      double * tran = (double *)pll_aligned_alloc(400*sizeof(double),
                                                  PLL_ALIGNMENT_AVX);
      if (!tran)
      {
        pll_aligned_free(expd);
        pll_aligned_free(temp);
        free(transposed);
        for (i = 0; i < n; ++i)
          if (tran_evecs[i]) pll_aligned_free(tran_evecs[i]);
        free(tran_evecs);

        pll_errno = PLL_ERROR_MEM_ALLOC;
        snprintf(pll_errmsg, 200, "Unable to allocate enough memory.");
        return PLL_FAILURE;
      }

      /* transpose eigen vectors */
      evecs = eigenvecs[index];
      for (i = 0; i < 20; ++i)
      {
        for (j = 0; j < 20; ++j)
          tran[i*20+j] = evecs[j*20+i];
      }
      
      /* update pointers and indicate that the eigen vector for the current
         rate matrix with index was updated */
      tran_evecs[index] = tran;
      transposed[index] = 1;
    }
  }
  free(transposed);

  __m256d xmm0,xmm1,xmm2,xmm3,xmm4,xmm5,xmm6,xmm7,xmm8;
  __m256d ymm0,ymm1,ymm2,ymm3,ymm4,ymm5,ymm6,ymm7,ymm8,ymm9;
  __m256d zmm0,zmm1,zmm2,zmm3;

  double * tran = NULL;
  for (i = 0; i < count; ++i)
  {
    assert(branch_lengths[i] >= 0);
    
    xmm3 = _mm256_set1_pd(branch_lengths[i]);
    pmat = pmatrix[matrix_indices[i]];

    /* compute effective pmatrix location */
    for (n = 0; n < rate_cats; ++n)
    {
      pinvar = prop_invar[params_indices[n]];
      tran = tran_evecs[params_indices[n]];
      inv_evecs = inv_eigenvecs[params_indices[n]];
      evals = eigenvals[params_indices[n]];

      /* if branch length is zero then set the p-matrix to identity matrix */
      if (!branch_lengths[i])
      {
        xmm0 = _mm256_setzero_pd();
        for (j = 0; j < 20; ++j)
        {
          _mm256_store_pd(pmat+0,xmm0);
          _mm256_store_pd(pmat+4,xmm0);
          _mm256_store_pd(pmat+8,xmm0);
          _mm256_store_pd(pmat+12,xmm0);
          _mm256_store_pd(pmat+16,xmm0);
          pmat[j] = 1;
          pmat += 20;
        }
        continue;
      }

      /* exponentiate eigenvalues */
      xmm2 = _mm256_set1_pd(rates[n]);

      if (pinvar > PLL_MISC_EPSILON)
        xmm6 = _mm256_set1_pd(1.0 - pinvar);

      for (k = 0; k < 5; ++k)
      {
        xmm1 = _mm256_load_pd(evals+k*4);

        /* scalar multiplication with rates */
        xmm4 = _mm256_mul_pd(xmm1,xmm2);

        /* scalar multiplication with branch lengths */
        xmm5 = _mm256_mul_pd(xmm4,xmm3);

        if (pinvar > PLL_MISC_EPSILON)
        {
          xmm5 = _mm256_div_pd(xmm5,xmm6);
        }

        _mm256_store_pd(expd+k*4,xmm5);
      }

      for (k = 0; k < 20; ++k)
        expd[k] = exp(expd[k]);
        
      /* load expd */
      xmm4 = _mm256_load_pd(expd+0);
      xmm5 = _mm256_load_pd(expd+4);
      xmm6 = _mm256_load_pd(expd+8);
      xmm7 = _mm256_load_pd(expd+12);
      xmm8 = _mm256_load_pd(expd+16);

      /* compute temp matrix */
      for (k = 0; k < 400; k += 20)
      {
        ymm0 = _mm256_load_pd(inv_evecs+k+0);
        ymm1 = _mm256_load_pd(inv_evecs+k+4);
        ymm2 = _mm256_load_pd(inv_evecs+k+8);
        ymm3 = _mm256_load_pd(inv_evecs+k+12);
        ymm4 = _mm256_load_pd(inv_evecs+k+16);

        ymm5 = _mm256_mul_pd(xmm4,ymm0);
        ymm6 = _mm256_mul_pd(xmm5,ymm1);
        ymm7 = _mm256_mul_pd(xmm6,ymm2);
        ymm8 = _mm256_mul_pd(xmm7,ymm3);
        ymm9 = _mm256_mul_pd(xmm8,ymm4);

        _mm256_store_pd(temp+k+0,ymm5);
        _mm256_store_pd(temp+k+4,ymm6);
        _mm256_store_pd(temp+k+8,ymm7);
        _mm256_store_pd(temp+k+12,ymm8);
        _mm256_store_pd(temp+k+16,ymm9);
      }

      for (j = 0; j < 400; j += 20)
      {
        xmm4 = _mm256_load_pd(temp+j+0);
        xmm5 = _mm256_load_pd(temp+j+4);
        xmm6 = _mm256_load_pd(temp+j+8);
        xmm7 = _mm256_load_pd(temp+j+12);
        xmm8 = _mm256_load_pd(temp+j+16);

        /* process four rows at a time */
        for (k = 0; k < 400; k += 80)
        {
          /* row 0 */
          ONESTEP(zmm0,tran+k+0);

          /* row 1 */
          ONESTEP(zmm1,tran+k+20);

          /* row 2 */
          ONESTEP(zmm2,tran+k+40);

          /* row 3 */
          ONESTEP(zmm3,tran+k+60);

          /* create a vector with the sums of zmm0, zmm1, zmm2, zmm3 */
          ymm4 = _mm256_unpackhi_pd(zmm0,zmm1);
          ymm5 = _mm256_unpacklo_pd(zmm0,zmm1);

          ymm6 = _mm256_unpackhi_pd(zmm2,zmm3);
          ymm7 = _mm256_unpacklo_pd(zmm2,zmm3);

          ymm0 = _mm256_add_pd(ymm4,ymm5);
          ymm1 = _mm256_add_pd(ymm6,ymm7);

          ymm2 = _mm256_permute2f128_pd(ymm0,ymm1, _MM_SHUFFLE(0,2,0,1));
          ymm3 = _mm256_blend_pd(ymm0,ymm1,12);
          ymm0 = _mm256_add_pd(ymm2,ymm3);

          _mm256_store_pd(pmat,ymm0);

          pmat += 4;
        }
      }
    }
  }

  pll_aligned_free(expd);
  pll_aligned_free(temp);

  for (i = 0; i < rate_cats; ++i)
    if (tran_evecs[i]) pll_aligned_free(tran_evecs[i]); 

  free(tran_evecs);
  return PLL_SUCCESS;
}
