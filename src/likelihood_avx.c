/*
    Copyright (C) 2015 Tomas Flouri, Kassian Kobert

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

void pll_update_partials_avx(pll_partition_t * partition,
                             const pll_operation_t * op)
{
  unsigned int i,j,k,n;
  unsigned int scaling;

  const double * left_clv = partition->clv[op->child1_clv_index];
  const double * right_clv = partition->clv[op->child2_clv_index];
  const double * left_matrix = partition->pmatrix[op->child1_matrix_index];
  const double * right_matrix = partition->pmatrix[op->child2_matrix_index];
  double * parent_clv = partition->clv[op->parent_clv_index];
  unsigned int * scaler = (op->parent_scaler_index == PLL_SCALE_BUFFER_NONE) ?
                        NULL : partition->scale_buffer[op->parent_scaler_index];

  const double * lmat;
  const double * rmat;

  unsigned int states = partition->states;
  unsigned int states_padded = partition->states_padded;
  unsigned int span = states_padded * partition->rate_cats;

  /* the same loop as update_partials() in likelihood.c but vectorized */
  for (n = 0; n < partition->sites; ++n)
  {
    lmat = left_matrix;
    rmat = right_matrix;
    scaling = (scaler) ? 1 : 0;
    for (k = 0; k < partition->rate_cats; ++k)
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

        for (j = 0; j < states; ++j)
          scaling = scaling && (parent_clv[j] < PLL_SCALE_THRESHOLD);
      }
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

      parent_clv -= span;
      for (i = 0; i < span; i += 4)
      {
        __m256d v_prod = _mm256_load_pd(parent_clv + i);
        v_prod = _mm256_mul_pd(v_prod,v_scale_factor);
        _mm256_store_pd(parent_clv + i, v_prod);
      }
      parent_clv += span;
      scaler[n] += 1;
    }
  }
}
