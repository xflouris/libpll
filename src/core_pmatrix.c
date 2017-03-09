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

PLL_EXPORT int pll_core_update_pmatrix(double ** pmatrix,
                                       unsigned int states,
                                       unsigned int rate_cats,
                                       double * rates,
                                       const double * branch_lengths,
                                       const unsigned int * matrix_indices,
                                       const unsigned int * params_indices,
                                       double * prop_invar,
                                       double ** eigenvals,
                                       double ** eigenvecs,
                                       double ** inv_eigenvecs,
                                       unsigned int count,
                                       unsigned int attrib)
{
  unsigned int i,n,j,k,m;
  unsigned int states_padded = states;
  double * expd;
  double * temp;

  double pinvar;
  double * evecs;
  double * inv_evecs;
  double * evals;
  double * pmat;


  #ifdef HAVE_SSE3
  if (attrib & PLL_ATTRIB_ARCH_SSE)
  {
    if (states == 4)
    {
      return pll_core_update_pmatrix_4x4_sse(pmatrix,
                                             rate_cats,
                                             rates,
                                             branch_lengths,
                                             matrix_indices,
                                             params_indices,
                                             prop_invar,
                                             eigenvals,
                                             eigenvecs,
                                             inv_eigenvecs,
                                             count);
    }
    /* this line is never called, but should we disable the else case above,
       then states_padded must be set to this value */
    states_padded = (states+1) & 0xFFFFFFFE;
  }
  #endif
  #ifdef HAVE_AVX
  if (attrib & PLL_ATTRIB_ARCH_AVX)
  {
    if (states == 4)
    {
      return pll_core_update_pmatrix_4x4_avx(pmatrix,
                                             rate_cats,
                                             rates,
                                             branch_lengths,
                                             matrix_indices,
                                             params_indices,
                                             prop_invar,
                                             eigenvals,
                                             eigenvecs,
                                             inv_eigenvecs,
                                             count);
    }
    if (states == 20)
    {
      return pll_core_update_pmatrix_20x20_avx(pmatrix,
                                             rate_cats,
                                             rates,
                                             branch_lengths,
                                             matrix_indices,
                                             params_indices,
                                             prop_invar,
                                             eigenvals,
                                             eigenvecs,
                                             inv_eigenvecs,
                                             count);
    }
    /* this line is never called, but should we disable the else case above,
       then states_padded must be set to this value */
    states_padded = (states+3) & 0xFFFFFFFC;
  }
  #endif
  #ifdef HAVE_AVX2
  if (attrib & PLL_ATTRIB_ARCH_AVX2)
  {
    if (states == 4)
    {
      /* use AVX version here since FMA doesn't make much sense */
      return pll_core_update_pmatrix_4x4_avx(pmatrix,
                                             rate_cats,
                                             rates,
                                             branch_lengths,
                                             matrix_indices,
                                             params_indices,
                                             prop_invar,
                                             eigenvals,
                                             eigenvecs,
                                             inv_eigenvecs,
                                             count);
    }
    if (states == 20)
    {
      return pll_core_update_pmatrix_20x20_avx2(pmatrix,
                                             rate_cats,
                                             rates,
                                             branch_lengths,
                                             matrix_indices,
                                             params_indices,
                                             prop_invar,
                                             eigenvals,
                                             eigenvecs,
                                             inv_eigenvecs,
                                             count);
    }
    /* this line is never called, but should we disable the else case above,
       then states_padded must be set to this value */
    states_padded = (states+3) & 0xFFFFFFFC;
  }
  #endif

  expd = (double *)malloc(states * sizeof(double));
  temp = (double *)malloc(states*states*sizeof(double));

  if (!expd || !temp)
  {
    if (expd) free(expd);
    if (temp) free(temp);

    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf(pll_errmsg, 200, "Unable to allocate enough memory.");
    return PLL_FAILURE;
  }

  for (i = 0; i < count; ++i)
  {
    assert(branch_lengths[i] >= 0);

    /* compute effective pmatrix location */
    for (n = 0; n < rate_cats; ++n)
    {
      pmat = pmatrix[matrix_indices[i]] + n*states*states_padded;

      pinvar = prop_invar[params_indices[n]];
      evecs = eigenvecs[params_indices[n]];
      inv_evecs = inv_eigenvecs[params_indices[n]];
      evals = eigenvals[params_indices[n]];

      /* if branch length is zero then set the p-matrix to identity matrix */
      if (!branch_lengths[i])
      {
        for (j = 0; j < states; ++j)
          for (k = 0; k < states; ++k)
            pmat[j*states_padded + k] = (j == k) ? 1 : 0;
      }
      else
      {
        /* exponentiate eigenvalues */
        if (pinvar > PLL_MISC_EPSILON)
        {
          for (j = 0; j < states; ++j)
            expd[j] = exp(evals[j] * rates[n] * branch_lengths[i]
                                       / (1.0 - pinvar));
        }
        else
        {
          for (j = 0; j < states; ++j)
            expd[j] = exp(evals[j] * rates[n] * branch_lengths[i]);
        }

        /* check if all values of expd are approximately one */
        for (k=0, j=0; j < states; ++j)
          if ((expd[j] > PLL_ONE_MIN) && (expd[j] < PLL_ONE_MAX))
            ++k;

        /* if yes, it means we are multiplying the inverse eigenvectors matrix
           by the eigenvectors matrix, and essentially the resulting pmatrix is
           the identity matrix. This is done to prevent having numerical issues
           (negative entries in the pmatrix) which can occur due to the
           different floating point representations of one in expd */
        if (k == states)
        {
          /* set identity matrix */
          for (j = 0; j < states; ++j)
            for (k = 0; k < states; ++k)
              pmat[j*states_padded + k] = (j == k) ? 1 : 0;

          continue;
        }

        for (j = 0; j < states; ++j)
          for (k = 0; k < states; ++k)
            temp[j*states+k] = inv_evecs[j*states+k] * expd[k];

        for (j = 0; j < states; ++j)
        {
          for (k = 0; k < states; ++k)
          {
            pmat[j*states_padded+k] = 0;
            for (m = 0; m < states; ++m)
            {
              pmat[j*states_padded+k] +=
                  temp[j*states+m] * evecs[m*states+k];
            }
          }
        }
      }
      #ifdef DEBUG
      for (j = 0; j < states; ++j)
        for (k = 0; k < states; ++k)
          assert(pmat[j*states_padded+k] >= 0);
      #endif
    }
  }

  free(expd);
  free(temp);
  return PLL_SUCCESS;
}
