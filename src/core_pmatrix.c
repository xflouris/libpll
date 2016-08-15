#include "pll.h"

PLL_EXPORT int pll_core_update_pmatrix(double * pmatrix,
                                       unsigned int states,
                                       double rate,
                                       double prop_invar,
                                       double branch_length,
                                       double * eigenvals,
                                       double * eigenvecs,
                                       double * inv_eigenvecs,
                                       unsigned int attrib)
{
  unsigned int j,k,m;
  unsigned int states_padded = states;
  double * expd;
  double * temp;

  expd = (double *)malloc(states * sizeof(double));
  temp = (double *)malloc(states*states* sizeof(double));
  if (!expd || !temp)
  {
    if (expd) free(expd);
    if (temp) free(temp);

    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf(pll_errmsg, 200, "Unable to allocate enough memory.");
    return PLL_FAILURE;
  }


  #ifdef HAVE_AVX
  if (attrib & PLL_ATTRIB_ARCH_AVX)
  {
    states_padded = (states+3) & 0xFFFFFFFC;
  }
  #endif

  /* if branch length is zero then set the p-matrix to identity matrix */
  if (!branch_length)
  {
    for (j = 0; j < states; ++j)
      for (k = 0; k < states; ++k)
        pmatrix[j*states_padded + k] = (j == k) ? 1 : 0;
  }
  else
  {
    /* exponentiate eigenvalues */
    for (j = 0; j < states; ++j)
      expd[j] = exp(eigenvals[j] * rate * branch_length
                                 / (1.0 - prop_invar));

    for (j = 0; j < states; ++j)
      for (k = 0; k < states; ++k)
        temp[j*states+k] = inv_eigenvecs[j*states+k] * expd[k];

    for (j = 0; j < states; ++j)
      for (k = 0; k < states; ++k)
      {
        pmatrix[j*states_padded+k] = 0;
        for (m = 0; m < states; ++m)
        {
          pmatrix[j*states_padded+k] +=
              temp[j*states+m] * eigenvecs[m*states+k];
        }
      }
  }
  free(expd);
  free(temp);
  return PLL_SUCCESS;
}
