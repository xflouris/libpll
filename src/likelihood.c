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

static void update_partials(pll_partition_t * partition,
                            const pll_operation_t * op)
{
  unsigned int i,j,k,n;
  int scaling;

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
  unsigned int span = states * partition->rate_cats;

  for (n = 0; n < partition->sites; ++n)
  {
    lmat = left_matrix;
    rmat = right_matrix;

    scaling = (scaler) ? 1 : 0;

    for (k = 0; k < partition->rate_cats; ++k)
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
      scaler[n] += 1;
    }
  }
}

PLL_EXPORT void pll_update_partials(pll_partition_t * partition,
                                    const pll_operation_t * operations,
                                    unsigned int count)
{
  unsigned int i,j;

  unsigned int * parent_scaler;
  unsigned int * c1_scaler;
  unsigned int * c2_scaler;

  /* if the parent is supposed to have a scaler then initialized it as the
     element-wise sum of child1 and child2 scalers */
  for (i = 0; i < count; ++i)
  {
    if (operations[i].parent_scaler_index == PLL_SCALE_BUFFER_NONE)
    {
      parent_scaler = NULL;
    }
    else
    {
      /* get scaler for parent */
      parent_scaler = partition->scale_buffer[operations[i].parent_scaler_index];
      memset(parent_scaler, 0, sizeof(unsigned int) * partition->sites);

      /* if child1 has a scaler copy it to the parent */
      if (operations[i].child1_scaler_index != PLL_SCALE_BUFFER_NONE)
      {
        c1_scaler = partition->scale_buffer[operations[i].child1_scaler_index];
        memcpy(parent_scaler, c1_scaler, sizeof(unsigned int) * partition->sites);
      }

      /* if child2 has a scaler add its values to the parent scaler */
      if (operations[i].child2_scaler_index != PLL_SCALE_BUFFER_NONE)
      {
        c2_scaler = partition->scale_buffer[operations[i].child2_scaler_index];
        for (j = 0; j < partition->sites; ++j)
          parent_scaler[j] += c2_scaler[j];
      }
    }
    /* select vectorization method */
    #ifdef HAVE_AVX
    if (partition->attributes & PLL_ATTRIB_ARCH_AVX)
      pll_update_partials_avx(partition, &(operations[i]));
    else
    #endif
    #ifdef HAVE_SSE
    if (partition->attributes & PLL_ATTRIB_ARCH_SSE)
      assert(0);
    else
    #endif
      update_partials(partition, &(operations[i]));
  }

}

PLL_EXPORT double pll_compute_root_loglikelihood(pll_partition_t * partition,
                                                 unsigned int clv_index,
                                                 int scaler_index,
                                                 unsigned int freqs_index)
{
  unsigned int i,j,k;

  double logl = 0;
  double term, term_r;
  const double * clv = partition->clv[clv_index];
  unsigned int rates = partition->rate_cats;
  double * weights = partition->rate_weights;
  unsigned int states = partition->states;
  unsigned int states_padded = partition->states_padded;
  const double * freqs = NULL;
  double prop_invar = partition->prop_invar[freqs_index];
  double site_lk, inv_site_lk;
  unsigned int * scaler;

  /* get scaler array if specified */
  if (scaler_index == PLL_SCALE_BUFFER_NONE)
    scaler = NULL;
  else
    scaler = partition->scale_buffer[scaler_index];

  /* iterate through sites */
  for (i = 0; i < partition->sites; ++i)
  {
    freqs = partition->frequencies[freqs_index];
    term = 0;
    for (j = 0; j < rates; ++j)
    {
      term_r = 0;
      for (k = 0; k < states; ++k)
      {
        term_r += clv[k] * freqs[k];
      }
      term += term_r * weights[j];
      clv += states;
      if (partition->mixture > 1)
        freqs += states_padded;
    }

    site_lk = term;

    /* account for invariant sites */
    if (prop_invar > 0)
    {
      inv_site_lk = (partition->invariant[i] == -1) ?
                         0 : freqs[partition->invariant[i]];
      site_lk = site_lk * (1. - prop_invar)
          + inv_site_lk * prop_invar;
    }

    logl += log (site_lk) * partition->pattern_weights[i];

    /* scale log-likelihood of site if needed */
    if (scaler && scaler[i])
      logl += scaler[i] * log(PLL_SCALE_THRESHOLD);
  }

  return logl;
}

PLL_EXPORT double pll_compute_edge_persite_loglikelihood(
                                                 pll_partition_t * partition,
                                                 unsigned int parent_clv_index,
                                                 int parent_scaler_index,
                                                 unsigned int child_clv_index,
                                                 int child_scaler_index,
                                                 unsigned int matrix_index,
                                                 unsigned int freqs_index,
                                                 double * persite_lk)
{
  unsigned int n,i,j,k;
  double logl = 0;
  double terma_r, termb;
  double site_lk, inv_site_lk;

  const double * clvp = partition->clv[parent_clv_index];
  const double * clvc = partition->clv[child_clv_index];
  const double * freqs = NULL;
  const double * pmatrix = partition->pmatrix[matrix_index];
  double prop_invar = partition->prop_invar[freqs_index];
  unsigned int states = partition->states;
  unsigned int states_padded = partition->states_padded;
  double * weights = partition->rate_weights;
  unsigned int scale_factors;

  unsigned int * parent_scaler;
  unsigned int * child_scaler;

  assert(persite_lk);

  if (child_scaler_index == PLL_SCALE_BUFFER_NONE)
    child_scaler = NULL;
  else
    child_scaler = partition->scale_buffer[child_scaler_index];

  if (parent_scaler_index == PLL_SCALE_BUFFER_NONE)
    parent_scaler = NULL;
  else
    parent_scaler = partition->scale_buffer[parent_scaler_index];

  for (n = 0; n < partition->sites; ++n)
  {
    double sitecat_lk;
    site_lk = 0;
    pmatrix = partition->pmatrix[matrix_index];
    freqs = partition->frequencies[freqs_index];

    inv_site_lk = (!partition->invariant || partition->invariant[n] == -1) ?
                              0 : freqs[partition->invariant[n]];

    scale_factors = (parent_scaler) ? parent_scaler[n] : 0;
    scale_factors += (child_scaler) ? child_scaler[n] : 0;

    for (i = 0; i < partition->rate_cats; ++i)
    {
      terma_r = 0;
      for (j = 0; j < partition->states; ++j)
      {
        termb = 0;
        for (k = 0; k < partition->states; ++k)
          termb += pmatrix[k] * clvc[k];
        terma_r += clvp[j] * freqs[j] * termb;
        pmatrix += states;
      }
      sitecat_lk = terma_r * weights[i];

      /* account for invariant sites */
      sitecat_lk = sitecat_lk * (1. - prop_invar) +
                   inv_site_lk * prop_invar;

      persite_lk[n*partition->rate_cats + i] = sitecat_lk ;
      if (scale_factors)
        for (j=0; j<scale_factors; ++j)
          persite_lk[n*partition->rate_cats + i] *= PLL_SCALE_THRESHOLD;

      clvp += states;
      clvc += states;
      if (partition->mixture > 1)
        freqs += states_padded;

      site_lk += sitecat_lk;
    }

    logl += log(site_lk) * partition->pattern_weights[n];
    if (scale_factors)
      logl += scale_factors * log(PLL_SCALE_THRESHOLD);
  }

  return logl;
}

PLL_EXPORT double pll_compute_edge_loglikelihood(pll_partition_t * partition,
                                                 unsigned int parent_clv_index,
                                                 int parent_scaler_index,
                                                 unsigned int child_clv_index,
                                                 int child_scaler_index,
                                                 unsigned int matrix_index,
                                                 unsigned int freqs_index)
{
  unsigned int n,i,j,k;
  double logl = 0;
  double terma, terma_r, termb;
  double site_lk, inv_site_lk;

  const double * clvp = partition->clv[parent_clv_index];
  const double * clvc = partition->clv[child_clv_index];
  const double * freqs = NULL;
  const double * pmatrix = partition->pmatrix[matrix_index];
  double prop_invar = partition->prop_invar[freqs_index];
  unsigned int states = partition->states;
  unsigned int states_padded = partition->states_padded;
  double * weights = partition->rate_weights;
  unsigned int scale_factors;

  unsigned int * parent_scaler;
  unsigned int * child_scaler;

  if (child_scaler_index == PLL_SCALE_BUFFER_NONE)
    child_scaler = NULL;
  else
    child_scaler = partition->scale_buffer[child_scaler_index];

  if (parent_scaler_index == PLL_SCALE_BUFFER_NONE)
    parent_scaler = NULL;
  else
    parent_scaler = partition->scale_buffer[parent_scaler_index];

  for (n = 0; n < partition->sites; ++n)
  {
    pmatrix = partition->pmatrix[matrix_index];
    freqs = partition->frequencies[freqs_index];
    terma = 0;
    for (i = 0; i < partition->rate_cats; ++i)
    {
      terma_r = 0;
      for (j = 0; j < partition->states; ++j)
      {
        termb = 0;
        for (k = 0; k < partition->states; ++k)
          termb += pmatrix[k] * clvc[k];
        terma_r += clvp[j] * freqs[j] * termb;
        pmatrix += states;
      }
      terma += terma_r * weights[i];
      clvp += states;
      clvc += states;
      if (partition->mixture > 1)
        freqs += states_padded;
    }

    site_lk = terma;

    /* account for invariant sites */
    if (prop_invar > 0)
    {
      inv_site_lk = (partition->invariant[n] == -1) ? 
                        0 : freqs[partition->invariant[n]];

      site_lk = site_lk * (1. - prop_invar) +
                inv_site_lk * prop_invar;
    }

    scale_factors = (parent_scaler) ? parent_scaler[n] : 0;
    scale_factors += (child_scaler) ? child_scaler[n] : 0;

    logl += log(site_lk) * partition->pattern_weights[n];
    if (scale_factors)
      logl += scale_factors * log(PLL_SCALE_THRESHOLD);
  }

  return logl;
}

PLL_EXPORT int pll_update_sumtable(pll_partition_t * partition,
                                      unsigned int parent_clv_index,
                                      unsigned int child_clv_index,
                                      unsigned int params_index,
                                      unsigned int freqs_index,
                                      double *sumtable)
{
  unsigned int i,j,k,n;
  const double * clvp = partition->clv[parent_clv_index];
  const double * clvc = partition->clv[child_clv_index];
  const double * eigenvecs = partition->eigenvecs[params_index];
  const double * inv_eigenvecs = partition->inv_eigenvecs[params_index];
  const double * freqs = partition->frequencies[freqs_index];
  unsigned int states = partition->states;
  unsigned int n_rates = partition->rate_cats;
  unsigned int sites = partition->sites;

  /* so far not available for mixture models */
  assert(partition->mixture == 1);

  /* build sumtable */
  double * sum = sumtable;
  const double * t_clvp = clvp;
  const double * t_clvc = clvc;
  for (n = 0; n < sites; n++)
  {
    for (i = 0; i < n_rates; ++i)
    {
      for (j = 0; j < states; ++j)
      {
        double lefterm = 0;
        double righterm = 0;
        for (k = 0; k < states; ++k)
        {
          lefterm += t_clvp[k] * freqs[k] * inv_eigenvecs[k * states + j];
          righterm += eigenvecs[j * states + k] * t_clvc[k];
        }
        sum[j] = lefterm * righterm;
      }
      t_clvc += states;
      t_clvp += states;
      sum += states;
    }
  }

  return PLL_SUCCESS;
}

PLL_EXPORT double pll_compute_likelihood_derivatives(pll_partition_t * partition,
                                                 unsigned int parent_clv_index,
                                                 int parent_scaler_index,
                                                 unsigned int child_clv_index,
                                                 int child_scaler_index,
                                                 double branch_length,
                                                 unsigned int params_index,
                                                 unsigned int freqs_index,
                                                 double * sumtable,
                                                 double * d_f,
                                                 double * dd_f)
{
  unsigned int n, i, j;
  unsigned int sites = partition->sites;
  double site_lk[3];
  unsigned int states = partition->states;
  unsigned int n_rates = partition->rate_cats;
  const double * eigenvals = partition->eigenvals[params_index];
  const double * rates = partition->rates;
  const double * freqs = partition->frequencies[freqs_index];
  double prop_invar = partition->prop_invar[params_index];
  const double * sum;
  double logLK = 0.0;

  unsigned int * parent_scaler;
  unsigned int * child_scaler;
  double deriv1, deriv2;

  /* so far not available for mixture models */
  assert(partition->mixture == 1);

  if (child_scaler_index == PLL_SCALE_BUFFER_NONE)
     child_scaler = NULL;
   else
     child_scaler = partition->scale_buffer[child_scaler_index];

   if (parent_scaler_index == PLL_SCALE_BUFFER_NONE)
     parent_scaler = NULL;
   else
     parent_scaler = partition->scale_buffer[parent_scaler_index];

  *d_f = 0.0;
  *dd_f = 0.0;

  double
    diagptable[1024] = {0}, /* TODO make this dynamic */
    *diagp,
    ki;

  branch_length /= (1.0 - prop_invar);

  /* pre-compute the derivatives of the P matrix for all discrete GAMMA rates */
  diagp = diagptable;
  for(i = 0; i < n_rates; ++i)
  {
    ki = rates[i];
    for(j = 0; j < states; ++j)
    {
      diagp[0] = exp(eigenvals[j] * ki * branch_length);
      diagp[1] = eigenvals[j] * ki * diagp[0];
      diagp[2] = eigenvals[j] * ki * eigenvals[j] * ki * diagp[0];
      diagp += 3;
    }
  }

  sum = sumtable;
  for (n = 0; n < sites; ++n)
  {
    site_lk[0] = site_lk[1] = site_lk[2] = 0;
    diagp = diagptable;
    for (i = 0; i < n_rates; ++i)
    {
      ki = rates[i];
      for (j = 0; j < states; ++j)
      {

        site_lk[0] += sum[j] * diagp[0];
        site_lk[1] += sum[j] * diagp[1];
        site_lk[2] += sum[j] * diagp[2];
        diagp += 3;
      }
      sum += states;
    }
    site_lk[0] /= n_rates;
    site_lk[1] /= n_rates;
    site_lk[2] /= n_rates;

    double inv_site_lk = 0.0;

    /* account for invariant sites */
    if (prop_invar > 0)
    {
      inv_site_lk =
          (partition->invariant[n] == -1) ? 0 : freqs[partition->invariant[n]] * prop_invar;
      site_lk[0] = site_lk[0] * (1. - prop_invar) + inv_site_lk;
      site_lk[1] = site_lk[1] * (1. - prop_invar);
      site_lk[2] = site_lk[2] * (1. - prop_invar);
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

    logLK += log (site_lk[0]) * partition->pattern_weights[n];
    if (scale_factors)
    {
      logLK += scale_factors * log (PLL_SCALE_THRESHOLD);
    }

    /* build derivatives */
    deriv1 = (-site_lk[1] / site_lk[0]);
    deriv2 = (deriv1 * deriv1 - (site_lk[2] / site_lk[0]));
    *d_f += partition->pattern_weights[n] * deriv1;
    *dd_f += partition->pattern_weights[n] * deriv2;
  }

  return logLK;
}
