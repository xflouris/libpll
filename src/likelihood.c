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

static void case_tiptip(pll_partition_t * partition,
                        const pll_operation_t * op)
{
  const double * left_matrix = partition->pmatrix[op->child1_matrix_index];
  const double * right_matrix = partition->pmatrix[op->child2_matrix_index];
  double * parent_clv = partition->clv[op->parent_clv_index];
  unsigned int * parent_scaler;

  /* get parent scaler */
  if (op->parent_scaler_index == PLL_SCALE_BUFFER_NONE)
    parent_scaler = NULL;
  else
    parent_scaler = partition->scale_buffer[op->parent_scaler_index];

  /* precompute lookup table */
  pll_core_create_lookup(partition->states,
                         partition->rate_cats,
                         partition->ttlookup,
                         left_matrix,
                         right_matrix,
                         partition->tipmap,
                         partition->maxstates,
                         partition->attributes);
                         

  /* and update CLV at inner node */
  pll_core_update_partial_tt(partition->states,
                             partition->sites,
                             partition->rate_cats,
                             parent_clv,
                             parent_scaler,
                             partition->tipchars[op->child1_clv_index],
                             partition->tipchars[op->child2_clv_index],
                             partition->tipmap,
                             partition->maxstates,
                             partition->ttlookup,
                             partition->attributes);
}

static void case_tipinner(pll_partition_t * partition,
                          const pll_operation_t * op)
{
  double * parent_clv = partition->clv[op->parent_clv_index];
  unsigned int tip_clv_index;
  unsigned int inner_clv_index;
  unsigned int tip_matrix_index;
  unsigned int inner_matrix_index;
  unsigned int * right_scaler;
  unsigned int * parent_scaler;

  /* get parent scaler */
  if (op->parent_scaler_index == PLL_SCALE_BUFFER_NONE)
    parent_scaler = NULL;
  else
    parent_scaler = partition->scale_buffer[op->parent_scaler_index];

  /* find which of the two child nodes is the tip */
  if (op->child1_clv_index < partition->tips)
  {
    tip_clv_index = op->child1_clv_index;
    tip_matrix_index = op->child1_matrix_index;
    inner_clv_index = op->child2_clv_index;
    inner_matrix_index = op->child2_matrix_index;
    if (op->child2_scaler_index == PLL_SCALE_BUFFER_NONE)
      right_scaler = NULL;
    else
      right_scaler = partition->scale_buffer[op->child2_scaler_index];
  }
  else
  {
    tip_clv_index = op->child2_clv_index;
    tip_matrix_index = op->child2_matrix_index;
    inner_clv_index = op->child1_clv_index;
    inner_matrix_index = op->child1_matrix_index;
    if (op->child1_scaler_index == PLL_SCALE_BUFFER_NONE)
      right_scaler = NULL;
    else
      right_scaler = partition->scale_buffer[op->child1_scaler_index];
  }

  pll_core_update_partial_ti(partition->states,
                             partition->sites,
                             partition->rate_cats,
                             parent_clv,
                             parent_scaler,
                             partition->tipchars[tip_clv_index],
                             partition->clv[inner_clv_index],
                             partition->pmatrix[tip_matrix_index],
                             partition->pmatrix[inner_matrix_index],
                             right_scaler,
                             partition->tipmap,
                             partition->attributes);
}

static void case_innerinner(pll_partition_t * partition,
                            const pll_operation_t * op)
{
  const double * left_matrix = partition->pmatrix[op->child1_matrix_index];
  const double * right_matrix = partition->pmatrix[op->child2_matrix_index];
  double * parent_clv = partition->clv[op->parent_clv_index];
  double * left_clv = partition->clv[op->child1_clv_index];
  double * right_clv = partition->clv[op->child2_clv_index];
  unsigned int * parent_scaler;
  unsigned int * left_scaler;
  unsigned int * right_scaler;

  /* get parent scaler */
  if (op->parent_scaler_index == PLL_SCALE_BUFFER_NONE)
    parent_scaler = NULL;
  else
    parent_scaler = partition->scale_buffer[op->parent_scaler_index];

  if (op->child1_scaler_index != PLL_SCALE_BUFFER_NONE)
    left_scaler = partition->scale_buffer[op->child1_scaler_index];
  else
    left_scaler = NULL;

  /* if child2 has a scaler add its values to the parent scaler */
  if (op->child2_scaler_index != PLL_SCALE_BUFFER_NONE)
    right_scaler = partition->scale_buffer[op->child2_scaler_index];
  else
    right_scaler = NULL;

  pll_core_update_partial_ii(partition->states,
                             partition->sites,
                             partition->rate_cats,
                             parent_clv,
                             parent_scaler,
                             left_clv,
                             right_clv,
                             left_matrix,
                             right_matrix,
                             left_scaler,
                             right_scaler,
                             partition->attributes);
}

PLL_EXPORT void pll_update_partials(pll_partition_t * partition,
                                    const pll_operation_t * operations,
                                    unsigned int count)
{
  unsigned int i;
  const pll_operation_t * op;

  for (i = 0; i < count; ++i)
  { 
    op = &(operations[i]);

    if (partition->attributes & PLL_ATTRIB_PATTERN_TIP)
    {
      if ((op->child1_clv_index < partition->tips) &&
          (op->child2_clv_index < partition->tips))
      {
        /* tip-tip case */
        case_tiptip(partition,op);
      }
      else if ((operations[i].child1_clv_index < partition->tips) ||
               (operations[i].child2_clv_index < partition->tips))
      {
        /* tip-inner */
        case_tipinner(partition,op);
      }
      else
      {
        /* inner-inner */
        case_innerinner(partition,op);
      }
    }
    else
    {
      /* inner-inner */
      case_innerinner(partition,op);
    }
  }
}

PLL_EXPORT double pll_compute_root_loglikelihood(pll_partition_t * partition,
                                                 unsigned int clv_index,
                                                 int scaler_index,
                                                 unsigned int * freqs_index)
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
  double prop_invar = 0;
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
    term = 0;
    for (j = 0; j < rates; ++j)
    {
      freqs = partition->frequencies[freqs_index[j]];
      prop_invar = partition->prop_invar[freqs_index[j]];
      term_r = 0;
      for (k = 0; k < states; ++k)
      {
        term_r += clv[k] * freqs[k];
      }
      term += term_r * weights[j];
      clv += states_padded;
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
                                                 unsigned int * freqs_index,
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
  double prop_invar = 0;
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

    scale_factors = (parent_scaler) ? parent_scaler[n] : 0;
    scale_factors += (child_scaler) ? child_scaler[n] : 0;

    for (i = 0; i < partition->rate_cats; ++i)
    {
      prop_invar = partition->prop_invar[freqs_index[i]];
      freqs = partition->frequencies[freqs_index[i]];

      inv_site_lk = (!partition->invariant || partition->invariant[n] == -1) ?
                              0 : freqs[partition->invariant[n]];
      terma_r = 0;
      for (j = 0; j < states; ++j)
      {
        termb = 0;
        for (k = 0; k < states; ++k)
          termb += pmatrix[k] * clvc[k];
        terma_r += clvp[j] * freqs[j] * termb;
        pmatrix += states_padded;
      }
      sitecat_lk = terma_r * weights[i];

      /* account for invariant sites */
      sitecat_lk = sitecat_lk * (1. - prop_invar) +
                   inv_site_lk * prop_invar;

      persite_lk[n*partition->rate_cats + i] = sitecat_lk ;
      if (scale_factors)
        for (j=0; j<scale_factors; ++j)
          persite_lk[n*partition->rate_cats + i] *= PLL_SCALE_THRESHOLD;

      clvp += states_padded;
      clvc += states_padded;

      site_lk += sitecat_lk;
    }

    logl += log(site_lk) * partition->pattern_weights[n];
    if (scale_factors)
      logl += scale_factors * log(PLL_SCALE_THRESHOLD);
  }

  return logl;
}

static double edge_loglikelihood_tipinner(pll_partition_t * partition,
                                          unsigned int parent_clv_index,
                                          int parent_scaler_index,
                                          unsigned int child_clv_index,
                                          unsigned int matrix_index,
                                          unsigned int * freqs_index)
{
  unsigned int n,i,j,k;
  double logl = 0;
  double terma, terma_r, termb;
  double site_lk, inv_site_lk;

  const double * clvp = partition->clv[parent_clv_index];
  const double * freqs = NULL;
  const double * pmatrix = partition->pmatrix[matrix_index];
  double prop_invar = 0;
  unsigned int states = partition->states;
  unsigned int states_padded = partition->states_padded;
  double * weights = partition->rate_weights;
  unsigned int scale_factors;

  unsigned int * parent_scaler;

  /* child is the tip sequence, gets its tipchar */
  char * tipchar = partition->tipchars[child_clv_index];
  unsigned int * tipmap= partition->tipmap;
  unsigned int cstate;

  if (parent_scaler_index == PLL_SCALE_BUFFER_NONE)
    parent_scaler = NULL;
  else
    parent_scaler = partition->scale_buffer[parent_scaler_index];

  if (states == 4)
  {
    /* 4x4 case */
    for (n = 0; n < partition->sites; ++n)
    {
      pmatrix = partition->pmatrix[matrix_index];
      terma = 0;
      for (i = 0; i < partition->rate_cats; ++i)
      {
        freqs = partition->frequencies[freqs_index[i]];
        prop_invar = partition->prop_invar[freqs_index[i]];
        terma_r = 0;
        for (j = 0; j < states; ++j)
        {
          termb = 0;
          cstate = (unsigned int) (*tipchar);
          for (k = 0; k < states; ++k)
          {
            if (cstate & 1)
              termb += pmatrix[k];
            cstate >>= 1;
          }
          terma_r += clvp[j] * freqs[j] * termb;
          pmatrix += states_padded;
        }
      
        /* account for invariant sites */
        if (prop_invar > 0)
        {
          inv_site_lk = (partition->invariant[n] == -1) ? 
                            0 : freqs[partition->invariant[n]];
          terma += weights[i] * ( terma_r * (1. - prop_invar) +
                  inv_site_lk * prop_invar);
        }
        else
        {
          terma += terma_r * weights[i];
        }

        clvp += states_padded;
      }
    
      site_lk = terma;

      scale_factors = (parent_scaler) ? parent_scaler[n] : 0;

      logl += log(site_lk) * partition->pattern_weights[n];
      if (scale_factors)
        logl += scale_factors * log(PLL_SCALE_THRESHOLD);

      tipchar++;
    }
  }
  else
  {
    for (n = 0; n < partition->sites; ++n)
    {
      pmatrix = partition->pmatrix[matrix_index];
      terma = 0;
      for (i = 0; i < partition->rate_cats; ++i)
      {
        freqs = partition->frequencies[freqs_index[i]];
        prop_invar = partition->prop_invar[freqs_index[i]];
        terma_r = 0;
        for (j = 0; j < states; ++j)
        {
          termb = 0;
          cstate = tipmap[(int)(*tipchar)];
          for (k = 0; k < states; ++k)
          {
            if (cstate & 1)
              termb += pmatrix[k];
            cstate >>= 1;
          }
          terma_r += clvp[j] * freqs[j] * termb;
          pmatrix += states_padded;
        }

        /* account for invariant sites */
        if (prop_invar > 0)
        {
          inv_site_lk = (partition->invariant[n] == -1) ? 
                            0 : freqs[partition->invariant[n]];
          terma += weights[i] * ( terma_r * (1. - prop_invar) +
                  inv_site_lk * prop_invar);
        }
        else
        {
          terma += terma_r * weights[i];
        }

        clvp += states_padded;
      }
    
      site_lk = terma;

      scale_factors = (parent_scaler) ? parent_scaler[n] : 0;

      logl += log(site_lk) * partition->pattern_weights[n];
      if (scale_factors)
        logl += scale_factors * log(PLL_SCALE_THRESHOLD);

      tipchar++;
    }
  }

  return logl;
}

static double edge_loglikelihood(pll_partition_t * partition,
                                 unsigned int parent_clv_index,
                                 int parent_scaler_index,
                                 unsigned int child_clv_index,
                                 int child_scaler_index,
                                 unsigned int matrix_index,
                                 unsigned int * freqs_index)
{
  unsigned int n,i,j,k;
  double logl = 0;
  double terma, terma_r, termb;
  double site_lk, inv_site_lk;

  const double * clvp = partition->clv[parent_clv_index];
  const double * clvc = partition->clv[child_clv_index];
  const double * freqs = NULL;
  const double * pmatrix = partition->pmatrix[matrix_index];
  double prop_invar = 0;
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
    terma = 0;
    for (i = 0; i < partition->rate_cats; ++i)
    {
      freqs = partition->frequencies[freqs_index[i]];
      prop_invar = partition->prop_invar[freqs_index[i]];
      terma_r = 0;
      for (j = 0; j < states; ++j)
      {
        termb = 0;
        for (k = 0; k < states; ++k)
        {
          termb += pmatrix[k] * clvc[k];
        }
        terma_r += clvp[j] * freqs[j] * termb;
        pmatrix += states_padded;
      }

      /* account for invariant sites */
      if (prop_invar > 0)
      {
        inv_site_lk = (partition->invariant[n] == -1) ? 
                          0 : freqs[partition->invariant[n]];
        terma += weights[i] * ( terma_r * (1. - prop_invar) +
                inv_site_lk * prop_invar);
      }
      else
      {
        terma += terma_r * weights[i];
      }

      clvp += states_padded;
      clvc += states_padded;
    }

    site_lk = terma;

    scale_factors = (parent_scaler) ? parent_scaler[n] : 0;
    scale_factors += (child_scaler) ? child_scaler[n] : 0;

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
                                                 unsigned int * freqs_index)
{
  double logl;

  if (partition->attributes & PLL_ATTRIB_PATTERN_TIP)
    if ((parent_clv_index < partition->tips) ||
        (child_clv_index < partition->tips))
      return edge_loglikelihood_tipinner(partition,
                                        (parent_clv_index < partition->tips) ? 
                                            child_clv_index : parent_clv_index,
                                        (parent_clv_index < partition->tips) ? 
                                            child_scaler_index : parent_scaler_index,
                                        (parent_clv_index < partition->tips) ? 
                                            parent_clv_index : child_clv_index,
                                        matrix_index,
                                        freqs_index);

  logl = edge_loglikelihood(partition,
                            parent_clv_index,
                            parent_scaler_index,
                            child_clv_index,
                            child_scaler_index,
                            matrix_index,
                            freqs_index);
  return logl; 
}

static int update_sumtable_tipinner(const char * t_tipchars,
                                    const double * clv_inner,
                                    const double * eigenvecs,
                                    const double * inv_eigenvecs,
                                    const double * freqs,
                                    unsigned int states,
                                    unsigned int states_padded,
                                    unsigned int n_rates,
                                    unsigned int sites,
                                    unsigned int * tipmap,
                                    double *sumtable)
{
  unsigned int i, j, k, n;
  unsigned int tipstate;
  double lefterm  = 0;
  double righterm = 0;

  double * sum          = sumtable;
  const double * t_clvc = clv_inner;

  /* build sumtable */
  for (n = 0; n < sites; n++)
  {
    for (i = 0; i < n_rates; ++i)
    {
      for (j = 0; j < states; ++j)
      {
        tipstate = tipmap[(int)t_tipchars[n]];
        lefterm = 0;
        righterm = 0;
        for (k = 0; k < states; ++k)
        {
          lefterm += (tipstate & 1) * freqs[k] * inv_eigenvecs[k * states + j];
          righterm += eigenvecs[j * states + k] * t_clvc[k];
          tipstate >>= 1;
        }
        sum[j] = lefterm * righterm;
      }
      t_clvc += states_padded;
      sum += states_padded;
    }
  }

  return PLL_SUCCESS;
}

static int update_sumtable(const double * clvp,
                           const double * clvc,
                           const double * eigenvecs,
                           const double * inv_eigenvecs,
                           const double * freqs,
                           unsigned int states,
                           unsigned int states_padded,
                           unsigned int n_rates,
                           unsigned int sites,
                           double *sumtable)
{
  unsigned int i, j, k, n;
  double lefterm  = 0;
  double righterm = 0;

  double * sum          = sumtable;
  const double * t_clvp = clvp;
  const double * t_clvc = clvc;

  /* build sumtable */
  for (n = 0; n < sites; n++)
  {
    for (i = 0; i < n_rates; ++i)
    {
      for (j = 0; j < states; ++j)
      {
        lefterm = 0;
        righterm = 0;
        for (k = 0; k < states; ++k)
        {
          lefterm += t_clvp[k] * freqs[k] * inv_eigenvecs[k * states + j];
          righterm += eigenvecs[j * states + k] * t_clvc[k];
        }
        sum[j] = lefterm * righterm;
      }
      t_clvc += states_padded;
      t_clvp += states_padded;
      sum += states_padded;
    }
  }

  return PLL_SUCCESS;
}

/* computes the table containing the constant parts of the likelihood function
 * partial derivatives on the branch lengths.
 * sumtable: [output] must be allocated for storing (rates x states) values */
PLL_EXPORT int pll_update_sumtable(pll_partition_t * partition,
                                      unsigned int parent_clv_index,
                                      unsigned int child_clv_index,
                                      unsigned int params_index,
                                      unsigned int freqs_index,
                                      double *sumtable)
{
  int retval;
  const double * clvp          = partition->clv[parent_clv_index];
  const double * clvc          = partition->clv[child_clv_index];
  const double * eigenvecs     = partition->eigenvecs[params_index];
  const double * inv_eigenvecs = partition->inv_eigenvecs[params_index];
  const double * freqs         = partition->frequencies[params_index];
  unsigned int states          = partition->states;
  unsigned int states_padded   = partition->states_padded;
  unsigned int n_rates         = partition->rate_cats;
  unsigned int sites           = partition->sites;

  if ((partition->attributes & PLL_ATTRIB_PATTERN_TIP) &&
      ((parent_clv_index < partition->tips) ||
       (child_clv_index < partition->tips)))
  {
    if (parent_clv_index < partition->tips)
      retval = update_sumtable_tipinner(partition->tipchars[parent_clv_index],
                                        clvc, eigenvecs, inv_eigenvecs,
                                        freqs, states, states_padded,
                                        n_rates, sites, partition->tipmap,
                                        sumtable);
    else
      retval = update_sumtable_tipinner(partition->tipchars[child_clv_index],
                                        clvp, eigenvecs, inv_eigenvecs,
                                        freqs, states, states_padded,
                                        n_rates, sites, partition->tipmap,
                                        sumtable);
  }
  else
    retval = update_sumtable(clvp, clvc, eigenvecs, inv_eigenvecs,
                             freqs, states, states_padded, n_rates, sites,
                             sumtable);

  return retval;
}

/* Computes partial derivatives on the branch lengths.
 * branch_length: [input] value where the derivative is computed
 * sumtable: [input] must be computed at the edge where the derivatives will
 *                   be computed
 * d_f:  [output] first derivative
 * dd_f: [output] second derivative
 */
PLL_EXPORT double pll_compute_likelihood_derivatives(pll_partition_t * partition,
                                                 int parent_scaler_index,
                                                 int child_scaler_index,
                                                 double branch_length,
                                                 unsigned int params_index,
                                                 unsigned int freqs_index,
                                                 double * sumtable,
                                                 double * d_f,
                                                 double * dd_f)
{
  unsigned int n, i, j;
  double site_lk[3];
  unsigned int sites = partition->sites;
  unsigned int states = partition->states;
  unsigned int states_padded = partition->states_padded;
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
    *diagptable = (double *) calloc (n_rates * states * 3, sizeof(double)),
    *diagp,
    ki;

  if (!diagptable)
  {
    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf (pll_errmsg, 200, "Cannot allocate memory for diagptable");
    return -INFINITY;
  }

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
      sum += states_padded;
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

  free(diagptable);
  return logLK;
}
