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
                            double * parent_clv, 
                            double * left_clv, 
                            double * right_clv, 
                            double * left_matrix, 
                            double * right_matrix)
{
  int i,j,k,n;

  double * lmat;
  double * rmat;

  int states = partition->states;

  for (n = 0; n < partition->sites; ++n)
  {
    lmat = left_matrix;
    rmat = right_matrix;
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
      }
      parent_clv += states;
      left_clv   += states;
      right_clv  += states;
    }
  }
}

void pll_update_partials(pll_partition_t * partition, 
                         pll_operation_t * operations, 
                         int count)
{
  int i;

  
  for (i = 0; i < count; ++i)
  {
    update_partials(partition,
                    partition->clv[operations[i].parent_clv_index],
                    partition->clv[operations[i].child1_clv_index],
                    partition->clv[operations[i].child2_clv_index],
                    partition->pmatrix[operations[i].child1_matrix_index],
                    partition->pmatrix[operations[i].child2_matrix_index]);
  }
}

static double compute_invariant_likelihood (pll_partition_t * partition,
                                            int freqs_index,
                                            int site)
{
  double inv_site_lk = 0.0;

  if ((partition->prop_invar > 0.0)
      && (partition->invariant[site] != PLL_INVALID_STATE))
  {
      inv_site_lk =
        partition->
          frequencies[freqs_index][(size_t)partition->invariant[site]];
  }
  return inv_site_lk;
}

double pll_compute_root_loglikelihood(pll_partition_t * partition, 
                                      int clv_index, 
                                      int freqs_index)
{
  int i,j,k;

  double logl = 0;
  double term;
  double * clv = partition->clv[clv_index];
  int rates = partition->rate_cats;
  int states = partition->states;
  double * freqs = partition->frequencies[freqs_index];
  double site_lk, inv_site_lk;

  for (i = 0; i < partition->sites; ++i)
  {
    term = 0;
    for (j = 0; j < rates; ++j)
    {
      for (k = 0; k < states; ++k)
      {
        term += clv[k] * freqs[k];
      }
      clv += states;
    }

    site_lk = term / rates;

    if (partition->prop_invar > 0.0)
    {
      inv_site_lk = compute_invariant_likelihood (partition,
                                                  freqs_index,
                                                  i);
      site_lk = site_lk * (1. - partition->prop_invar)
          + inv_site_lk * partition->prop_invar;
    }

    logl += log (site_lk);
  }


  return logl;
}

double pll_compute_edge_loglikelihood(pll_partition_t * partition, 
                                      int parent_clv_index, 
                                      int child_clv_index, 
                                      int matrix_index,
                                      int freqs_index)
{
  int n,i,j,k;
  double logl = 0;
  double terma, termb;
  double site_lk, inv_site_lk;

  double * clvp = partition->clv[parent_clv_index];
  double * clvc = partition->clv[child_clv_index];
  double * freqs = partition->frequencies[freqs_index];
  double * pmatrix = partition->pmatrix[matrix_index];
  int states = partition->states;
  int rates = partition->rate_cats;

  for (n = 0; n < partition->sites; ++n)
  {
    pmatrix = partition->pmatrix[matrix_index];
    terma = 0;
    for (i = 0; i < partition->rate_cats; ++i)
    {
      for (j = 0; j < partition->states; ++j)
      {
        termb = 0;
        for (k = 0; k < partition->states; ++k)
          termb += pmatrix[k] * clvc[k];
        terma += clvp[j] * freqs[j] * termb;
        pmatrix += states;
      }
      clvp += states;
      clvc += states;
    }

    site_lk = terma / rates;

    if (partition->prop_invar > 0.0) {
      inv_site_lk = compute_invariant_likelihood(partition,
                                                 freqs_index,
                                                 n);
      site_lk = site_lk * (1. - partition->prop_invar) + inv_site_lk * partition->prop_invar;
    }

    logl += log(site_lk);
  }

  return logl;
}
