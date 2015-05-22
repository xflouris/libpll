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
                            const double * left_clv,
                            const double * right_clv,
                            const double * left_matrix,
                            const double * right_matrix)
{
  int i,j,k,n;

  const double * lmat;
  const double * rmat;

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
                         const pll_operation_t * operations,
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

double pll_compute_root_loglikelihood(pll_partition_t * partition, 
                                      int clv_index, 
                                      int freqs_index)
{
  int i,j,k;

  double logl = 0;
  double term;
  const double * clv = partition->clv[clv_index];
  int rates = partition->rate_cats;
  int states = partition->states;
  const double * freqs = partition->frequencies[freqs_index];
  double prop_invar = partition->prop_invar[freqs_index];
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

    /* account for invariant sites */
    /* @MTH: This is not quite right, at least it does not cover
        all of the corner cases, if the library is intending
        to deal with partial ambiguity.
      It is possible for a DNA site to be all N's and Y's.
      This could be dealt with by expanding partition->frequencies
        to deal with ambiguity codes, and then allowing 
        partition->invariant to index that large state space.
      Not sure if it is worth it...
      Depending on the preprocessing, it may be necessary to deal with
        the all missing case too (those could be removed before scoring
        as their likelihood contribution is always 1.0)
    */
    if (prop_invar > 0)
    {
      inv_site_lk = (partition->invariant[i] == -1) ?
                         0 : freqs[partition->invariant[i]];
      site_lk = site_lk * (1. - prop_invar)
          + inv_site_lk * prop_invar;
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

  const double * clvp = partition->clv[parent_clv_index];
  const double * clvc = partition->clv[child_clv_index];
  const double * freqs = partition->frequencies[freqs_index];
  const double * pmatrix = partition->pmatrix[matrix_index];
  double prop_invar = partition->prop_invar[freqs_index];
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

    /* account for invariant sites */
    /* See @MTH comment above */
    if (prop_invar > 0)
    {
      inv_site_lk = (partition->invariant[n] == -1) ? 
                        0 : freqs[partition->invariant[n]];

      site_lk = site_lk * (1. - prop_invar) +
                inv_site_lk * prop_invar;
    }

    logl += log(site_lk);
  }

  return logl;
}
