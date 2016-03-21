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

static void fill_parent_scaler(unsigned int sites,
                               unsigned int * parent_scaler,
                               const unsigned int * left_scaler,
                               const unsigned int * right_scaler)
{
  unsigned int i;

  if (!left_scaler && !right_scaler)
    memset(parent_scaler, 0, sizeof(unsigned int) * sites);
  else if (left_scaler && right_scaler)
  {
    memcpy(parent_scaler, left_scaler, sizeof(unsigned int) * sites);
    for (i = 0; i < sites; ++i)
      parent_scaler[i] += right_scaler[i];
  }
  else
  {
    if (left_scaler)
      memcpy(parent_scaler, left_scaler, sizeof(unsigned int) * sites);
    else
      memcpy(parent_scaler, right_scaler, sizeof(unsigned int) * sites);
  }
}

void pll_update_partials_tiptip_avx(pll_partition_t * partition,
                                    const pll_operation_t * op)
{
  unsigned int j,k,n;
  unsigned int log2_maxstates = partition->log2_maxstates;
  unsigned int log2_states = partition->log2_states;
  unsigned int log2_rates = partition->log2_rates;
  unsigned int * parent_scaler;

  const char * left_tipchars = partition->tipchars[op->child1_clv_index];
  const char * right_tipchars = partition->tipchars[op->child2_clv_index];
  double * lh_statepair;

  double * parent_clv = partition->clv[op->parent_clv_index];

  unsigned int states = partition->states;
  unsigned int states_padded = partition->states_padded;
  unsigned int span = states * partition->rate_cats;
  unsigned int span_padded = states_padded * partition->rate_cats;

  /* check for dedicated functions */
  if (states == 4)
  {
    pll_core_update_partials_tt_4x4_avx(partition->sites,
                                        partition->rate_cats,
                                        partition->clv[op->parent_clv_index],
                                        (op->parent_scaler_index == PLL_SCALE_BUFFER_NONE) ?
                                          NULL : partition->scale_buffer[op->parent_scaler_index],
                                        left_tipchars,
                                        right_tipchars,
                                        partition->lh_statepair);
    return;

                        
  }

  if (op->parent_scaler_index != PLL_SCALE_BUFFER_NONE)
  {
    parent_scaler = partition->scale_buffer[op->parent_scaler_index];
    memset(parent_scaler, 0, sizeof(unsigned int) * partition->sites);
  }

  for (n = 0; n < partition->sites; ++n)
  {
    j = (unsigned int)(left_tipchars[n]);
    k = (unsigned int)(right_tipchars[n]);

    lh_statepair = partition->lh_statepair;
    lh_statepair += ((j << log2_maxstates) + k) << (log2_states+log2_rates);

    memcpy(parent_clv, lh_statepair, span*sizeof(double));

    parent_clv += span_padded;
  }
}

void pll_update_partials_tipinner_avx(pll_partition_t * partition,
                                      const pll_operation_t * op)
{
  int scaling;
  unsigned int i,j,k,n;
  unsigned int tip_clv_index;
  unsigned int inner_clv_index;
  unsigned int tip_matrix_index;
  unsigned int inner_matrix_index;
  unsigned int * inner_scaler;

  /* determine which child is the tip and mark it as the left (first) child */
  if (op->child1_clv_index < partition->tips)
  {
    tip_clv_index = op->child1_clv_index;
    tip_matrix_index = op->child1_matrix_index;
    inner_clv_index = op->child2_clv_index;
    inner_matrix_index = op->child2_matrix_index;
    if (op->child2_scaler_index == PLL_SCALE_BUFFER_NONE)
      inner_scaler = NULL;
    else
      inner_scaler = partition->scale_buffer[op->child2_scaler_index];
  }
  else
  {
    tip_clv_index = op->child2_clv_index;
    tip_matrix_index = op->child2_matrix_index;
    inner_clv_index = op->child1_clv_index;
    inner_matrix_index = op->child1_matrix_index;
    if (op->child1_scaler_index == PLL_SCALE_BUFFER_NONE)
      inner_scaler = NULL;
    else
      inner_scaler = partition->scale_buffer[op->child1_scaler_index];
  }

  double * parent_clv = partition->clv[op->parent_clv_index];
  unsigned int * scaler = (op->parent_scaler_index == PLL_SCALE_BUFFER_NONE) ?
                        NULL : partition->scale_buffer[op->parent_scaler_index];

  /* check for dedicated functions */
  if (partition->states == 4)
  {
    pll_core_update_partial_ti_4x4_avx(partition->sites,
                                       partition->rate_cats,
                                       parent_clv,
                                       scaler,
                                       partition->tipchars[tip_clv_index],
                                       partition->clv[inner_clv_index],
                                       partition->pmatrix[tip_matrix_index],
                                       partition->pmatrix[inner_matrix_index],
                                       inner_scaler);
    return;
  }

  const char * left_tipchars = partition->tipchars[tip_clv_index];
  const double * right_clv = partition->clv[inner_clv_index];
  const double * tip_matrix = partition->pmatrix[tip_matrix_index];
  const double * inner_matrix = partition->pmatrix[inner_matrix_index];

  
  const double * lmat;
  const double * rmat;

  unsigned int states = partition->states;
  unsigned int states_padded = partition->states_padded;
  unsigned int span = states * partition->rate_cats;
  unsigned int span_padded = states_padded * partition->rate_cats;

  unsigned int * revmap = partition->revmap;

  if (scaler)
    fill_parent_scaler(partition->sites, scaler, NULL, inner_scaler);

  for (n = 0; n < partition->sites; ++n)
  {
    lmat = tip_matrix;
    rmat = inner_matrix;

    scaling = (scaler) ? 1 : 0;

    for (k = 0; k < partition->rate_cats; ++k)
    {
      for (i = 0; i < states; ++i)
      {
        double terma = 0;
        double termb = 0;
        unsigned int lstate = revmap[(int)left_tipchars[n]];
        for (j = 0; j < states; ++j)
        {
          if (lstate & 1)
            terma += lmat[j];

          termb += rmat[j] * right_clv[j];

          lstate >>= 1;
        }
        parent_clv[i] = terma*termb;
        lmat += states_padded;
        rmat += states_padded;

        scaling = scaling && (parent_clv[i] < PLL_SCALE_THRESHOLD);
      }
      parent_clv += states_padded;
      right_clv  += states_padded;
    }
    /* if *all* entries of the site CLV were below the threshold then scale
       (all) entries by PLL_SCALE_FACTOR */
    if (scaling)
    {
      parent_clv -= span_padded;
      for (i = 0; i < span; ++i)
        parent_clv[i] *= PLL_SCALE_FACTOR;
      parent_clv += span_padded;
      scaler[n] += 1;
    }
  }
}

void pll_update_partials_avx(pll_partition_t * partition,
                             const pll_operation_t * op)
{
  const double * left_clv = partition->clv[op->child1_clv_index];
  const double * right_clv = partition->clv[op->child2_clv_index];
  const double * left_matrix = partition->pmatrix[op->child1_matrix_index];
  const double * right_matrix = partition->pmatrix[op->child2_matrix_index];
  double * parent_clv = partition->clv[op->parent_clv_index];
  unsigned int * scaler = (op->parent_scaler_index == PLL_SCALE_BUFFER_NONE) ?
                        NULL : partition->scale_buffer[op->parent_scaler_index];

  unsigned int * left_scaler;
  unsigned int * right_scaler;
  
  left_scaler = (op->child1_scaler_index == PLL_SCALE_BUFFER_NONE) ?
                  NULL : partition->scale_buffer[op->child1_scaler_index];
  right_scaler = (op->child2_scaler_index == PLL_SCALE_BUFFER_NONE) ?
                   NULL : partition->scale_buffer[op->child2_scaler_index];

  pll_core_update_partial_avx(partition->states,
                                partition->sites,
                                partition->rate_cats,
                                parent_clv,
                                scaler,
                                left_clv,
                                right_clv,
                                left_matrix,
                                right_matrix,
                                left_scaler,
                                right_scaler,
                                partition->attributes);
}
