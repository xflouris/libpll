/* Copyright (C) 2015 Diego Darriba, Tomas Flouri
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * Contact: Diego Darriba <Diego.Darriba@h-its.org>,
 * Heidelberg Institute for Theoretical Studies,
 * Schloss-Wolfsbrunnenweg 35, D-69118 Heidelberg, Germany
 */

/**
 * Example of heterotachous models
 * We define a sample tree with 4 tips (4 DNA sequences / 6 sites)
 *
 * Most of the stuff here was already commented in the "unrooted" example
 * Only the relevant parts are pointed out
 *
 * We use 3 different models. One for root branch, and one for each subtree.
 */
#include "pll.h"

static const unsigned int rmatrix_count = 3;

static void update_pmatrices(pll_partition_t * partition,
                             unsigned int * matrix_indices,
                             unsigned int * matrix_start,
                             unsigned int * matrix_count,
                             double * branch_lengths)
{
  unsigned int i,j;
  unsigned int params_indices[4];

  for (i=0; i<rmatrix_count; i++)
  {
    for (j=0; j < 4; ++j)
      params_indices[j] = i;
    pll_update_prob_matrices(partition,
                             params_indices,
                             matrix_indices + matrix_start[i],
                             branch_lengths + matrix_start[i],
                             matrix_count[i]);
  }
}

int main(int argc, char * argv[])
{
  unsigned int i, j;
  pll_partition_t * partition;
  pll_operation_t * operations;
  double alpha = 1.0;

  /* create the PLL partition instance */
  partition = pll_partition_create(4,             /* tips */
                                   2,             /* number of clv buffers */
                                   4,             /* states */
                                   6,             /* sites */
                                   rmatrix_count, /* number of different rate matrices */
                                   5,             /* number of p-matrices */
                                   4,             /* number of rate categories */
                                   2,             /* number of scale buffers */
                                   PLL_ATTRIB_ARCH_CPU);

  pll_set_tip_states(partition, 0, pll_map_nt, "WAAAAB");
  pll_set_tip_states(partition, 1, pll_map_nt, "CACACD");
  pll_set_tip_states(partition, 2, pll_map_nt, "AGGACA");
  pll_set_tip_states(partition, 3, pll_map_nt, "CGTAGT");

  double branch_lengths[5] = { 0.2, 0.4, 0.3, 0.5, 0.6};

  /* 3 sets of frequencies */
  double frequencies[3][4] = {
                               { 0.17, 0.19, 0.25, 0.39 },
                               { 0.17, 0.19, 0.25, 0.39 },
                               { 0.17, 0.19, 0.25, 0.39 }
                             };

  /* matrix indices and assignments. The plan is to create the first two
     transition probability matrices (p-matrices) using the first rate matrix,
     p-matrices 3 and 4 using the second rate matrix, and p-matrix 5 using
     the third rate matrix */
  unsigned int matrix_indices[5] = {0,1,2,3,4};
  unsigned int matrix_count[3]   = {2,2,1};
  unsigned int matrix_start[3]   = {0,2,4};

  /* 3 sets of substitution parameters */
  double subst_params[3][6] = {
                                {  1, 1, 1,  1, 1, 1},
                                {  1, 2, 1,  1, 2, 1},
                                { .5, 2, 1, .5, 2, 1}
                              };

  /* fixed rate categories */
  double rate_cats[4];
  pll_compute_gamma_cats(alpha, 4, rate_cats, PLL_GAMMA_RATES_MEAN);

  /* set frequencies to the parameter sets */
  pll_set_frequencies(partition, 0, frequencies[0]);
  pll_set_frequencies(partition, 1, frequencies[1]);
  pll_set_frequencies(partition, 2, frequencies[2]);

  /* set substitution parameters */
  pll_set_subst_params(partition, 0, subst_params[0]);
  pll_set_subst_params(partition, 1, subst_params[1]);
  pll_set_subst_params(partition, 2, subst_params[2]);

  /* set (single set of) rate categories */
  pll_set_category_rates(partition, rate_cats);

  /* update probability matrices for the corresponding branch lengths */
  update_pmatrices(partition,
                   matrix_indices,
                   matrix_start,
                   matrix_count,
                   branch_lengths);

  /* output the two probability matrices (for each rate category) on screen */

  for (i = 0; i < rmatrix_count; ++i)
  {
    for (j = matrix_start[i]; j < matrix_start[i]+matrix_count[i]; ++j)
    {
      printf ("P-matrix for model %d and branch length %.1f\n", i, branch_lengths[j]);
      pll_show_pmatrix(partition, j, 4);
      printf ("\n");
    }
  }


  /* create an operations array for specifying the traversal
     descriptor when computing the CLVs */
  operations = (pll_operation_t *)malloc(4 * sizeof(pll_operation_t));

  operations[0].parent_clv_index    = 4;
  operations[0].child1_clv_index    = 0;
  operations[0].child2_clv_index    = 1;
  operations[0].child1_matrix_index = 0;
  operations[0].child2_matrix_index = 1;
  operations[0].parent_scaler_index = 0;
  operations[0].child1_scaler_index = PLL_SCALE_BUFFER_NONE;
  operations[0].child2_scaler_index = PLL_SCALE_BUFFER_NONE;

  operations[1].parent_clv_index    = 5;
  operations[1].child1_clv_index    = 2;
  operations[1].child2_clv_index    = 3;
  operations[1].child1_matrix_index = 2;
  operations[1].child2_matrix_index = 3;
  operations[1].parent_scaler_index = 1;
  operations[1].child1_scaler_index = PLL_SCALE_BUFFER_NONE;
  operations[1].child2_scaler_index = PLL_SCALE_BUFFER_NONE;

  /* use the operations array to compute 4 CLVs. Operations will be carried out
     starting from operation 0 to 3 */
  pll_update_partials(partition, operations, 2);

  /* print out the CLVs at inner nodes*/
  printf ("CLV 4: ");
  pll_show_clv(partition,4,0,7);
  printf ("CLV 5: ");
  pll_show_clv(partition,5,1,7);

  unsigned int freqs_indices[4] = {2,2,2,2};
  double logl = pll_compute_edge_loglikelihood(partition,
                                               4,0, /* parent */
                                               5,1, /* child */
                                               4,   /* P-matrix */
                                               freqs_indices,   /* frequencies */
                                               NULL);

  printf("Log-L: %f\n", logl);

  pll_update_invariant_sites(partition);
  pll_update_invariant_sites_proportion(partition, 0, 0.5);

  /* we need to update the probability matrices after stating that we want
     to use invariant sites */
  update_pmatrices(partition,
                   matrix_indices,
                   matrix_start,
                   matrix_count,
                   branch_lengths);

  pll_update_partials(partition, operations, 2);
  logl = pll_compute_edge_loglikelihood(partition,4,0,5,1,4,freqs_indices, NULL);

  printf("Log-L (Inv+Gamma 0.5): %f\n", logl);

  pll_update_invariant_sites_proportion(partition, 0, 0.75);

  update_pmatrices(partition,
                   matrix_indices,
                   matrix_start,
                   matrix_count,
                   branch_lengths);

  pll_update_partials(partition, operations, 2);
  logl = pll_compute_edge_loglikelihood(partition,4,0,5,1,4,freqs_indices, NULL);

  printf("Log-L (Inv+Gamma 0.75): %f\n", logl);

  free(operations);
  pll_partition_destroy(partition);

  return (0);
}
