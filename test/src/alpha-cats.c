/*
    Copyright (C) 2015 Diego Darriba

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

    Contact: Diego Darriba <Diego.Darriba@h-its.org>,
    Exelixis Lab, Heidelberg Instutute for Theoretical Studies
    Schloss-Wolfsbrunnenweg 35, D-69118 Heidelberg, Germany
*/
#include "pll.h"

#define NUM_ALPHAS 9
#define NUM_CATS   5
#define N_STATES_NT 4

#define FLOAT_PRECISION 5

double titv = 2.5;

  double alpha[NUM_ALPHAS] = {0.1, 0.5, 0.75, 1, 1.5, 5, 10, 50, 100};
  int n_cat_gamma[NUM_CATS] = {1, 2, 4, 8, 16};

int main(int argc, char * argv[])
{
  int i,j, k;
  double lk_scores[NUM_ALPHAS * NUM_CATS];
  int n_sites     = 20;
  int n_tips      = 5;
  pll_operation_t * operations;

  operations = (pll_operation_t *)malloc(4* sizeof(pll_operation_t));

  operations[0].parent_clv_index    = 5;
  operations[0].child1_clv_index    = 0;
  operations[0].child2_clv_index    = 1;
  operations[0].child1_matrix_index = 1;
  operations[0].child2_matrix_index = 1;

  operations[1].parent_clv_index    = 6;
  operations[1].child1_clv_index    = 5;
  operations[1].child2_clv_index    = 2;
  operations[1].child1_matrix_index = 0;
  operations[1].child2_matrix_index = 1;

  operations[2].parent_clv_index    = 7;
  operations[2].child1_clv_index    = 3;
  operations[2].child2_clv_index    = 4;
  operations[2].child1_matrix_index = 1;
  operations[2].child2_matrix_index = 1;

  for (k = 0; k < NUM_CATS; ++k) {
    pll_partition_t * partition;
    partition = pll_create_partition(n_tips,
                               4, /* clv buffers */
                               N_STATES_NT, /* number of states */
                               n_sites,     /* sequence length */
                               1,           /* different rate parameters */
                               2*n_tips-3,  /* probability matrices */
                               n_cat_gamma[k], /* gamma categories */
                               1,           /* scale buffers */
                               1);          /* attributes */

    if (partition == PLL_FAILURE) 
    {
      printf("Fail creating partition");
      continue;
    }

    double branch_lengths[4] = { 0.1, 0.2, 1, 1};
    double frequencies[4] = { 0.3, 0.4, 0.1, 0.2 };
    int matrix_indices[4] = { 0, 1, 2, 3 };
    double subst_params[6] = {1,titv,1,1,titv,1};

    pll_set_frequencies(partition, 0, frequencies);
    pll_set_subst_params(partition, 0, subst_params, 6);

    pll_set_tip_states(partition, 0, pll_map_nt, "WAACTCGCTA--ATTCTAAT");
    pll_set_tip_states(partition, 1, pll_map_nt, "CACCATGCTA--ATTGTCTT");
    pll_set_tip_states(partition, 2, pll_map_nt, "AG-C-TGCAG--CTTCTACT");
    pll_set_tip_states(partition, 3, pll_map_nt, "CGTCTTGCAA--AT-C-AAG");
    pll_set_tip_states(partition, 4, pll_map_nt, "CGACTTGCCA--AT-T-AAG");

    for (i = 0; i < NUM_ALPHAS; ++i) {

      printf("\n\n TEST alpha(ncats) = %6.2f(%2d)\n\n", alpha[i], n_cat_gamma[k]);

      double * rate_cats = (double *) malloc(n_cat_gamma[k] * sizeof(double));

      if (pll_compute_gamma_cats(alpha[i], n_cat_gamma[k], rate_cats) == PLL_FAILURE)
      {
        printf("Fail computing the gamma rates\n");
        continue;
      }

	for (j=0; j<n_cat_gamma[k]; j++) {
		printf("%f ", rate_cats[j]);
	}
	printf("\n");
      pll_set_category_rates(partition, rate_cats);
      free(rate_cats);

      pll_update_prob_matrices(partition, 0, matrix_indices, branch_lengths, 4);
      pll_update_partials(partition, operations, 3);

      for (j = 0; j < 4; ++j)
      {
        printf ("[%d] P-matrix for branch length %f\n", i, branch_lengths[j]);
        pll_show_pmatrix(partition, j, FLOAT_PRECISION);
        printf ("\n");
      }

      printf ("[%d] Tip 0: ", i);
      pll_show_clv(partition,0,FLOAT_PRECISION+1);
      printf ("[%d] Tip 1: ", i);
      pll_show_clv(partition,1,FLOAT_PRECISION+1);
      printf ("[%d] Tip 2: ", i);
      pll_show_clv(partition,2,FLOAT_PRECISION+1);
      printf ("[%d] Tip 3: ", i);
      pll_show_clv(partition,3,FLOAT_PRECISION+1);
      printf ("[%d] Tip 4: ", i);
      pll_show_clv(partition,4,FLOAT_PRECISION+1);
      printf ("[%d] CLV 5: ", i);
      pll_show_clv(partition,5,FLOAT_PRECISION+1);
      printf ("[%d] CLV 6: ", i);
      pll_show_clv(partition,6,FLOAT_PRECISION+1);
      printf ("[%d] CLV 7: ", i);
      pll_show_clv(partition,7,FLOAT_PRECISION+1);

      lk_scores[k*NUM_ALPHAS + i] = pll_compute_edge_loglikelihood(partition,6,7,0,0);
    }

    pll_destroy_partition(partition);
  }

    printf("\n");
    for (k = 0; k < NUM_CATS; ++k) 
    {
      for (i = 0; i < NUM_ALPHAS; ++i) 
      {
        printf("ti/tv:alpha(ncats) = %6.2f(%2d)   logL: %17.12f\n", alpha[i], n_cat_gamma[k], lk_scores[k*NUM_ALPHAS + i]);
      }
    }

    free(operations);

    return (0);
}
