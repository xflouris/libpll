/*
    Copyright (C) 2015 Diego Darriba, Tomas Flouri

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

#define NUM_TESTS 10
#define N_STATES_NT 4
#define N_CAT_GAMMA 4
#define FLOAT_PRECISION 4

static double titv[NUM_TESTS] = { 
    0.175, 1, 1.5, 2.25, 2.725, 4, 7.125, 8.19283745, 9.73647382, 10 
};

int main(int argc, char * argv[])
{
  unsigned int i,j;
  double lk_scores[NUM_TESTS];

  double alpha    = 1.0;

  unsigned int n_sites     = 20;
  unsigned int n_tips      = 5;
  unsigned int params_indices[N_CAT_GAMMA] = {0,0,0,0};

  pll_partition_t * partition;
  pll_operation_t * operations;
  partition = pll_partition_create(n_tips,
                                   4,           /* clv buffers */
                                   N_STATES_NT, /* number of states */
                                   n_sites,     /* sequence length */
                                   1,           /* different rate parameters */
                                   2*n_tips-3,  /* probability matrices */
                                   N_CAT_GAMMA, /* gamma categories */
                                   0,           /* scale buffers */
                                   PLL_ATTRIB_ARCH_AVX //| PLL_ATTRIB_PATTERN_TIP
                                   );          /* attributes */
  double branch_lengths[4] = { 0.1, 0.2, 1, 1};
  double frequencies[4] = { 0.3, 0.4, 0.1, 0.2 };
  unsigned int matrix_indices[4] = { 0, 1, 2, 3 };
  double subst_params[6] = {1,1,1,1,1,1};
  double rate_cats[4];

  pll_compute_gamma_cats(alpha, N_CAT_GAMMA, rate_cats, PLL_GAMMA_RATES_MEAN);

  pll_set_frequencies(partition, 0, frequencies);

  pll_set_category_rates(partition, rate_cats);

  pll_set_tip_states(partition, 0, pll_map_nt, "WAACTCGCTA--ATTCTAAT");
  pll_set_tip_states(partition, 1, pll_map_nt, "CACCATGCTA--ATTGTCTT");
  pll_set_tip_states(partition, 2, pll_map_nt, "AG-C-TGCAG--CTTCTACT");
  pll_set_tip_states(partition, 3, pll_map_nt, "CGTCTTGCAA--AT-C-AAG");
  pll_set_tip_states(partition, 4, pll_map_nt, "CGACTTGCCA--AT-T-AAG");

  operations = (pll_operation_t *)malloc(4* sizeof(pll_operation_t));

  operations[0].parent_clv_index    = 5;
  operations[0].child1_clv_index    = 0;
  operations[0].child2_clv_index    = 1;
  operations[0].child1_matrix_index = 1;
  operations[0].child2_matrix_index = 1;
  operations[0].parent_scaler_index = PLL_SCALE_BUFFER_NONE;
  operations[0].child1_scaler_index = PLL_SCALE_BUFFER_NONE;
  operations[0].child2_scaler_index = PLL_SCALE_BUFFER_NONE;

  operations[1].parent_clv_index    = 6;
  operations[1].child1_clv_index    = 5;
  operations[1].child2_clv_index    = 2;
  operations[1].child1_matrix_index = 0;
  operations[1].child2_matrix_index = 1;
  operations[1].parent_scaler_index = PLL_SCALE_BUFFER_NONE;
  operations[1].child1_scaler_index = PLL_SCALE_BUFFER_NONE;
  operations[1].child2_scaler_index = PLL_SCALE_BUFFER_NONE;

  operations[2].parent_clv_index    = 7;
  operations[2].child1_clv_index    = 3;
  operations[2].child2_clv_index    = 4;
  operations[2].child1_matrix_index = 1;
  operations[2].child2_matrix_index = 1;
  operations[2].parent_scaler_index = PLL_SCALE_BUFFER_NONE;
  operations[2].child1_scaler_index = PLL_SCALE_BUFFER_NONE;
  operations[2].child2_scaler_index = PLL_SCALE_BUFFER_NONE;

  for (i = 0; i < NUM_TESTS; ++i) 
  {
    subst_params[1] = subst_params[4] = titv[i];
    pll_set_subst_params(partition, 0, subst_params);

    pll_update_prob_matrices(partition,
                             params_indices,
                             matrix_indices,
                             branch_lengths,
                             4);

    pll_update_partials(partition, operations, 3);

    printf("\n\n TEST ti/tv = %.4f\n\n", titv[i]);

    for (j = 0; j < 4; ++j)
    {
      printf ("[%d] P-matrix for branch length %.4f\n", i, branch_lengths[j]);
      pll_show_pmatrix(partition, j, FLOAT_PRECISION);
      printf ("\n");
    }

    printf ("[%d] Tip 0: ", i);
    pll_show_clv(partition,0,PLL_SCALE_BUFFER_NONE,FLOAT_PRECISION+1);
    printf ("[%d] Tip 1: ", i);
    pll_show_clv(partition,1,PLL_SCALE_BUFFER_NONE,FLOAT_PRECISION+1);
    printf ("[%d] Tip 2: ", i);
    pll_show_clv(partition,2,PLL_SCALE_BUFFER_NONE,FLOAT_PRECISION+1);
    printf ("[%d] Tip 3: ", i);
    pll_show_clv(partition,3,PLL_SCALE_BUFFER_NONE,FLOAT_PRECISION+1);
    printf ("[%d] Tip 4: ", i);
    pll_show_clv(partition,4,PLL_SCALE_BUFFER_NONE,FLOAT_PRECISION+1);
    printf ("[%d] CLV 5: ", i);
    pll_show_clv(partition,5,PLL_SCALE_BUFFER_NONE,FLOAT_PRECISION+1);
    printf ("[%d] CLV 6: ", i);
    pll_show_clv(partition,6,PLL_SCALE_BUFFER_NONE,FLOAT_PRECISION+1);
    printf ("[%d] CLV 7: ", i);
    pll_show_clv(partition,7,PLL_SCALE_BUFFER_NONE,FLOAT_PRECISION+1);

    lk_scores[i] = pll_compute_edge_loglikelihood(partition,
                                                  6,
                                                  PLL_SCALE_BUFFER_NONE,
                                                  7,
                                                  PLL_SCALE_BUFFER_NONE,
                                                  0,
                                                  params_indices,
                                                  NULL);
  }

  printf("\n");
  for (i = 0; i < NUM_TESTS; ++i) 
  {
    printf("ti/tv: %14.4f      logL: %17.4f\n", titv[i], lk_scores[i]);
  }

  free(operations);
  pll_partition_destroy(partition);

  return (0);
}
