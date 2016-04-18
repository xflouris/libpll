/*
    Copyright (C) 2016 Diego Darriba

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

#include "pll_optimize.h"
#include "common.h"

int main(int argc, char * argv[])
{
  unsigned int i;
  pll_partition_t * partition;
  pll_operation_t * operations;
  double alpha = 0.841;
  double logl;

  unsigned int attributes = get_attributes(argc, argv);

  partition = pll_partition_create(3,       /* Tip CLVs */
                                   1,       /* Inner CLVs */
                                   4,       /* States */
                                   4,       /* Sequence length */
                                   1,       /* Models (sets of subst. params)*/
                                   3,       /* P matrices */
                                   4,       /* Rate categories */
                                   0,       /* Scale buffers */
                                   pll_map_nt,
                                   attributes);
  
  double branch_lengths[3] = { 0.105361, 0.166920, 0.166920 };
  double frequencies[4] = { 0.25, 0.25, 0.25, 0.25 };
  unsigned int matrix_indices[3] = { 0, 1, 2 };
  double subst_params[6] = {1.452176, 0.937951, 0.462880, 0.617729, 1.745312, 1.000000};
  double rate_cats[4];
  pll_compute_gamma_cats(alpha, 4, rate_cats);

  /* set */
  pll_set_frequencies(partition, 0, frequencies);
  pll_set_subst_params(partition, 0, subst_params);
  pll_set_category_rates(partition, rate_cats);

  if (attributes & PLL_ATTRIB_PATTERN_TIP)
  {
    skip_test();
  }
  else
  {
    double tip1[64] = {0.0000000000,0.0000000000,0.0000000000,1.0000000000,0.0000000000,0.0000000000,0.0000000000,1.0000000000,
                       0.0000000000,0.0000000000,0.0000000000,1.0000000000,0.0000000000,0.0000000000,0.0000000000,1.0000000000,
                       0.0000000000,0.0000000000,0.0000000000,1.0000000000,0.0000000000,0.0000000000,0.0000000000,1.0000000000,
                       0.0000000000,0.0000000000,0.0000000000,1.0000000000,0.0000000000,0.0000000000,0.0000000000,1.0000000000,
                       0.0000000000,0.0000000000,1.0000000000,0.0000000000,0.0000000000,0.0000000000,1.0000000000,0.0000000000,
                       0.0000000000,0.0000000000,1.0000000000,0.0000000000,0.0000000000,0.0000000000,1.0000000000,0.0000000000,
                       0.0000000000,0.0000000000,1.0000000000,0.0000000000,0.0000000000,0.0000000000,1.0000000000,0.0000000000,
                       0.0000000000,0.0000000000,1.0000000000,0.0000000000,0.0000000000,0.0000000000,1.0000000000,0.0000000000};
    double tip2[64] = {1.0000000000,0.0000000000,0.0000000000,0.0000000000,1.0000000000,0.0000000000,0.0000000000,0.0000000000,
                       1.0000000000,0.0000000000,0.0000000000,0.0000000000,1.0000000000,0.0000000000,0.0000000000,0.0000000000,
                       0.0000000000,0.0000000000,0.0000000000,1.0000000000,0.0000000000,0.0000000000,0.0000000000,1.0000000000,
                       0.0000000000,0.0000000000,0.0000000000,1.0000000000,0.0000000000,0.0000000000,0.0000000000,1.0000000000,
                       1.0000000000,0.0000000000,0.0000000000,0.0000000000,1.0000000000,0.0000000000,0.0000000000,0.0000000000,
                       1.0000000000,0.0000000000,0.0000000000,0.0000000000,1.0000000000,0.0000000000,0.0000000000,0.0000000000,
                       0.0000000000,0.0000000000,1.0000000000,0.0000000000,0.0000000000,0.0000000000,1.0000000000,0.0000000000,
                       0.0000000000,0.0000000000,1.0000000000,0.0000000000,0.0000000000,0.0000000000,1.0000000000,0.0000000000};
    double tip3[64] = {0.0187458510,0.0000024231,0.0000002543,0.0000000729,0.0182452872,0.0000026797,0.0000002672,0.0000000766,
                       0.0178965003,0.0000028695,0.0000002763,0.0000000793,0.0173815188,0.0000031672,0.0000002902,0.0000000832,
                       0.0000000285,0.0000040274,0.0000000481,0.0017162677,0.0000000329,0.0000045174,0.0000000543,0.0017287813,
                       0.0000000363,0.0000048804,0.0000000590,0.0017359452,0.0000000418,0.0000054510,0.0000000666,0.0017442835,
                       0.0000002234,0.0000001064,0.0171331364,0.0000005364,0.0000002357,0.0000001125,0.0167737124,0.0000005903,
                       0.0000002445,0.0000001169,0.0165200071,0.0000006303,0.0000002578,0.0000001236,0.0161407319,0.0000006936,
                       0.0000006741,0.0000006317,0.0001509298,0.0000005609,0.0000007244,0.0000006695,0.0001579430,0.0000006157,
                       0.0000007612,0.0000006968,0.0001627960,0.0000006564,0.0000008184,0.0000007388,0.0001699152,0.0000007205};
    pll_set_tip_clv(partition, 0, tip1);
    pll_set_tip_clv(partition, 1, tip2);
    pll_set_tip_clv(partition, 2, tip3);
  }

  pll_update_prob_matrices(partition, 0, matrix_indices, branch_lengths, 3);

  for (i = 0; i < 3; ++i)
  {
    printf ("P-matrix for branch length %f\n", branch_lengths[i]);
    pll_show_pmatrix(partition, i, 4);
    printf ("\n");
  }

  operations = (pll_operation_t *)malloc(1 * sizeof(pll_operation_t));

  operations[0].parent_clv_index    = 3;
  operations[0].child1_clv_index    = 0;
  operations[0].child2_clv_index    = 1;
  operations[0].child1_matrix_index = 0;
  operations[0].child2_matrix_index = 1;
  operations[0].parent_scaler_index = PLL_SCALE_BUFFER_NONE;
  operations[0].child1_scaler_index = PLL_SCALE_BUFFER_NONE;
  operations[0].child2_scaler_index = PLL_SCALE_BUFFER_NONE;

  pll_update_partials(partition, operations, 1);

  unsigned int params_indices[4] = {0,0,0,0};
  logl = pll_compute_edge_loglikelihood(partition,
                                               3,PLL_SCALE_BUFFER_NONE, /* parent clv/scaler */
                                               2,PLL_SCALE_BUFFER_NONE, /* child clv/scaler  */
                                               2,   /* P-matrix          */
                                               params_indices);

  printf("Initial Log-L: %.10f\n", logl);

  pll_utree_t * tree = (pll_utree_t *) malloc (6 * sizeof(pll_utree_t));
  for (i=0; i<6; ++i) tree[i].scaler_index = PLL_SCALE_BUFFER_NONE;
  tree[0].next = tree[1].next = tree[2].next = NULL;
  tree[3].next = tree + 4; tree[4].next = tree + 5; tree[5].next = tree + 3;
  tree[0].back = tree + 3; tree[3].back = tree;
  tree[1].back = tree + 4; tree[4].back = tree + 1;
  tree[2].back = tree + 5; tree[5].back = tree + 2;
  tree[0].clv_index = 0;tree[1].clv_index = 1;tree[2].clv_index = 2;
  tree[3].clv_index = tree[4].clv_index = tree[5].clv_index = 3;
  tree[0].pmatrix_index = tree[3].pmatrix_index = 0;
  tree[1].pmatrix_index = tree[4].pmatrix_index = 1;
  tree[2].pmatrix_index = tree[5].pmatrix_index = 2;
  tree[0].length = tree[3].length = branch_lengths[0];
  tree[1].length = tree[4].length = branch_lengths[1];
  tree[2].length = tree[5].length = branch_lengths[2];
  tree[0].label = "TIP 1";
  tree[1].label = "TIP 2";
  tree[2].label = "INNER";
  tree[3].label = tree[4].label = tree[5].label = NULL;

  char * newick = pll_utree_export_newick(tree);
  printf("Tree (reference): %s\n", newick);
  free(newick);

  double test_logl = pll_optimize_branch_lengths_local (partition,
                                     tree[2].back,
                                     params_indices,    /* params index */
                                     1e-4, /* tolerance    */
                                     1,    /* smoothings   */
                                     1,    /* radius       */
                                     1);   /* keep update  */

   logl = pll_compute_edge_loglikelihood(partition,
                                         3,PLL_SCALE_BUFFER_NONE, /* parent clv/scaler */
                                         2,PLL_SCALE_BUFFER_NONE, /* child clv/scaler  */
                                         2,   /* P-matrix          */
                                         params_indices);

   printf("\n");
   printf("-Log-L returned by BL-opt:       %.10f\n", test_logl);
   printf(" Log-L recomputed after BL-opt: %.10f\n", logl);

   newick = pll_utree_export_newick(tree);
   printf("Tree (optimized): %s\n", newick);
   free(newick);

   /* clean */
  free(operations);
  free(tree);
  pll_partition_destroy(partition);

  return (0);
}
