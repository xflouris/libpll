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

#define RATE_CATS 4

/* odd map with 5 states (a,b,c,d,e) */
const unsigned int odd_map[256] =
  { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0x1f, 0, 0, 0x1f, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0x1f, 0, 0x01, 0x02, 0x04,
    0x08, 0x0c, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0x01, 0x02, 0x04, 0x08, 0x0c, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

int main(int argc, char * argv[])
{
  unsigned int i;
  pll_partition_t * partition;
  pll_operation_t * operations;
  double rate_cats[RATE_CATS];
  double alpha = 0.841;
  double logl;
  double test_logl;
  char * newick;

  unsigned int params_indices[4] = {0, 0, 0, 0};
  unsigned int attributes = get_attributes(argc, argv);

  partition = pll_partition_create(3,         /* Tip CLVs */
                                   1,         /* Inner CLVs */
                                   5,         /* States */
                                   4,         /* Sequence length */
                                   1,         /* Models */
                                   3,         /* P matrices */
                                   RATE_CATS, /* Rate categories */
                                   0,         /* Scale buffers */
                                   odd_map,
                                   attributes);
  
  double branch_lengths[3] = { 0.105361, 0.166920, 0.166920 };
  double frequencies[5] = { 0.2, 0.2, 0.2, 0.2, 0.2 };
  unsigned int matrix_indices[3] = { 0, 1, 2 };
  double subst_params[10] = {1.452176, 0.937951, 0.462880, 0.617729, 1.745312, 0.937951, 0.462880, 0.617729, 1.745312, 1.000000};


  pll_compute_gamma_cats(alpha, RATE_CATS, rate_cats);

  /* set */
  pll_set_frequencies(partition, 0, frequencies);
  pll_set_subst_params(partition, 0, subst_params);
  pll_set_category_rates(partition, rate_cats);

  pll_set_tip_states (partition, 0, odd_map, "DABC");
  pll_set_tip_states (partition, 1, odd_map, "DAEC");
  pll_set_tip_states (partition, 2, odd_map, "DEEC");

  pll_update_prob_matrices(partition, 0, matrix_indices, branch_lengths, 3);

  for (i = 0; i < 3; ++i)
  {
    printf ("P-matrix for branch length %f\n", branch_lengths[i]);
    pll_show_pmatrix(partition, i, 4);
    printf ("\n");
  }

  printf("Update ops\n");
  operations = (pll_operation_t *)malloc(1 * sizeof(pll_operation_t));

  operations[0].parent_clv_index    = 3;
  operations[0].child1_clv_index    = 0;
  operations[0].child2_clv_index    = 1;
  operations[0].child1_matrix_index = 0;
  operations[0].child2_matrix_index = 1;
  operations[0].parent_scaler_index = PLL_SCALE_BUFFER_NONE;
  operations[0].child1_scaler_index = PLL_SCALE_BUFFER_NONE;
  operations[0].child2_scaler_index = PLL_SCALE_BUFFER_NONE;

  printf("Update partials\n");
  pll_update_partials(partition, operations, 1);

//  pll_show_clv(partition, 0, PLL_SCALE_BUFFER_NONE, 3);
//  pll_show_clv(partition, 1, PLL_SCALE_BUFFER_NONE, 3);
//  pll_show_clv(partition, 2, PLL_SCALE_BUFFER_NONE, 3);
//  for (i=0;i<2;++i)
//    printf("%p %f\n", partition->clv[3]+i, partition->clv[3][i]);
//  pll_show_clv(partition, 3, PLL_SCALE_BUFFER_NONE, 3);

  printf("Compute lnL\n");
  logl = pll_compute_edge_loglikelihood(partition,
                                               3,PLL_SCALE_BUFFER_NONE, /* parent clv/scaler */
                                               2,PLL_SCALE_BUFFER_NONE, /* child clv/scaler  */
                                               2,   /* P-matrix          */
                                               params_indices,
                                               NULL);

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
  tree[2].label = "TIP 3";
  tree[3].label = tree[4].label = tree[5].label = NULL;

  newick = pll_utree_export_newick(tree);
  printf("Tree (reference): %s\n", newick);
  free(newick);

  test_logl = pll_optimize_branch_lengths_local (partition,
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
                                         params_indices,
                                         NULL);  /* freqs index       */

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
