#include "pll.h"

int main(int argc, char * argv[])
{
  unsigned int i;
  pll_partition_t * partition;
  pll_operation_t * operations;

  /* create the PLL partition instance */
  partition = pll_partition_create(5,       /* How many tip sequences do we have */
                                   4,       /* How many extra CLV buffers (apart from the tip sequences) should we allocate */
                                   4,       /* How many states do our data have */
                                   6,       /* How long are the tip sequences (number of sites) */
                                   1,       /* How many different substitution models (or eigen decompositions) do we want to use concurrently (i.e. 4 for LG4) */
                                   5,       /* How many probability matrices should we allocate */
                                   4,       /* Number of rate categories */
                                   4,       /* How many scale buffers do we want */
                                   PLL_ATTRIB_ARCH_AVX);        /* various attributes */

  /* initialize an array of two different branch lengths */
  double branch_lengths[5] = { 0.36, 0.722, 0.985, 0.718, 1.44};

  /* initialize an array of frequencies */
  double frequencies[4] = { 0.17, 0.19, 0.25, 0.39 };

  /* To be used together with branch_lengths to map branch lengths to
     probability matrices */
  unsigned int matrix_indices[5] = { 0, 1, 2, 3, 4};

  /* substitution rates for the GTR model */
  double subst_params[6] = {1,1,1,1,1,1};

  /* discretized category rates from a gamma distribution with alpha shape 1 */
  double rate_cats[4] = {0.13695378267140107,
                         0.47675185617665189,
                         0.99999999997958422,
                         2.38629436117236260};

  /* set frequencies */
  pll_set_frequencies(partition, 0, frequencies);

  /* set substitution parameters */
  pll_set_subst_params(partition, 0, subst_params);

  /* set rate categories */
  pll_set_category_rates(partition, rate_cats);

  /* set the 5 tip CLVs, and use the pll_map_nt map for converting
     the sequences to CLVs */
  pll_set_tip_states(partition, 0, pll_map_nt, "WAAAAB");
  pll_set_tip_states(partition, 1, pll_map_nt, "CACACD");
  pll_set_tip_states(partition, 2, pll_map_nt, "AGGACA");
  pll_set_tip_states(partition, 3, pll_map_nt, "CGTAGT");
  pll_set_tip_states(partition, 4, pll_map_nt, "CGAATT");

  /* update five probability matrices using the rate matrix with
     index 0. The i-th matrix (i ranges from 0 to matrix_count - 1) is
     generated using branch length branch_lengths[i] and rate matrix
     (substitution rates + frequencies) params_indices[i], and can be refered
     to with index matrix_indices[i] */

  unsigned int params_indices[4] = {0,0,0,0};

  pll_update_prob_matrices(partition,
                           params_indices,
                           matrix_indices,
                           branch_lengths,
                           5);

  /* output the two probability matrices (for each rate category) on screen */
  for (i = 0; i < 5; ++i)
  {
    printf ("P-matrix for branch length %f\n", branch_lengths[i]);
    pll_show_pmatrix(partition, i, 7);
    printf ("\n");
  }


  /* create an operations array for specifying the traversal
     descriptor when computing the CLVs */
  operations = (pll_operation_t *)malloc(4 * sizeof(pll_operation_t));

  operations[0].parent_clv_index    = 5;
  operations[0].parent_scaler_index = 0;
  operations[0].child1_clv_index    = 0;
  operations[0].child2_clv_index    = 1;
  operations[0].child1_matrix_index = 0;
  operations[0].child2_matrix_index = 0;
  operations[0].child1_scaler_index = PLL_SCALE_BUFFER_NONE;
  operations[0].child2_scaler_index = PLL_SCALE_BUFFER_NONE;

  operations[1].parent_clv_index    = 6;
  operations[1].parent_scaler_index = 1;
  operations[1].child1_clv_index    = 5;
  operations[1].child2_clv_index    = 2;
  operations[1].child1_matrix_index = 1;
  operations[1].child2_matrix_index = 2;
  operations[1].child1_scaler_index = 0;
  operations[1].child2_scaler_index = PLL_SCALE_BUFFER_NONE;

  operations[2].parent_clv_index    = 7;
  operations[2].parent_scaler_index = 2;
  operations[2].child1_clv_index    = 3;
  operations[2].child2_clv_index    = 4;
  operations[2].child1_matrix_index = 0;
  operations[2].child2_matrix_index = 0;
  operations[2].child1_scaler_index = PLL_SCALE_BUFFER_NONE;
  operations[2].child2_scaler_index = PLL_SCALE_BUFFER_NONE;

  operations[3].parent_clv_index    = 8;
  operations[3].parent_scaler_index = 3;
  operations[3].child1_clv_index    = 6;
  operations[3].child2_clv_index    = 7;
  operations[3].child1_matrix_index = 3;
  operations[3].child2_matrix_index = 4;
  operations[3].child1_scaler_index = 1;
  operations[3].child2_scaler_index = 2;

  /* use the operations array to compute 4 CLVs. Operations will be carried out
     starting from operation 0 to 3 */
  pll_update_partials(partition, operations, 4);

  /* print out the CLVs at tip and inner nodes*/
  if (!(partition->attributes & PLL_ATTRIB_PATTERN_TIP))
  {
    printf ("Tip 0: ");
    pll_show_clv(partition,0,PLL_SCALE_BUFFER_NONE,7);
    printf ("Tip 1: ");
    pll_show_clv(partition,1,PLL_SCALE_BUFFER_NONE,7);
    printf ("Tip 2: ");
    pll_show_clv(partition,2,PLL_SCALE_BUFFER_NONE,7);
    printf ("Tip 3: ");
    pll_show_clv(partition,3,PLL_SCALE_BUFFER_NONE,7);
    printf ("Tip 4: ");
    pll_show_clv(partition,4,PLL_SCALE_BUFFER_NONE,7);
  }
  printf ("CLV 5: ");
  pll_show_clv(partition,5,0,7);
  printf ("CLV 6: ");
  pll_show_clv(partition,6,1,7);
  printf ("CLV 7: ");
  pll_show_clv(partition,7,2,7);
  printf ("CLV 8: ");
  pll_show_clv(partition,8,3,7);

  /* compute the likelihood at the root of the rooted tree by specifying the CLV
     index of the root CLV and the index of the frequency vector to be used */
  double logl = pll_compute_root_loglikelihood(partition,8,3,params_indices,NULL);

  printf("Log-L: %f\n", logl);

  /* What if we want to consider invariant sites? Let's first update the
     partition with information about which sites are invariant */
  pll_update_invariant_sites(partition);

  /* Now let's set the log-likelihood proportion that
     invariant sites affect to 0.5 */
  pll_update_invariant_sites_proportion(partition, 0, 0.5);

  /* we need to update the probability matrices after stating that we want
     to use invariant sites */
  pll_update_prob_matrices(partition,
                           params_indices,
                           matrix_indices,
                           branch_lengths,
                           5);

  /* recompute the CLVs using the same traversal */
  pll_update_partials(partition, operations, 4);

  /* re-evaluate the log-likelihood */
  logl = pll_compute_root_loglikelihood(partition,8,3,params_indices,NULL);

  printf("Log-L (Inv+Gamma 0.5): %f\n", logl);

  /* Let's assume now we want to use a proportion of 0.75 for invariants. Since
     tip states haven't changed, we should only update the proportion and
     then update the probability matrices */
  pll_update_invariant_sites_proportion(partition, 0, 0.75);
  pll_update_prob_matrices(partition,
                           params_indices,
                           matrix_indices,
                           branch_lengths,
                           5);

  /* recompute the CLVs using the same traversal */
  pll_update_partials(partition, operations, 4);

  /* re-evaluate the log-likelihood */
  logl = pll_compute_root_loglikelihood(partition,8,3,params_indices,NULL);

  printf("Log-L (Inv+Gamma 0.75): %f\n", logl);

  /* we may now free the operations structure */
  free(operations);

  /* destroy all structures allocated for the concrete PLL partition instance */
  pll_partition_destroy(partition);

  return (0);
}
