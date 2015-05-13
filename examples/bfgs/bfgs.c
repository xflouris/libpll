#include "optimize.h"

int main(int argc, char * argv[])
{
  pll_partition_t * partition;
  pll_operation_t * operations;

  /* create the PLL partition instance */
  partition = pll_create_partition(5,       /* How many tip sequences do we have */
                                   4,       /* How many extra CLV buffers (apart from the tip sequences) should we allocate */
                                   4,       /* How many states do our data have */
                                   4*491,   /* How long are the tip sequences (number of sites) */
                                   1,       /* How many different substitution models (or eigen decompositions) do we want to use concurrently (i.e. 4 for LG4) */
                                   5,       /* How many probability matrices should we allocate */
                                   4,       /* Number of rate categories */
                                   1,       /* How many scale buffers do we want (not implemented currently) */
                                   PLL_ATTRIB_ARCH_SSE);        /* various attributes (not yet implemented) */
  
  /* initialize an array of two different branch lengths */
  double branch_lengths[5] = { 0.36, 0.722, 0.985, 0.718, 1.44};

  /* initialize an array of frequencies */
  double frequencies[4] = { 0.17, 0.19, 0.25, 0.39 };

  /* To be used together with branch_lengths to map branch lengths to 
     probability matrices */
  int matrix_indices[5] = { 0, 1, 2, 3, 4};

  /* substitution rates for the GTR model */
  double subst_params[6] = {1,1,1,1,1,1};
//      5.52541,
//      3.44044,
//      4.09842,
//      0.01000,
//     200.00000,
//      1.00000
//  };

  /* discretized category rates from a gamma distribution with alpha shape 1 */
  double rate_cats[4] = {0.13695378267140107,  
                         0.47675185617665189,  
                         0.99999999997958422,  
                         2.38629436117236260};

  /* set frequencies */
  pll_set_frequencies(partition, 0, frequencies);

  /* set substitution parameters */
  pll_set_subst_params(partition, 0, subst_params, 6);

  /* set rate categories */
  pll_set_category_rates(partition, rate_cats);

  /* set the 5 tip CLVs, and use the pll_map_nt map for converting
     the sequences to CLVs */
  pll_fasta_t * fp;
  int i;
  char * seq, *header;
  long seq_len, header_len, seqno;
  fp = pll_fasta_open ("small.fas", pll_map_fasta);
  i = 0;
  while (pll_fasta_getnext (fp, &header, &header_len, &seq, &seq_len, &seqno))
  {
    if (!pll_set_tip_states (partition, i, pll_map_nt, seq))
      exit (PLL_FAILURE);
    free (header);
    free (seq);
    ++i;
  }
  pll_fasta_close (fp);

  /* update two probability matrices for the corresponding branch lengths */
  pll_update_prob_matrices(partition, 0, matrix_indices, branch_lengths, 5);

  /* output the two probability matrices (for each rate category) on screen */
//  for (i = 0; i < 5; ++i)
//  {
//    printf ("P-matrix for branch length %f\n", branch_lengths[i]);
//    pll_show_pmatrix(partition, i);
//    printf ("\n");
//  }


  /* create an operations array for specifying the traversal
     descriptor when computing the CLVs */
  operations = (pll_operation_t *)malloc(4 * sizeof(pll_operation_t));

  operations[0].parent_clv_index    = 5;
  operations[0].child1_clv_index    = 0;
  operations[0].child2_clv_index    = 1;
  operations[0].child1_matrix_index = 0;
  operations[0].child2_matrix_index = 0;

  operations[1].parent_clv_index    = 6;
  operations[1].child1_clv_index    = 5;
  operations[1].child2_clv_index    = 2;
  operations[1].child1_matrix_index = 1;
  operations[1].child2_matrix_index = 2;

  operations[2].parent_clv_index    = 7;
  operations[2].child1_clv_index    = 3;
  operations[2].child2_clv_index    = 4;
  operations[2].child1_matrix_index = 0;
  operations[2].child2_matrix_index = 0;

  operations[3].parent_clv_index    = 8;
  operations[3].child1_clv_index    = 6;
  operations[3].child2_clv_index    = 7;
  operations[3].child1_matrix_index = 3;
  operations[3].child2_matrix_index = 4;

  /* use the operations array to compute 4 CLVs. Operations will be carried out
     starting from operation 0 to 3 */
  pll_update_partials(partition, operations, 4);

  /* compute the likelihood at the root of the rooted tree by specifying the CLV
     index of the root CLV and the index of the frequency vector to be used */
  double logl = pll_compute_root_loglikelihood(partition,8,0);

  printf("Log-L: %f\n", logl);


  /* DO BFGS OPTIMIZATION */

  pll_update_invariant_sites(partition);

  lk_stuff tree;
  tree.partition = partition;
  tree.params_index = 0;
  tree.operations = operations;
  tree.branch_lengths = branch_lengths;
  tree.matrix_indices = matrix_indices;
  tree.n_prob_matrices = 5;
  tree.n_operations = 4;
  tree.n_cat_gamma = 4;
  tree.alpha = 1;

  tree.is_rooted = 1;
  tree.clv_index1 = 8;
  tree.clv_index2 = -1;
  tree.branch_matrix_index = -1;
  tree.freqs_index = 0;


  logl *= -1;
  double newlnl = logl+1;

  while (fabs(newlnl-logl) > 1e-2)
  {
    logl = newlnl;

    tree.which_parameter = PARAM_SUBST_RATES;
    newlnl = optimize (&tree, 0.0001);

    printf("        %f\n", newlnl);

    tree.which_parameter = PARAM_ALPHA; //PARAM_SUBST_RATES;
    newlnl = optimize (&tree, 0.0001);

    printf("        %f\n", newlnl);

    tree.which_parameter = PARAM_PINV; //PARAM_SUBST_RATES;
    newlnl = optimize (&tree, 1.0);

    printf("-logL = %f\n", newlnl);
  //  break;
    if (newlnl > logl) break;
  }

  /* we may now free the operations structure */
  free(operations);

  /* destroy all structures allocated for the concrete PLL partition instance */
  pll_destroy_partition(partition);

  return (0);
}
