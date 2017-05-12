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

/*
    derivatives.c

    This test evaluates the computation of the likelihood derivatives at
    different branch lengths on a small tree and msa.

    The derivatives are computed twice at an inner edge and at a tip edge
    using 3 different alphas, 4 proportion of invariant sites, 3 sets of
    rate categories and 9 branches ranging from 0.1 to 90.
 */
#include "common.h"

#define NUM_ALPHAS   3
#define NUM_BRANCHES 9
#define NUM_CATS     3
#define NUM_PINV     4
#define N_STATES_NT  4

#define FLOAT_PRECISION 4

static double titv = 2.5;

static double alpha[NUM_ALPHAS] = {0.1, 0.75, 1.5};
static double pinvar[NUM_PINV] = {0.0, 0.3, 0.6, 0.9};
static unsigned int n_cat_gamma[NUM_CATS] = {1, 2, 4};
unsigned int params_indices[16] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
static double testbranches[NUM_BRANCHES] = {0.1, 0.2, 0.5, 0.9, 1.5, 5, 10, 50, 90};

int main(int argc, char * argv[])
{
  unsigned int i,j, k, b, p;
  double lk_scores[NUM_ALPHAS * NUM_CATS];
  double f, d_f, dd_f;
  unsigned int n_sites = 20;
  unsigned int n_tips = 5;
  pll_operation_t * operations;
  double * sumtable;

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

  /* additional operation for moving the root into the tip branch */

  operations[3].parent_clv_index    = 7;
  operations[3].child1_clv_index    = 6;
  operations[3].child2_clv_index    = 3;
  operations[3].child1_matrix_index = 0;
  operations[3].child2_matrix_index = 0;
  operations[3].parent_scaler_index = PLL_SCALE_BUFFER_NONE;
  operations[3].child1_scaler_index = PLL_SCALE_BUFFER_NONE;
  operations[3].child2_scaler_index = PLL_SCALE_BUFFER_NONE;

  /* check attributes */
  unsigned int attributes = get_attributes(argc, argv);

  for (k = 0; k < NUM_CATS; ++k) {
    pll_partition_t * partition;
    printf("FREE CREATE\n");
    partition = pll_partition_create(
                                n_tips,      /* numer of tips */
                                4,           /* clv buffers */
                                N_STATES_NT, /* number of states */
                                n_sites,     /* sequence length */
                                1,           /* different rate parameters */
                                2*n_tips-3,  /* probability matrices */
                                n_cat_gamma[k], /* gamma categories */
                                0,           /* scale buffers */
                                attributes
                                );          /* attributes */

    if (!partition)
    {
      printf("Fail creating partition");
      return(-1);
    }

    sumtable = pll_aligned_alloc(
      partition->sites * partition->rate_cats * partition->states_padded *
      sizeof(double), partition->alignment);

    if (!sumtable)
    {
      printf("Fail creating sumtable");
      pll_partition_destroy(partition);
      return(-1);
    }

    double branch_lengths[4] = { 0.1, 0.2, 0.3, 0.4};
    double frequencies[4] = { 0.3, 0.4, 0.1, 0.2 };
    unsigned int matrix_indices[4] = { 0, 1, 2, 3 };
    double subst_params[6] = {1,titv,1,1,titv,1};

    pll_set_frequencies(partition, 0, frequencies);
    pll_set_subst_params(partition, 0, subst_params);

    pll_set_tip_states(partition, 0, pll_map_nt, "WAACTCGCTA--ATTCTAAT");
    pll_set_tip_states(partition, 1, pll_map_nt, "CACCATGCTA--ATTGTCTT");
    pll_set_tip_states(partition, 2, pll_map_nt, "AG-C-TGCAG--CTTCTACT");
    pll_set_tip_states(partition, 3, pll_map_nt, "CGTCTTGCAA--AT-C-AAG");
    pll_set_tip_states(partition, 4, pll_map_nt, "CGACTTGCCA--AT-T-AAG");

    for (i = 0; i < NUM_ALPHAS; ++i) {
      for (p = 0; p < NUM_PINV; ++p) {

        printf("\n\n TEST alpha(ncats) = %6.2f(%2d) ; pinv = %.2f\n\n", alpha[i], n_cat_gamma[k], pinvar[p]);

        double * rate_cats = (double *) malloc(n_cat_gamma[k] * sizeof(double));

        if (pll_compute_gamma_cats(alpha[i], n_cat_gamma[k], rate_cats, PLL_GAMMA_RATES_MEAN) == PLL_FAILURE)
        {
          printf("Fail computing the gamma rates\n");
          continue;
        }

      	for (j=0; j<n_cat_gamma[k]; j++)
        {
      		printf("%f ", rate_cats[j]);
      	}
      	printf("\n");

        pll_set_category_rates(partition, rate_cats);
        free(rate_cats);

        for (j=0; j<partition->rate_matrices; ++j)
        {
          pll_update_invariant_sites_proportion(partition, j, pinvar[p]);
        }

        pll_update_prob_matrices(partition, params_indices, matrix_indices, branch_lengths, 4);
        pll_update_partials(partition, operations, 3);

        lk_scores[k*NUM_ALPHAS + i] = pll_compute_edge_loglikelihood(partition,
                                                           6,
                                                           PLL_SCALE_BUFFER_NONE,
                                                           7,
                                                           PLL_SCALE_BUFFER_NONE,
                                                           0,
                                                           params_indices,
                                                           NULL);

        pll_update_sumtable(partition, 6, 7,
                            PLL_SCALE_BUFFER_NONE, PLL_SCALE_BUFFER_NONE,
                            params_indices, sumtable);

        for (b = 0; b < NUM_BRANCHES; ++b) {
          if (!pll_compute_likelihood_derivatives(partition,
                                             PLL_SCALE_BUFFER_NONE,
                                             PLL_SCALE_BUFFER_NONE,
                                             testbranches[b],
                                             params_indices,
                                             sumtable,
                                             &d_f, &dd_f))
          {
            printf("Error computing likelihood derivatives\n");
            exit(1);
          }

          /* update logLikelihood */
          unsigned int pmatrix_index = 0;
          pll_update_prob_matrices(partition,
                                   params_indices,
                                   &pmatrix_index,
                                   &testbranches[b],
                                   1);
          f = pll_compute_edge_loglikelihood(partition,
                                             6,
                                             PLL_SCALE_BUFFER_NONE,
                                             7,
                                             PLL_SCALE_BUFFER_NONE,
                                             0,
                                             params_indices,
                                             NULL);

          printf("Branch %6.1f : %10.6f %12.4e %12.4e\n", testbranches[b], f, d_f, dd_f);
        }

        /* test original branch length */
        if (!pll_compute_likelihood_derivatives(partition,
                                           PLL_SCALE_BUFFER_NONE,
                                           PLL_SCALE_BUFFER_NONE,
                                           branch_lengths[0],
                                           params_indices,
                                           sumtable,
                                           &d_f, &dd_f))
        {
          printf("Error computing likelihood derivatives\n");
          exit(1);
        }

        /* update logLikelihood */
        unsigned int pmatrix_index = 0;
        pll_update_prob_matrices(partition,
                                 params_indices,
                                 &pmatrix_index,
                                 branch_lengths,
                                 1);
        f = pll_compute_edge_loglikelihood(partition,
                                           6,
                                           PLL_SCALE_BUFFER_NONE,
                                           7,
                                           PLL_SCALE_BUFFER_NONE,
                                           0,
                                           params_indices,
                                           NULL);

        printf("Test %10.6f = %10.6f\n", f, lk_scores[k*NUM_ALPHAS + i]);
        assert(fabs(f - lk_scores[k*NUM_ALPHAS + i]) < 1e-7);


        /* move to a tip branch */
        pll_update_partials(partition, operations + 3, 1);

        lk_scores[k*NUM_ALPHAS + i] = pll_compute_edge_loglikelihood(partition,
                                                           4,
                                                           PLL_SCALE_BUFFER_NONE,
                                                           7,
                                                           PLL_SCALE_BUFFER_NONE,
                                                           1,
                                                           params_indices,
                                                           NULL);

        pll_update_sumtable(partition, 4, 7,
                            PLL_SCALE_BUFFER_NONE, PLL_SCALE_BUFFER_NONE,
                            params_indices, sumtable);

        for (b = 0; b < NUM_BRANCHES; ++b) {
          if (!pll_compute_likelihood_derivatives(partition,
                                             PLL_SCALE_BUFFER_NONE,
                                             PLL_SCALE_BUFFER_NONE,
                                             testbranches[b],
                                             params_indices,
                                             sumtable,
                                             &d_f, &dd_f))
          {
            printf("Error computing likelihood derivatives\n");
            exit(1);
          }

          /* update logLikelihood */
          pmatrix_index = 1;
          pll_update_prob_matrices(partition,
                                   params_indices,
                                   &pmatrix_index,
                                   &testbranches[b],
                                   1);
          f = pll_compute_edge_loglikelihood(partition,
                                             4,
                                             PLL_SCALE_BUFFER_NONE,
                                             7,
                                             PLL_SCALE_BUFFER_NONE,
                                             1,
                                             params_indices,
                                             NULL);
          printf("Branch(Tip) %6.1f : %10.6f %12.4e %12.4e\n", testbranches[b], f, d_f, dd_f);
        }

        /* test original branch length */
        if (!pll_compute_likelihood_derivatives(partition,
                                           PLL_SCALE_BUFFER_NONE,
                                           PLL_SCALE_BUFFER_NONE,
                                           branch_lengths[1],
                                           params_indices,
                                           sumtable,
                                           &d_f, &dd_f))
        {
          printf("Error computing likelihood derivatives\n");
          exit(1);
        }

        /* update logLikelihood */
        pmatrix_index = 1;
        pll_update_prob_matrices(partition,
                                 params_indices,
                                 &pmatrix_index,
                                 &branch_lengths[1],
                                 1);
        f = pll_compute_edge_loglikelihood(partition,
                                           4,
                                           PLL_SCALE_BUFFER_NONE,
                                           7,
                                           PLL_SCALE_BUFFER_NONE,
                                           1,
                                           params_indices,
                                           NULL);

        printf("Test %10.6f = %10.6f\n", f, lk_scores[k*NUM_ALPHAS + i]);
        assert(fabs(f - lk_scores[k*NUM_ALPHAS + i]) < 1e-7);

      }
    }

    pll_aligned_free(sumtable);
    pll_partition_destroy(partition);
  }

    printf("\n");
    for (k = 0; k < NUM_CATS; ++k)
    {
      for (i = 0; i < NUM_ALPHAS; ++i)
      {
        printf("ti/tv:alpha(ncats) = %6.2f(%2d)   logL: %17.6f\n",
            alpha[i], n_cat_gamma[k], lk_scores[k*NUM_ALPHAS + i]);
      }
    }

    free(operations);

    return (0);
}
