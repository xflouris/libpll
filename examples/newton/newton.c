/*
    Copyright (C) 2015 Tomas Flouri

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

#define RATES 4         /* number of rate categories */
#define STATES 4        /* number of states */
#define MAX_ITER 32     /* max iterations when optimizing branch lengths */
#define EPSILON 1e-5    /* threshold for detecting zero */


/* optimize the length of branch having parent and child as end-points */
static double newton(pll_partition_t * partition,
                     unsigned int parent_clv_index,
                     int parent_scaler_index,
                     unsigned int child_clv_index,
                     int child_scaler_index,
                     const unsigned int * params_indices,
                     double initial_length)
{
  int i;
  double d1, d2;
  double len;

  /* allocate space for the sumtable, for storing the constant part of the
     conditional probability vector computation that remains constant even if we change the
     branch length. This is done to avoid unnecessary computations when
     optimizing the branch length */
  double * sumtable = (double *)pll_aligned_alloc(partition->sites *
                                                  partition->rate_cats *
                                                  partition->states_padded *
                                                  sizeof(double),
                                                  PLL_ALIGNMENT_CPU);

  /* compute the sumtable for the particular branch once before proceeding with the
     optimization */
  pll_update_sumtable(partition,
                      parent_clv_index,
                      child_clv_index,
                      parent_scaler_index,
                      child_scaler_index,
                      params_indices,
                      sumtable);


  /* compute the derivatives of the likelihood function at most MAX_ITER times, */
  len = initial_length;
  for (i=0; i<MAX_ITER; ++i)
  {
    /* for the given branch length (len) compute the first (d1) and second (d2)
       derivatives of the likelihood function */
    double opt_logl = pll_compute_likelihood_derivatives(partition,
                                                         parent_scaler_index,
                                                         child_scaler_index,
                                                         len,
                                                         params_indices,
                                                         sumtable,
                                                         &d1,
                                                         &d2);

    printf("Branch length: %f log-L: %f Derivative: %f\n", len, opt_logl, d1);

    /* if the derivative is zero it means we reached a maximum and hence stop
       the computation */
    if (fabs(d1) < EPSILON) break;

    /* Newton's method for finding the optimum of a function. The iteration to
       reach the optimum is

       x_{i+1} = x_i - f'(x_i) / f''(x_i)

       where x_i is the current branch, f'(x_i) the first derivative and f''(x_i)
       the second derivative of the likelihood function */
    len -= d1/d2;
  }

  /* deallocate sumtable */
  pll_aligned_free(sumtable);

  /* return computed branch length */
  return len;
}

int main(int argc, char * argv[])
{
  unsigned int i;
  pll_partition_t * partition;
  pll_operation_t * operations;
  double alpha = 1.0;

  /* create the PLL partition instance */
  partition = pll_partition_create(4,       /* How many tip sequences do we have */
                                   2,       /* How many extra CLV buffers (apart from the tip sequences) should we allocate */
                                   STATES,  /* How many states do our data have */
                                   6,       /* How long are the tip sequences (number of sites) */
                                   1,       /* Number of different substitution models (or eigen decompositions) to use (i.e. 4 for LG4) */
                                   5,       /* How many probability matrices should we allocate */
                                   RATES,   /* Number of rate categories */
                                   2,       /* How many scale buffers do we want */
                                   PLL_ATTRIB_ARCH_CPU);        /* do not use vectorizations */

  /* initialize an array of two different branch lengths */
  double branch_lengths[5] = { 0.2, 0.4, 0.3, 0.5, 0.6};

  /* initialize an array of frequencies */
  double frequencies[4] = { 0.17, 0.19, 0.25, 0.39 };

  /* To be used together with branch_lengths to map branch lengths to
     probability matrices */
  unsigned int matrix_indices[5] = { 0, 1, 2, 3, 4};

  /* substitution rates for the GTR model (states*states-1)/2 */
  double subst_params[6] = {1,1,1,1,1,1};

  /* discretized category rates from a gamma distribution with alpha shape 1 */
  double rate_cats[RATES];
  pll_compute_gamma_cats(alpha, RATES, rate_cats, PLL_GAMMA_RATES_MEAN);

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

  /* update five probability matrices using the rate matrix with
     index 0. The i-th matrix (i ranges from 0 to matrix_count - 1) is
     generated using branch length branch_lengths[i] and rate matrix
     (substitution rates + frequencies) params_indices[i], and can be refered
     to with index matrix_indices[i] */

  unsigned int params_indices[RATES] = {0,0,0,0};

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
  operations = (pll_operation_t *)malloc(2 * sizeof(pll_operation_t));

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

  /* use the operations array to compute 2 CLVs. Operations will be carried out
     going from operation 0 to 1 */
  pll_update_partials(partition, operations, 2);

  /* print out the CLVs at tip and inner nodes*/
  printf ("Tip 0: ");
  pll_show_clv(partition,0,PLL_SCALE_BUFFER_NONE,7);
  printf ("Tip 1: ");
  pll_show_clv(partition,1,PLL_SCALE_BUFFER_NONE,7);
  printf ("Tip 2: ");
  pll_show_clv(partition,2,PLL_SCALE_BUFFER_NONE,7);
  printf ("Tip 3: ");
  pll_show_clv(partition,3,PLL_SCALE_BUFFER_NONE,7);

  printf ("CLV 4: ");
  pll_show_clv(partition,4,0,7);
  printf ("CLV 5: ");
  pll_show_clv(partition,5,1,7);

  /* compute the likelihood at the root of the rooted tree by specifying the CLV
     index of the root CLV and the index of the frequency vector to be used */

  double logl = pll_compute_edge_loglikelihood(partition,
                                               4,
                                               0,
                                               5,
                                               1,
                                               4,
                                               params_indices,
                                               NULL);

  printf("Log-L: %f\n\n", logl);

  printf("-*- Optimizing branch length -*-\n\n");

  /* optimize the branch length between nodes 4 and 5 */
  newton(partition, 4, 0, 5, 1, params_indices, branch_lengths[4]);

  /* we may now free the operations structure */
  free(operations);

  /* destroy all structures allocated for the concrete PLL partition instance */
  pll_partition_destroy(partition);

  return (0);
}
