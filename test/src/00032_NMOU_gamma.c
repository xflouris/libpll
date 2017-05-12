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
#include "common.h"

#define N_ALPHAS  3
#define N_CATS    3
#define N_MODES   2
#define N_STATES  7
#define N_SUBST_PARAMS 21 // N_STATES*(N_STATES-1)/2
#define FLOAT_PRECISION 4

#define MODENAME(m) (m == PLL_GAMMA_RATES_MEAN ? "MEAN" : "MEDIAN")

static double alpha[N_ALPHAS] =
  { 0.1, 1.25, 100 };
static unsigned int n_cat_gamma[N_CATS] =
  { 1, 4, 6 };

static int modes[N_MODES] =
  { PLL_GAMMA_RATES_MEDIAN, PLL_GAMMA_RATES_MEAN };


/* odd map with 7 states (a,b,c,d,e,f,g) */
const unsigned int odd_map[256] =
  { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0x3f, 0, 0, 0x3f, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0x3f, 0, 0x01, 0x02, 0x04,
      0x08, 0x0c, 0x10, 0x20, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0x01, 0x02, 0x04, 0x08, 0x0c, 0x10, 0x20, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

int main (int argc, char * argv[])
{
  unsigned int i, j, k, m;
  double lk_scores[N_ALPHAS * N_CATS * N_MODES];
  unsigned int n_sites = 20;
  unsigned int n_tips = 5;
  pll_operation_t * operations;

  unsigned int attributes = get_attributes (argc, argv);

  operations = (pll_operation_t *) malloc (4 * sizeof(pll_operation_t));
  
  unsigned int params_indices[6] = {0,0,0,0,0,0};

  operations[0].parent_clv_index = 5;
  operations[0].child1_clv_index = 0;
  operations[0].child2_clv_index = 1;
  operations[0].child1_matrix_index = 1;
  operations[0].child2_matrix_index = 1;
  operations[0].parent_scaler_index = PLL_SCALE_BUFFER_NONE;
  operations[0].child1_scaler_index = PLL_SCALE_BUFFER_NONE;
  operations[0].child2_scaler_index = PLL_SCALE_BUFFER_NONE;

  operations[1].parent_clv_index = 6;
  operations[1].child1_clv_index = 5;
  operations[1].child2_clv_index = 2;
  operations[1].child1_matrix_index = 0;
  operations[1].child2_matrix_index = 1;
  operations[1].parent_scaler_index = PLL_SCALE_BUFFER_NONE;
  operations[1].child1_scaler_index = PLL_SCALE_BUFFER_NONE;
  operations[1].child2_scaler_index = PLL_SCALE_BUFFER_NONE;

  operations[2].parent_clv_index = 7;
  operations[2].child1_clv_index = 3;
  operations[2].child2_clv_index = 4;
  operations[2].child1_matrix_index = 1;
  operations[2].child2_matrix_index = 1;
  operations[2].parent_scaler_index = PLL_SCALE_BUFFER_NONE;
  operations[2].child1_scaler_index = PLL_SCALE_BUFFER_NONE;
  operations[2].child2_scaler_index = PLL_SCALE_BUFFER_NONE;

  /* test illegal alpha value */
  double invalid_alpha = 0;
  double * rate_cats = (double *) malloc (4 * sizeof(double));
  if (pll_compute_gamma_cats (invalid_alpha, 4, rate_cats, PLL_GAMMA_RATES_MEAN) == PLL_FAILURE)
  {
    if (pll_errno != PLL_ERROR_PARAM_INVALID)
      printf ("Error is %d instead of %d\n", pll_errno,
      PLL_ERROR_PARAM_INVALID);
  }
  else
  {
    printf ("Computing gamma rates for alpha = %f should have failed\n",
            invalid_alpha);
  }
  free (rate_cats);

  for (k = 0; k < N_CATS; ++k)
  {
    pll_partition_t * partition;
    partition = pll_partition_create (n_tips, /* numer of tips */
                                      4, /* clv buffers */
                                      N_STATES, /* number of states */
                                      n_sites, /* sequence length */
                                      1, /* different rate parameters */
                                      2 * n_tips - 3, /* probability matrices */
                                      n_cat_gamma[k], /* gamma categories */
                                      0, /* scale buffers */
                                      attributes); /* attributes */

    if (!partition)
    {
      printf ("Fail creating partition");
      continue;
    }

    double branch_lengths[4] =
      { 0.1, 0.2, 1, 1 };
    double frequencies[N_STATES] =
      { 0.12, 0.14, 0.13, 0.11, 0.15, 0.13, 0.12 };
    unsigned int matrix_indices[4] =
      { 0, 1, 2, 3 };
    double subst_params[N_SUBST_PARAMS] =
      { 0.5, 2, 3, 4, 5, 1.1, 1.2, 1.3, 1.4, 1.5, 2.1, 2.2, 2.3, 2.4, 2.5, 3.1,
          3.2, 3.3, 3.4, 3.5, 1 };

    printf ("Subst params: ");
    for (j = 0; j < N_SUBST_PARAMS; j++)
    {
      printf ("%8.5f ", subst_params[j]);
      if ((j % 7) == 6)
        printf ("\n              ");
    }
    printf ("\n");

    pll_set_frequencies (partition, 0, frequencies);
    pll_set_subst_params (partition, 0, subst_params);

    int check_states = 1;
    check_states &= pll_set_tip_states (partition, 0, odd_map,
                                        "BAACDCGCDA--AEECFAAD");
    check_states &= pll_set_tip_states (partition, 1, odd_map,
                                        "CACCABGCBA--BDDGFCDA");
    check_states &= pll_set_tip_states (partition, 2, odd_map,
                                        "AG-C-CGCAG--CGFCFACC");
    check_states &= pll_set_tip_states (partition, 3, odd_map,
                                        "CGDCBDGCAA--AB-C-AAG");
    check_states &= pll_set_tip_states (partition, 4, odd_map,
                                        "CGACFFGCCA--AF-D-AAG");

    if (!check_states)
    {
      printf ("Error: Wrong states\n");
      return (1);
    }

    for (i = 0; i < N_ALPHAS; ++i)
    {
      for (m = 0; m < N_MODES; ++m)
      {
        printf ("\n\n TEST alpha(ncats) = %6.2f(%2d), mode = %s\n\n",
                alpha[i], n_cat_gamma[k], MODENAME(m));

        double * rate_cats = (double *) malloc (n_cat_gamma[k] * sizeof(double));

        if (pll_compute_gamma_cats (alpha[i], n_cat_gamma[k],
                                    rate_cats, modes[m]) == PLL_FAILURE)
        {
          printf ("Fail computing the gamma rates\n");
          continue;
        }

        printf ("Rates: ");
        for (j = 0; j < n_cat_gamma[k]; j++)
        {
          printf ("%8.5f ", rate_cats[j]);
          if ((j % 4) == 3)
            printf ("\n       ");
        }
        printf ("\n");

        pll_set_category_rates (partition, rate_cats);
        free (rate_cats);

        pll_update_prob_matrices(partition,
                                 params_indices,
                                 matrix_indices,
                                 branch_lengths,
                                 4);

        pll_update_partials (partition, operations, 3);

        for (j = 0; j < 4; ++j)
        {
          printf ("[%d] P-matrix for branch length %f\n", i, branch_lengths[j]);
          pll_show_pmatrix (partition, j, FLOAT_PRECISION);
          printf ("\n");
        }

        printf ("[%d] CLV 5: ", i);
        pll_show_clv (partition, 5, PLL_SCALE_BUFFER_NONE, FLOAT_PRECISION + 1);
        printf ("[%d] CLV 6: ", i);
        pll_show_clv (partition, 6, PLL_SCALE_BUFFER_NONE, FLOAT_PRECISION + 1);
        printf ("[%d] CLV 7: ", i);
        pll_show_clv (partition, 7, PLL_SCALE_BUFFER_NONE, FLOAT_PRECISION);

        lk_scores[k * N_ALPHAS * N_MODES + i*N_MODES + m] = pll_compute_edge_loglikelihood (
            partition, 6,
            PLL_SCALE_BUFFER_NONE,
            7,
            PLL_SCALE_BUFFER_NONE,
            0, params_indices, NULL);
      }
    }
    pll_partition_destroy (partition);
  }

  printf ("\n");
  for (k = 0; k < N_CATS; ++k)
  {
    for (i = 0; i < N_ALPHAS; ++i)
    {
      for (m = 0; m < N_MODES; ++m)
      {
        printf ("ti/tv:alpha(ncats) = %6.2f(%2d), mode = %6s    logL: %17.6f\n",
                alpha[i], n_cat_gamma[k], MODENAME(modes[m]),
                lk_scores[k * N_ALPHAS * N_MODES + i * N_MODES + m]);
      }
    }
  }

  free (operations);

  return (0);
}
