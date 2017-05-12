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
#include <string.h>
#include <assert.h>

#define ALPHA        0.5
#define N_STATES     4
#define N_RATE_CATS  4

#define N_TAXA_BIG    43
#define N_TAXA_SMALL   5
#define N_SITES      491

static int failtest(unsigned int attributes)
{
  pll_fasta_t * fp;
  fp = pll_fasta_open ("unexistent-file", pll_map_fasta);
  assert(!fp && pll_errno == PLL_ERROR_FILE_OPEN);

  return PLL_SUCCESS;
}

static int bigtest(unsigned int attributes)
{
  unsigned int i;
  char * seq, *header;
  long seq_len, header_len, seqno;
  pll_fasta_t * fp;
  pll_partition_t * partition;

  printf ("Creating PLL partition\n");

  partition = pll_partition_create(N_TAXA_BIG, /* tips */
                                    4, /* clv buffers */
                                    N_STATES, /* states */
                                    N_SITES, /* sites */
                                    1, /* different rate parameters */
                                    8, /* probability matrices */
                                    N_RATE_CATS, /* rate categories */
                                    1, 
                                    attributes
                                    );

  fp = pll_fasta_open ("testdata/worms16s.fas", pll_map_fasta);
  if (!fp)
  {
    printf (" ERROR opening file (%d): %s\n", pll_errno, pll_errmsg);
    exit (PLL_FAILURE);
  }

  i = 0;
  while (pll_fasta_getnext (fp, &header, &header_len, &seq, &seq_len, &seqno))
  {
    if (seq_len != (N_SITES + 1))
    {
      printf (
          " ERROR: Mismatching sequence length for sequence %d (%ld, and it should be %d)\n",
          i, seq_len - 1, N_SITES);
      exit (PLL_FAILURE);
    }
    if (!pll_set_tip_states (partition, i, pll_map_nt, seq))
    {
      printf (" ERROR setting states (%d): %s\n", pll_errno, pll_errmsg);
      exit (PLL_FAILURE);
    }
    printf ("Header of sequence %d(%ld) %s (%ld sites)\n", i, seqno, header,
            seq_len);
    free (header);
    free (seq);
    ++i;
  }

  if (pll_errno != PLL_ERROR_FILE_EOF)
  {
    printf (" ERROR at the end (%d): %s\n", pll_errno, pll_errmsg);
    exit (PLL_FAILURE);
  }

  assert(i == N_TAXA_BIG);

  pll_fasta_close (fp);
  pll_partition_destroy(partition);

  return PLL_SUCCESS;
}

static int smalltest (unsigned int attributes)
{
  unsigned int i;
  char * seq, *header;
  long seq_len, header_len, seqno;
  pll_fasta_t * fp;
  pll_partition_t * partition;
  pll_operation_t * operations;
  double rate_cats[4];
  unsigned int params_indices[N_RATE_CATS] = {0,0,0,0};

  unsigned int num_sites = 4 * N_SITES;

  double branch_lengths[4] =
    { 0.1, 0.2, 1, 1 };
  double frequencies[4] =
    { 0.1, 0.2, 0.3, 0.4 };
  unsigned int matrix_indices[4] = { 0, 1, 2, 3 };
  double subst_params[6] = { 1, 5, 1, 1, 5, 1 };

  partition = pll_partition_create(N_TAXA_SMALL, 4,
                                   N_STATES,
                                   num_sites, 1, 2 * N_TAXA_SMALL - 3,
                                   N_RATE_CATS,
                                   1, 
                                   attributes);

  fp = pll_fasta_open ("testdata/small.fas", pll_map_fasta);
  i = 0;
  while (pll_fasta_getnext (fp, &header, &header_len, &seq, &seq_len, &seqno))
  {
    if (!pll_set_tip_states (partition, i, pll_map_nt, seq))
      exit (PLL_FAILURE);
    free (header);
    free (seq);
    ++i;
  }
  assert(i == (N_TAXA_SMALL));

  operations = (pll_operation_t *) malloc (4 * sizeof(pll_operation_t));

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

  pll_compute_gamma_cats (ALPHA, N_RATE_CATS, rate_cats, PLL_GAMMA_RATES_MEAN);
  pll_set_subst_params (partition, 0, subst_params);
  pll_set_frequencies (partition, 0, frequencies);
  pll_set_category_rates (partition, rate_cats);
  pll_update_prob_matrices (partition, params_indices, matrix_indices, branch_lengths, 4);
  pll_update_partials (partition, operations, 3);

  printf ("logL: %17.6f\n", 
  pll_compute_edge_loglikelihood(partition,
                                 6,
                                 PLL_SCALE_BUFFER_NONE,
                                 7,
                                 PLL_SCALE_BUFFER_NONE,
                                 0,
                                 params_indices,
                                 NULL));

  free (operations);
  pll_fasta_close (fp);
  pll_partition_destroy(partition);

  return PLL_SUCCESS;
}

int main (int argc, char * argv[])
{
  unsigned int attributes = get_attributes(argc, argv);

  if (bigtest (attributes))
    printf ("Big test OK\n\n");

  if (smalltest (attributes))
    printf ("Small test OK\n\n");

  if (failtest (attributes))
    printf ("Fail test OK\n");

  return PLL_SUCCESS;
}
