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
#include "optimize.h"

#include <stdarg.h>
#include <search.h>
#include <time.h>

#define STATES    4
#define RATE_CATS 4

#define OPT_EPSILON 1e-2

static void set_missing_branch_length_recursive (pll_utree_t * tree,
                                                 double length)
{
  if (tree)
  {
    /* set branch length to default if not set */
    if (!tree->length)
      tree->length = length;

    if (tree->next)
    {
      if (!tree->next->length)
        tree->next->length = length;

      if (!tree->next->next->length)
        tree->next->next->length = length;

      set_missing_branch_length_recursive (tree->next->back, length);
      set_missing_branch_length_recursive (tree->next->next->back, length);
    }
  }
}

/* branch lengths not present in the newick file get a value of 0.000001 */
static void set_missing_branch_length (pll_utree_t * tree, double length)
{
  set_missing_branch_length_recursive (tree, 0.000001);
  set_missing_branch_length_recursive (tree->back, 0.000001);
}

static void fatal (const char * format, ...)
{
  va_list argptr;
  va_start(argptr, format);
  vfprintf (stderr, format, argptr);
  va_end(argptr);
  fprintf (stderr, "\n");
  exit (EXIT_FAILURE);
}

int main (int argc, char * argv[])
{

  int i, tip_count;
  pll_partition_t * partition;
  pll_operation_t * operations = NULL;
  double * branch_lengths = NULL;
  int * matrix_indices = NULL;
  opt_params params;
  time_t start_time, end_time;
  int parameters_to_optimize;

  pll_utree_t * tree;
  char ** tipnames;
  int * data;

  if (argc != 3)
    fatal (" syntax: %s [newick] [fasta]", argv[0]);

  tree = pll_parse_newick_utree (argv[1], &tip_count);

  /* fix all missing branch lengths (i.e. those that did not appear in the
   newick) to 0.00001 */
  set_missing_branch_length (tree, 0.00001);

  /*  obtain an array of pointers to tip names */
  tipnames = pll_query_utree_tipnames (tree, tip_count);

  /* create a libc hash table of size tip_count */
  hcreate (tip_count);

  /* populate a libc hash table with tree tip labels */
  data = (int *) malloc (tip_count * sizeof(int));
  for (i = 0; i < tip_count; ++i)
  {
    data[i] = i;
    ENTRY entry;
    entry.key = tipnames[i];
    entry.data = (void *) (data + i);
    hsearch (entry, ENTER);
  }

  /* open FASTA file */
  pll_fasta_t * fp = pll_fasta_open (argv[2], pll_map_fasta);
  if (!fp)
    fatal ("Error opening file %s", argv[2]);

  {
    char * seq = NULL;
    char * hdr = NULL;
    long seqlen;
    long hdrlen;
    long seqno;

    /* allocate arrays to store FASTA headers and sequences */
    char ** headers = (char **) calloc (tip_count, sizeof(char *));
    char ** seqdata = (char **) calloc (tip_count, sizeof(char *));

    /* read FASTA sequences and make sure they are all of the same length */
    int sites = -1;
    for (i = 0; pll_fasta_getnext (fp, &hdr, &hdrlen, &seq, &seqlen, &seqno);
        ++i)
    {
      if (i >= tip_count)
        fatal ("FASTA file contains more sequences than expected");

      if (sites != -1 && sites != seqlen)
      {
        fatal ("FASTA file does not contain equal size sequences\n");
      }

      if (sites == -1)
      {
        sites = seqlen;
      }

      headers[i] = hdr;
      seqdata[i] = seq;
    }

    /* did we stop reading the file because we reached EOF? */
    if (pll_errno != PLL_ERROR_FILE_EOF)
      fatal ("Error while reading file %s", argv[2]);

    /* close FASTA file */
    pll_fasta_close (fp);

    if (sites == -1)
      fatal ("Unable to read alignment");

    if (i != tip_count)
      fatal ("Some taxa are missing from FASTA file");

    /* create the PLL partition instance

     tip_count : the number of tip sequences we want to have
     tip_count-2 : the number of CLV buffers to be allocated for inner nodes
     STATES : the number of states that our data have
     1 : number of different substitution models (or eigen decomposition)
     to use concurrently (i.e. 4 for LG4)
     2*tip_count - 3: number of probability matrices to be allocated
     RATE_CATS : number of rate categories we will use
     tip_count-2 : how many scale buffers to use (not yet implemented)
     PLL_ATTRIB_ARCH_SSE : list of flags for hardware acceleration

     */

    partition = pll_create_partition (tip_count, tip_count - 2,
    STATES,
                                      sites, 1, 2 * tip_count - 3,
                                      RATE_CATS,
                                      tip_count - 2,
                                      PLL_ATTRIB_ARCH_SSE);

    /* find sequences in hash table and link them with the corresponding taxa */
    for (i = 0; i < tip_count; ++i)
    {
      ENTRY query;
      query.key = headers[i];
      ENTRY * found = NULL;

      found = hsearch (query, FIND);

      if (!found)
        fatal ("Sequence with header %s does not appear in the tree", hdr);

      int tip_clv_index = *((int *) (found->data));

      pll_set_tip_states (partition, tip_clv_index, pll_map_nt, seqdata[i]);
    }

    /* destroy hash table */
    hdestroy ();

    /* we no longer need these two arrays (keys and values of hash table... */
    free (data);
    free (tipnames);

    /* ...neither the sequences and the headers as they are already
     present in the form of probabilities in the tip CLVs */
    for (i = 0; i < tip_count; ++i)
    {
      free (seqdata[i]);
      free (headers[i]);
    }
    free (seqdata);
    free (headers);
  }
  int edge_pmatrix_index;
  int clv1;
  int clv2;

  pll_traverse_utree (tree, tip_count, &branch_lengths, &matrix_indices,
                      &operations, &edge_pmatrix_index, &clv1, &clv2);

  /* we will no longer need the tree structure */
  pll_destroy_utree (tree);

  /* initialize the array of base frequencies */
  double frequencies[4] =
    { 0.17, 0.19, 0.25, 0.39 };

  /* substitution rates for the 4x4 GTR model. This means we need exactly
   (4*4-4)/2 = 6 values, i.e. the number of elements above the diagonal */
  double subst_params[6] =
    { 1, 1, 1, 1, 1, 1 };

  /* we'll use 4 rate categories, and currently initialize them to 0 */
  double rate_cats[4] =
    { 0 };

  /* compute the discretized category rates from a gamma distribution
   with alpha shape 1 and store them in rate_cats  */
  pll_compute_gamma_cats (0.43, 4, rate_cats);

  /* set frequencies at model with index 0 (we currently have only one model) */
  pll_set_frequencies (partition, 0, frequencies);

  /* set 6 substitution parameters at model with index 0 */
  pll_set_subst_params (partition, 0, subst_params);

  /* set rate categories */
  pll_set_category_rates (partition, rate_cats);

  pll_update_prob_matrices (partition, 0, matrix_indices, branch_lengths,
                            2 * tip_count - 3);

  pll_update_partials (partition, operations, tip_count - 2);

  double logl = pll_compute_edge_loglikelihood (partition, clv1, clv2,
                                                edge_pmatrix_index, 0);
  printf ("Log-L: %f\n", logl);

  /* pll stuff */
  params.partition = partition;
  params.operations = operations;
  params.branch_lengths = branch_lengths;
  params.matrix_indices = matrix_indices;
  params.clv1 = clv1;
  params.clv2 = clv2;
  params.edge_pmatrix_index = edge_pmatrix_index;
  params.num_gamma_cats = RATE_CATS;
  params.params_index = 0;
  params.alpha_value = 1.0;

  /* optimization parameters */
  params.factr = 1e8;
  params.pgtol = 0.1;

  parameters_to_optimize = PARAM_SUBST_RATES | PARAM_ALPHA; // | PARAM_PINV; // | PARAM_BRANCH_LENGTHS

  start_time = time (NULL);
  logl *= -1;
  double cur_logl = logl + 10;
  while (fabs (cur_logl - logl) > OPT_EPSILON)
  {
    logl = cur_logl;

    if (parameters_to_optimize & PARAM_BRANCH_LENGTHS)
    {
      params.which_parameters = PARAM_BRANCH_LENGTHS;
      optimize_parameters (&params);
    }

    if (parameters_to_optimize & PARAM_SUBST_RATES)
    {
      params.which_parameters = PARAM_SUBST_RATES; //| PARAM_ALPHA;// | PARAM_PINV;
      optimize_parameters (&params);
    }

    if (parameters_to_optimize & PARAM_ALPHA)
    {
      params.which_parameters = PARAM_ALPHA;
      cur_logl = optimize_parameters (&params);
    }

    if (parameters_to_optimize & PARAM_PINV)
    {
      params.which_parameters = PARAM_PINV;
      cur_logl = optimize_parameters (&params);
    }

    printf ("  iter: %ld s. : %f\n", time (NULL) - start_time, cur_logl);
  }
  end_time = time (NULL);
  cur_logl *= -1;

  printf ("Final Log-L: %f\n", cur_logl);
  printf ("Time:  %ld s.\n", end_time - start_time);

  printf ("Alpha: %f\n", params.alpha_value);
  printf ("P-inv: %f\n", partition->prop_invar[0]);
  printf ("Rates: %f %f %f %f %f %f\n", partition->subst_params[0][0],
          partition->subst_params[0][1], partition->subst_params[0][2],
          partition->subst_params[0][3], partition->subst_params[0][4],
          partition->subst_params[0][5]);

  printf ("Final Log-L: %f\n", logl);

  /* CLEAN */
  /* destroy all structures allocated for the concrete PLL partition instance */
  pll_destroy_partition (partition);

  free (branch_lengths);
  free (matrix_indices);
  free (operations);

  return (0);
}
