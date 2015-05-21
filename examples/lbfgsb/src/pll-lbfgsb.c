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
#include "pll_optimize.h"

#include <assert.h>
#include <stdarg.h>
#include <search.h>
#include <time.h>

#define STATES    4
#define RATE_CATS 4
#define SUBST_PARAMS (STATES*(STATES-1)/2)

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

int * build_model_simmetries (const char * modelmatrix)
{
  int i, next_index;
  assert(strlen(modelmatrix) == SUBST_PARAMS);

  next_index = 0;
  int * subst_matrix_symmetries = (int *) calloc (SUBST_PARAMS, sizeof(int));
  int * t = (int *) alloca(10 * sizeof(int));
  for (i = 0; i < 10; i++)
    t[i] = -1;
  for (i = 0; i < SUBST_PARAMS; i++)
  {
    int v = (int) modelmatrix[i] - '0';
    if (v < 0 || v > 9)
      fatal ("Error in the model symmetries matrix");
    if (t[v] == -1)
    {
      t[v] = next_index;
      next_index++;
    }
    subst_matrix_symmetries[i] = t[v];
  }
  return subst_matrix_symmetries;
}

int read_fasta (pll_partition_t ** partition_ptr, char *file, int tip_count)
{

  int i;
  pll_partition_t * partition;

  /* open FASTA file */
  pll_fasta_t * fp = pll_fasta_open (file, pll_map_fasta);
  if (!fp)
    return PLL_FAILURE;
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
      {
        pll_errno = PLL_ERROR_TAXA_MISMATCH;
      }
      snprintf (pll_errmsg, 200,
                "FASTA file contains more sequences than expected");

      if (sites != -1 && sites != seqlen)
      {
        pll_errno = PLL_ERROR_SEQLEN_MISMATCH;
        snprintf (pll_errmsg, 200,
                  "FASTA file does not contain equal size sequences");
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
      return PLL_FAILURE;

    /* close FASTA file */
    pll_fasta_close (fp);

    if (sites == -1)
    {
      snprintf (pll_errmsg, 200, "Unable to read alignment");
      pll_errno = PLL_ERROR_ALIGN_UNREADABLE;
      return PLL_FAILURE;
    }

    if (i != tip_count)
    {
      snprintf (pll_errmsg, 200, "Some taxa are missing from FASTA file");
      pll_errno = PLL_ERROR_TAXA_MISMATCH;
      return PLL_FAILURE;
    }

    (*partition_ptr) = pll_create_partition (tip_count, tip_count - 2,
    STATES,
                                             sites, 1, 2 * tip_count - 3,
                                             RATE_CATS,
                                             tip_count - 2,
                                             PLL_ATTRIB_ARCH_SSE);
    partition = *partition_ptr;

    /* find sequences in hash table and link them with the corresponding taxa */
    for (i = 0; i < tip_count; ++i)
    {
      ENTRY query;
      query.key = headers[i];
      ENTRY * found = NULL;

      found = hsearch (query, FIND);

      if (!found)
      {
        snprintf (pll_errmsg, 200,
                  "Sequence with header %s does not appear in the tree", hdr);
        pll_errno = PLL_ERROR_TAXA_MISMATCH;
        return PLL_FAILURE;
      }

      int tip_clv_index = *((int *) (found->data));

      pll_set_tip_states (partition, tip_clv_index, pll_map_nt, seqdata[i]);
    }

    /* destroy hash table */
    hdestroy ();

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

  return PLL_SUCCESS;
}

int main (int argc, char * argv[])
{

  int i, tip_count;
  pll_partition_t * partition;
  pll_operation_t * operations = NULL;
  double * branch_lengths = NULL;
  int * matrix_indices = NULL;
  pll_optimize_options_t params;
  time_t start_time, end_time;
  int parameters_to_optimize;
  int edge_pmatrix_index;
  int clv1;
  int clv2;

  pll_utree_t * tree;
  char ** tipnames;
  int * data;
  int * subst_params_symmetries;

  if (argc != 4)
    fatal (" syntax: %s [newick] [fasta] [model]", argv[0]);

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

  partition = pll_create_partition_fasta(argv[2], STATES, 1, RATE_CATS,
                                         PLL_ATTRIB_ARCH_SSE, 0,
                                         tip_count, 0);
  if (!partition)
    fatal("Error %d: %s\n", pll_errno, pll_errmsg);

  /* destroy hash table */
  hdestroy ();

  /* we no longer need these two arrays (keys and values of hash table... */
  free (data);
  free (tipnames);

  subst_params_symmetries = build_model_simmetries(argv[3]);
  printf("Model: ");
  for (i=0; i<SUBST_PARAMS; i++)
    printf("%d", subst_params_symmetries[i]);
  printf("\n");

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
  params.lk_params.partition = partition;
  params.lk_params.operations = operations;
  params.lk_params.branch_lengths = branch_lengths;
  params.lk_params.matrix_indices = matrix_indices;
  params.lk_params.alpha_value = 1.0;
  params.lk_params.freqs_index = 0;
  params.lk_params.rooted = 0;
  params.lk_params.where.unrooted_t.parent_clv_index = clv1;
  params.lk_params.where.unrooted_t.child_clv_index = clv2;
  params.lk_params.where.unrooted_t.edge_pmatrix_index = edge_pmatrix_index;

  /* optimization parameters */
  params.params_index = 0;
  params.subst_params_symmetries = subst_params_symmetries;
  params.factr = 1e8;
  params.pgtol = 0.01;

  parameters_to_optimize =
      PLL_PARAMETER_SUBST_RATES |
      PLL_PARAMETER_ALPHA;
  // | PLL_PARAMETER_PINV;
  // | PLL_PARAMETER_BRANCH_LENGTHS

  start_time = time (NULL);
  logl *= -1;
  double cur_logl = logl + 10;
  while (fabs (cur_logl - logl) > OPT_EPSILON)
  {
    logl = cur_logl;

    if (parameters_to_optimize & PLL_PARAMETER_BRANCH_LENGTHS)
    {
      params.which_parameters = PLL_PARAMETER_BRANCH_LENGTHS;
      cur_logl = pll_optimize_parameters_lbfgsb (&params);
      printf ("  [branches]: %ld s. : %f\n", time (NULL) - start_time, cur_logl);
    }

    if (parameters_to_optimize & PLL_PARAMETER_SUBST_RATES)
    {
      params.which_parameters = PLL_PARAMETER_SUBST_RATES;
      cur_logl = pll_optimize_parameters_lbfgsb (&params);
      printf ("  [s_rates]: %ld s. : %f\n", time (NULL) - start_time, cur_logl);
      printf ("             %f %f %f %f %f %f\n", partition->subst_params[0][0],
                partition->subst_params[0][1], partition->subst_params[0][2],
                partition->subst_params[0][3], partition->subst_params[0][4],
                partition->subst_params[0][5]);
    }

    if (parameters_to_optimize & PLL_PARAMETER_ALPHA)
    {
      params.which_parameters = PLL_PARAMETER_ALPHA;
      cur_logl = pll_optimize_parameters_lbfgsb (&params);
      printf ("  [alpha]: %ld s. : %f\n", time (NULL) - start_time,
              cur_logl);
      printf ("             %f\n", params.lk_params.alpha_value);
    }

    if (parameters_to_optimize & PLL_PARAMETER_PINV)
    {
      params.which_parameters = PLL_PARAMETER_PINV;
      cur_logl = pll_optimize_parameters_lbfgsb (&params);
      printf ("  [p-inv]: %ld s. :%f\n", time (NULL) - start_time,
              cur_logl);
      printf ("             %f\n", partition->prop_invar[0]);
    }

    printf ("Iteration: %ld s. : %f\n", time (NULL) - start_time, cur_logl);
  }
  end_time = time (NULL);
  cur_logl *= -1;

  printf ("Final Log-L: %f\n", cur_logl);
  printf ("Time:  %ld s.\n", end_time - start_time);

  printf ("Alpha: %f\n", params.lk_params.alpha_value);
  printf ("P-inv: %f\n", partition->prop_invar[0]);
  printf ("Rates: %f %f %f %f %f %f\n", partition->subst_params[0][0],
          partition->subst_params[0][1], partition->subst_params[0][2],
          partition->subst_params[0][3], partition->subst_params[0][4],
          partition->subst_params[0][5]);

  printf ("Final Log-L: %f\n", logl);

  /* CLEAN */
  /* destroy all structures allocated for the concrete PLL partition instance */
  pll_destroy_partition (partition);

  free(subst_params_symmetries);
  free (branch_lengths);
  free (matrix_indices);
  free (operations);

  return (0);
}
