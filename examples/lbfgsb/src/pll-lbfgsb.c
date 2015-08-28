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

#define OPTIMIZE_BRANCHES     1
#define OPTIMIZE_SUBST_PARAMS 1
#define OPTIMIZE_ALPHA        1
#define OPTIMIZE_FREQS        1
#define OPTIMIZE_PINV         0

/* tolerances */
#define OPT_EPSILON       1e-1
#define OPT_PARAM_EPSILON 1e-1

static void fatal (const char * format, ...) __attribute__ ((noreturn));

/* a callback function for performing a full traversal */
static int cb_full_traversal(pll_utree_t * node)
{
  return 1;
}

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
  set_missing_branch_length_recursive (tree, length);
  set_missing_branch_length_recursive (tree->back, length);
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

static int * build_model_symmetries (const char * modelmatrix)
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

/**
 * if tip_count == 0, this function reads the file twice. The first
 *           time for getting the number of tips.
 * taxa_list allows mapping the sequences to their respective tips.
 *           if taxa_list == NULL, we assume that a hashing table
 *           has been already created.
 */
static pll_partition_t * partition_fasta_create (const char *file,
                                                         unsigned int states,
                                                         unsigned int n_rate_matrices,
                                                         unsigned int n_rate_cats,
                                                         int attributes,
                                                         int rooted,
                                                         unsigned int tip_count,
                                                         const char **tipnames)
{

  int unsigned i, j;
  pll_partition_t * partition;
  pll_fasta_t * fp;

  pll_errno = 0;

  /* open FASTA file */
  fp = pll_fasta_open (file, pll_map_fasta);
  if (!fp)
    return PLL_FAILURE;
  {
    char * seq = NULL;
    char * hdr = NULL;
    long seqlen;
    long hdrlen;
    long seqno;

    if (!tip_count)
    {
      /* get the number of tips */
      while (pll_fasta_getnext (fp, &hdr, &hdrlen, &seq, &seqlen, &seqno))
      {
        free (seq);
        free (hdr);
        ++tip_count;
      }

      /* reset the fasta pointer */

      if (!pll_fasta_rewind (fp))
      {
        pll_fasta_close (fp);
        return PLL_FAILURE;
      }
    }

    /* allocate arrays to store FASTA headers and sequences */
    char ** headers = (char **) calloc ((size_t)tip_count, sizeof(char *));
    char ** seqdata = (char **) calloc ((size_t)tip_count, sizeof(char *));

    /* read FASTA sequences and make sure they are all of the same length */
    long sites = -1;
    for (i = 0; pll_fasta_getnext (fp, &hdr, &hdrlen, &seq, &seqlen, &seqno);
        ++i)
    {
      if (i >= tip_count)
      {
        snprintf (pll_errmsg, 200,
                  "FASTA file contains more sequences than expected");
        pll_errno = PLL_ERROR_TAXA_MISMATCH;
        return PLL_FAILURE;
      }

      if (sites != -1 && sites != seqlen)
      {
        snprintf (pll_errmsg, 200,
                  "FASTA file does not contain equal size sequences");
        pll_errno = PLL_ERROR_SEQLEN_MISMATCH;
        return PLL_FAILURE;
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
    pll_errno = 0;

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

    partition = pll_partition_create (
        tip_count, rooted ? (tip_count - 1) : (tip_count - 2),
        states,
        (unsigned int) sites,
        n_rate_matrices,
        rooted ? (2 * tip_count - 2) : (2 * tip_count - 3),
        n_rate_cats,
        rooted ? (tip_count - 1) : (tip_count - 2),
        attributes);
    if (!partition)
    {
      return PLL_FAILURE;
    }

    /* find sequences and link them with the corresponding taxa */
    for (i = 0; i < tip_count; ++i)
    {
      int tip_clv_index = -1;

      for (j = 0; j < tip_count; ++j)
      {
        if (!strcmp (tipnames[j], headers[i]))
        {
          tip_clv_index = (int) j;
          break;
        }
      }

      if (tip_clv_index == -1)
      {
        snprintf (pll_errmsg, 200,
                  "Sequence with header %s does not appear in the tree",
                  headers[i]);
        pll_errno = PLL_ERROR_TAXA_MISMATCH;
        return PLL_FAILURE;
      }

      const unsigned int * pll_map = 0;
      switch (states)
      {
        case 4:
          pll_map = pll_map_nt;
          break;
        case 20:
          pll_map = pll_map_aa;
          break;
        default:
        {
            /* clean and return */
            for (j = 0; j < tip_count; ++j)
            {
              free (seqdata[j]);
              free (headers[j]);
            }
            free (seqdata);
            free (headers);
          return PLL_FAILURE;
        }
      }
      int set_states = pll_set_tip_states (partition,
                                           (unsigned int) tip_clv_index,
                                           pll_map,
                                           seqdata[i]);

      if (set_states == PLL_FAILURE)
      {
        pll_partition_destroy(partition);
        partition = 0;
      }
    }

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

  return partition;
} /* create_partition_fasta */

int main (int argc, char * argv[])
{

  unsigned int i, tip_count, nodes_count, branch_count, inner_nodes_count;
  pll_partition_t * partition;
  pll_operation_t * operations = NULL;
  double * branch_lengths = NULL;
  unsigned int * matrix_indices = NULL;
  pll_optimize_options_t params;
  time_t start_time, end_time;
  int parameters_to_optimize;

  pll_utree_t * tree;
  pll_utree_t ** travbuffer;
  unsigned int * data;
  int * subst_params_symmetries;

  unsigned int matrix_count, ops_count;

  if (argc != 4)
    fatal (" syntax: %s [newick] [fasta] [model]", argv[0]);

  tree = pll_utree_parse_newick(argv[1], &tip_count);
  nodes_count = 2*tip_count - 2;
  branch_count = 2*tip_count - 3;
  inner_nodes_count = tip_count - 2;

  /* fix all missing branch lengths to 0.00001 */
  set_missing_branch_length (tree, 0.00001);

  /*  obtain an array of pointers to tip nodes */
    pll_utree_t ** tipnodes = (pll_utree_t  **)calloc((size_t)tip_count,
                                                      sizeof(pll_utree_t *));
    pll_utree_query_tipnodes(tree, tipnodes);

  /* create and populate the list of tipnames */
  char **tipnames = (char **) malloc (tip_count * sizeof(char *));
  data = (unsigned int *) malloc ((size_t)tip_count * sizeof(unsigned int));
  for (i = 0; i < tip_count; ++i)
    tipnames[i] = tipnodes[i]->label;

  partition = partition_fasta_create (
      argv[2],
      STATES,
      1,
      RATE_CATS,
      PLL_ATTRIB_ARCH_SSE,
      0,
      tip_count,
      (const char **)tipnames);

  if (!partition)
    fatal ("Error %d: %s\n", pll_errno, pll_errmsg);

  /* we no longer need these arrays */
  free (data);
  free (tipnodes);
  free (tipnames);

  subst_params_symmetries = build_model_symmetries (argv[3]);
  printf ("Model: ");
  for (i = 0; i < SUBST_PARAMS; i++)
    printf ("%d", subst_params_symmetries[i]);
  printf ("\n");

  travbuffer = (pll_utree_t **)malloc((size_t)nodes_count * sizeof(pll_utree_t *));
  branch_lengths = (double *)malloc((size_t)branch_count * sizeof(double));
  matrix_indices = (unsigned int *)malloc((size_t)branch_count * sizeof(unsigned int));
  operations = (pll_operation_t *)malloc((size_t)inner_nodes_count *
                                                sizeof(pll_operation_t));
  /* perform a postorder traversal of the unrooted tree */
  unsigned int traversal_size;
  if (!pll_utree_traverse (tree,
                           cb_full_traversal,
                           travbuffer, &
                           traversal_size))
    fatal("Function pll_utree_traverse() requires inner nodes as parameters");

  /* given the computed traversal descriptor, generate the operations
     structure, and the corresponding probability matrix indices that
     may need recomputing */
  pll_utree_create_operations(travbuffer,
                              traversal_size,
                              branch_lengths,
                              matrix_indices,
                              operations,
                              &matrix_count,
                              &ops_count);

  /* initialize the array of base frequencies */
  double frequencies[4] =
    { 0.25, 0.25, 0.25, 0.25 };

  /* substitution rates for the 4x4 GTR model. This means we need exactly
   (4*4-4)/2 = 6 values, i.e. the number of elements above the diagonal */
  double subst_params[6] =
    { 1, 1, 1, 1, 1, 1 };

  /* we'll use 4 rate categories, and currently initialize them to 0 */
  double rate_cats[RATE_CATS] =
    { 0 };

  /* compute the discretized category rates from a gamma distribution
   with alpha shape 1 and store them in rate_cats  */
  pll_compute_gamma_cats (0.1, RATE_CATS, rate_cats);

  /* set frequencies at model with index 0 (we currently have only one model) */
  pll_set_frequencies (partition, 0, frequencies);

  /* set 6 substitution parameters at model with index 0 */
  pll_set_subst_params (partition, 0, subst_params);

  /* set rate categories */
  pll_set_category_rates (partition, rate_cats);

  pll_update_prob_matrices (partition, 0, matrix_indices, branch_lengths,
                            2 * tip_count - 3);

  pll_update_partials (partition, operations, tip_count - 2);

  double logl = pll_compute_edge_loglikelihood (partition,
                                                tree->clv_index,
                                                tree->scaler_index,
                                                tree->back->clv_index,
                                                tree->back->scaler_index,
                                                tree->pmatrix_index,
                                                0);

  char * newick = pll_utree_export_newick(tree);
  printf("Starting tree: %s\n", newick);
  free(newick);
  printf ("Log-L: %f\n", logl);

  /* pll stuff */
  params.lk_params.partition = partition;
  params.lk_params.operations = operations;
  params.lk_params.branch_lengths = branch_lengths;
  params.lk_params.matrix_indices = matrix_indices;
  params.lk_params.alpha_value = 0.1;
  params.lk_params.freqs_index = 0;
  params.lk_params.rooted = 0;
  params.lk_params.where.unrooted_t.parent_clv_index = tree->clv_index;
  params.lk_params.where.unrooted_t.parent_scaler_index = tree->scaler_index;
  params.lk_params.where.unrooted_t.child_clv_index = tree->back->clv_index;
  params.lk_params.where.unrooted_t.child_scaler_index = tree->back->scaler_index;
  params.lk_params.where.unrooted_t.edge_pmatrix_index = tree->pmatrix_index;

  /* optimization parameters */
  params.params_index = 0;
  params.subst_params_symmetries = subst_params_symmetries;
  params.factr = 1e8;
  params.pgtol = OPT_PARAM_EPSILON;

  parameters_to_optimize =
      (OPTIMIZE_SUBST_PARAMS * PLL_PARAMETER_SUBST_RATES) |
      (OPTIMIZE_ALPHA * PLL_PARAMETER_ALPHA) |
      (OPTIMIZE_BRANCHES * PLL_PARAMETER_BRANCHES_ALL) |
      (OPTIMIZE_PINV * PLL_PARAMETER_PINV) |
      (OPTIMIZE_FREQS * PLL_PARAMETER_FREQUENCIES);

  start_time = time (NULL);
  logl *= -1;
  double cur_logl = logl + 10;
  while (fabs (cur_logl - logl) > OPT_EPSILON)
  {
    logl = cur_logl;

    if (parameters_to_optimize & PLL_PARAMETER_FREQUENCIES)
    {
      params.freq_ratios = calloc ((size_t) partition->states - 1,
                                   sizeof(double));
      for (i = 0; i < (partition->states - 1); i++)
      {
        params.freq_ratios[i] =
          partition->frequencies[params.lk_params.freqs_index][i] /
          partition->frequencies[params.lk_params.freqs_index]
                                 [partition->states - 1];
      }
      params.which_parameters = PLL_PARAMETER_FREQUENCIES;
      cur_logl = pll_optimize_parameters_lbfgsb (&params);
      printf ("  %5ld s [freqs]: %f\n", time (NULL) - start_time, cur_logl);
      printf ("             ");
      for (i = 0; i < partition->states; i++)
        printf ("%f ", partition->frequencies[0][i]);
      printf ("\n");
      free (params.freq_ratios);
    }

    if (parameters_to_optimize & PLL_PARAMETER_BRANCHES_ALL)
    {
      params.which_parameters = PLL_PARAMETER_BRANCHES_SINGLE;
      cur_logl = pll_optimize_branch_lengths_iterative (&params, tree, 1);

      /* reset variables */
      //params.lk_params.branch_lengths[0] = branch_lengths[0];
      pll_utree_create_operations(travbuffer,
                                  traversal_size,
                                  branch_lengths,
                                  matrix_indices,
                                  operations,
                                  &matrix_count,
                                  &ops_count);
      params.lk_params.where.unrooted_t.parent_clv_index = tree->clv_index;
      params.lk_params.where.unrooted_t.parent_scaler_index = tree->scaler_index;
      params.lk_params.where.unrooted_t.child_clv_index = tree->back->clv_index;
      params.lk_params.where.unrooted_t.child_scaler_index = tree->back->scaler_index;
      params.lk_params.where.unrooted_t.edge_pmatrix_index = tree->pmatrix_index;

      printf ("  %5ld s [branches]: %f\n", time (NULL) - start_time,
              cur_logl);
    }

    if (parameters_to_optimize & PLL_PARAMETER_SUBST_RATES)
    {
      params.which_parameters = PLL_PARAMETER_SUBST_RATES;
      cur_logl = pll_optimize_parameters_lbfgsb (&params);
      printf ("  %5ld s [s_rates]: %f\n", time (NULL) - start_time, cur_logl);
      printf ("             %f %f %f %f %f %f\n", partition->subst_params[0][0],
              partition->subst_params[0][1], partition->subst_params[0][2],
              partition->subst_params[0][3], partition->subst_params[0][4],
              partition->subst_params[0][5]);
    }

    if (parameters_to_optimize & PLL_PARAMETER_ALPHA)
    {
      params.which_parameters = PLL_PARAMETER_ALPHA;
      cur_logl = pll_optimize_parameters_lbfgsb (&params);
      printf ("  %5ld s [alpha]: %f\n", time (NULL) - start_time, cur_logl);
      printf ("             %f\n", params.lk_params.alpha_value);
    }

    if (parameters_to_optimize & PLL_PARAMETER_PINV)
    {
      params.which_parameters = PLL_PARAMETER_PINV;
      cur_logl = pll_optimize_parameters_lbfgsb (&params);
      printf ("  %5ld s [p-inv]: %f\n", time (NULL) - start_time, cur_logl);
      printf ("             %f\n", partition->prop_invar[0]);
    }

    printf ("Iteration: %5ld s. : %f\n", time (NULL) - start_time, cur_logl);
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

  newick = pll_utree_export_newick(tree);
  printf("Final tree: %s\n", newick);
  free(newick);
  printf ("Final Log-L: %f\n", logl);

  /* CLEAN */

  pll_utree_destroy (tree);
  pll_partition_destroy (partition);

  free (travbuffer);
  free (subst_params_symmetries);
  free (branch_lengths);
  free (matrix_indices);
  free (operations);

  return (0);
}
