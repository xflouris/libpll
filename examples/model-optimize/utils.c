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
#include "utils.h"

int cb_full_traversal(pll_utree_t * node)
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
void set_missing_branch_length (pll_utree_t * tree, double length)
{
  set_missing_branch_length_recursive (tree, length);
  set_missing_branch_length_recursive (tree->back, length);
}

void fatal (const char * format, ...)
{
  va_list argptr;
  va_start(argptr, format);
  vfprintf (stderr, format, argptr);
  va_end(argptr);
  fprintf (stderr, "\n");
  exit (EXIT_FAILURE);
}

int * build_model_symmetries (const char * modelmatrix)
{
  int i, next_index;
  if(strlen(modelmatrix) != SUBST_PARAMS) return PLL_FAILURE;

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
pll_partition_t * partition_fasta_create (const char *file,
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

    partition = pll_partition_create (
        tip_count, rooted ? (tip_count - 1) : (tip_count - 2),
        states,
        (unsigned int) sites,
        n_rate_matrices,
        rooted ? (2 * tip_count - 2) : (2 * tip_count - 3),
        n_rate_cats,
        rooted ? (tip_count - 1) : (tip_count - 2),
        pll_map,
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

static void update_clvs(pll_partition_t * partition,
                        unsigned int params_index,
                        unsigned int * matrix_indices,
                        double * branch_lengths,
                        pll_operation_t * operations)
{
  unsigned int n_branches = 2*partition->tips - 3;
  unsigned int n_inner    = partition->tips - 2;

  if (matrix_indices)
    pll_update_prob_matrices (partition, 0,
                              matrix_indices,
                              branch_lengths,
                              n_branches);
  if (operations)
    pll_update_partials (partition,
                         operations,
                         n_inner);
}

double target_rates_opt (void * p, double * x)
{
  my_params_t * params = (my_params_t *) p;
  pll_partition_t * partition = params->partition;
  double score;

  /* set x to partition */
  int * symm;
  int n_subst_rates;
  double * subst_rates;

  symm = params->symmetries;
  n_subst_rates = partition->states * (partition->states - 1) / 2;
  subst_rates = (double *) malloc ((size_t) n_subst_rates * sizeof(double));

  assert(subst_rates);

  if (symm)
  {
    int i, j, k;
    int n_subst_free_params = params->n_subst_params;

    /* assign values to the substitution rates */
    k = 0;
    for (i = 0; i <= n_subst_free_params; i++)
    {
      double next_value = (i == symm[n_subst_rates - 1]) ? 1.0 : x[k++];
      for (j = 0; j < n_subst_rates; j++)
        if (symm[j] == i)
        {
          subst_rates[j] = next_value;
        }
    }
  }
  else
  {
    memcpy (subst_rates, x, ((size_t) n_subst_rates - 1) * sizeof(double));
    subst_rates[n_subst_rates - 1] = 1.0;
  }
  pll_set_subst_params (partition, params->params_index, subst_rates);
  free (subst_rates);

  update_clvs(partition,
              params->params_index,
              params->matrix_indices,
              params->branch_lengths,
              params->operations);

  score = -1
      * pll_compute_edge_loglikelihood (partition, params->parent_clv_index,
                                        params->parent_scaler_index,
                                        params->child_clv_index,
                                        params->child_scaler_index,
                                        params->edge_pmatrix_index,
                                        params->freqs_indices,
                                        NULL);

  return score;
}

double target_freqs_opt (void * p, double * x)
{
  my_params_t * params = (my_params_t *) p;
  pll_partition_t * partition = params->partition;
  double score;
  unsigned int i;

  unsigned int n_states = partition->states;
  unsigned int cur_index;
  double sum_ratios = 1.0;
  double *freqs = (double *) malloc ((size_t) n_states * sizeof(double));
  for (i = 0; i < (n_states - 1); ++i)
  {
    sum_ratios += x[i];
  }
  cur_index = 0;
  for (i = 0; i < (n_states); ++i)
  {
    if (i != params->highest_freq_state)
    {
      freqs[i] = x[cur_index] / sum_ratios;
      cur_index++;
    }
  }
  freqs[params->highest_freq_state] = 1.0 / sum_ratios;
  pll_set_frequencies (partition, params->params_index, freqs);
  free (freqs);

  update_clvs (partition, params->params_index, params->matrix_indices,
               params->branch_lengths, params->operations);

  score = -1
      * pll_compute_edge_loglikelihood (partition, params->parent_clv_index,
                                        params->parent_scaler_index,
                                        params->child_clv_index,
                                        params->child_scaler_index,
                                        params->edge_pmatrix_index,
                                        params->freqs_indices,
                                        NULL);

  return score;
}
