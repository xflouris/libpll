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
#include "lbfgsb/lbfgsb.h"

typedef struct
{
  int initialized;
  int clv_index;
  int matrix_index;
  int scaler_index;
  int branch_len_index;
  double * branch_length;
} node_info_t;

static int v_int_max (int * v, int n)
{
  int i, max = v[0];
  for (i = 1; i < n; i++)
    if (v[i] > max)
      max = v[i];
  return max;
}

static void inittravtree (pll_utree_t *p, pll_operation_t *ops,
                          int n_operations, int depth, int clv_index,
                          int matrix_index, int scaler_index,
                          node_info_t *node_infos, int *next_node_info)
{
  int i;

  /* traverse tree to set initialized and v to initial values */
  pll_operation_t *op = ops;
  int children_clv_index[2] =
    { -1, -1 };
  int children_matrix_index[2] =
    { -1, -1 };
  int children_scaler_index[2] =
    {
    PLL_SCALE_BUFFER_NONE,
    PLL_SCALE_BUFFER_NONE };
  node_info_t *node_info = &(node_infos[(*next_node_info)++]);
  node_info->initialized = 0;
  node_info->clv_index = clv_index;
  node_info->matrix_index = matrix_index;
  node_info->scaler_index = scaler_index;
  node_info->branch_length = &p->length;

  p->data = (void *) node_info;

  if (p->next != NULL)
  {
    for (i = 0; i < n_operations; i++, op++)
    {
      if (op->parent_clv_index == clv_index)
      {
        children_clv_index[0] = op->child1_clv_index;
        children_clv_index[1] = op->child2_clv_index;
        children_matrix_index[0] = op->child1_matrix_index;
        children_matrix_index[1] = op->child2_matrix_index;
        children_scaler_index[0] = op->child1_scaler_index;
        children_scaler_index[1] = op->child2_scaler_index;
      }
    }
  }

  if (p->length < PLL_OPT_MIN_BRANCH_LEN)
  {
    p->length = PLL_OPT_DEFAULT_BRANCH_LEN;
    p->back->length = PLL_OPT_DEFAULT_BRANCH_LEN;
  }

  if (p->next != NULL)
  {
    /* is inner */
    inittravtree (p->next->back, ops, n_operations, depth + 1,
                  children_clv_index[0], children_matrix_index[0],
                  children_scaler_index[0], node_infos, next_node_info);
    inittravtree (p->next->next->back, ops, n_operations, depth + 1,
                  children_clv_index[1], children_matrix_index[1],
                  children_scaler_index[1], node_infos, next_node_info);

    node_info = &(node_infos[(*next_node_info)++]);
    node_info->initialized = 0;
    node_info->clv_index = clv_index;
    node_info->matrix_index =
        ((node_info_t *) p->next->back->data)->matrix_index;
    node_info->scaler_index = scaler_index;
    node_info->branch_length = &(p->next->length);
    p->next->data = (void *) node_info;

    node_info = &(node_infos[(*next_node_info)++]);
    node_info->initialized = 0;
    node_info->clv_index = clv_index;
    node_info->matrix_index =
        ((node_info_t *) p->next->next->back->data)->matrix_index;
    node_info->scaler_index = scaler_index;
    node_info->branch_length = &(p->next->next->length);
    p->next->next->data = (void *) node_info;
  }
} /* inittravtree */

static double recomp_iterative (pll_optimize_options_t * params,
                              pll_utree_t * tree, double prev_lnl)
{

  /* evaluate at edge */
  node_info_t *info1 = (node_info_t *) tree->data;
  node_info_t *info2 = (node_info_t *) tree->back->data;
  double lnl, new_lnl;

  lnl = prev_lnl;

  /* set Branch Length */
  params->lk_params.branch_lengths[0] = *(info1->branch_length);
  params->lk_params.where.unrooted_t.child_clv_index = info2->clv_index;
  params->lk_params.where.unrooted_t.child_scaler_index = info2->scaler_index;
  params->lk_params.where.unrooted_t.parent_clv_index = info1->clv_index;
  params->lk_params.where.unrooted_t.parent_scaler_index = info1->scaler_index;
  params->lk_params.where.unrooted_t.edge_pmatrix_index = info1->matrix_index;

  new_lnl = -1 * pll_optimize_parameters_lbfgsb(params);
  if (new_lnl < lnl)
  {
    pll_update_prob_matrices(params->lk_params.partition,
                             params->params_index,
                             &(info1->matrix_index),
                             info1->branch_length,
                             1);
    lnl = pll_compute_edge_loglikelihood (params->lk_params.partition,
                                            info2->clv_index, info2->scaler_index,
                                            info1->clv_index, info1->scaler_index,
                                            info1->matrix_index,
                                            params->lk_params.freqs_index);
  }
  else
  {
    lnl = new_lnl;
    *(info1->branch_length) = params->lk_params.branch_lengths[0];
    *(info2->branch_length) = params->lk_params.branch_lengths[0];
  }

  DBG("forward lnL: %f (%f)\n", lnl, info1->branch_length);

  if (tree->next)
  {
    node_info_t *info_child1 = (node_info_t *) tree->next->back->data;
    node_info_t *info_child2 = (node_info_t *) tree->next->next->back->data;
    pll_operation_t new_op;

    /* set CLV */
    new_op.parent_clv_index = info1->clv_index;
    new_op.parent_scaler_index = info1->scaler_index;
    new_op.child1_clv_index = info2->clv_index;
    new_op.child1_matrix_index = info2->matrix_index;
    new_op.child1_scaler_index = info2->scaler_index;
    new_op.child2_clv_index = info_child2->clv_index;
    new_op.child2_matrix_index = info_child2->matrix_index;
    new_op.child2_scaler_index = info_child2->scaler_index;
    pll_update_partials (params->lk_params.partition, &new_op, 1);
    /* eval */
    recomp_iterative (params, tree->next->back, lnl);

    /* set CLV */
    new_op.parent_clv_index = info1->clv_index;
    new_op.parent_scaler_index = info1->scaler_index;
    new_op.child1_clv_index = info2->clv_index;
    new_op.child1_matrix_index = info2->matrix_index;
    new_op.child1_scaler_index = info2->scaler_index;
    new_op.child2_clv_index = info_child1->clv_index;
    new_op.child2_matrix_index = info_child1->matrix_index;
    new_op.child2_scaler_index = info_child1->scaler_index;
    pll_update_partials (params->lk_params.partition, &new_op, 1);

    /* eval */
    recomp_iterative (params, tree->next->next->back, lnl);

    /* reset CLV */
    new_op.parent_clv_index = info1->clv_index;
    new_op.parent_scaler_index = info1->scaler_index;
    new_op.child1_clv_index = info_child1->clv_index;
    new_op.child1_matrix_index = info_child1->matrix_index;
    new_op.child1_scaler_index = info_child1->scaler_index;
    new_op.child2_clv_index = info_child2->clv_index;
    new_op.child2_matrix_index = info_child2->matrix_index;
    new_op.child2_scaler_index = info_child2->scaler_index;
    pll_update_partials (params->lk_params.partition, &new_op, 1);
  }

  return lnl;
} /* recomp_iterative */

static int set_x_to_parameters(pll_optimize_options_t * params,
                               double *x)
{
  pll_partition_t * partition = params->lk_params.partition;
  pll_operation_t * operations = params->lk_params.operations;
  double * branch_lengths = params->lk_params.branch_lengths;
  int * matrix_indices = params->lk_params.matrix_indices;
  int params_index = params->params_index;
  int n_branches, n_inner_nodes;
  double score;
  double * xptr = x;

  if (params->lk_params.rooted)
  {
    n_branches = 2 * partition->tips - 2;
    n_inner_nodes = partition->tips - 1;
  }
  else
  {
    n_branches = 2 * partition->tips - 3;
    n_inner_nodes = partition->tips - 2;
  }

  /* update substitution rate parameters */
  if (params->which_parameters & PLL_PARAMETER_SUBST_RATES)
  {
    int * symm;
    int n_subst_rates;
    double * subst_rates;

    symm = params->subst_params_symmetries;
    n_subst_rates = partition->states * (partition->states - 1) / 2;
    subst_rates = (double *) malloc ((size_t) n_subst_rates * sizeof(double));

    assert(subst_rates);

    if (symm)
    {
      int i, j, k;
      int n_subst_free_params = 0;

      /* compute the number of free parameters */
      n_subst_free_params = v_int_max (symm, n_subst_rates);

      /* assign values to the substitution rates */
      k = 0;
      for (i = 0; i <= n_subst_free_params; i++)
      {
        double next_value = (i == symm[n_subst_rates - 1]) ? 1.0 : xptr[k++];
        for (j = 0; j < n_subst_rates; j++)
          if (symm[j] == i)
          {
            subst_rates[j] = next_value;
          }
      }
      pll_set_subst_params (partition, 0, subst_rates);
      xptr += n_subst_free_params;
    }
    else
    {
      memcpy (subst_rates, xptr, (n_subst_rates - 1) * sizeof(double));
      subst_rates[n_subst_rates - 1] = 1.0;
    }
    free (subst_rates);
  }
  /* update stationary frequencies */
  if (params->which_parameters & PLL_PARAMETER_FREQUENCIES)
  {
    int i;
    double sum_ratios = 1.0;
    for (i = 0; i < (partition->states - 1); ++i)
      sum_ratios += xptr[i];
    for (i = 0; i < (partition->states - 1); ++i)
      partition->frequencies[params->lk_params.freqs_index][i] = xptr[i]
          / sum_ratios;
    partition->frequencies[params->lk_params.freqs_index][partition->states - 1] =
        1.0 / sum_ratios;
    xptr += (partition->states - 1);
  }
  /* update proportion of invariant parameters */
  if (params->which_parameters & PLL_PARAMETER_PINV)
  {
    if (!pll_update_invariant_sites_proportion (partition,
                                                params_index,
                                                xptr[0]))
    {
      pll_errno = PLL_ERROR_INVALID_PINV;
      snprintf (pll_errmsg, 200,
                        "Invalid proportion of invariant sites: %f", xptr[0]);
      return PLL_FAILURE;
    }
    xptr++;
  }
  if (params->which_parameters & PLL_PARAMETER_ALPHA)
  {
    /* assign discrete rates */
    double * rate_cats;
    rate_cats = malloc((size_t)partition->rate_cats * sizeof(double));
    params->lk_params.alpha_value = xptr[0];
    if (!pll_compute_gamma_cats (xptr[0], partition->rate_cats, rate_cats))
    {
      pll_errno = PLL_ERROR_INVALID_PINV;
      snprintf (pll_errmsg, 200, "Invalid alpha shape: %f", xptr[0]);
      return PLL_FAILURE;
    }
    pll_set_category_rates (partition, rate_cats);

    free(rate_cats);
    xptr++;
  }
  if (params->which_parameters & PLL_PARAMETER_BRANCH_LENGTHS)
  {
    /* assign branch lengths */
    memcpy (branch_lengths, xptr, n_branches * sizeof(double));
    xptr += n_branches;
  }
  if (params->which_parameters & PLL_PARAMETER_SINGLE_BRANCH)
   {
     /* assign branch lengths */
     *branch_lengths = *xptr;
     pll_update_prob_matrices (partition,
                         0,
                         &params->lk_params.where.unrooted_t.edge_pmatrix_index,
                         xptr,
                         1);
     xptr += 1;
   }
  else
   {
       pll_update_prob_matrices (partition, 0, matrix_indices, branch_lengths,
                                 n_branches);

       pll_update_partials (partition, operations, n_inner_nodes);
   }
  return PLL_SUCCESS;
}

static double compute_lnl_unrooted (pll_optimize_options_t * params, double *x)
{
  pll_partition_t * partition = params->lk_params.partition;
  pll_operation_t * operations = params->lk_params.operations;
  double * branch_lengths = params->lk_params.branch_lengths;
  int * matrix_indices = params->lk_params.matrix_indices;
  int params_index = params->params_index;
  int num_branch_lengths;
  int num_inner_nodes;
  double score;
  double * xptr = x;

  if (!set_x_to_parameters(params, x))
  {
    return -INFINITY;
  }

  if (params->lk_params.rooted)
  {
    score = -1
        * pll_compute_root_loglikelihood (
            partition, params->lk_params.where.rooted_t.root_clv_index,
            params->lk_params.where.rooted_t.scaler_index,
            params->lk_params.freqs_index);
  }
  else
  {
    score = -1
        * pll_compute_edge_loglikelihood (
            partition, params->lk_params.where.unrooted_t.parent_clv_index,
            params->lk_params.where.unrooted_t.parent_scaler_index,
            params->lk_params.where.unrooted_t.child_clv_index,
            params->lk_params.where.unrooted_t.child_scaler_index,
            params->lk_params.where.unrooted_t.edge_pmatrix_index,
            params->lk_params.freqs_index);
  }

  return score;
} /* compute_lnl_unrooted */

static int count_n_free_variables (pll_optimize_options_t * params)
{
  int num_variables = 0;
  pll_partition_t * partition = params->lk_params.partition;

  /* count number of variables for dynamic allocation */
  if (params->which_parameters & PLL_PARAMETER_SUBST_RATES)
  {
    int n_subst_rates = partition->states * (partition->states - 1) / 2;
    num_variables +=
        params->subst_params_symmetries ?
            v_int_max (params->subst_params_symmetries, n_subst_rates) :
            n_subst_rates - 1;
  }
  if (params->which_parameters & PLL_PARAMETER_FREQUENCIES)
    num_variables += partition->states - 1;
  num_variables += (params->which_parameters & PLL_PARAMETER_PINV) != 0;
  num_variables += (params->which_parameters & PLL_PARAMETER_ALPHA) != 0;
  num_variables += (params->which_parameters & PLL_PARAMETER_SINGLE_BRANCH)
      != 0;
  if (params->which_parameters & PLL_PARAMETER_BRANCH_LENGTHS)
  {
    int num_branch_lengths =
        params->lk_params.rooted ?
            (2 * partition->tips - 3) : (2 * partition->tips - 2);
    num_variables += num_branch_lengths;
  }
  return num_variables;
} /* count_n_free_variables */

/**
 * if tip_count == 0, this function reads the file twice. The first
 *           time for getting the number of tips.
 * taxa_list allows mapping the sequences to their respective tips.
 *           if taxa_list == NULL, we assume that a hashing table
 *           has been already created.
 */
PLL_EXPORT pll_partition_t * pll_create_partition_fasta (char *file, int states,
                                                         int n_rate_matrices,
                                                         int n_rate_cats,
                                                         int attributes,
                                                         int rooted,
                                                         int tip_count,
                                                         char **tipnames)
{

  int i, j;
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
    char ** headers = (char **) calloc (tip_count, sizeof(char *));
    char ** seqdata = (char **) calloc (tip_count, sizeof(char *));

    /* read FASTA sequences and make sure they are all of the same length */
    int sites = -1;
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

    partition = pll_create_partition (
        tip_count, rooted ? (tip_count - 1) : (tip_count - 2), states, sites,
        n_rate_matrices, rooted ? (2 * tip_count - 2) : (2 * tip_count - 3),
        n_rate_cats, rooted ? (tip_count - 1) : (tip_count - 2), attributes);
    if (!partition)
    {
      return PLL_FAILURE;
    }

    /* find sequences and link them with the corresponding taxa */
    for (i = 0; i < tip_count; ++i)
    {
      int tip_clv_index = -1;

      if (tipnames)
      {
        for (j = 0; j < tip_count; ++j)
        {
          if (!strcmp (tipnames[j], headers[i]))
          {
            tip_clv_index = j;
            break;
          }
        }
      }
      else
      {
        ENTRY query;
        query.key = headers[i];
        ENTRY * found = NULL;

        found = hsearch (query, FIND);

        if (found)
        {
          tip_clv_index = *((int *) (found->data));
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

      pll_set_tip_states (partition, tip_clv_index, pll_map_nt, seqdata[i]);
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
} /* pll_create_partition_fasta */

PLL_EXPORT double pll_optimize_parameters_lbfgsb (
    pll_optimize_options_t * params)
{
  int i;
  pll_partition_t * partition = params->lk_params.partition;

  /* L-BFGS-B */
  int max_corrections, num_variables;
  double score;
  double *x, *g, *lower_bounds, *upper_bounds, *wa;
  int *bound_type, *iwa;

  /*     static char task[60]; */
  int taskValue;
  int *task = &taskValue; /* must initialize !! */

  /*     static char csave[60]; */
  int csaveValue;
  int *csave = &csaveValue;
  double dsave[29];
  int isave[44];
  logical lsave[4];

  int iprint = -1;

  max_corrections = 5;
  num_variables = count_n_free_variables (params);

  x = (double *) calloc ((size_t) num_variables, sizeof(double));
  g = (double *) calloc ((size_t) num_variables, sizeof(double));
  lower_bounds = (double *) calloc ((size_t) num_variables, sizeof(double));
  upper_bounds = (double *) calloc ((size_t) num_variables, sizeof(double));
  bound_type = (int *) calloc ((size_t) num_variables, sizeof(int));

  {
    int * nbd_ptr = bound_type;
    double * l_ptr = lower_bounds, *u_ptr = upper_bounds;
    int check_n = 0;
    if (params->which_parameters & PLL_PARAMETER_SUBST_RATES)
    {
      int n_subst_rates;
      int n_subst_free_params;

      n_subst_rates = partition->states * (partition->states - 1) / 2;
      if (params->subst_params_symmetries)
        n_subst_free_params = v_int_max (params->subst_params_symmetries,
                                         n_subst_rates);
      else
        n_subst_free_params = n_subst_rates;

      for (i = 0; i < n_subst_free_params; i++)
      {
        nbd_ptr[i] = PLL_LBFGSB_BOUND_BOTH;
        l_ptr[i] = 0.001;
        u_ptr[i] = 1000.;
        x[i] = partition->subst_params[params->params_index][i];
      }
      nbd_ptr += n_subst_free_params;
      l_ptr += n_subst_free_params;
      u_ptr += n_subst_free_params;
      check_n += n_subst_free_params;
    }
    if (params->which_parameters & PLL_PARAMETER_FREQUENCIES)
    {
      int n_freqs_free_params;

      n_freqs_free_params = params->lk_params.partition->states - 1;
      for (i = 0; i < n_freqs_free_params; i++)
      {
        nbd_ptr[i] = PLL_LBFGSB_BOUND_BOTH;
        l_ptr[i] = 0.00001;
        u_ptr[i] = 1000;
        x[check_n + i] = params->freq_ratios[i];
      }
      check_n += n_freqs_free_params;
      nbd_ptr += n_freqs_free_params;
      l_ptr += n_freqs_free_params;
      u_ptr += n_freqs_free_params;
    }
    if (params->which_parameters & PLL_PARAMETER_PINV)
    {
      *nbd_ptr = PLL_LBFGSB_BOUND_BOTH;
      *l_ptr = 0.0;
      *u_ptr = 0.99;
      x[check_n] = partition->prop_invar[params->params_index];
      check_n += 1;
      nbd_ptr++;
      l_ptr++;
      u_ptr++;
    }
    if (params->which_parameters & PLL_PARAMETER_ALPHA)
    {
      *nbd_ptr = PLL_LBFGSB_BOUND_BOTH;
      *l_ptr = 0.020001;
      *u_ptr = 100.0;
      x[check_n] = params->lk_params.alpha_value;
      check_n += 1;
      nbd_ptr++;
      l_ptr++;
      u_ptr++;
    }
    if (params->which_parameters & PLL_PARAMETER_TOPOLOGY)
    {
      return PLL_FAILURE;
    }
    if (params->which_parameters & PLL_PARAMETER_SINGLE_BRANCH)
    {
        *nbd_ptr = PLL_LBFGSB_BOUND_LOWER;
        *l_ptr = PLL_OPT_MIN_BRANCH_LEN;
        x[check_n] = params->lk_params.branch_lengths[0];
        check_n += 1;
        nbd_ptr++;
        l_ptr++;
        u_ptr++;
    }
    if (params->which_parameters & PLL_PARAMETER_BRANCH_LENGTHS)
    {
      int num_branch_lengths =
          params->lk_params.rooted ?
              (2 * partition->tips - 3) : (2 * partition->tips - 2);
      for (i = 0; i < num_branch_lengths; i++)
      {
        nbd_ptr[i] = PLL_LBFGSB_BOUND_LOWER;
        l_ptr[i] = 0.00001;
        x[check_n + i] = params->lk_params.branch_lengths[i];
      }
      check_n += num_branch_lengths;
      nbd_ptr += num_branch_lengths;
      l_ptr += num_branch_lengths;
      u_ptr += num_branch_lengths;
    }
    assert(check_n == num_variables);
  }

  /*     We start the iteration by initializing task. */
  *task = (int) START;

  iwa = (int *) calloc ((size_t) 3 * num_variables, sizeof(int));
  wa = (double *) calloc (
      (size_t) (2 * max_corrections + 5) * num_variables
          + 12 * max_corrections * (max_corrections + 1),
      sizeof(double));

  int continue_opt = 1;
  while (continue_opt)
  {
    /*     This is the call to the L-BFGS-B code. */
    setulb (&num_variables, &max_corrections, x, lower_bounds, upper_bounds,
            bound_type, &score, g, &(params->factr), &(params->pgtol), wa, iwa,
            task, &iprint, csave, lsave, isave, dsave);
    if (IS_FG(*task))
    {
      /*
       * the minimization routine has returned to request the
       * function f and gradient g values at the current x.
       * Compute function value f for the sample problem.
       */

      score = compute_lnl_unrooted (params, x);

      double h, temp;
      for (i = 0; i < num_variables; i++)
      {
        double ERROR_X = 1.0e-4;
        temp = x[i];
        h = ERROR_X * fabs (temp);
        if (h < 1e-12)
          h = ERROR_X;

        x[i] = temp + h;
        h = x[i] - temp;
        double lnderiv = compute_lnl_unrooted (params, x);

        g[i] = (lnderiv - score) / h;

        /* reset variable */
        x[i] = temp;
      }
      if (!set_x_to_parameters (params, x))
      {
        return -INFINITY;
      }
    }
    else if (*task != NEW_X)
    {
      continue_opt = 0;
    }
  }

  free (iwa);
  free (wa);
  free (x);
  free (g);
  free (lower_bounds);
  free (upper_bounds);
  free (bound_type);

  return score;
} /* pll_optimize_parameters_lbfgsb */

PLL_EXPORT double pll_optimize_branch_lengths_iterative (
    pll_optimize_options_t * params, pll_utree_t * tree, int smoothings)
{
  int i, j;
  double lnl = 0.0;

  params->which_parameters = PLL_PARAMETER_SINGLE_BRANCH;

  int n_node_infos = params->lk_params.partition->tips * 4 - 6;
  node_info_t * node_infos = calloc(
      (size_t) n_node_infos,
          sizeof(node_info_t));
  int next_node_info = 0;

  inittravtree (tree, params->lk_params.operations,
                params->lk_params.partition->tips - 2, 0,
                params->lk_params.where.unrooted_t.child_clv_index,
                params->lk_params.where.unrooted_t.edge_pmatrix_index,
                params->lk_params.where.unrooted_t.child_scaler_index,
                node_infos, &next_node_info);
  inittravtree (tree->back, params->lk_params.operations,
                params->lk_params.partition->tips - 2, 0,
                params->lk_params.where.unrooted_t.parent_clv_index,
                params->lk_params.where.unrooted_t.edge_pmatrix_index,
                params->lk_params.where.unrooted_t.parent_scaler_index,
                node_infos, &next_node_info);
  assert(next_node_info == n_node_infos);

  /* set the branch length indices */
  for (i=0; i<n_node_infos; i++)
  {
    int m_index = node_infos[i].matrix_index;
    for (j=0; j<(2*params->lk_params.partition->tips - 2); j++)
    {
      if (m_index == params->lk_params.matrix_indices[j])
      {
        node_infos[i].branch_len_index = j;
        break;
      }
    }
  }

  for (i=0; i<smoothings; i++)
  {
    lnl = recomp_iterative (params, tree, PLL_OPT_LNL_UNLIKELY);
    lnl = recomp_iterative (params, tree->back, lnl);
  }

  for (i=0; i<n_node_infos; i++)
  {
    params->lk_params.branch_lengths[node_infos[i].branch_len_index]
      = *(node_infos[i].branch_length);
  }

  free(node_infos);

  return -1 * lnl;
} /* pll_optimize_branch_lengths_iterative */
