/*
    Copyright (C) 2015 Tomas Flouri, Diego Darriba

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
#include <stdarg.h>
#include <search.h>


#define STATES    4
#define RATE_CATS 4

static void fatal(const char * format, ...) __attribute__ ((noreturn));

static void * xmalloc(size_t size)
{ 
  void * t;
  t = malloc(size);
  if (!t)
    fatal("Unable to allocate enough memory.");
  
  return t;
} 
  
static char * xstrdup(const char * s)
{ 
  size_t len = strlen(s);
  char * p = (char *)xmalloc(len+1);
  return strcpy(p,s);
} 

/* a callback function for performing a full traversal */
static int cb_full_traversal(pll_rnode_t * node)
{
  return 1;
}

/* branch lengths not present in the newick file get a value of 0.000001 */
static void set_missing_branch_length(pll_rtree_t * tree, double length)
{
  unsigned int i;

  for (i = 0; i < tree->tip_count + tree->inner_count; ++i)
    if (!tree->nodes[i]->length)
      tree->nodes[i]->length = length;
}

static void fatal(const char * format, ...)
{
  va_list argptr;
  va_start(argptr, format);
  vfprintf(stderr, format, argptr);
  va_end(argptr);
  fprintf(stderr, "\n");
  exit(EXIT_FAILURE);
}

int main(int argc, char * argv[])
{
  unsigned int i;
  unsigned int tip_nodes_count, inner_nodes_count, nodes_count, branch_count;
  unsigned int matrix_count, ops_count;
  unsigned int * matrix_indices;
  double * branch_lengths;
  pll_partition_t * partition;
  pll_operation_t * operations;
  pll_rnode_t ** travbuffer;

  if (argc != 3)
    fatal(" syntax: %s [newick] [fasta]", argv[0]);

  pll_rtree_t * tree = pll_rtree_parse_newick(argv[1]);
  if (!tree)
    fatal("Tree must be an rooted binary tree");
  
  /* fix all missing branch lengths (i.e. those that did not appear in the
     newick) to 0.000001 */
  set_missing_branch_length(tree, 0.000001);

  /* compute and show node count information */
  tip_nodes_count = tree->tip_count;
  inner_nodes_count = tip_nodes_count - 1;
  nodes_count = inner_nodes_count + tip_nodes_count;
  branch_count = nodes_count - 1;

  printf("Number of tip/leaf nodes in tree: %d\n", tip_nodes_count);
  printf("Number of inner nodes in tree: %d\n", inner_nodes_count);
  printf("Total number of nodes in tree: %d\n", nodes_count);
  printf("Number of branches in tree: %d\n", branch_count);

  /* Uncomment to display ASCII tree and newick format

  printf("Number of tips in tree: %d\n", tip_nodes_count);
  pll_rtree_show_ascii(tree->nodes[nodes_count-1],
                       PLL_UTREE_SHOW_LABEL |
                       PLL_UTREE_SHOW_BRANCH_LENGTH |
                       PLL_UTREE_SHOW_CLV_INDEX);
  char * newick = pll_rtree_export_newick(tree->nodes[nodes_count-1], NULL);
  printf("%s\n", newick);
  free(newick);

  */

  /* create a libc hash table of size tip_count */
  hcreate(tip_nodes_count);

  /* populate a libc hash table with tree tip labels */
  unsigned int * data = (unsigned int *)xmalloc(tip_nodes_count *
                                                sizeof(unsigned int));
  for (i = 0; i < tip_nodes_count; ++i)
  {
    data[i] = tree->nodes[i]->clv_index;
    ENTRY entry;
#ifdef __APPLE__
    entry.key = xstrdup(tree->nodes[i]->label);
#else
    entry.key = tree->nodes[i]->label;
#endif
    entry.data = (void *)(data+i);
    hsearch(entry, ENTER);
  }

  /* open FASTA file */
  pll_fasta_t * fp = pll_fasta_open(argv[2], pll_map_fasta);
  if (!fp)
    fatal("Error opening file %s", argv[2]);

  char * seq = NULL;
  char * hdr = NULL;
  long seqlen;
  long hdrlen;
  long seqno;

  /* allocate arrays to store FASTA headers and sequences */
  char ** headers = (char **)calloc(tip_nodes_count, sizeof(char *));
  char ** seqdata = (char **)calloc(tip_nodes_count, sizeof(char *));

  /* read FASTA sequences and make sure they are all of the same length */
  int sites = -1;
  for (i = 0; pll_fasta_getnext(fp,&hdr,&hdrlen,&seq,&seqlen,&seqno); ++i)
  {
    if (i >= tip_nodes_count)
      fatal("FASTA file contains more sequences than expected");

    /* assert that sequence length is not longer than the number unsigned int
       can store */
    assert(!(seqlen & ((unsigned long)(~0) ^ (unsigned int)(~0))));

    if (sites != -1 && sites != seqlen)
      fatal("FASTA file does not contain equal size sequences\n");

    if (sites == -1) sites = seqlen;

    headers[i] = hdr;
    seqdata[i] = seq;
  }

  /* did we stop reading the file because we reached EOF? */
  if (pll_errno != PLL_ERROR_FILE_EOF)
    fatal("Error while reading file %s", argv[2]);

  /* close FASTA file */
  pll_fasta_close(fp);

  if (sites == -1)
    fatal("Unable to read alignment");

  if (i != tip_nodes_count)
    fatal("Some taxa are missing from FASTA file");

  /* create the PLL partition instance

  tip_nodes_count : the number of tip sequences we want to have
  inner_nodes_count : the number of CLV buffers to be allocated for inner nodes
  STATES : the number of states that our data have
  1 : number of different substitution models (or eigen decomposition)
      to use concurrently (i.e. 4 for LG4)
  branch_count: number of probability matrices to be allocated
  RATE_CATS : number of rate categories we will use
  inner_nodes_count : how many scale buffers to use
  PLL_ATTRIB_ARCH_SSE : list of flags for hardware acceleration

  */

  partition = pll_partition_create(tip_nodes_count,
                                   inner_nodes_count,
                                   STATES,
                                   (unsigned int)sites,
                                   1,
                                   branch_count,
                                   RATE_CATS,
                                   inner_nodes_count,
                                   PLL_ATTRIB_ARCH_CPU);

  /* find sequences in hash table and link them with the corresponding taxa */
  for (i = 0; i < tip_nodes_count; ++i)
  {
    ENTRY query;
    query.key = headers[i];
    ENTRY * found = NULL;

    found = hsearch(query,FIND);

    if (!found)
      fatal("Sequence with header %s does not appear in the tree", headers[i]);

    unsigned int tip_clv_index = *((unsigned int *)(found->data));

    pll_set_tip_states(partition, tip_clv_index, pll_map_nt, seqdata[i]);
  }

  /* destroy hash table */
  hdestroy();

  /* we no longer need these two arrays (keys and values of hash table... */
  free(data);

  /* ...neither the sequences and the headers as they are already
     present in the form of probabilities in the tip CLVs */
  for(i = 0; i < tip_nodes_count; ++i)
  {
    free(seqdata[i]);
    free(headers[i]);
  }
  free(seqdata);
  free(headers);


  /* allocate a buffer for storing pointers to nodes of the tree in postorder
     traversal */
  travbuffer = (pll_rnode_t **)xmalloc(nodes_count * sizeof(pll_rnode_t *));

  branch_lengths = (double *)xmalloc(branch_count * sizeof(double));
  matrix_indices = (unsigned int *)xmalloc(branch_count * sizeof(unsigned int));
  operations = (pll_operation_t *)xmalloc(inner_nodes_count *
                                          sizeof(pll_operation_t));

  /* perform a postorder traversal of the rooted tree */
  unsigned int traversal_size;
  if (!pll_rtree_traverse(tree->root,
                          PLL_TREE_TRAVERSE_POSTORDER,
                          cb_full_traversal,
                          travbuffer,
                          &traversal_size))
    fatal("Function pll_rtree_traverse() requires inner nodes as parameters");

  /* given the computed traversal descriptor, generate the operations
     structure, and the corresponding probability matrix indices that
     may need recomputing */
  pll_rtree_create_operations(travbuffer,
                              traversal_size,
                              branch_lengths,
                              matrix_indices,
                              operations,
                              &matrix_count,
                              &ops_count);


  printf ("Traversal size: %d\n", traversal_size);
  printf ("Operations: %d\n", ops_count);
  printf ("Probability Matrices: %d\n", matrix_count);

  /* initialize the array of base frequencies */
  double frequencies[4] = { 0.17, 0.19, 0.25, 0.39 };

  /* substitution rates for the 4x4 GTR model. This means we need exactly
     (4*4-4)/2 = 6 values, i.e. the number of elements above the diagonal */
  double subst_params[6] = {1,1,1,1,1,1};

  /* we'll use 4 rate categories, and currently initialize them to 0 */
  double rate_cats[4] = {0};

  /* compute the discretized category rates from a gamma distribution
     with alpha shape 1 and store them in rate_cats  */
  pll_compute_gamma_cats(1, 4, rate_cats, PLL_GAMMA_RATES_MEAN);

  /* set frequencies at model with index 0 (we currently have only one model) */
  pll_set_frequencies(partition, 0, frequencies);

  /* set 6 substitution parameters at model with index 0 */
  pll_set_subst_params(partition, 0, subst_params);

  /* set rate categories */
  pll_set_category_rates(partition, rate_cats);

  /* update matrix_count probability matrices for model with index 0. The i-th
     matrix (i ranges from 0 to matrix_count-1) is generated using branch length
     branch_lengths[i] and can be refered to with index matrix_indices[i]. Note
     that for our case, branch_count = matrix_count since we have one matrix
     per branch */
  unsigned int params_indices[4] = {0,0,0,0};

  pll_update_prob_matrices(partition,
                           params_indices,
                           matrix_indices,
                           branch_lengths,
                           matrix_count);

  /* Uncomment to output the probability matrices (for each branch and each rate
     category) on screen
  for (i = 0; i < branch_count; ++i)
  {
    printf ("P-matrix (%d) for branch length %f\n", i, branch_lengths[i]);
    pll_show_pmatrix(partition, i,17);
    printf ("\n");
  }

  */

  /* use the operations array to compute all inner_nodes_count inner CLVs.
     Operations will be carried out sequentially starting from operation 0 and
     upwards */
  pll_update_partials(partition, operations, inner_nodes_count);

  /* Uncomment to print on screen the CLVs at tip and inner nodes. From 0 to
     0 to tip_nodes_count-1 are tip CLVs, the rest are inner node CLVs.

  for (i = 0; i < nodes_count; ++i)
  {
    printf ("CLV %d: ", i);
    pll_show_clv(partition,
                 tree->nodes[i]->clv_index,
                 tree->nodes[i]->scaler_index,
                 17);
  }

  */

  /* compute the likelihood on an edge of the unrooted tree by specifying
     the CLV indices at the two end-point of the branch, the probability matrix
     index for the concrete branch length, and the index of the model of whose
     frequency vector is to be used */

  double logl = pll_compute_root_loglikelihood(partition,
                                               tree->root->clv_index,
                                               tree->root->scaler_index,
                                               params_indices,
                                               NULL);

  printf("Log-L: %f\n", logl);

  /* destroy all structures allocated for the concrete PLL partition instance */
  pll_partition_destroy(partition);

  /* deallocate traversal buffer, branch lengths array, matrix indices
     array and operations */
  free(travbuffer);
  free(branch_lengths);
  free(matrix_indices);
  free(operations);

  /* we will no longer need the tree structure */
  pll_rtree_destroy(tree,NULL);

  return (EXIT_SUCCESS);
}
