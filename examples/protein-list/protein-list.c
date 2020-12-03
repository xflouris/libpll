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


#define STATES    20
#define RATE_CATS 4
#define PROT_MODELS_COUNT 19

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

static const double * protein_models_rates_list[PROT_MODELS_COUNT] =
 {
   pll_aa_rates_dayhoff,
   pll_aa_rates_lg,
   pll_aa_rates_dcmut,
   pll_aa_rates_jtt,
   pll_aa_rates_mtrev,
   pll_aa_rates_wag,
   pll_aa_rates_rtrev,
   pll_aa_rates_cprev,
   pll_aa_rates_vt,
   pll_aa_rates_blosum62,
   pll_aa_rates_mtmam,
   pll_aa_rates_mtart,
   pll_aa_rates_mtzoa,
   pll_aa_rates_pmb,
   pll_aa_rates_hivb,
   pll_aa_rates_hivw,
   pll_aa_rates_jttdcmut,
   pll_aa_rates_flu,
   pll_aa_rates_stmtrev
 };

static const double * protein_models_freqs_list[PROT_MODELS_COUNT] =
 {
   pll_aa_freqs_dayhoff,
   pll_aa_freqs_lg,
   pll_aa_freqs_dcmut,
   pll_aa_freqs_jtt,
   pll_aa_freqs_mtrev,
   pll_aa_freqs_wag,
   pll_aa_freqs_rtrev,
   pll_aa_freqs_cprev,
   pll_aa_freqs_vt,
   pll_aa_freqs_blosum62,
   pll_aa_freqs_mtmam,
   pll_aa_freqs_mtart,
   pll_aa_freqs_mtzoa,
   pll_aa_freqs_pmb,
   pll_aa_freqs_hivb,
   pll_aa_freqs_hivw,
   pll_aa_freqs_jttdcmut,
   pll_aa_freqs_flu,
   pll_aa_freqs_stmtrev
 };

static const char * protein_models_names_list[PROT_MODELS_COUNT] =
{
   "DAYHOFF",
   "LG",
   "DCMUT",
   "JTT",
   "MTREV",
   "WAG",
   "RTREV",
   "CPREV",
   "VT",
   "BLOSUM62",
   "MTMAM",
   "MTART",
   "MTZOA",
   "PMB",
   "HIVB",
   "HIVW",
   "JTTDCMUT",
   "FLU",
   "STMTREV"
};

/* a callback function for performing a full traversal */
static int cb_full_traversal(pll_unode_t * node)
{
  return 1;
}

/* branch lengths not present in the newick file get a value of 0.000001 */
static void set_missing_branch_length(pll_utree_t * tree, double length)
{
  unsigned int i;

  for (i = 0; i < tree->tip_count; ++i)
    if (!tree->nodes[i]->length)
      tree->nodes[i]->length = length;

  for (i = tree->tip_count; i < tree->tip_count + tree->inner_count; ++i)
  {
    if (!tree->nodes[i]->length)
      tree->nodes[i]->length = length;
    if (!tree->nodes[i]->next->length)
      tree->nodes[i]->next->length = length;
    if (!tree->nodes[i]->next->next->length)
      tree->nodes[i]->next->next->length = length;
  } 
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
  pll_unode_t ** travbuffer;

  /* we accept only two arguments - the newick tree (unrooted binary) and the
     alignment in the form of FASTA reads */
  if (argc != 3)
    fatal(" syntax: %s [newick] [fasta]", argv[0]);

  /* parse the unrooted binary tree in newick format, and store the number
     of tip nodes in tip_nodes_count */
  pll_utree_t * tree = pll_utree_parse_newick(argv[1]);
  if (!tree)
    fatal("Tree must be an unrooted binary tree");

  tip_nodes_count = tree->tip_count;

  /* fix all missing branch lengths (i.e. those that did not appear in the
     newick) to 0.000001 */
  set_missing_branch_length(tree, 0.000001);


  /* compute and show node count information */
  inner_nodes_count = tip_nodes_count - 2;
  nodes_count = inner_nodes_count + tip_nodes_count;
  branch_count = nodes_count - 1;

  printf("Number of tip/leaf nodes in tree: %d\n", tip_nodes_count);
  printf("Number of inner nodes in tree: %d\n", inner_nodes_count);
  printf("Total number of nodes in tree: %d\n", nodes_count);
  printf("Number of branches in tree: %d\n", branch_count);

  /* Uncomment to display the parsed tree ASCII tree together with information
     as to which CLV index, branch length and label is associated with each
     node. The code will also write (and print on screen) the newick format
     of the tree.

  pll_utree_show_ascii(tree->nodes[nodes_count-1],
                       PLL_UTREE_SHOW_LABEL |
                       PLL_UTREE_SHOW_BRANCH_LENGTH |
                       PLL_UTREE_SHOW_CLV_INDEX);
  char * newick = pll_utree_export_newick(tree->nodes[nodes_count-1],NULL);
  printf("%s\n", newick);
  free(newick);

  */

  /* create a libc hash table of size tip_nodes_count */
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
  PLL_ATTRIB_ARCH_SSE : list of flags for hardware acceleration (not yet implemented)

  */

  partition = pll_partition_create(tip_nodes_count,
                                   inner_nodes_count,
                                   STATES,
                                   (unsigned int)sites,
                                   1,
                                   branch_count,
                                   RATE_CATS,
                                   inner_nodes_count,
                                   PLL_ATTRIB_ARCH_AVX);

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

    pll_set_tip_states(partition, tip_clv_index, pll_map_aa, seqdata[i]);
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
  travbuffer = (pll_unode_t **)xmalloc(nodes_count * sizeof(pll_unode_t *));

  branch_lengths = (double *)xmalloc(branch_count * sizeof(double));
  matrix_indices = (unsigned int *)xmalloc(branch_count * sizeof(unsigned int));
  operations = (pll_operation_t *)xmalloc(inner_nodes_count *
                                          sizeof(pll_operation_t));

  /* compute a partial traversal starting from the randomly selected
     inner node */
  pll_unode_t * root = tree->nodes[tip_nodes_count+inner_nodes_count-1];
  unsigned int traversal_size;

  if (!pll_utree_traverse(root,
                          PLL_TREE_TRAVERSE_POSTORDER,
                          cb_full_traversal,
                          travbuffer,
                          &traversal_size))
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

  /* we'll use 4 rate categories, and currently initialize them to 0 */
  double rate_cats[4] = {0};

  /* compute the discretized category rates from a gamma distribution
     with alpha shape 1 and store them in rate_cats  */
  pll_compute_gamma_cats(1, 4, rate_cats, PLL_GAMMA_RATES_MEAN);

  /* set rate categories */
  pll_set_category_rates(partition, rate_cats);


  /* iterate through all protein models and compute the log-likelihood */
  for (i = 0; i < PROT_MODELS_COUNT; ++i)
  {

    /* set frequencies at model with index 0 (we currently have only one model) */
    pll_set_frequencies(partition, 0, protein_models_freqs_list[i]);

    /* set 6 substitution parameters at model with index 0 */
    pll_set_subst_params(partition, 0, protein_models_rates_list[i]);

    /* update matrix_count probability matrices using the rate matrix with
       index 0. The i-th matrix (i ranges from 0 to matrix_count - 1) is
       generated using branch length branch_lengths[i] and rate matrix
       (substitution rates + frequencies) params_indices[i], and can be refered
       to with index matrix_indices[i] */
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

    /* use the operations array to compute all tip_count-2 inner CLVs. Operations
       will be carried out sequentially starting from operation 0 and upwards */
    pll_update_partials(partition, operations, ops_count);

    /* Uncomment to print on screen the CLVs at tip and inner nodes. From 0 to
       tip_count-1 are tip CLVs, the rest are inner node CLVs.

    for (i = 0; i < 2*tip_count-2; ++i)
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
     index for the concrete branch length, and the array of indices of rate matrix
     whose frequency vector is to be used for each rate category */

  double logl = pll_compute_edge_loglikelihood(partition,
                                               root->clv_index,
                                               root->scaler_index,
                                               root->back->clv_index,
                                               root->back->scaler_index,
                                               root->pmatrix_index,
                                               params_indices,
                                               NULL);

    printf("Log-L (%s): %f\n", protein_models_names_list[i], logl);
  }

  /* destroy all structures allocated for the concrete PLL partition instance */
  pll_partition_destroy(partition);

  free(travbuffer);
  free(branch_lengths);
  free(matrix_indices);
  free(operations);

  /* we will no longer need the tree structure */
  pll_utree_destroy(tree,NULL);

  return (EXIT_SUCCESS);
}
