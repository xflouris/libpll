#include "pll.h"
#include <stdarg.h>
#include <search.h>


#define STATES    20
#define MIXTURE 4
#define RATE_CATS MIXTURE
#define PROT_MODELS_COUNT 19

static void fatal(const char * format, ...) __attribute__ ((noreturn));

/* a callback function for performing a full traversal */
static int cb_full_traversal(pll_utree_t * node)
{
  return 1;
}

static void set_missing_branch_length_recursive(pll_utree_t * tree, 
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

      set_missing_branch_length_recursive(tree->next->back, length);
      set_missing_branch_length_recursive(tree->next->next->back, length);
    }
  }
}

/* branch lengths not present in the newick file get a value of 0.000001 */
static void set_missing_branch_length(pll_utree_t * tree, double length)
{
  set_missing_branch_length_recursive(tree, length);
  set_missing_branch_length_recursive(tree->back, length);
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
  pll_utree_t ** travbuffer;

  /* we accept only two arguments - the newick tree (unrooted binary) and the
     alignment in the form of FASTA reads */
  if (argc != 3)
    fatal(" syntax: %s [newick] [fasta]", argv[0]);

  /* parse the unrooted binary tree in newick format, and store the number
     of tip nodes in tip_nodes_count */
  pll_utree_t * tree = pll_utree_parse_newick(argv[1], &tip_nodes_count);
  
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

  /*
  pll_utree_show_ascii(tree, PLL_UTREE_SHOW_LABEL |
                             PLL_UTREE_SHOW_BRANCH_LENGTH |
                             PLL_UTREE_SHOW_CLV_INDEX);
  char * newick = pll_utree_export_newick(tree);
  printf("%s\n", newick);
  free(newick);
 // */

  /*  obtain an array of pointers to tip nodes */
  pll_utree_t ** tipnodes = (pll_utree_t  **)calloc(tip_nodes_count, 
                                                    sizeof(pll_utree_t *));
  pll_utree_query_tipnodes(tree, tipnodes);

  /* create a libc hash table of size tip_nodes_count */
  hcreate(tip_nodes_count);

  /* populate a libc hash table with tree tip labels */
  unsigned int * data = (unsigned int *)malloc(tip_nodes_count * 
                                               sizeof(unsigned int));
  for (i = 0; i < tip_nodes_count; ++i)
  {
    data[i] = i;
    ENTRY entry;
    entry.key = tipnodes[i]->label;
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
  4 : mixture size
  4 : number of different substitution models (or eigen decomposition)
      to use concurrently (multiple of mixture size)
  branch_count: number of probability matrices to be allocated
  RATE_CATS : number of rate categories we will use
  inner_nodes_count : how many scale buffers to use
  PLL_ATTRIB_ARCH_SSE : list of flags for hardware acceleration
  PLL_ATTRIB_MIXT_LINKED : flag for linking the matrices to the gamma rates
  */

  partition = pll_partition_create(tip_nodes_count,
                                   inner_nodes_count,
                                   STATES,
                                   (unsigned int)sites,
                                   MIXTURE,
                                   1,
                                   branch_count,
                                   RATE_CATS,
                                   inner_nodes_count,
                                   PLL_ATTRIB_ARCH_CPU | PLL_ATTRIB_MIXT_LINKED);

  /* find sequences in hash table and link them with the corresponding taxa */
  for (i = 0; i < tip_nodes_count; ++i)
  {
    ENTRY query;
    query.key = headers[i];
    ENTRY * found = NULL;

    found = hsearch(query,FIND);

    if (!found)
      fatal("Sequence with header %s does not appear in the tree", hdr);
        
    unsigned int tip_clv_index = *((unsigned int *)(found->data));

    pll_set_tip_states(partition, tip_clv_index, pll_map_aa, seqdata[i]);
  }

  /* destroy hash table */
  hdestroy();

  /* we no longer need these two arrays (keys and values of hash table... */
  free(data);
  free(tipnodes);

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
  travbuffer = (pll_utree_t **)malloc(nodes_count * sizeof(pll_utree_t *));

  branch_lengths = (double *)malloc(branch_count * sizeof(double));
  matrix_indices = (unsigned int *)malloc(branch_count * sizeof(unsigned int));
  operations = (pll_operation_t *)malloc(inner_nodes_count * 
                                                sizeof(pll_operation_t));

  /* compute a partial traversal starting from the randomly selected 
     inner node */
  unsigned int traversal_size;
  if (!pll_utree_traverse(tree, 
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
  pll_compute_gamma_cats(1.00000, 4, rate_cats);

  /* set rate categories */
  pll_set_category_rates(partition, rate_cats);

  for (i=0; i<partition->mixture; i++)
  {
    /* set frequencies at model with index i */
    pll_set_frequencies(partition, 0, i, pll_aa_freqs_lg4m[i]);

    /* set substitution parameters at model with index i */
    pll_set_subst_params(partition, 0, i, pll_aa_rates_lg4m[i]);
  }

    /* update 2*tip_count-3 probability matrices for model with index 0. The i-th
       matrix (i ranges from 0 to 2*tip_count-4) is generated using branch length
       branch_lengths[i] and can be refered to with index matrix_indices[i] */
    pll_update_prob_matrices(partition, 
                             0, 
                             matrix_indices, 
                             branch_lengths, 
                             matrix_count);

    /*
    for (i = 0; i < branch_count; ++i)
    {
      printf ("P-matrices (%d) for branch length %f\n", i, branch_lengths[i]);
      pll_show_pmatrix(partition, i,4);
      printf ("\n");
    }
    */

    /* use the operations array to compute all tip_count-2 inner CLVs. Operations
       will be carried out sequentially starting from operation 0 and upwards */
    pll_update_partials(partition, operations, ops_count);

  /* compute the likelihood on an edge of the unrooted tree by specifying
     the CLV indices at the two end-point of the branch, the probability matrix
     index for the concrete branch length, and the index of the model of whose
     frequency vector is to be used */
  double logl = pll_compute_edge_loglikelihood(partition,
                                               tree->clv_index,
                                               tree->scaler_index,
                                               tree->back->clv_index,
                                               tree->back->scaler_index,
                                               tree->pmatrix_index,
                                               0);

  printf("\nRates:   ");
      for (i = 0; i < partition->mixture; i++)
        printf("%.6f ", partition->rates[i]);
  printf("\nLog-L (LG4M): %f\n", logl);

  /* now let's switch into LG4X model with fixed rates and weights */

  for (i = 0; i < partition->mixture; i++)
  {
    /* set frequencies at model with index i */
    pll_set_frequencies (partition, 0, i, pll_aa_freqs_lg4x[i]);

    /* set substitution parameters at model with index i */
    pll_set_subst_params (partition, 0, i, pll_aa_rates_lg4x[i]);
  }

  /* set free rates and weights */
  double fixedweights[] =
    { 0.209224645, 0.224707726, 0.277599198, 0.288468431 };
  double fixedrates[] =
    { 0.498991136, 0.563680734, 0.808264032, 1.887769458 };
  pll_set_category_rates (partition, fixedrates);
  pll_set_category_weights (partition, fixedweights);

  pll_update_prob_matrices(partition,
                           0,
                           matrix_indices,
                           branch_lengths,
                           matrix_count);

  /* Uncomment to display the P matrices

  for (i = 0; i < branch_count; ++i)
  {
      printf ("P-matrices (%d) for branch length %f\n", i, branch_lengths[i]);
      pll_show_pmatrix (partition, i, 4);
      printf ("\n");
  }
  */

  pll_update_partials(partition, operations, ops_count);

  logl = pll_compute_edge_loglikelihood(partition,
                                        tree->clv_index,
                                        tree->scaler_index,
                                        tree->back->clv_index,
                                        tree->back->scaler_index,
                                        tree->pmatrix_index,
                                        0);

  printf("\nWeights: ");
  for (i = 0; i < partition->mixture; i++)
    printf("%.6f ", partition->rate_weights[i]);
  printf("\nRates:   ");
    for (i = 0; i < partition->mixture; i++)
      printf("%.6f ", partition->rates[i]);
  printf("\nLog-L (LG4X): %f\n", logl);

  /* destroy all structures allocated for the concrete PLL partition instance */
  pll_partition_destroy(partition);

  free(travbuffer);
  free(branch_lengths);
  free(matrix_indices);
  free(operations);

  /* we will no longer need the tree structure */
  pll_utree_destroy(tree);

  return (EXIT_SUCCESS);
}
