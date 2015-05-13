#include "pll.h"
#include <stdarg.h>
#include <search.h>


#define STATES    4
#define RATE_CATS 4


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
  set_missing_branch_length_recursive(tree, 0.000001);
  set_missing_branch_length_recursive(tree->back, 0.000001);
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
  int i, tip_count;
  pll_partition_t * partition;
  pll_operation_t * operations = NULL;
  double * branch_lengths= NULL;
  int * matrix_indices = NULL;

  if (argc != 3)
    fatal(" syntax: %s [newick] [fasta]", argv[0]);

  pll_utree_t * tree = pll_parse_newick_utree(argv[1], &tip_count);
  
  /* fix all missing branch lengths (i.e. those that did not appear in the 
     newick) to 0.000001 */
  set_missing_branch_length(tree, 0.000001);


  /* Uncomment to display ASCII tree and newick format

  printf("Number of tips in tree: %d\n", tip_count);
  pll_show_ascii_utree(tree);
  char * newick = pll_write_newick_utree(tree);
  printf("%s\n", newick);
  free(newick);

  */

  /*  obtain an array of pointers to tip names */
  char ** tipnames = pll_query_utree_tipnames(tree, tip_count);

  /* create a libc hash table of size tip_count */
  hcreate(tip_count);

  /* populate a libc hash table with tree tip labels */
  int * data = (int *)malloc(tip_count * sizeof(int));
  for (i = 0; i < tip_count; ++i)
  {
    data[i] = i;
    ENTRY entry;
    entry.key = tipnames[i];
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
  char ** headers = (char **)calloc(tip_count, sizeof(char *));
  char ** seqdata = (char **)calloc(tip_count, sizeof(char *));

  /* read FASTA sequences and make sure they are all of the same length */
  int sites = -1;
  for (i = 0; pll_fasta_getnext(fp,&hdr,&hdrlen,&seq,&seqlen,&seqno); ++i)
  {
    if (i >= tip_count)
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

  if (i != tip_count)
    fatal("Some taxa are missing from FASTA file");

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

  partition = pll_create_partition(tip_count,
                                   tip_count-2,
                                   STATES,
                                   sites,
                                   1,
                                   2*tip_count-3,
                                   RATE_CATS,
                                   tip_count-2,
                                   PLL_ATTRIB_ARCH_SSE);

  /* find sequences in hash table and link them with the corresponding taxa */
  for (i = 0; i < tip_count; ++i)
  {
    ENTRY query;
    query.key = headers[i];
    ENTRY * found = NULL;

    found = hsearch(query,FIND);

    if (!found)
      fatal("Sequence with header %s does not appear in the tree", hdr);
        
    int tip_clv_index = *((int *)(found->data));

    pll_set_tip_states(partition, tip_clv_index, pll_map_nt, seqdata[i]);
  }

  /* destroy hash table */
  hdestroy();

  /* we no longer need these two arrays (keys and values of hash table... */
  free(data);
  free(tipnames);

  /* ...neither the sequences and the headers as they are already
     present in the form of probabilities in the tip CLVs */
  for(i = 0; i < tip_count; ++i)
  {
    free(seqdata[i]);
    free(headers[i]);
  }
  free(seqdata);
  free(headers);



  int edge_pmatrix_index;
  int clv1;
  int clv2;

  /* We perform a simple traversal on the unrooted tree topology. The following
     function allocates branch_lengths, matrix_indices and operations in case
     they were set to NULL, and fills them up.  branch_lengths will contain the
     branch lengths of the postorder traversal starting from 'tree' and
     visiting its three children in order tree->back, tree->next->back and
     tree->next->next->back. Note that tree must be an inner ternary node.
     matrix_indices will assign a number to each branch length, (0 to
     tip_count-1) for branches leading to tips, and (tip_count to
     2*tip_count-4) for the other branches. These numbers are used to assign
     slots for the probability matrix of each branch.  operations will be
     filled to compute all inner CLVs. Finally, edge_matrix_index will point to
     the edge between tree and tree->back, clv1 is set to the index of
     tree->back and clv2 to the index of tree. These last three parameters can
     then be used to evaluate the log-likelihood using the
     pll_compute_edge_loglikelihood function 
  */
  pll_traverse_utree(tree, 
                     tip_count, 
                     &branch_lengths, 
                     &matrix_indices, 
                     &operations, 
                     &edge_pmatrix_index, 
                     &clv1, 
                     &clv2);

  /* we will no longer need the tree structure */
  pll_destroy_utree(tree);

  /* initialize the array of base frequencies */
  double frequencies[4] = { 0.17, 0.19, 0.25, 0.39 };

  /* substitution rates for the 4x4 GTR model. This means we need exactly
     (4*4-4)/2 = 6 values, i.e. the number of elements above the diagonal */
  double subst_params[6] = {1,1,1,1,1,1};

  /* we'll use 4 rate categories, and currently initialize them to 0 */
  double rate_cats[4] = {0};

  /* compute the discretized category rates from a gamma distribution
     with alpha shape 1 and store them in rate_cats  */
  pll_compute_gamma_cats(1, 4, rate_cats);

  /* set frequencies at model with index 0 (we currently have only one model) */
  pll_set_frequencies(partition, 0, frequencies);

  /* set 6 substitution parameters at model with index 0 */
  pll_set_subst_params(partition, 0, subst_params, 6);

  /* set rate categories */
  pll_set_category_rates(partition, rate_cats);

  /* update 2*tip_count-3 probability matrices for model with index 0. The i-th
     matrix (i ranges from 0 to 2*tip_count-4) is generated using branch length
     branch_lengths[i] and can be refered to with index matrix_indices[i] */
  pll_update_prob_matrices(partition, 
                           0, 
                           matrix_indices, 
                           branch_lengths, 
                           2*tip_count-3);

  /* Uncomment to output the probability matrices (for each branch and each rate
     category) on screen 
  for (i = 0; i < 2*tip_count-3; ++i)
  {
    printf ("P-matrix (%d) for branch length %f\n", i, branch_lengths[i]);
    pll_show_pmatrix(partition, i,17);
    printf ("\n");
  }
  
  */

  /* use the operations array to compute all tip_count-2 inner CLVs. Operations
     will be carried out sequentially starting from operation 0 and upwards */
  pll_update_partials(partition, operations, tip_count-2);

  /* Uncomment to print on screen the CLVs at tip and inner nodes. From 0 to
     tip_count-1 are tip CLVs, the rest are inner node CLVs.

  for (i = 0; i < 2*tip_count-2; ++i)
  {
    printf ("CLV %d: ", i);
    pll_show_clv(partition,i,17);
  }

  */

  /* compute the likelihood on an edge of the unrooted tree by specifying
     the CLV indices at the two end-point of the branch, the probability matrix
     index for the concrete branch length, and the index of the model of whose
     frequency vector is to be used */
  double logl = pll_compute_edge_loglikelihood(partition,
                                               clv1,
                                               clv2,
                                               edge_pmatrix_index,
                                               0);

  printf("Log-L: %f\n", logl);

  /* destroy all structures allocated for the concrete PLL partition instance */
  pll_destroy_partition(partition);

  free(branch_lengths);
  free(matrix_indices);
  free(operations);

  return (EXIT_SUCCESS);
}
