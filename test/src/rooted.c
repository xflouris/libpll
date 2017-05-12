#include "common.h"

#include <stdarg.h>
#include <search.h>

#define STATES      4
#define N_RATE_CATS 4

#define FASTAFILE     "testdata/small.fas"
#define TREEFILE      "testdata/small.rooted.tree"
#define TREEFILE_TIP  "testdata/small.rooted.tip.tree"

static double prop_invar_list[4] = {0.0, 0.1, 0.5, 0.9};

typedef struct
{
  int clv_valid;
} node_info_t;

int main(int argc, char * argv[])
{
  unsigned int i, j;
  unsigned int tip_nodes_count, inner_nodes_count, nodes_count, branch_count;
  unsigned int matrix_count, ops_count;
  unsigned int * matrix_indices;
  double * branch_lengths;
  pll_partition_t * partition;
  pll_operation_t * operations;
  pll_rnode_t ** travbuffer;
  pll_rnode_t ** inner_nodes_list;
  unsigned int params_indices[N_RATE_CATS] = {0,0,0,0};

  /* parse the unrooted binary tree in newick format, and store the number
     of tip nodes in tip_nodes_count */
  pll_rtree_t * tree = pll_rtree_parse_newick(TREEFILE);
  if (!tree)
  {
    printf("Error reading tree\n");
    exit(1);
  }

  tip_nodes_count = tree->tip_count;

  unsigned int attributes = get_attributes(argc, argv);

  /* compute and show node count information */
  inner_nodes_count = tip_nodes_count - 1;
  nodes_count = inner_nodes_count + tip_nodes_count;
  branch_count = nodes_count - 1;

  printf("Number of tip/leaf nodes in tree: %d\n", tip_nodes_count);
  printf("Number of inner nodes in tree: %d\n", inner_nodes_count);
  printf("Total number of nodes in tree: %d\n", nodes_count);
  printf("Number of branches in tree: %d\n", branch_count);

  pll_rtree_show_ascii(tree->root,
                       PLL_UTREE_SHOW_LABEL |
                       PLL_UTREE_SHOW_BRANCH_LENGTH |
                       PLL_UTREE_SHOW_CLV_INDEX);
  char * newick = pll_rtree_export_newick(tree->root,NULL);
  printf("%s\n", newick);
  free(newick);

  /* create a libc hash table of size tip_nodes_count */
  hcreate(tip_nodes_count);

  /* populate a libc hash table with tree tip labels */
  unsigned int * data = (unsigned int *)malloc(tip_nodes_count *
                                               sizeof(unsigned int));
  for (i = 0; i < tip_nodes_count; ++i)
  {
    data[i] = i;
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
  pll_fasta_t * fp = pll_fasta_open(FASTAFILE, pll_map_fasta);
  if (!fp)
    fatal("Error opening file");

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
    fatal("Error while reading file");

  /* close FASTA file */
  pll_fasta_close(fp);

  if (sites == -1)
    fatal("Unable to read alignment");

  if (i != tip_nodes_count)
    fatal("Some taxa are missing from FASTA file");

  partition = pll_partition_create(tip_nodes_count,
                                   inner_nodes_count,
                                   STATES,
                                   (unsigned int)sites,
                                   1,
                                   branch_count,
                                   N_RATE_CATS,
                                   inner_nodes_count,
                                   attributes
                                   );

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
  travbuffer = (pll_rnode_t **)malloc(nodes_count * sizeof(pll_rnode_t *));

  branch_lengths = (double *)malloc(branch_count * sizeof(double));
  matrix_indices = (unsigned int *)malloc(branch_count * sizeof(int));
  operations = (pll_operation_t *)malloc(inner_nodes_count *
                                                sizeof(pll_operation_t));

  /* get inner nodes */
  inner_nodes_list = (pll_rnode_t **)malloc(inner_nodes_count *
                                                sizeof(pll_rnode_t *));
  memcpy(inner_nodes_list,
         tree->nodes+tip_nodes_count,
         inner_nodes_count*sizeof(pll_rnode_t *));

  unsigned int traversal_size;

  /* compute a partial traversal starting from the randomly selected
     inner node */

  if (!pll_rtree_traverse(tree->root,
                          PLL_TREE_TRAVERSE_POSTORDER,
                          cb_rfull_traversal,
                          travbuffer,
                          &traversal_size))
    fatal("Function pll_rtree_traverse() root node as parameter");

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
  printf ("Matrices: %d\n", matrix_count);

  for (j=0; j<4; ++j)
  {
    pll_update_invariant_sites_proportion(partition,
                                          0,
                                          prop_invar_list[j]);
    /* update matrix_count probability matrices for model with index 0. The i-th
       matrix (i ranges from 0 to matrix_count - 1) is generated using branch
       length branch_lengths[i] and can be refered to with index
       matrix_indices[i] */
    pll_update_prob_matrices(partition,
                             params_indices,
                             matrix_indices,
                             branch_lengths,
                             matrix_count);

    for (i = 0; i < branch_count; ++i)
    {
      printf ("P-matrix (%d) for branch length %f\n", i, branch_lengths[i]);
      pll_show_pmatrix(partition, i,6);
      printf ("\n");
    }

    /* use the operations array to compute all ops_count inner CLVs. Operations
       will be carried out sequentially starting from operation 0 towrds ops_count-1 */
    pll_update_partials(partition, operations, ops_count);

    /* compute the likelihood on an edge of the unrooted tree by specifying
       the CLV indices at the two end-point of the branch, the probability matrix
       index for the concrete branch length, and the index of the model of whose
       frequency vector is to be used */
    double logl = pll_compute_root_loglikelihood(partition,
                                                 tree->root->clv_index,
                                                 tree->root->scaler_index,
                                                 params_indices,
                                                 NULL);
    printf("Log-L: %f (pinv = %f)\n", logl, prop_invar_list[j]);
  }



  /* deallocate the inner nodes list */
  free(inner_nodes_list);

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
