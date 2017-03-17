#include "common.h"
#include "rng.h"

#include <search.h>

#define STATES      4
#define N_RATE_CATS 4

#define FASTAFILE "testdata/246x4465.fas"
#define TREEFILE  "testdata/246x4465.tree"

typedef struct
{
  int clv_valid;
} node_info_t;

/* a callback function for performing a partial traversal */
static int cb_partial_traversal(pll_utree_t * node)
{
  node_info_t * node_info;

  /* if we don't want tips in the traversal we must return 0 here. For now,
     allow tips */
  if (!node->next) return 1;

  /* get the data element from the node and check if the CLV vector is
     oriented in the direction that we want to traverse. If the data
     element is not yet allocated then we allocate it, set the direction
     and instruct the traversal routine to place the node in the traversal array
     by returning 1 */
  node_info = (node_info_t *)(node->data);
  if (!node_info)
  {
    /* allocate data element */
    node->data             = (node_info_t *)calloc(1,sizeof(node_info_t));
    node->next->data       = (node_info_t *)calloc(1,sizeof(node_info_t));
    node->next->next->data = (node_info_t *)calloc(1,sizeof(node_info_t));

    /* set orientation on selected direction and traverse the subtree */
    node_info = node->data;
    node_info->clv_valid = 1;
    return 1;
  }

  /* if the data element was already there and the CLV on this direction is
     set, i.e. the CLV is valid, we instruct the traversal routine not to
     traverse the subtree rooted in this node/direction by returning 0 */
  if (node_info->clv_valid) return 0;

  /* otherwise, set orientation on selected direction */
  node_info->clv_valid = 1;

  /* reset orientation on the other two directions and return 1 (i.e. traverse
     the subtree */
  node_info = node->next->data;
  node_info->clv_valid = 0;
  node_info = node->next->next->data;
  node_info->clv_valid = 0;

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

int main(int argc, char * argv[])
{
  unsigned int i,j,r;
  unsigned int tip_nodes_count, inner_nodes_count, nodes_count, branch_count;
  unsigned int matrix_count, ops_count;
  unsigned int * matrix_indices;
  double * branch_lengths;
  pll_partition_t * partition;
  pll_operation_t * operations;
  pll_utree_t ** travbuffer;
  pll_utree_t ** inner_nodes_list;
  unsigned int params_indices[N_RATE_CATS] = {0,0,0,0};

  /* parse the unrooted binary tree in newick format, and store the number
     of tip nodes in tip_nodes_count */
  pll_utree_t * tree = pll_utree_parse_newick(TREEFILE, &tip_nodes_count);

  unsigned int attributes = get_attributes(argc, argv);

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

  /* create the PLL partition instance */
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
  pll_compute_gamma_cats(1, 4, rate_cats);

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
  matrix_indices = (unsigned int *)malloc(branch_count * sizeof(int));
  operations = (pll_operation_t *)malloc(inner_nodes_count *
                                                sizeof(pll_operation_t));

  /* get inner nodes */
  inner_nodes_list = (pll_utree_t **)malloc(inner_nodes_count *
                                                sizeof(pll_utree_t *));
  pll_utree_query_innernodes(tree, inner_nodes_list);

  /* get random directions for each inner node */
  for (i = 0; i < inner_nodes_count; ++i)
  {
    r = RAND % 3;
    for (j = 0; j < r; j++)
      inner_nodes_list[i] = inner_nodes_list[i]->next;
  }

  double cmplogl = 0.0;
  for (i = 0; i < 20; ++i)
  {
    unsigned int traversal_size;

    /* randomly select an inner node */
    r = RAND % inner_nodes_count;
    pll_utree_t * node = inner_nodes_list[r];

    /* compute a partial traversal starting from the randomly selected
       inner node */

    if (!pll_utree_traverse(node,
                            cb_partial_traversal,
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




    printf("\nComputing logL between CLV %d and %d - "
           "(pmatrix %d with branch length %f)\n",
            node->clv_index,
            node->back->clv_index,
            node->pmatrix_index,
            node->length);

    printf ("Traversal size: %d\n", traversal_size);
    printf ("Operations: %d\n", ops_count);
    printf ("Matrices: %d\n", matrix_count);

    /* update matrix_count probability matrices for model with index 0. The i-th
       matrix (i ranges from 0 to matrix_count - 1) is generated using branch
       length branch_lengths[i] and can be refered to with index
       matrix_indices[i] */
    pll_update_prob_matrices(partition,
                             params_indices,
                             matrix_indices,
                             branch_lengths,
                             matrix_count);

    /* use the operations array to compute all ops_count inner CLVs. Operations
       will be carried out sequentially starting from operation 0 towrds ops_count-1 */
    pll_update_partials(partition, operations, ops_count);

    /* compute the likelihood on an edge of the unrooted tree by specifying
       the CLV indices at the two end-point of the branch, the probability matrix
       index for the concrete branch length, and the index of the model of whose
       frequency vector is to be used */
    double logl = pll_compute_edge_loglikelihood(partition,
                                                 node->clv_index,
                                                 node->scaler_index,
                                                 node->back->clv_index,
                                                 node->back->scaler_index,
                                                 node->pmatrix_index,
                                                 params_indices,
                                                 NULL);

    if (cmplogl >= -1)
      cmplogl = logl;
    else
      if (fabs(cmplogl - logl) > 1e-5)
        printf("Log-L differs!\n");

    printf("Log-L: %f\n", logl);
  }

  /* deallocate the data elements at inner nodes */
  for (i = 0; i < inner_nodes_count; ++i)
  {
    free(inner_nodes_list[i]->data);
    free(inner_nodes_list[i]->next->data);
    free(inner_nodes_list[i]->next->next->data);
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
  pll_utree_destroy(tree,NULL);

  return (EXIT_SUCCESS);
}
