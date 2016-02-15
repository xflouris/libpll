#include "pll_optimize.h"
#include "pll_tree.h"
#include <stdarg.h>
#include <search.h>
#include <time.h>

#define STATES    4
#define RATE_CATS 4

#define SHOW_ASCII_TREE 0

static void fatal (const char * format, ...) __attribute__ ((noreturn));

const char * msa_filename;
const char * tree_filename;

const static float branch_opt_epsilon = 1e-2;
const static int branch_opt_smoothings = 2;

typedef struct
{
  int clv_valid;
} node_info_t;

/* a callback function for performing a full traversal */
static int cb_full_traversal (pll_utree_t * node)
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

static void show_tree (pll_utree_t * tree)
{
#if(SHOW_ASCII_TREE)
  printf ("\n");
  pll_utree_show_ascii (tree, PLL_UTREE_SHOW_LABEL |
  PLL_UTREE_SHOW_BRANCH_LENGTH |
  PLL_UTREE_SHOW_CLV_INDEX);
  char * newick = pll_utree_export_newick (tree);
  printf ("%s\n\n", newick);
  free (newick);
#else
  printf ("ASCII tree not shown (SHOW_ASCII_TREE flag)\n");
  return;
#endif
}

int main (int argc, char * argv[])
{
  unsigned int i;
  unsigned int tip_nodes_count, inner_nodes_count, nodes_count, branch_count;
  unsigned int matrix_count, ops_count;
  unsigned int * matrix_indices;
  double * branch_lengths;
  pll_partition_t * partition;
  pll_operation_t * operations;
  pll_utree_t ** travbuffer;

  if (argc != 3)
    fatal (" syntax: %s [newick] [fasta]", argv[0]);
  tree_filename = argv[1];
  msa_filename = argv[2];

  /* parse the unrooted binary tree in newick format, and store the number
   of tip nodes in tip_nodes_count */
  printf("Parsing tree: %s\n", tree_filename);
  pll_utree_t * tree = pll_utree_parse_newick (tree_filename, &tip_nodes_count);
  if (!tree)
    fatal ("%s does not exist", tree_filename);

  set_missing_branch_length (tree, 0.000001);

  /* compute and show node count information */
  inner_nodes_count = tip_nodes_count - 2;
  nodes_count = inner_nodes_count + tip_nodes_count;
  branch_count = nodes_count - 1;

  printf ("  Number of tip/leaf nodes in tree: %d\n", tip_nodes_count);
  printf ("  Number of inner nodes in tree: %d\n", inner_nodes_count);
  printf ("  Total number of nodes in tree: %d\n", nodes_count);
  printf ("  Number of branches in tree: %d\n", branch_count);

  /*  obtain an array of pointers to tip and inner nodes */
  pll_utree_t ** tipnodes = (pll_utree_t **) calloc (tip_nodes_count,
                                                     sizeof(pll_utree_t *));
  pll_utree_t ** innernodes = (pll_utree_t **) calloc (inner_nodes_count,
                                                       sizeof(pll_utree_t *));
  pll_utree_query_tipnodes (tree, tipnodes);
  pll_utree_query_innernodes (tree, innernodes);

  /* place the virtual root at a random inner node */
  tree = innernodes[rand() % inner_nodes_count];

  show_tree (tree);

  /* create a libc hash table of size tip_nodes_count */
  hcreate (tip_nodes_count);

  /* populate a libc hash table with tree tip labels */
  unsigned int * data = (unsigned int *) malloc (
      tip_nodes_count * sizeof(unsigned int));
  for (i = 0; i < tip_nodes_count; ++i)
  {
    data[i] = i;
    ENTRY entry;
    entry.key = tipnodes[i]->label;
    entry.data = (void *) (data + i);
    hsearch (entry, ENTER);
  }

  /* open FASTA file */
  printf("Reading FASTA file: %s\n", msa_filename);
  pll_fasta_t * fp = pll_fasta_open (msa_filename, pll_map_fasta);
  if (!fp)
    fatal ("%s does not exist", msa_filename);

  char * seq = NULL;
  char * hdr = NULL;
  long seqlen;
  long hdrlen;
  long seqno;

  /* allocate arrays to store FASTA headers and sequences */
  char ** headers = (char **) calloc (tip_nodes_count, sizeof(char *));
  char ** seqdata = (char **) calloc (tip_nodes_count, sizeof(char *));

  /* read FASTA sequences and make sure they are all of the same length */
  int sites = -1;
  for (i = 0; pll_fasta_getnext (fp, &hdr, &hdrlen, &seq, &seqlen, &seqno); ++i)
  {
    if (i >= tip_nodes_count)
      fatal ("FASTA file contains more sequences than expected");

    if (sites != -1 && sites != seqlen)
      fatal ("FASTA file does not contain equal size sequences\n");

    if (sites == -1)
      sites = seqlen;

    headers[i] = hdr;
    seqdata[i] = seq;
  }

  /* did we stop reading the file because we reached EOF? */
  if (pll_errno != PLL_ERROR_FILE_EOF)
    fatal ("Error in %s", msa_filename);

  /* close FASTA file */
  pll_fasta_close (fp);

  if (sites == -1)
    fatal ("Unable to read alignment");

  if (i != tip_nodes_count)
    fatal ("Some taxa are missing from FASTA file");

  printf ("  Length of sequences: %d\n", sites);

  /* create the PLL partition instance */
  partition = pll_partition_create (tip_nodes_count,
                                    inner_nodes_count,
                                    STATES,
                                    (unsigned int) sites,
                                    0,
                                    1,
                                    branch_count,
                                    RATE_CATS,
                                    inner_nodes_count,
                                    pll_map_nt,
                                    PLL_ATTRIB_ARCH_CPU);

  /* find sequences in hash table and link them with the corresponding taxa */
  for (i = 0; i < tip_nodes_count; ++i)
  {
    ENTRY query;
    query.key = headers[i];
    ENTRY * found = NULL;

    found = hsearch (query, FIND);

    if (!found)
      fatal ("Sequence with header %s does not appear in the tree", hdr);

    unsigned int tip_clv_index = *((unsigned int *) (found->data));

    pll_set_tip_states (partition, tip_clv_index, pll_map_nt, seqdata[i]);
  }

  /* destroy hash table */
  hdestroy ();

  /* we no longer need these two arrays (keys and values of hash table... */
  free (data);
  free (tipnodes);

  /* ...neither the sequences and the headers as they are already
   present in the form of probabilities in the tip CLVs */
  for (i = 0; i < tip_nodes_count; ++i)
  {
    free (seqdata[i]);
    free (headers[i]);
  }
  free (seqdata);
  free (headers);

  /* initialize base frequencies */
  double * frequencies = pll_compute_empirical_frequencies(partition);
  pll_set_frequencies (partition, 0, 0, frequencies);
  free(frequencies);

  /* initialize substitution rates */
  double * subst_params = pll_compute_empirical_subst_rates(partition);
  pll_set_subst_params (partition, 0, 0, subst_params);
  free(subst_params);

  /* compute the discretized category rates from a gamma distribution
     with alpha shape 1 and store them in rate_cats  */
  double rate_cats[RATE_CATS] = { 0 };
  pll_compute_gamma_cats (1, RATE_CATS, rate_cats);
  pll_set_category_rates (partition, rate_cats);

  printf ("Model paramters:\n");
  printf ("  Frequencies:  ");
  for (i=0; i<STATES; i++)
    printf("%.4f ", partition->frequencies[0][i]);
  printf("\n");
  printf ("  Subst. rates: ");
  for (i=0; i<(STATES*(STATES-1))/2; i++)
    printf("%.4f ", partition->subst_params[0][i]);
  printf("\n");
  printf ("  Gamma rates:  ");
  for (i=0; i<RATE_CATS; i++)
    printf("%.4f ", partition->rates[i]);
  printf("\n");

  /* allocate a buffer for storing pointers to nodes of the tree in postorder
   traversal */
  travbuffer = (pll_utree_t **) malloc (nodes_count * sizeof(pll_utree_t *));

  branch_lengths = (double *) malloc (branch_count * sizeof(double));
  matrix_indices = (unsigned int *) malloc (
      branch_count * sizeof(unsigned int));
  operations = (pll_operation_t *) malloc (
      inner_nodes_count * sizeof(pll_operation_t));

  /* perform a postorder traversal of the unrooted tree */
  unsigned int traversal_size;
  if (!pll_utree_traverse (tree, cb_full_traversal, travbuffer,
                           &traversal_size))
    fatal ("Function pll_utree_traverse() requires inner nodes as parameters");

  /* given the computed traversal descriptor, generate the operations
   structure, and the corresponding probability matrix indices that
   may need recomputing */
  pll_utree_create_operations (travbuffer, traversal_size, branch_lengths,
                               matrix_indices, operations, &matrix_count,
                               &ops_count);

  printf ("Traversal size: %d\n", traversal_size);
  printf ("Operations: %d\n", ops_count);
  printf ("Probability Matrices: %d\n", matrix_count);

  /* update matrix_count probability matrices for model with index 0. The i-th
   matrix (i ranges from 0 to matrix_count - 1) is generated using branch
   length branch_lengths[i] and can be refered to with index
   matrix_indices[i] */
  pll_update_prob_matrices (partition, 0, matrix_indices, branch_lengths,
                            matrix_count);

  /* use the operations array to compute all ops_count inner CLVs. Operations
   will be carried out sequentially starting from operation 0 towards ops_count-1 */
  pll_update_partials (partition, operations, ops_count);

  /* compute the likelihood on an edge of the unrooted tree by specifying
   the CLV indices at the two end-point of the branch, the probability matrix
   index for the concrete branch length, and the index of the model of whose
   frequency vector is to be used */
  double logl = pll_compute_edge_loglikelihood (partition,
                                                tree->clv_index,
                                                tree->scaler_index,
                                                tree->back->clv_index,
                                                tree->back->scaler_index,
                                                tree->pmatrix_index,
                                                0);

  printf ("Log-L at %s-%s: %f\n", tree->label, tree->back->label, logl);

  /* optimize branch lengths */
  printf ("Optimizing branch length parameters\n");
  pll_optimize_branch_lengths_iterative (partition, tree, 0, 0,
                                         branch_opt_epsilon,
                                         branch_opt_smoothings);

  logl = pll_compute_edge_loglikelihood (partition,
                                         tree->clv_index,
                                         tree->scaler_index,
                                         tree->back->clv_index,
                                         tree->back->scaler_index,
                                         tree->pmatrix_index,
                                         0);

  printf ("Log-L* at %s-%s: %f\n", tree->label, tree->back->label, logl);

  /* Test TBR */

  unsigned int distance = 3;
  unsigned int n_nodes_at_dist;
  pll_utree_t ** nodes_at_dist = (pll_utree_t **) malloc (
      sizeof(pll_utree_t *) * pow (2, distance));
  pll_edge_t reconnect;
  pll_utree_t * bisect_edge;

  printf("\n\n");
  printf("Obtaining random bisection edge\n");
  printf("Obtaining random reconnection edges at a distance of %u\n", distance);
  int max_tests = 100;
  while (max_tests > 0)
  {
    max_tests--;
    bisect_edge = innernodes[rand () % inner_nodes_count];

    /* the bisection edge should not be a tip branch */
    if (!(bisect_edge->next && bisect_edge->back->next))
      continue;

    /* find nodes at a certain distance */
    pll_utree_nodes_at_node_dist (bisect_edge, nodes_at_dist, &n_nodes_at_dist,
                                  distance, 1);
    if (!n_nodes_at_dist)
      continue;
    reconnect.edge.utree.parent = nodes_at_dist[rand () % n_nodes_at_dist];
    pll_utree_nodes_at_node_dist (bisect_edge->back, nodes_at_dist,
                                  &n_nodes_at_dist, distance, 1);
    if (!n_nodes_at_dist)
          continue;
    else
      max_tests = -1;
    reconnect.edge.utree.child = nodes_at_dist[rand () % n_nodes_at_dist];
  }

  if (!max_tests)
    fatal ("ERROR: Could not find valid bisection and reconnection points");

  printf ("Tree bisect at %s-%s\n", bisect_edge->label, bisect_edge->back->label);
  printf ("Tree reconnect %s-%s and %s-%s\n",
          reconnect.edge.utree.parent->label,
          reconnect.edge.utree.parent->back->label,
          reconnect.edge.utree.child->label,
          reconnect.edge.utree.child->back->label);
  reconnect.length = 0.555;

  if (!pll_utree_TBR (bisect_edge, &reconnect))
    fatal ("TBR move cannot be applied");

  tree = reconnect.edge.utree.parent;
  free(nodes_at_dist);

  /* if tree is a tip node, move to its neighbor inner node */
  if (!tree->next) tree = tree->back;

  /* validate tree integrity */
  printf ("Integrity check %s... ", tree->label);
  fflush(stdout);
  if (!pll_utree_check_integrity (tree))
    fatal ("Tree is not consistent");
  printf ("OK\n");

  /* traversing the new tree */
  if (!pll_utree_traverse (tree, cb_full_traversal, travbuffer,
                           &traversal_size))
    fatal ("Function pll_utree_traverse() requires inner nodes as parameters");

  pll_utree_create_operations (travbuffer, traversal_size, branch_lengths,
                               matrix_indices, operations, &matrix_count,
                               &ops_count);
  show_tree (reconnect.edge.utree.child);
  pll_update_prob_matrices (partition, 0, matrix_indices, branch_lengths,
                            matrix_count);
  pll_update_partials (partition, operations, ops_count);

  /* compute marginal likelihoods */
  printf("\nMarginal likelihoods:\n");
  logl = pll_compute_root_loglikelihood (partition,
                                         tree->clv_index,
                                         tree->scaler_index,
                                         0);
  printf ("  Log-L Partial at %s: %f\n", tree->label, logl);
  logl = pll_compute_root_loglikelihood (partition,
                                         tree->back->clv_index,
                                         tree->back->scaler_index,
                                         0);
  printf ("  Log-L Partial at %s: %f\n", tree->back->label, logl);

  /* compute global likelihood */
  logl = pll_compute_edge_loglikelihood (partition,
                                         tree->clv_index,
                                         tree->scaler_index,
                                         tree->back->clv_index,
                                         tree->back->scaler_index,
                                         tree->pmatrix_index,
                                         0);

  printf ("  Log-L at %s-%s: %f\n", tree->label, tree->back->label, logl);

  printf ("\nOptimizing branch length parameters locally (radius = 2)\n");
  pll_optimize_branch_lengths_local (partition, tree, 0, 0, 1e-2, 4, 2);
  logl = pll_compute_edge_loglikelihood (partition,
                                         tree->clv_index,
                                         tree->scaler_index,
                                         tree->back->clv_index,
                                         tree->back->scaler_index,
                                         tree->pmatrix_index,
                                         0);

  printf ("Log-L* at %s-%s: %f\n", tree->label, tree->back->label, logl);

  /* destroy all structures allocated for the concrete PLL partition instance */
  pll_partition_destroy (partition);

  /* deallocate traversal buffer, branch lengths array, matrix indices
   array and operations */
  printf ("\nDestroy buffers\n");
  free (innernodes);
  free (travbuffer);
  free (branch_lengths);
  free (matrix_indices);
  free (operations);

  printf ("Destroy tree\n");
  /* we will no longer need the tree structure */
  pll_utree_destroy (tree);

  return (EXIT_SUCCESS);
}
