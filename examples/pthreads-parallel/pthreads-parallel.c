#include "pll.h"

#include <stdarg.h>
#include <search.h>
#include <time.h>
#include <pthread.h>
#include <unistd.h>

#define STATES    4
#define RATE_CATS 4

static void fatal(const char * format, ...) __attribute__ ((noreturn));

typedef struct
{
  int clv_valid;
} node_info_t;

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

typedef struct
{
  long thread_id;
  long num_threads;
  pll_partition_t * partition;
  unsigned int * matrix_indices;
  unsigned int matrix_count;
  pll_operation_t * operations;
  unsigned int ops_count;
  double * branch_lengths;
  pll_utree_t * vroot;
  pthread_barrier_t * barrier_buf;
  double * result_buf;
  int trap;
} thread_data_t;

#define JOB_WAIT             0
#define JOB_UPDATE_MATRICES  1
#define JOB_FULL_LK          2
#define JOB_UPDATE_PARTIALS  3
#define JOB_REDUCE_LK_EDGE   4
#define JOB_FINALIZE        -1

volatile int thread_job = 0;
volatile double global_lnl = 0;

static int barrier(thread_data_t * data)
{
  return pthread_barrier_wait (data->barrier_buf);
}

static void dealloc_partition_local(pll_partition_t * partition)
{
  if (partition->clv)
    free(partition->clv);
  if (partition->scale_buffer)
    free(partition->scale_buffer);
  free(partition);
}

void * worker(void * void_data)
{
  thread_data_t * thread_data;
  unsigned int i;
  unsigned int id, mat_count, mat_offset, mat_start;
  unsigned int * matrix_indices;
  double * branch_lengths;

  thread_data = (thread_data_t *) void_data;
  id = thread_data->thread_id;
  mat_count  = (thread_data->matrix_count / thread_data->num_threads);
  mat_offset = (thread_data->matrix_count % thread_data->num_threads);
  mat_start = id * mat_count + ((id<mat_offset)?id:mat_offset);
  mat_count += (id < mat_offset)?1:0;
  matrix_indices = thread_data->matrix_indices + mat_start;
  branch_lengths = thread_data->branch_lengths + mat_start;

  do
  {
    /* wait */
    switch (thread_job)
      {
      case JOB_FINALIZE:
        {
          /* finalize */
          thread_data->trap = 0;
          break;
        }
      case JOB_UPDATE_MATRICES:
      {
        barrier (thread_data);

          /* check eigen */
          if (!thread_data->partition->eigen_decomp_valid[0])
          {
            barrier (thread_data);
            if (!id)
              pll_update_eigen (thread_data->partition, 0);
            barrier (thread_data);
            if (!id)
              thread_data->partition->eigen_decomp_valid[0] = 1;
          }

          barrier (thread_data);

          pll_update_prob_matrices (thread_data->partition, 0,
                                    matrix_indices,
                                    branch_lengths,
                                    mat_count);

          if (!id)
            thread_job = JOB_WAIT;
          barrier (thread_data);
          break;
        }
      case JOB_UPDATE_PARTIALS:
        {
          barrier (thread_data);

          pll_update_partials (thread_data->partition,
                               thread_data->operations,
                               thread_data->ops_count);
          if (!id)
            thread_job = JOB_WAIT;
          barrier (thread_data);
          break;
        }
      case JOB_REDUCE_LK_EDGE:
        {
          if (!id)
            global_lnl = 0;
          barrier (thread_data);

          thread_data->result_buf[id] =
              pll_compute_edge_loglikelihood (
                thread_data->partition,
                thread_data->vroot->clv_index,
                thread_data->vroot->scaler_index,
                thread_data->vroot->back->clv_index,
                thread_data->vroot->back->scaler_index,
                thread_data->vroot->pmatrix_index, 0);
#if(0)
          printf ("Thread %ld LogLK %f\n", id, logl);
#endif

          /* reduce */
          barrier (thread_data);
          if (!id)
            for (i=0; i<thread_data->num_threads; i++)
              global_lnl += thread_data->result_buf[i];
          barrier (thread_data);

          if (!id)
            thread_job = JOB_WAIT;
          barrier (thread_data);
          break;
        }
      case JOB_FULL_LK:
        {
          /* execute */
          if (!id)
            global_lnl = 0;
          barrier (thread_data);

          pll_update_partials (thread_data->partition, thread_data->operations,
                               thread_data->ops_count);

          thread_data->result_buf[id] =
              pll_compute_edge_loglikelihood (
                thread_data->partition,
                thread_data->vroot->clv_index,
                thread_data->vroot->scaler_index,
                thread_data->vroot->back->clv_index,
                thread_data->vroot->back->scaler_index,
                thread_data->vroot->pmatrix_index, 0);
#if(0)
          printf ("Thread %ld LogLK %f\n", id, logl);
#endif

          /* reduce */
          barrier (thread_data);
          if (!id)
            for (i = 0; i < thread_data->num_threads; i++)
              global_lnl += thread_data->result_buf[i];
          barrier (thread_data);

          if (!id)
            thread_job = JOB_WAIT;
          barrier (thread_data);
          break;
        }
      }
  } while (thread_data->trap);

  return 0;
}

static pll_partition_t * pll_partition_clone_partial(
                             pll_partition_t * partition,
                             unsigned int id,
                             unsigned int count)
{
  unsigned int i;
  unsigned int start, sites, offset;
  pll_partition_t * new_partition;

  sites = (partition->sites / count);
  offset = (partition->sites % count);
  start = sites * id + ((id<offset)?id:offset);
  sites += (id < offset)?1:0;

  new_partition = (pll_partition_t *) malloc(sizeof(pll_partition_t));
  new_partition->tips = partition->tips;
  new_partition->clv_buffers = partition->clv_buffers;
  new_partition->states = partition->states;
  new_partition->sites = sites;
  new_partition->rate_matrices = partition->rate_matrices;
  new_partition->prob_matrices = partition->prob_matrices;
  new_partition->rate_cats = partition->rate_cats;
  new_partition->scale_buffers = partition->scale_buffers;
  new_partition->attributes = partition->attributes;

    /* vectorization options */
  new_partition->alignment = partition->alignment;
  new_partition->states_padded = partition->states_padded;

  new_partition->mixture = partition->mixture;
  new_partition->clv = (double **)calloc(partition->tips + partition->clv_buffers,
                                            sizeof(double *));
  if (!new_partition->clv)
  {
    dealloc_partition_local (new_partition);
    return PLL_FAILURE;
  }
  for (i = 0; i < new_partition->tips + new_partition->clv_buffers; ++i)
    new_partition->clv[i] = partition->clv[i]
        + (start * partition->states_padded * partition->rate_cats);

  new_partition->scale_buffer = (unsigned int **) calloc (
      new_partition->scale_buffers, sizeof(unsigned int *));
  if (!new_partition->scale_buffer)
  {
    dealloc_partition_local (new_partition);
    return PLL_FAILURE;
  }
  for (i = 0; i < new_partition->scale_buffers; ++i)
    new_partition->scale_buffer[i] = &partition->scale_buffer[i][start];

  if (partition->invariant)
    new_partition->invariant = &partition->invariant[start];
  else
    new_partition->invariant = NULL;
  new_partition->pattern_weights = &partition->pattern_weights[start];

  new_partition->pmatrix = partition->pmatrix;
  new_partition->rates = partition->rates;
  new_partition->rate_weights = partition->rate_weights;
  new_partition->subst_params = partition->subst_params;
  new_partition->frequencies = partition->frequencies;
  new_partition->prop_invar = partition->prop_invar;

  new_partition->eigen_decomp_valid = partition->eigen_decomp_valid;
  new_partition->eigenvecs = partition->eigenvecs;
  new_partition->inv_eigenvecs = partition->inv_eigenvecs;
  new_partition->eigenvals = partition->eigenvals;

#if(0)
      pll_partition_create(partition->tips,
                                       partition->clv_buffers,
                                       partition->states,
                                       (unsigned int)sites,
                                       partition->mixture,
                                       partition->rate_matrices,
                                       partition->prob_matrices,
                                       partition->rate_cats,
                                       partition->scale_buffers,
                                       partition->attributes);
#endif

  return new_partition;
}

static unsigned long get_millis()
{
  struct timespec cur_time;
  clock_gettime(CLOCK_REALTIME, &cur_time);
  return cur_time.tv_sec*1000 + cur_time.tv_nsec/1e+6;
}

static void start_job_sync(int JOB)
{
  thread_job = JOB;
  while (thread_job != JOB_WAIT);
}

int main(int argc, char * argv[])
{
  unsigned int i;
  unsigned int tip_nodes_count, inner_nodes_count, nodes_count, branch_count;
  unsigned int matrix_count, ops_count;
  unsigned int * matrix_indices;
  unsigned int num_threads;
  double * branch_lengths;
  pll_partition_t * partition;
  pll_operation_t * operations;
  pll_utree_t ** travbuffer;
  unsigned long start_time, end_time, time_ms;

  if (argc != 4)
    fatal(" syntax: %s [newick] [fasta] [num_threads]", argv[0]);
  pll_utree_t * tree = pll_utree_parse_newick(argv[1], &tip_nodes_count);
  if (!tree)
    fatal("Tree must be an unrooted binary tree");
  set_missing_branch_length(tree, 0.000001);

  inner_nodes_count = tip_nodes_count - 2;
  nodes_count = inner_nodes_count + tip_nodes_count;
  branch_count = nodes_count - 1;

  printf("Number of tip/leaf nodes in tree: %d\n", tip_nodes_count);
  printf("Number of inner nodes in tree: %d\n", inner_nodes_count);
  printf("Total number of nodes in tree: %d\n", nodes_count);
  printf("Number of branches in tree: %d\n", branch_count);
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

  num_threads = atoi(argv[3]);
  if (atoi <= 0)
    fatal("Number of threads must be a positive integer");

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

  printf("Number of sites in MSA: %d\n", sites);
  partition = pll_partition_create(tip_nodes_count,
                                   inner_nodes_count,
                                   STATES,
                                   (unsigned int)sites,
                                   0,
                                   1,
                                   branch_count,
                                   RATE_CATS,
                                   inner_nodes_count,
                                   PLL_ATTRIB_ARCH_CPU);

  double frequencies[4] = { 0.17, 0.19, 0.25, 0.39 };
  double subst_params[6] = {1,1,1,1,1,1};
  double rate_cats[4] = {0};

  pll_compute_gamma_cats(1, 4, rate_cats);
  pll_set_frequencies(partition, 0, 0, frequencies);
  pll_set_subst_params(partition, 0, 0, subst_params);
  pll_set_category_rates(partition, rate_cats);

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
  matrix_indices = (unsigned int *)malloc(branch_count * sizeof(unsigned int));
  operations = (pll_operation_t *)malloc(inner_nodes_count *
                                                sizeof(pll_operation_t));

  /* perform a postorder traversal of the unrooted tree */
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

  pthread_t * threads = (pthread_t *) malloc(num_threads * sizeof(pthread_t));
  thread_data_t * thread_data = (thread_data_t *) malloc(num_threads * sizeof(thread_data_t));
  pll_partition_t ** partition_local = (pll_partition_t **) malloc(num_threads * sizeof(pll_partition_t *));
  double * result_buf = (double *) malloc(num_threads * sizeof(double));
  thread_job = JOB_WAIT;
  pthread_barrier_t barrier_buf;
  pthread_barrier_init (&barrier_buf,
                        NULL,
                        num_threads);
  for (i=0; i<num_threads; i++)
  {
    partition_local[i] = pll_partition_clone_partial(partition, i, num_threads);

    printf("Start thread %d\n", i);
    thread_data[i].thread_id = i;
    thread_data[i].num_threads = (long) num_threads;
    thread_data[i].partition = partition_local[i];
    thread_data[i].matrix_indices = matrix_indices;
    thread_data[i].matrix_count = matrix_count;
    thread_data[i].operations = operations;
    thread_data[i].ops_count = ops_count;
    thread_data[i].branch_lengths = branch_lengths;
    thread_data[i].vroot = tree;
    thread_data[i].barrier_buf = &barrier_buf;
    thread_data[i].trap = 1;
    thread_data[i].result_buf = result_buf;
    if(pthread_create(&(threads[i]), NULL, worker, &(thread_data[i]))) {
      fprintf(stderr, "Error creating thread\n");
      return 1;
    }
  }

  start_time = get_millis();
  unsigned int reps = 10;

  for (i=0; i<reps; i++)
    {
      unsigned long local_start = get_millis();
      start_job_sync(JOB_UPDATE_MATRICES);
      printf("[%3d] Done pmatrices %ld\n", i, get_millis() - local_start);
    }
    end_time = get_millis();
    time_ms = end_time - start_time;
    printf("Time: %ld\nAverage: %ld\n", time_ms, time_ms/reps);

  for (i=0; i<reps; i++)
  {
    unsigned long local_start = get_millis();
    start_job_sync(JOB_UPDATE_PARTIALS);
    printf("[%3d] Done partials %ld\n", i, get_millis() - local_start);
  }
  end_time = get_millis();
  time_ms = end_time - start_time;
  printf("Time: %ld\nAverage: %ld\n", time_ms, time_ms/reps);

  start_time = get_millis();
  for (i=0; i<reps; i++)
  {
    unsigned long local_start = get_millis();
    start_job_sync(JOB_FULL_LK);
    printf("[%3d] DONE job LK = %f %ld\n", i, global_lnl, get_millis() - local_start);
  }
  end_time = get_millis();
  time_ms = end_time - start_time;
  printf("Time: %ld\nAverage: %ld\n", time_ms, time_ms/reps);

  start_job_sync(JOB_FULL_LK);
  thread_job = JOB_FINALIZE;
  printf("[%3d] DONE job LK = %f\n", i, global_lnl);

  for (i=0; i<num_threads; i++)
    pthread_join (threads[i], NULL);

  /* clean */
  for (i=0; i<num_threads; i++)
    dealloc_partition_local(thread_data[i].partition);

  free(result_buf);
  free(threads);
  free(thread_data);
  free(partition_local);

  printf ("Traversal size: %d\n", traversal_size);
  printf ("Operations: %d\n", ops_count);
  printf ("Probability Matrices: %d\n", matrix_count);

  pll_update_prob_matrices(partition, 
                           0, 
                           matrix_indices, 
                           branch_lengths, 
                           matrix_count);
  pll_update_partials(partition, operations, ops_count);
  double logl = pll_compute_edge_loglikelihood(partition,
                                               tree->clv_index,
                                               tree->scaler_index,
                                               tree->back->clv_index,
                                               tree->back->scaler_index,
                                               tree->pmatrix_index,
                                               0);

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
  pll_utree_destroy(tree);

  return (EXIT_SUCCESS);
}
