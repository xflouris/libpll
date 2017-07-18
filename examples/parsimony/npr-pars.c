/*
    Copyright (C) 2015 Tomas Flouri

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
#include <time.h>

#define STATES 20
#define MAP pll_map_aa

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
  unsigned int i,j,n;
  unsigned int tip_nodes_count, inner_nodes_count, nodes_count, branch_count;
  unsigned int ops_count;
  pll_parsimony_t * pars;
  pll_pars_buildop_t * operations;
  pll_pars_recop_t * recops;
  pll_rnode_t ** travbuffer;

  /* we accept only two arguments - the newick tree (unrooted binary) and the
     alignment in the form of PHYLIP reads */
  if (argc != 3)
    fatal(" syntax: %s [newick] [phylip]", argv[0]);

  /* parse the unrooted binary tree in newick format, and store the number
     of tip nodes in tip_nodes_count */
  pll_rtree_t * tree = pll_rtree_parse_newick(argv[1]);
  if (!tree)
    fatal("Tree must be a rooted binary tree");

  /* compute and show node count information */
  tip_nodes_count = tree->tip_count;
  inner_nodes_count = tip_nodes_count - 1;
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

  /* read PHYLIP alignment */
  pll_phylip_t * fd = pll_phylip_open(argv[2], pll_map_phylip);
  if (!fd)
    fatal(pll_errmsg);

  pll_msa_t * msa = pll_phylip_parse_interleaved(fd);
  if (!msa)
    fatal(pll_errmsg);

  pll_phylip_close(fd);

  if ((unsigned int)(msa->count) != tip_nodes_count)
    fatal("Number of sequences does not match number of leaves in tree");

  /* create the PLL parsimony instance

     parameters:

     tips : number of tip sequences we want to have
     states : number of states that our data have
     sites : number of sites
     score_matrix : score matrix to use
     score_buffers: number of score buffers to allocate (typically one per
                    inner node)
     ancestral_buffers: number of ancestral state buffers to allocate
                        (typically one per inner node)
  */

  /* set a matrix where mutations are penalized equally by 1 */
  double score_matrix[STATES*STATES];
  for (i = 0; i < STATES*STATES; ++i)
    score_matrix[i] = 1;
  for (i = 0; i < STATES; ++i)
    score_matrix[i*STATES+i] = 0;

  pars = pll_parsimony_create(tip_nodes_count,
                              STATES,
                              (unsigned int)(msa->length),
                              score_matrix,
                              inner_nodes_count,
                              inner_nodes_count);

  /* find sequences in hash table and link them with the corresponding taxa */
  for (i = 0; i < tip_nodes_count; ++i)
  {
    ENTRY query;
    query.key = msa->label[i];
    ENTRY * found = NULL;

    found = hsearch(query,FIND);

    if (!found)
      fatal("Sequence with header %s does not appear in the tree", msa->label[i]);

    unsigned int tip_clv_index = *((unsigned int *)(found->data));

    printf("Setting sequence: %s\n", msa->sequence[i]);
    pll_set_parsimony_sequence(pars, tip_clv_index, MAP, msa->sequence[i]);
  }

  pll_msa_destroy(msa);

  /* destroy hash table */
  hdestroy();

  /* we no longer need these two arrays (keys and values of hash table... */
  free(data);

  /* allocate a buffer for storing pointers to nodes of the tree in postorder
     traversal */
  travbuffer = (pll_rnode_t **)xmalloc(nodes_count * sizeof(pll_rnode_t *));


  operations = (pll_pars_buildop_t *)xmalloc(inner_nodes_count *
                                             sizeof(pll_pars_buildop_t));

  /* perform a postorder traversal of the rooted tree */
  unsigned int traversal_size;
  if (!pll_rtree_traverse(tree->root,
                          PLL_TREE_TRAVERSE_POSTORDER,
                          cb_full_traversal,
                          travbuffer,
                          &traversal_size))
    fatal("Function pll_rtree_traverse() requires inner nodes as parameters");

  /* given the computed traversal descriptor, generate the build operations
     structure */
  pll_rtree_create_pars_buildops(travbuffer,
                                 traversal_size,
                                 operations,
                                 &ops_count);


  printf ("Traversal size: %d\n", traversal_size);
  printf ("Operations: %d\n", ops_count);

  double minscore = pll_parsimony_build(pars,
                                        operations,
                                        ops_count);
  printf("Minimum parsimony score: %f\n", minscore);


  /* print score buffer for each inner node */
  for (i = 0; i < traversal_size; ++i)
  {
    unsigned int id = travbuffer[i]->clv_index;
    if (id >= tip_nodes_count)
    {
      printf ("%s : ",travbuffer[i]->label);
      for (n = 0; n < pars->sites; ++n)
      {
        for (j = 0; j < STATES; ++j)
          printf("%.0f ", pars->sbuffer[id][n*pars->states+j]);
        printf("+ ");
       }
      printf("\n");
    }
  }

  /* perform a preorder traversal of the rooted tree */
  if (!pll_rtree_traverse(tree->root,
                          PLL_TREE_TRAVERSE_PREORDER,
                          cb_full_traversal,
                          travbuffer,
                          &traversal_size))
    fatal("Function pll_rtree_traverse() requires inner nodes as parameters");

  /* create the reconstruct operations */
  recops = (pll_pars_recop_t *)xmalloc(inner_nodes_count *
                                       sizeof(pll_pars_recop_t));
  pll_rtree_create_pars_recops(travbuffer,
                               traversal_size,
                               recops,
                               &ops_count);


  /* reconstruct ancestral sequences */
  pll_parsimony_reconstruct(pars, MAP, recops, ops_count);

  printf("\n\nReconstruction:\n");

  /* print reconstructed sequences at inner nodes */
  for (i = 0; i < traversal_size; ++i)
  {
    unsigned int id = travbuffer[i]->clv_index;
    if (id >= tip_nodes_count)
    {
      printf ("%s : ",travbuffer[i]->label);
      for (j = 0; j < pars->sites; ++j)
        printf("%c", pars->anc_states[id][j]);
      printf("\n");
    }
  }


  /* free traversal, build operations and reconstruct operations */
  free(travbuffer);
  free(operations);
  free(recops);

  /* we will no longer need the tree structure */
  pll_rtree_destroy(tree,NULL);
  
  /* destroy the parsimony structure */
  pll_parsimony_destroy(pars);

  return (EXIT_SUCCESS);
}
