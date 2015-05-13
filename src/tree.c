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
    Heidelberg Institute for Theoretical Studies,
    Schloss-Wolfsbrunnenweg 35, D-69118 Heidelberg, Germany
*/

#include "pll.h"

static int indend_space = 4;

static void print_tree_recurse(pll_utree_t * tree, 
                               int indend_level, 
                               int * active_node_order)
{
  int i,j;

  if (!tree) return;

  for (i = 0; i < indend_level; ++i)
  {
    if (active_node_order[i])
      printf("|");
    else
      printf(" ");

    for (j = 0; j < indend_space-1; ++j)
      printf(" ");
  }
  printf("\n");

  for (i = 0; i < indend_level-1; ++i)
  {
    if (active_node_order[i])
      printf("|");
    else
      printf(" ");

    for (j = 0; j < indend_space-1; ++j)
      printf(" ");
  }

  printf("+");
  for (j = 0; j < indend_space-1; ++j)
    printf ("-");
  if (tree->next) printf("+");

  printf (" %s:%f\n", tree->label, tree->length);

  if (active_node_order[indend_level-1] == 2) 
    active_node_order[indend_level-1] = 0;

  if (tree->next)
  {
    active_node_order[indend_level] = 1;
    print_tree_recurse(tree->next->back,
                       indend_level+1,
                       active_node_order);
    active_node_order[indend_level] = 2;
    print_tree_recurse(tree->next->next->back, 
                       indend_level+1,
                       active_node_order);
  }

}

static int tree_indend_level(pll_utree_t * tree, int indend)
{
  if (!tree->next) return indend+1;

  int a = tree_indend_level(tree->next->back,       indend+1);
  int b = tree_indend_level(tree->next->next->back, indend+1);

  return MAX(a,b);
}

void pll_show_ascii_utree(pll_utree_t * tree)
{
  int a, b;
  
  a = tree_indend_level(tree->back,1);
  b = tree_indend_level(tree,0);
  int max_indend_level = MAX(a,b);


  int * active_node_order = (int *)malloc((max_indend_level+1) * sizeof(int));
  active_node_order[0] = 1;
  active_node_order[1] = 1;

  print_tree_recurse(tree->back,             1, active_node_order);
  print_tree_recurse(tree->next->back,       1, active_node_order);
  active_node_order[0] = 2;
  print_tree_recurse(tree->next->next->back, 1, active_node_order);
  free(active_node_order);
}

static char * newick_utree_recurse(pll_utree_t * root)
{
  char * newick;

  if (!root->next)
    asprintf(&newick, "%s:%f", root->label, root->length);
  else
  {
    char * subtree1 = newick_utree_recurse(root->next->back);
    char * subtree2 = newick_utree_recurse(root->next->next->back);

    asprintf(&newick, "(%s,%s)%s:%f", subtree1,
                                      subtree2,
                                      root->label ? root->label : "",
                                      root->length);
    free(subtree1);
    free(subtree2);
  }

  return newick;
}

PLL_EXPORT char * pll_write_newick_utree(pll_utree_t * root)
{
  char * newick;

  if (!root) return NULL;

  char * subtree1 = newick_utree_recurse(root->back);
  char * subtree2 = newick_utree_recurse(root->next->back);
  char * subtree3 = newick_utree_recurse(root->next->next->back);

  asprintf(&newick, "(%s,%s,%s)%s:%f;", subtree1,
                                        subtree2,
                                        subtree3,
                                        root->label ? root->label : "",
                                        root->length);
  free(subtree1);
  free(subtree2);
  free(subtree3);

  return (newick);

}

static void traverse_utree(pll_utree_t * tree, 
                    double * branch_lengths, 
                    int * indices,
                    int * index,
                    int * tip_count,
                    int * inner_count,
                    pll_operation_t * ops,
                    int * ops_index)
{
  /* is it a tip? */
  if (!tree->next)      /* tip */
  {
    branch_lengths[*index] = tree->length;
    indices[*index] = *tip_count;

    *index = *index + 1;
    *tip_count = *tip_count + 1;
  }
  else  /* inner */
  {
    traverse_utree(tree->next->back, 
                   branch_lengths, 
                   indices, 
                   index, 
                   tip_count, 
                   inner_count, 
                   ops, 
                   ops_index);

    int child1_index = indices[*index - 1];

    traverse_utree(tree->next->next->back, 
                   branch_lengths, 
                   indices, 
                   index, 
                   tip_count, 
                   inner_count, 
                   ops, 
                   ops_index);

    int child2_index = indices[*index - 1];

    ops[*ops_index].parent_clv_index = *inner_count;

    ops[*ops_index].child1_clv_index = child1_index;
    ops[*ops_index].child1_matrix_index = child1_index;

    ops[*ops_index].child2_clv_index = child2_index;
    ops[*ops_index].child2_matrix_index = child2_index;

    branch_lengths[*index] = tree->length;
    indices[*index] = *inner_count;
    *index = *index + 1;
    *inner_count = *inner_count + 1;
    *ops_index = *ops_index + 1;
  }
}

void pll_traverse_utree(pll_utree_t * tree, 
                        int tips, 
                        double ** branch_lengths, 
                        int ** indices,
                        pll_operation_t ** ops,
                        int * edge_pmatrix_index,
                        int * edge_node1_clv_index,
                        int * edge_node2_clv_index)
{
  int all_nodes = tips*2 - 2;

  int tip_count = 0;
  int inner_count = tips;
  int ops_index = 0;
  int index = 0;
  
  if (!*branch_lengths)
    *branch_lengths = (double *)calloc(all_nodes-1, sizeof(double));
  if (!*indices)
    *indices = (int *)calloc(all_nodes-1, sizeof(int));
  if (!*ops)
    *ops = (pll_operation_t *)calloc(all_nodes - tips, sizeof(pll_operation_t));


  /* we will traverse an unrooted tree in the following way
      
              2
            /
      1  --*
            \
              3

     the last operation in the ops structure will drive the computation of the
     CLV of the parent of 2 and 3 (let's call this node as *), edge_pmatrix_index
     will point to edge 1-*, edge_node1_clv_index will be the CLV index of 1, and
     edge_node2_clv_index the CLV index of *. */


  /* traverse first subtree */

  traverse_utree(tree->back, 
                 *branch_lengths, 
                 *indices, &index, 
                 &tip_count, 
                 &inner_count, 
                 *ops, 
                 &ops_index);

  *edge_node1_clv_index = (*indices)[index-1];
  *edge_pmatrix_index = (*indices)[index-1];

  /* traverse second subtree */
  traverse_utree(tree->next->back, 
                 *branch_lengths, 
                 *indices, 
                 &index, 
                 &tip_count, 
                 &inner_count, 
                 *ops, 
                 &ops_index);
 
  int child1_index = (*indices)[index - 1];
  
  /* traverse third subtree */
  traverse_utree(tree->next->next->back, 
                 *branch_lengths, 
                 *indices, &index, 
                 &tip_count, 
                 &inner_count, 
                 *ops, 
                 &ops_index);

  int child2_index = (*indices)[index - 1];

  (*ops)[ops_index].parent_clv_index = inner_count;

  (*ops)[ops_index].child1_clv_index = child1_index;
  (*ops)[ops_index].child1_matrix_index = child1_index;

  (*ops)[ops_index].child2_clv_index = child2_index;
  (*ops)[ops_index].child2_matrix_index = child2_index;

  *edge_node2_clv_index = inner_count;
}

static void query_utree_tipnames_recursive(pll_utree_t * tree,
                                           char ** tipnames,
                                           int * index)
{
  if (!tree->next)
  {
    tipnames[*index] = tree->label;
    *index = *index + 1;
    return;
  }

  query_utree_tipnames_recursive(tree->next->back, tipnames, index);
  query_utree_tipnames_recursive(tree->next->next->back, tipnames, index);
}

PLL_EXPORT char ** pll_query_utree_tipnames(pll_utree_t * tree,
                                                    int tips)
{
  char ** tipnames = (char **)calloc(tips, sizeof(char *)); 
  int index = 0;

  query_utree_tipnames_recursive(tree->back, tipnames, &index);

  query_utree_tipnames_recursive(tree->next->back, tipnames, &index);
  query_utree_tipnames_recursive(tree->next->next->back, tipnames, &index);

  return tipnames;
}
