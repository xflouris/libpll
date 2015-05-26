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

static void print_tree_recurse(pll_rtree_t * tree, 
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
  if (tree->left || tree->right) printf("+");

  printf (" %s:%f\n", tree->label, tree->length);

  if (active_node_order[indend_level-1] == 2) 
    active_node_order[indend_level-1] = 0;

  active_node_order[indend_level] = 1;
  print_tree_recurse(tree->left,
                     indend_level+1,
                     active_node_order);
  active_node_order[indend_level] = 2;
  print_tree_recurse(tree->right,
                     indend_level+1,
                     active_node_order);

}

static int tree_indend_level(pll_rtree_t * tree, int indend)
{
  if (!tree) return indend;

  int a = tree_indend_level(tree->left,  indend+1);
  int b = tree_indend_level(tree->right, indend+1);

  return (a > b ? a : b);
}

void pll_show_ascii_rtree(pll_rtree_t * tree)
{
  
  int indend_max = tree_indend_level(tree,0);

  int * active_node_order = (int *)malloc((indend_max+1) * sizeof(int));
  active_node_order[0] = 1;
  active_node_order[1] = 1;

  printf(" %s:%f\n", tree->label, tree->length);
  print_tree_recurse(tree->left,  1, active_node_order);
  print_tree_recurse(tree->right, 1, active_node_order);
  free(active_node_order);
}

PLL_EXPORT char * pll_write_newick_rtree(pll_rtree_t * root)
{
  char * newick;

  if (!root) return NULL;

  if (!(root->left) || !(root->right))
    asprintf(&newick, "%s:%f", root->label, root->length);
  else
  {
    char * subtree1 = pll_write_newick_rtree(root->left);
    char * subtree2 = pll_write_newick_rtree(root->right);

    asprintf(&newick, "(%s,%s)%s:%f", subtree1,
                                      subtree2,
                                      root->label ? root->label : "",
                                      root->length);
    free(subtree1);
    free(subtree2);
  }

  return newick;
}

static void traverse_rtree(pll_rtree_t * tree,
                           int tips,
                           double * branch_lengths,
                           int * indices,
                           int * index,
                           int * tip_count,
                           int * inner_count,
                           pll_operation_t * ops,
                           int * ops_index)
{
  if (!tree->left)
  {
    branch_lengths[*index] = tree->length;
    indices[*index] = *tip_count;

    *index = *index + 1;
    *tip_count = *tip_count + 1;
  }
  else
  {
    traverse_rtree(tree->left,
                   tips,
                   branch_lengths,
                   indices,
                   index,
                   tip_count,
                   inner_count,
                   ops,
                   ops_index);

    int left_index = indices[*index-1];

    traverse_rtree(tree->right,
                   tips,
                   branch_lengths,
                   indices,
                   index,
                   tip_count,
                   inner_count,
                   ops,
                   ops_index);         

    int right_index = indices[*index - 1];

    ops[*ops_index].parent_clv_index = *inner_count;

    ops[*ops_index].child1_clv_index = left_index;
    ops[*ops_index].child1_matrix_index = left_index;

    ops[*ops_index].child2_clv_index = right_index;
    ops[*ops_index].child2_matrix_index = right_index;

    ops[*ops_index].parent_scaler_index = *inner_count - tips;
    ops[*ops_index].child1_scaler_index = (left_index >= tips)
                                  ? left_index - tips : PLL_SCALE_BUFFER_NONE;
    ops[*ops_index].child2_scaler_index = (right_index >= tips)
                                  ? right_index - tips : PLL_SCALE_BUFFER_NONE;

    branch_lengths[*index] = tree->length;
    indices[*index] = *inner_count;
    *index = *index + 1;
    *inner_count = *inner_count + 1;
    *ops_index = *ops_index + 1;

  }
}

void pll_traverse_rtree(pll_rtree_t * tree,
                        int tips,
                        double ** branch_lengths,
                        int ** indices,
                        pll_operation_t ** ops,
                        int * root_clv_index,
                        int * root_scaler_index)
{
  int all_nodes = tips*2 - 1;

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

  traverse_rtree(tree->left,
                 tips,
                 *branch_lengths,
                 *indices, 
                 &index,
                 &tip_count,
                 &inner_count,
                 *ops,
                 &ops_index);

  int left_index = (*indices)[index - 1];

  traverse_rtree(tree->right,
                 tips,
                 *branch_lengths,
                 *indices, 
                 &index,
                 &tip_count,
                 &inner_count,
                 *ops,
                 &ops_index);

  int right_index = (*indices)[index - 1];

  /* set the last record of operations */
  (*ops)[ops_index].parent_clv_index    = inner_count;
  (*ops)[ops_index].parent_scaler_index = inner_count - tips;

  (*ops)[ops_index].child1_clv_index    = left_index;
  (*ops)[ops_index].child1_matrix_index = left_index;
  (*ops)[ops_index].child1_scaler_index = (left_index >= tips)
                                  ? left_index - tips : PLL_SCALE_BUFFER_NONE;

  (*ops)[ops_index].child2_clv_index    = right_index;
  (*ops)[ops_index].child2_matrix_index = right_index;
  (*ops)[ops_index].child2_scaler_index = (right_index >= tips)
                                  ? right_index - tips : PLL_SCALE_BUFFER_NONE;


  *root_clv_index = inner_count;
  *root_scaler_index = inner_count - tips;

}

static void query_utree_tipnames_recursive(pll_rtree_t * tree,
                                           char ** tipnames,
                                           int * index)
{
  if (!tree) return;

  if (!tree->left)
  {
    tipnames[*index] = tree->label;
    *index = *index + 1;
    return;
  }

  query_utree_tipnames_recursive(tree->left,  tipnames, index);
  query_utree_tipnames_recursive(tree->right, tipnames, index);
}

PLL_EXPORT char ** pll_query_rtree_tipnames(pll_rtree_t * tree,
                                            int tips)
{
  char ** tipnames = (char **)calloc(tips, sizeof(char *)); 
  int index = 0;

  query_utree_tipnames_recursive(tree->left,  tipnames, &index);
  query_utree_tipnames_recursive(tree->right, tipnames, &index);

  return tipnames;
}
