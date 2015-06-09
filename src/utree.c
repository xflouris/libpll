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

static void print_node_info(pll_utree_t * tree, int options)
{
  if (options & PLL_UTREE_SHOW_LABEL)
    printf (" %s", tree->label);
  if (options & PLL_UTREE_SHOW_BRANCH_LENGTH)
    printf (" %f", tree->length);
  if (options & PLL_UTREE_SHOW_CLV_INDEX)
    printf (" %d", tree->clv_index);
  if (options & PLL_UTREE_SHOW_SCALER_INDEX)
    printf (" %d", tree->scaler_index);
  if (options & PLL_UTREE_SHOW_PMATRIX_INDEX)
    printf (" %d", tree->pmatrix_index);
  printf("\n");
}

static void print_tree_recurse(pll_utree_t * tree, 
                               int indend_level, 
                               int * active_node_order,
                               int options)
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

  print_node_info(tree, options);

  if (active_node_order[indend_level-1] == 2) 
    active_node_order[indend_level-1] = 0;

  if (tree->next)
  {
    active_node_order[indend_level] = 1;
    print_tree_recurse(tree->next->back,
                       indend_level+1,
                       active_node_order,
                       options);
    active_node_order[indend_level] = 2;
    print_tree_recurse(tree->next->next->back, 
                       indend_level+1,
                       active_node_order,
                       options);
  }

}

static int tree_indend_level(pll_utree_t * tree, int indend)
{
  if (!tree->next) return indend+1;

  int a = tree_indend_level(tree->next->back,       indend+1);
  int b = tree_indend_level(tree->next->next->back, indend+1);

  return (a > b ? a : b);
}

void pll_utree_show_ascii(pll_utree_t * tree, int options)
{
  int a, b;
  
  a = tree_indend_level(tree->back,1);
  b = tree_indend_level(tree,0);
  int max_indend_level = (a > b ? a : b);


  int * active_node_order = (int *)malloc((max_indend_level+1) * sizeof(int));
  active_node_order[0] = 1;
  active_node_order[1] = 1;

  print_tree_recurse(tree->back,             1, active_node_order, options);
  print_tree_recurse(tree->next->back,       1, active_node_order, options);
  active_node_order[0] = 2;
  print_tree_recurse(tree->next->next->back, 1, active_node_order, options);
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

PLL_EXPORT char * pll_utree_export_newick(pll_utree_t * root)
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

PLL_EXPORT void pll_utree_create_operations(pll_utree_t ** trav_buffer,
                                            int trav_buffer_size,
                                            double * branches,
                                            int * pmatrix_indices,
                                            pll_operation_t * ops,
                                            int * matrix_count,
                                            int * ops_count)
{
  pll_utree_t * node;
  int i;

  *ops_count = 0;
  *matrix_count = 0;

  for (i = 0; i < trav_buffer_size; ++i)
  {
    node = trav_buffer[i];

    /* if the current node is the second end-point of the edge
    shared with the root node, then do not add the edge to the
    list as it will be added in the end (avoid duplicate edges
    in the list) */
    if (node != trav_buffer[trav_buffer_size - 1]->back)
    {
      *branches++ = node->length;
      *pmatrix_indices++ = node->pmatrix_index;
      *matrix_count = *matrix_count + 1;
    }

    if (node->next)
    {
      ops[*ops_count].parent_clv_index = node->clv_index;
      ops[*ops_count].parent_scaler_index = node->scaler_index;

      ops[*ops_count].child1_clv_index = node->next->back->clv_index;
      ops[*ops_count].child1_scaler_index = node->next->back->scaler_index;
      ops[*ops_count].child1_matrix_index = node->next->back->pmatrix_index;

      ops[*ops_count].child2_clv_index = node->next->next->back->clv_index;
      ops[*ops_count].child2_scaler_index = node->next->next->back->scaler_index;
      ops[*ops_count].child2_matrix_index = node->next->next->back->pmatrix_index;

      *ops_count = *ops_count + 1;
    }
  }
}

static void utree_traverse(pll_utree_t * node,
                           pll_utree_t * prev,
                           int (*cbtrav)(pll_utree_t *, pll_utree_t *),
                           int * index,
                           pll_utree_t ** outbuffer)
{
  if (!node->next)
  {
    if (cbtrav(node, prev))
    {
      outbuffer[*index] = node;
      *index = *index + 1;
    }
    return;
  }
  if (!cbtrav(node, prev))
    return;
  utree_traverse(node->next->back, node, cbtrav, index, outbuffer);

  utree_traverse(node->next->next->back, node, cbtrav, index, outbuffer);

  outbuffer[*index] = node;
  *index = *index + 1;
}

PLL_EXPORT int pll_utree_traverse(pll_utree_t * root,
                                  int (*cbtrav)(pll_utree_t *, pll_utree_t *),
                                  pll_utree_t ** outbuffer)
{
  int index = 0;

  if (!root->next) return -1;

  /* we will traverse an unrooted tree in the following way
      
              2
            /
      1  --*
            \
              3

     at each node the callback function is called to decide whether we
     are going to traversing the subtree rooted at the specific node */

  utree_traverse(root->back, root, cbtrav, &index, outbuffer);
  utree_traverse(root, NULL, cbtrav, &index, outbuffer);

  return index;
}


static void utree_query_tipnodes_recursive(pll_utree_t * node,
                                           pll_utree_t ** node_list,
                                           int * index)
{
  if (!node->next)
  {
    node_list[*index] = node;
    *index = *index + 1;
    return;
  }

  utree_query_tipnodes_recursive(node->next->back, node_list, index);
  utree_query_tipnodes_recursive(node->next->next->back, node_list, index);
}

PLL_EXPORT int pll_utree_query_tipnodes(pll_utree_t * root,
                                        pll_utree_t ** node_list)
{
  int index = 0;

  if (!root) return 0;

  if (!root->next) root = root->back;

  utree_query_tipnodes_recursive(root->back, node_list, &index);

  utree_query_tipnodes_recursive(root->next->back, node_list, &index);
  utree_query_tipnodes_recursive(root->next->next->back, node_list, &index);

  return index;
}

static void utree_query_innernodes_recursive(pll_utree_t * node,
                                           pll_utree_t ** node_list,
                                           int * index)
{
  if (!node->next) return;

  /* postorder traversal */

  utree_query_innernodes_recursive(node->next->back, node_list, index);
  utree_query_innernodes_recursive(node->next->next->back, node_list, index);

  node_list[*index] = node;
  *index = *index + 1;
  return;
}

PLL_EXPORT int pll_utree_query_innernodes(pll_utree_t * root,
                                          pll_utree_t ** node_list)
{
  int index = 0;

  if (!root) return 0;
  if (!root->next) root = root->back;

  utree_query_innernodes_recursive(root->back, node_list, &index);

  utree_query_innernodes_recursive(root->next->back, node_list, &index);
  utree_query_innernodes_recursive(root->next->next->back, node_list, &index);

  node_list[index++] = root;

  return index;
}

