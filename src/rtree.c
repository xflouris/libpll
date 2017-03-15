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

static int indent_space = 4;

static void print_node_info(const pll_rtree_t * tree, int options)
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

static void print_tree_recurse(const pll_rtree_t * tree,
                               int indent_level,
                               int * active_node_order,
                               int options)
{
  int i,j;

  if (!tree) return;

  for (i = 0; i < indent_level; ++i)
  {
    if (active_node_order[i])
      printf("|");
    else
      printf(" ");

    for (j = 0; j < indent_space-1; ++j)
      printf(" ");
  }
  printf("\n");

  for (i = 0; i < indent_level-1; ++i)
  {
    if (active_node_order[i])
      printf("|");
    else
      printf(" ");

    for (j = 0; j < indent_space-1; ++j)
      printf(" ");
  }

  printf("+");
  for (j = 0; j < indent_space-1; ++j)
    printf ("-");
  if (tree->left || tree->right) printf("+");

  print_node_info(tree, options);

  if (active_node_order[indent_level-1] == 2)
    active_node_order[indent_level-1] = 0;

  active_node_order[indent_level] = 1;
  print_tree_recurse(tree->left,
                     indent_level+1,
                     active_node_order,
                     options);
  active_node_order[indent_level] = 2;
  print_tree_recurse(tree->right,
                     indent_level+1,
                     active_node_order,
                     options);

}

static unsigned int tree_indent_level(const pll_rtree_t * tree, unsigned int indent)
{
  if (!tree) return indent;

  unsigned int a = tree_indent_level(tree->left,  indent+1);
  unsigned int b = tree_indent_level(tree->right, indent+1);

  return (a > b ? a : b);
}

void pll_rtree_show_ascii(const pll_rtree_t * tree, int options)
{

  unsigned int indent_max = tree_indent_level(tree,0);

  int * active_node_order = (int *)malloc((indent_max+1) * sizeof(int));
  if (!active_node_order)
  {
    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf(pll_errmsg, 200, "Unable to allocate enough memory.");
    return;
  }
  active_node_order[0] = 1;
  active_node_order[1] = 1;

  print_node_info(tree, options);
  print_tree_recurse(tree->left,  1, active_node_order, options);
  print_tree_recurse(tree->right, 1, active_node_order, options);
  free(active_node_order);
}

static char * rtree_export_newick_recursive(const pll_rtree_t * root)
{
  char * newick;
  int size_alloced;
  assert(root != NULL);

  if (!(root->left) || !(root->right))
    size_alloced = asprintf(&newick, "%s:%f", root->label, root->length);
  else
  {
    char * subtree1 = rtree_export_newick_recursive(root->left);
    if (subtree1 == NULL)
    {
      return NULL;
    }
    char * subtree2 = rtree_export_newick_recursive(root->right);
    if (subtree2 == NULL)
    {
      free(subtree1);
      return NULL;
    }

    size_alloced = asprintf(&newick,
                            "(%s,%s)%s:%f",
                            subtree1,
                            subtree2,
                            root->label ? root->label : "",
                            root->length);
    free(subtree1);
    free(subtree2);
  }
  if (size_alloced < 0)
  {
    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf(pll_errmsg, 200, "memory allocation during newick export failed.");
    return NULL;
  }

  return newick;
}

PLL_EXPORT char * pll_rtree_export_newick(const pll_rtree_t * root)
{
  char * newick;
  int size_alloced;
  if (!root) return NULL;

  if (!(root->left) || !(root->right))
    size_alloced = asprintf(&newick, "%s:%f", root->label, root->length);
  else
  {
    char * subtree1 = rtree_export_newick_recursive(root->left);
    if (subtree1 == NULL)
    {
      pll_errno = PLL_ERROR_MEM_ALLOC;
      snprintf(pll_errmsg, 200, "Unable to allocate enough memory.");
      return NULL;
    }
    char * subtree2 = rtree_export_newick_recursive(root->right);
    if (subtree2 == NULL)
    {
      free(subtree1);
      pll_errno = PLL_ERROR_MEM_ALLOC;
      snprintf(pll_errmsg, 200, "Unable to allocate enough memory.");
      return NULL;
    }

    size_alloced = asprintf(&newick,
                            "(%s,%s)%s:%f;",
                            subtree1,
                            subtree2,
                            root->label ? root->label : "",
                            root->length);
    free(subtree1);
    free(subtree2);
  }
  if (size_alloced < 0)
  {
    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf(pll_errmsg, 200, "memory allocation during newick export failed");
    return NULL;
  }

  return newick;
}


PLL_EXPORT void pll_rtree_create_operations(pll_rtree_t ** trav_buffer,
                                            unsigned int trav_buffer_size,
                                            double * branches,
                                            unsigned int * pmatrix_indices,
                                            pll_operation_t * ops,
                                            unsigned int * matrix_count,
                                            unsigned int * ops_count)
{
  pll_rtree_t * node;
  unsigned int i;

  *ops_count = 0;
  *matrix_count = 0;

  for (i = 0; i < trav_buffer_size; ++i)
  {
    node = trav_buffer[i];

    /* do not store the branch of the root, since it does not exist */
    if (i < trav_buffer_size-1)
    {
      *branches++ = node->length;
      *pmatrix_indices++ = node->pmatrix_index;
      *matrix_count = *matrix_count + 1;
    }

    if (node->left)
    {
      ops[*ops_count].parent_clv_index = node->clv_index;
      ops[*ops_count].parent_scaler_index = node->scaler_index;

      ops[*ops_count].child1_clv_index = node->left->clv_index;
      ops[*ops_count].child1_scaler_index = node->left->scaler_index;
      ops[*ops_count].child1_matrix_index = node->left->pmatrix_index;

      ops[*ops_count].child2_clv_index = node->right->clv_index;
      ops[*ops_count].child2_scaler_index = node->right->scaler_index;
      ops[*ops_count].child2_matrix_index = node->right->pmatrix_index;

      *ops_count = *ops_count + 1;
    }
  }
}

static void rtree_traverse(pll_rtree_t * node,
                           int (*cbtrav)(pll_rtree_t *),
                           unsigned int * index,
                           pll_rtree_t ** outbuffer)
{
  if (!node->left)
  {
    if (cbtrav(node))
    {
      outbuffer[*index] = node;
      *index = *index + 1;
    }
    return;
  }
  if (!cbtrav(node))
    return;
  rtree_traverse(node->left, cbtrav, index, outbuffer);
  rtree_traverse(node->right, cbtrav, index, outbuffer);

  outbuffer[*index] = node;
  *index = *index + 1;
}

int pll_rtree_traverse(pll_rtree_t * root,
                       int (*cbtrav)(pll_rtree_t *),
                       pll_rtree_t ** outbuffer,
                       unsigned int * trav_size)
{
  *trav_size = 0;
  if (!root->left) return PLL_FAILURE;

  /* we will traverse an unrooted tree in the following way

           root
            /\
           /  \
        left   right

     at each node the callback function is called to decide whether we
     are going to traversing the subtree rooted at the specific node */

  rtree_traverse(root, cbtrav, trav_size, outbuffer);
  return PLL_SUCCESS;
}

static void rtree_traverse_preorder(pll_rtree_t * node,
                                    int (*cbtrav)(pll_rtree_t *),
                                    unsigned int * index,
                                    pll_rtree_t ** outbuffer)
{
  if (!node->left)
  {
    if (cbtrav(node))
    {
      outbuffer[*index] = node;
      *index = *index + 1;
    }
    return;
  }
  if (!cbtrav(node))
    return;

  outbuffer[*index] = node;
  *index = *index + 1;

  rtree_traverse_preorder(node->left, cbtrav, index, outbuffer);
  rtree_traverse_preorder(node->right, cbtrav, index, outbuffer);

}

PLL_EXPORT int pll_rtree_traverse_preorder(pll_rtree_t * root,
                                           int (*cbtrav)(pll_rtree_t *),
                                           pll_rtree_t ** outbuffer,
                                           unsigned int * trav_size)
{
  *trav_size = 0;
  if (!root->left) return PLL_FAILURE;

  /* we will traverse an unrooted tree in the following way

           root
            /\
           /  \
        left   right

     at each node the callback function is called to decide whether we
     are going to traversing the subtree rooted at the specific node */

  rtree_traverse_preorder(root, cbtrav, trav_size, outbuffer);
  return PLL_SUCCESS;
}

static void rtree_query_tipnodes_recursive(pll_rtree_t * node,
                                           pll_rtree_t ** node_list,
                                           unsigned int * index)
{
  if (!node) return;

  if (!node->left)
  {
    node_list[*index] = node;
    *index = *index + 1;
    return;
  }

  rtree_query_tipnodes_recursive(node->left,  node_list, index);
  rtree_query_tipnodes_recursive(node->right, node_list, index);
}

PLL_EXPORT unsigned int pll_rtree_query_tipnodes(pll_rtree_t * root,
                                                 pll_rtree_t ** node_list)
{
  unsigned int index = 0;

  if (!root) return 0;
  if (!root->left)
  {
    node_list[index++] = root;
    return index;
  }

  rtree_query_tipnodes_recursive(root->left,  node_list, &index);
  rtree_query_tipnodes_recursive(root->right, node_list, &index);

  return index;
}

static void rtree_query_innernodes_recursive(pll_rtree_t * root,
                                             pll_rtree_t ** node_list,
                                             unsigned int * index)
{
  if (!root) return;
  if (!root->left) return;

  /* postorder traversal */

  rtree_query_innernodes_recursive(root->left,  node_list, index);
  rtree_query_innernodes_recursive(root->right, node_list, index);

  node_list[*index] = root;
  *index = *index + 1;
  return;
}

PLL_EXPORT unsigned int pll_rtree_query_innernodes(pll_rtree_t * root,
                                        pll_rtree_t ** node_list)
{
  unsigned int index = 0;

  if (!root) return 0;
  if (!root->left) return 0;

  rtree_query_innernodes_recursive(root->left,  node_list, &index);
  rtree_query_innernodes_recursive(root->right, node_list, &index);

  node_list[index++] = root;

  return index;
}

PLL_EXPORT void pll_rtree_create_pars_buildops(pll_rtree_t ** trav_buffer,
                                               unsigned int trav_buffer_size,
                                               pll_pars_buildop_t * ops,
                                               unsigned int * ops_count)
{
  pll_rtree_t * node;
  unsigned int i;

  *ops_count = 0;

  for (i = 0; i < trav_buffer_size; ++i)
  {
    node = trav_buffer[i];

    if (node->left)
    {
      ops[*ops_count].parent_score_index = node->clv_index;
      ops[*ops_count].child1_score_index = node->left->clv_index;
      ops[*ops_count].child2_score_index = node->right->clv_index;

      *ops_count = *ops_count + 1;
    }
  }
}

PLL_EXPORT void pll_rtree_create_pars_recops(pll_rtree_t ** trav_buffer,
                                             unsigned int trav_buffer_size,
                                             pll_pars_recop_t * ops,
                                             unsigned int * ops_count)
{
  pll_rtree_t * node;
  unsigned int i;

  *ops_count = 0;

  for (i = 0; i < trav_buffer_size; ++i)
  {
    node = trav_buffer[i];

    if (node->left)
    {
      ops[*ops_count].node_score_index = node->clv_index;
      ops[*ops_count].node_ancestral_index = node->clv_index;

      if (node->parent)
      {
        ops[*ops_count].parent_score_index = node->parent->clv_index;
        ops[*ops_count].parent_ancestral_index = node->parent->clv_index;
      }
      else
      {
        /* invalid entries for the root - they will never be used */
        ops[*ops_count].parent_score_index = 0;
        ops[*ops_count].parent_ancestral_index = 0;
      }

      *ops_count = *ops_count + 1;
    }
  }
}

