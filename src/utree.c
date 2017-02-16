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

static void print_node_info(const pll_utree_t * tree, int options)
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

static char * xstrdup(const char * s)
{
  size_t len = strlen(s);
  char * p = (char *)malloc(len+1);
  if (!p)
  {
    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf(pll_errmsg, 200, "Memory allocation failed");
    return NULL;
  }
  return strcpy(p,s);
}

static void print_tree_recurse(pll_utree_t * tree,
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
  if (tree->next) printf("+");

  print_node_info(tree, options);

  if (active_node_order[indent_level-1] == 2)
    active_node_order[indent_level-1] = 0;

  if (tree->next)
  {
    active_node_order[indent_level] = 1;
    print_tree_recurse(tree->next->back,
                       indent_level+1,
                       active_node_order,
                       options);
    active_node_order[indent_level] = 2;
    print_tree_recurse(tree->next->next->back,
                       indent_level+1,
                       active_node_order,
                       options);
  }

}

static unsigned int tree_indent_level(const pll_utree_t * tree, unsigned int indent)
{
  if (!tree->next) return indent+1;

  unsigned int a = tree_indent_level(tree->next->back,       indent+1);
  unsigned int b = tree_indent_level(tree->next->next->back, indent+1);

  return (a > b ? a : b);
}

PLL_EXPORT void pll_utree_show_ascii(const pll_utree_t * tree, int options)
{
  unsigned int a, b;

  if (!tree->next) tree=tree->back;

  a = tree_indent_level(tree->back,1);
  b = tree_indent_level(tree,0);
  unsigned int max_indent_level = (a > b ? a : b);

  int * active_node_order = (int *)malloc((max_indent_level+1) * sizeof(int));
  if (!active_node_order)
  {
    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf(pll_errmsg, 200, "Unable to allocate enough memory.");
    return;
  }
  active_node_order[0] = 1;
  active_node_order[1] = 1;

  print_tree_recurse(tree->back,             1, active_node_order, options);
  print_tree_recurse(tree->next->back,       1, active_node_order, options);
  active_node_order[0] = 2;
  print_tree_recurse(tree->next->next->back, 1, active_node_order, options);
  free(active_node_order);
}

static char * newick_utree_recurse(const pll_utree_t * root)
{
  char * newick;
  int size_alloced;
  assert(root != NULL);
  if (!root->next)
    size_alloced = asprintf(&newick, "%s:%f", root->label, root->length);
  else
  {
    char * subtree1 = newick_utree_recurse(root->next->back);
    if (subtree1 == NULL)
    {
      pll_errno = PLL_ERROR_MEM_ALLOC;
      snprintf(pll_errmsg, 200, "Unable to allocate enough memory.");
      return NULL;
    }
    char * subtree2 = newick_utree_recurse(root->next->next->back);
    if (subtree2 == NULL)
    {
      free(subtree1);
      pll_errno = PLL_ERROR_MEM_ALLOC;
      snprintf(pll_errmsg, 200, "Unable to allocate enough memory.");
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
    snprintf(pll_errmsg, 200, "memory allocation during newick export failed");
    return NULL;
  }

  return newick;
}

PLL_EXPORT char * pll_utree_export_newick(const pll_utree_t * root)
{
  char * newick;
  int size_alloced;
  if (!root) return NULL;

  if (!root->next) root=root->back;

  char * subtree1 = newick_utree_recurse(root->back);
  if (subtree1 == NULL)
  {
    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf(pll_errmsg, 200, "Unable to allocate enough memory.");
    return NULL;
  }
  char * subtree2 = newick_utree_recurse(root->next->back);
  if (subtree2 == NULL)
  {
    free(subtree1);
    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf(pll_errmsg, 200, "Unable to allocate enough memory.");
    return NULL;
  }
  char * subtree3 = newick_utree_recurse(root->next->next->back);
  if (subtree3 == NULL)
  {
    free(subtree1);
    free(subtree2);
    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf(pll_errmsg, 200, "Unable to allocate enough memory.");
    return NULL;
  }

  size_alloced = asprintf(&newick,
                          "(%s,%s,%s)%s:0.0;",
                          subtree1,
                          subtree2,
                          subtree3,
                          root->label ? root->label : "");
  free(subtree1);
  free(subtree2);
  free(subtree3);
  if (size_alloced < 0)
  {
    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf(pll_errmsg, 200, "memory allocation during newick export failed");
    return NULL;
  }

  return (newick);
}

PLL_EXPORT void pll_utree_create_operations(pll_utree_t ** trav_buffer,
                                            unsigned int trav_buffer_size,
                                            double * branches,
                                            unsigned int * pmatrix_indices,
                                            pll_operation_t * ops,
                                            unsigned int * matrix_count,
                                            unsigned int * ops_count)
{
  pll_utree_t * node;
  unsigned int i;

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

static int utree_every_recursive(pll_utree_t * node,
                                 int (*cb)(pll_utree_t *))
{
  if (!node->next)
    return cb(node);

  if (!cb(node))
    return 0;

  return (utree_every_recursive(node->next->back,cb) &&
          utree_every_recursive(node->next->next->back,cb));
}

PLL_EXPORT int pll_utree_every(pll_utree_t * node,
                               int (*cb)(pll_utree_t *))
{
  return (utree_every_recursive(node,cb) &&
          utree_every_recursive(node->back,cb));
}

static void utree_traverse(pll_utree_t * node,
                           int (*cbtrav)(pll_utree_t *),
                           unsigned int * index,
                           pll_utree_t ** outbuffer)
{
  if (!node->next)
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
  utree_traverse(node->next->back, cbtrav, index, outbuffer);

  utree_traverse(node->next->next->back, cbtrav, index, outbuffer);

  outbuffer[*index] = node;
  *index = *index + 1;
}

PLL_EXPORT int pll_utree_traverse(pll_utree_t * root,
                                  int (*cbtrav)(pll_utree_t *),
                                  pll_utree_t ** outbuffer,
                                  unsigned int * trav_size)
{
  *trav_size = 0;
  if (!root->next) return PLL_FAILURE;

  /* we will traverse an unrooted tree in the following way

              2
            /
      1  --*
            \
              3

     at each node the callback function is called to decide whether we
     are going to traversing the subtree rooted at the specific node */

  utree_traverse(root->back, cbtrav, trav_size, outbuffer);
  utree_traverse(root, cbtrav, trav_size, outbuffer);

  return PLL_SUCCESS;
}

static void utree_query_tipnodes_recursive(pll_utree_t * node,
                                           pll_utree_t ** node_list,
                                           unsigned int * index)
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

PLL_EXPORT unsigned int pll_utree_query_tipnodes(pll_utree_t * root,
                                                 pll_utree_t ** node_list)
{
  unsigned int index = 0;

  if (!root) return 0;

  if (!root->next) root = root->back;

  utree_query_tipnodes_recursive(root->back, node_list, &index);

  utree_query_tipnodes_recursive(root->next->back, node_list, &index);
  utree_query_tipnodes_recursive(root->next->next->back, node_list, &index);

  return index;
}

static void utree_query_innernodes_recursive(pll_utree_t * node,
                                             pll_utree_t ** node_list,
                                             unsigned int * index)
{
  if (!node->next) return;

  /* postorder traversal */

  utree_query_innernodes_recursive(node->next->back, node_list, index);
  utree_query_innernodes_recursive(node->next->next->back, node_list, index);

  node_list[*index] = node;
  *index = *index + 1;
  return;
}

PLL_EXPORT unsigned int pll_utree_query_innernodes(pll_utree_t * root,
                                                   pll_utree_t ** node_list)
{
  unsigned int index = 0;

  if (!root) return 0;
  if (!root->next) root = root->back;

  utree_query_innernodes_recursive(root->back, node_list, &index);

  utree_query_innernodes_recursive(root->next->back, node_list, &index);
  utree_query_innernodes_recursive(root->next->next->back, node_list, &index);

  node_list[index++] = root;

  return index;
}

/* a callback function for checking tree integrity */
static int cb_check_integrity(pll_utree_t * node)
{
  unsigned int clv_index = node->clv_index;
  int scaler_index = node->scaler_index;
  unsigned int pmatrix_index = node->pmatrix_index;
  char * label = node->label;
  double length = node->length;

  /* edge attributes */
  if (node->back->length != length
      || node->back->pmatrix_index != pmatrix_index)
    return 0;
  if (node->next)
  {
    /* node attributes */
    if (node->next->clv_index != clv_index ||
        node->next->next->clv_index != clv_index)
      return 0;
    if (node->next->scaler_index != scaler_index ||
            node->next->next->scaler_index != scaler_index)
          return 0;
    if (node->next->label != label ||
            node->next->next->label != label)
          return 0;
  }
  return 1;
}

PLL_EXPORT int pll_utree_check_integrity(pll_utree_t * root)
{
  pll_utree_t * start_node = root->next?root:root->back;

  return pll_utree_every(start_node->back, cb_check_integrity);
}

/* TODO: Memory allocation checks were not implemented in this function!!! */
static pll_utree_t * clone_node(pll_utree_t * node)
{
  pll_utree_t * new = (pll_utree_t *)malloc(sizeof(pll_utree_t));
  memcpy(new, node, sizeof(pll_utree_t));

  if (node->label)
  {
    new->label = (char *)malloc(strlen(node->label)+1);
    strcpy(new->label,node->label);
  }

  if (node->next)
  {
    new->next = (pll_utree_t *)malloc(sizeof(pll_utree_t));
    memcpy(new->next, node->next, sizeof(pll_utree_t));

    new->next->next = (pll_utree_t *)malloc(sizeof(pll_utree_t));
    memcpy(new->next->next, node->next->next, sizeof(pll_utree_t));

    new->next->next->next = new;
    new->next->label = new->next->next->label = new->label;
  }
  return new;
}

static void utree_recurse_clone(pll_utree_t * new, pll_utree_t * root)
{
  if (root->back)
  {
    new->back = clone_node(root->back);
    new->back->back = new;

    if (root->back->next)
    {
      utree_recurse_clone(new->back->next,       root->back->next);
      utree_recurse_clone(new->back->next->next, root->back->next->next);
    }
  }
}

PLL_EXPORT pll_utree_t * pll_utree_clone(pll_utree_t * root)
{
  pll_utree_t * new = clone_node(root);

  utree_recurse_clone(new, root);
  if (root->next)
  {
    utree_recurse_clone(new->next, root->next);
    utree_recurse_clone(new->next->next, root->next->next);
  }

  return new;
}

static pll_utree_t * rtree_unroot(pll_rtree_t * root, pll_utree_t * back)
{
  pll_utree_t * uroot = (void *)calloc(1,sizeof(pll_utree_t));
  if (!uroot)
  {
    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf(pll_errmsg, 200, "Unable to allocate enough memory.");
    return NULL;
  }

  uroot->back = back;
  uroot->label = (root->label) ? xstrdup(root->label) : NULL;
  uroot->length = uroot->back->length;

  if (!root->left)
  {
    uroot->next = NULL;
    return uroot;
  }

  uroot->next = (void *)calloc(1,sizeof(pll_utree_t));
  if (!uroot->next)
  {
    free(uroot);
    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf(pll_errmsg, 200, "Unable to allocate enough memory.");
    return NULL;
  }

  uroot->next->next = (void *)calloc(1,sizeof(pll_utree_t));
  if (!uroot->next->next)
  {
    free(uroot->next);
    free(uroot);
    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf(pll_errmsg, 200, "Unable to allocate enough memory.");
    return NULL;
  }

  uroot->next->next->next = uroot;

  uroot->next->length = root->left->length;
  uroot->next->back = rtree_unroot(root->left, uroot->next);
  uroot->next->next->length = root->right->length;
  uroot->next->next->back = rtree_unroot(root->right, uroot->next->next);

  return uroot;
}

PLL_EXPORT pll_utree_t * pll_rtree_unroot(pll_rtree_t * root)
{
  if (!root->left->left && !root->right->left)
  {
    pll_errno = PLL_ERROR_TREE_CONVERSION;
    snprintf(pll_errmsg,
             200,
             "Tree requires at least three tips to be converted to unrooted");
    return NULL;
  }

  pll_rtree_t * new_root;

  pll_utree_t * uroot = (void *)calloc(1,sizeof(pll_utree_t));
  if (!uroot)
  {
    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf(pll_errmsg, 200, "Unable to allocate enough memory.");
    return NULL;
  }

  uroot->next = (void *)calloc(1,sizeof(pll_utree_t));
  if (!uroot->next)
  {
    free(uroot);
    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf(pll_errmsg, 200, "Unable to allocate enough memory.");
    return NULL;
  }

  uroot->next->next = (void *)calloc(1,sizeof(pll_utree_t));
  if (!uroot->next->next)
  {
    free(uroot->next);
    free(uroot);
    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf(pll_errmsg, 200, "Unable to allocate enough memory.");
    return NULL;
  }

  uroot->next->next->next = uroot;
  uroot->length = root->left->length + root->right->length;

  /* get the first root child that has descendants and make  it the new root */
  if (root->left->left)
  {
    new_root = root->left;
    uroot->back = rtree_unroot(root->right,uroot);
    /* TODO: Need to clean uroot in case of error */
    if (!uroot->back) return NULL;
  }
  else
  {
    new_root = root->right;
    uroot->back = rtree_unroot(root->left,uroot);
    /* TODO: Need to clean uroot in case of error*/
    if (!uroot->back) return NULL;
  }

  uroot->label = (new_root->label) ? xstrdup(new_root->label) : NULL;

  uroot->next->label = uroot->label;
  uroot->next->length = new_root->left->length;
  uroot->next->back = rtree_unroot(new_root->left, uroot->next);
  /* TODO: Need to clean uroot in case of error*/
  if (!uroot->next->back) return NULL;

  uroot->next->next->label = uroot->label;
  uroot->next->next->length = new_root->right->length;
  uroot->next->next->back = rtree_unroot(new_root->right, uroot->next->next);
  /* TODO: Need to clean uroot in case of error*/
  if (!uroot->next->next->back) return NULL;

  return uroot;
}

PLL_EXPORT void pll_utree_create_pars_buildops(pll_utree_t ** trav_buffer,
                                               unsigned int trav_buffer_size,
                                               pll_pars_buildop_t * ops,
                                               unsigned int * ops_count)
{
  pll_utree_t * node;
  unsigned int i;

  *ops_count = 0;

  for (i = 0; i < trav_buffer_size; ++i)
  {
    node = trav_buffer[i];

    if (node->next)
    {
      ops[*ops_count].parent_score_index = node->clv_index;
      ops[*ops_count].child1_score_index = node->next->back->clv_index;
      ops[*ops_count].child2_score_index = node->next->next->back->clv_index;

      *ops_count = *ops_count + 1;
    }
  }
}
