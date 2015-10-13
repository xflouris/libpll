#include "pll_tree.h"

/**
 * @brief Creates a new circular node
 *
 *           n2
 *          / |
 *        n1  |
 *          \ |
 *           n3
 *
 * All parameters are shared among the nodes in the triplet
 *
 * @param clv_index    the clv_index
 * @param scaler_index the scaler index
 * @param label        the node label
 * @param data         the data pointer
 *
 * @returns the new node
 */
static pll_utree_t * utree_create_node(unsigned int clv_index,
                                       unsigned int scaler_index,
                                       char * label,
                                       void * data)
{
  pll_utree_t * new_node = (pll_utree_t *)calloc(1, sizeof(pll_utree_t));
  new_node->next         = (pll_utree_t *)calloc(1, sizeof(pll_utree_t));
  new_node->next->next   = (pll_utree_t *)calloc(1, sizeof(pll_utree_t));
  new_node->next->next->next = new_node;
  new_node->label = label;
  new_node->next->label = new_node->next->next->label = new_node->label;
  new_node->next->data = new_node->next->next->data = new_node->data = data;
  new_node->next->length = new_node->next->next->length = new_node->length = 0;
  new_node->next->clv_index = new_node->next->next->clv_index = new_node->clv_index = clv_index;
  new_node->next->scaler_index = new_node->next->next->scaler_index = new_node->scaler_index = scaler_index;
  new_node->back = new_node->next->back = new_node->next->next->back = NULL;
  return new_node;
}

/**
 * removes a circular node
 */
static int utree_destroy_node(pll_utree_t * node)
{
  if (node->next)
  {
    free(node->next->next);
    free(node->next);
  }
  if (node->label)
    free(node->label);
  free(node);
  return PLL_SUCCESS;
}

/**
 * connects 2 nodes
 */
static int utree_connect_nodes(pll_utree_t * parent,
                               pll_utree_t * child,
                               double length)
{
//  if (parent->back || child->back)
//    return PLL_FAILURE;
  parent->back = child;
  child->back = parent;
  parent->length = child->length = length;

  /* PMatrix index is set to parent node */
  child->pmatrix_index = parent->pmatrix_index;

  return PLL_SUCCESS;
}

PLL_EXPORT int pll_utree_bisect(pll_utree_t * edge,
                                pll_edge_t * parent_subtree,
                                pll_edge_t * child_subtree)
{
  assert(parent_subtree);
  assert(child_subtree);

  if (!edge->next)
    return PLL_FAILURE;

  pll_utree_t * c_edge = edge->back;

  /* connect parent subtree */
  parent_subtree->edge.utree.parent = edge->next->back;
  parent_subtree->edge.utree.child = edge->next->next->back;
  parent_subtree->length =
      parent_subtree->edge.utree.parent->length +
      parent_subtree->edge.utree.child->length;

  /* save removed pmatrix index */
  parent_subtree->additional_pmatrix_index = parent_subtree->edge.utree.child->pmatrix_index;

  utree_connect_nodes(parent_subtree->edge.utree.parent,
                      parent_subtree->edge.utree.child,
                      parent_subtree->length);

  /* connect child subtree */
  child_subtree->edge.utree.parent = c_edge->next->back;
  child_subtree->edge.utree.child = c_edge->next->next->back;
  child_subtree->length =
      child_subtree->edge.utree.parent->length +
      child_subtree->edge.utree.child->length;

  /* save removed pmatrix index */
  child_subtree->additional_pmatrix_index = child_subtree->edge.utree.child->pmatrix_index;

  utree_connect_nodes(child_subtree->edge.utree.parent,
                      child_subtree->edge.utree.child,
                      child_subtree->length);

  /* remove edges */
  utree_destroy_node(edge);
  utree_destroy_node(c_edge);
  return PLL_SUCCESS;
}

PLL_EXPORT pll_edge_t pll_utree_reconnect(pll_edge_t * edge,
                                          unsigned int parent_pmatrix_index,
                                          unsigned int parent_clv_index,
                                          unsigned int parent_scaler_index,
                                          unsigned int child_pmatrix_index,
                                          unsigned int child_clv_index,
                                          unsigned int child_scaler_index,
                                          unsigned int edge_pmatrix_index)
{
  /* create and connect 2 new nodes */
  char * parent_label = (char *) malloc(20);
  char * child_label = (char *) malloc(20);
  snprintf(parent_label, 6, "NEW_P");
  snprintf(child_label, 6, "NEW_C");
  pll_utree_t * parent_node = utree_create_node(parent_clv_index,
                                                parent_scaler_index,
                                                parent_label,
                                                NULL);
  pll_utree_t * child_node = utree_create_node(child_clv_index,
                                               child_scaler_index,
                                               child_label,
                                               NULL);
  pll_edge_t new_edge;
  new_edge.edge.utree.parent = parent_node;
  new_edge.edge.utree.child = child_node;
  new_edge.length = edge->length;
  new_edge.additional_pmatrix_index = -1;

  utree_connect_nodes(parent_node, child_node, edge->length);
  parent_node->pmatrix_index = child_node->pmatrix_index = edge_pmatrix_index;

  /* reconnect parent close to edge.parent */
  utree_connect_nodes(edge->edge.utree.parent->back,
                      parent_node->next->next,
                      edge->edge.utree.parent->back->length);
  parent_node->next->pmatrix_index = parent_pmatrix_index;
  utree_connect_nodes(parent_node->next,
                      edge->edge.utree.parent,
                      0);

  /* reconnect child close to edge.child */
  utree_connect_nodes(edge->edge.utree.child->back,
                      child_node->next->next,
                      edge->edge.utree.child->back->length);
  child_node->next->pmatrix_index = child_pmatrix_index;
  utree_connect_nodes(child_node->next,
                      edge->edge.utree.child,
                      0);

  return new_edge;
}

PLL_EXPORT static int utree_find_node_in_subtree(pll_utree_t * root,
                                                 pll_utree_t * node)
{
  if (root == node)
  {
    return PLL_SUCCESS;
  }

  if (root->next)
  {
    if (root->next == node || root->next->next == node)
    {
      return PLL_SUCCESS;
    }

    return utree_find_node_in_subtree(root->next->back, node)
        || utree_find_node_in_subtree(root->next->next->back, node);
  }

  return PLL_FAILURE;
}

PLL_EXPORT int pll_utree_TBR(pll_utree_t * b_edge, pll_edge_t * r_edge)
{
  pll_edge_t parent, child;
  unsigned int parent_clv_index, child_clv_index;
  unsigned int parent_scaler_index, child_scaler_index;
  unsigned int edge_pmatrix_index, parent_pmatrix_index, child_pmatrix_index;

  /* validate if the move can be applied */

  /* 1. bisection point must not be a leaf branch */
  if (!(b_edge->next && b_edge->back->next))
  {
    pll_errno = PLL_ERROR_TBR_LEAF_BISECTION;
    return PLL_FAILURE;
  }

  /* 2. reconnection edges are different from bisection point */
  if (b_edge == r_edge->edge.utree.parent ||
      b_edge == r_edge->edge.utree.parent->back ||
      b_edge == r_edge->edge.utree.child ||
      b_edge == r_edge->edge.utree.child->back ||
      b_edge->back == r_edge->edge.utree.parent ||
      b_edge->back == r_edge->edge.utree.parent->back ||
      b_edge->back == r_edge->edge.utree.child ||
      b_edge->back == r_edge->edge.utree.child->back)
  {
    pll_errno = PLL_ERROR_TBR_OVERLAPPED_NODES;
    return PLL_FAILURE;
  }

  /* 3. reconnection edges must belong to different subtrees rooted at b_edge
   *    and b_edge->back
   */
  if (!(utree_find_node_in_subtree(b_edge, r_edge->edge.utree.parent) &&
        utree_find_node_in_subtree(b_edge->back, r_edge->edge.utree.child)) &&
      !(utree_find_node_in_subtree(b_edge->back, r_edge->edge.utree.parent) &&
        utree_find_node_in_subtree(b_edge, r_edge->edge.utree.child)))
  {
    pll_errno = PLL_ERROR_TBR_SAME_SUBTREE;
    return PLL_FAILURE;
  }

  /* keep removed CLVs, scalers and pmatrix indices */
  parent_clv_index = b_edge->clv_index;
  child_clv_index  = b_edge->back->clv_index;
  parent_scaler_index = b_edge->scaler_index;
  child_scaler_index = b_edge->back->scaler_index;

  edge_pmatrix_index = b_edge->pmatrix_index;

  /* bisect at b_edge */
  pll_utree_bisect(b_edge, &parent, &child);
  parent_pmatrix_index = parent.additional_pmatrix_index;
  child_pmatrix_index  = child.additional_pmatrix_index;

  /* reconnect at r_edge */
  pll_utree_reconnect(r_edge,
                      parent_pmatrix_index,
                      parent_clv_index,
                      parent_scaler_index,
                      child_pmatrix_index,
                      child_clv_index,
                      child_scaler_index,
                      edge_pmatrix_index);

  /* set CLV scalers and pmatrix indices */
  pll_utree_t *p, *q;
  p = r_edge->edge.utree.parent->back;
  q = r_edge->edge.utree.child->back;

  p->clv_index =
      p->next->clv_index =
      p->next->next->clv_index = parent_clv_index;
  p->scaler_index =
      p->next->scaler_index =
      p->next->next->scaler_index = parent_scaler_index;
  p->pmatrix_index = p->back->pmatrix_index = parent_pmatrix_index;

  q->clv_index =
        q->next->clv_index =
        q->next->next->clv_index = child_clv_index;
  q->scaler_index =
        q->next->scaler_index =
        q->next->next->scaler_index = child_scaler_index;
  q->pmatrix_index = q->back->pmatrix_index = child_pmatrix_index;

  return PLL_SUCCESS;
}

static void utree_nodes_at_dist(pll_utree_t * node,
                                pll_utree_t ** outbuffer,
                                unsigned int * index,
                                unsigned int distance,
                                unsigned int depth,
                                int fixed)
{
  if (depth == distance || !fixed)
  {
    outbuffer[*index] = node;
    *index = *index + 1;
  }

  if (depth >= distance || !(node->next)) return;

  utree_nodes_at_dist(node->next->back, outbuffer, index,
                      distance, depth+1, fixed);
  utree_nodes_at_dist(node->next->next->back, outbuffer, index,
                      distance, depth+1, fixed);
}

int pll_utree_nodes_at_node_dist(pll_utree_t * node,
                                 pll_utree_t ** outbuffer,
                                 unsigned int * n_nodes,
                                 unsigned int distance,
                                 int fixed)
{
  unsigned int depth = 0;
  if (!node->next) return PLL_FAILURE;

  *n_nodes = 0;

  /* we will traverse an unrooted tree in the following way

               1
             /
          --*
             \
               2
    */

  utree_nodes_at_dist(node, outbuffer, n_nodes, distance, depth, fixed);

  return PLL_SUCCESS;
}

PLL_EXPORT int pll_utree_nodes_at_edge_dist(pll_utree_t * root,
                                            pll_utree_t ** outbuffer,
                                            unsigned int * n_nodes,
                                            unsigned int distance,
                                            int fixed)
{
  unsigned int depth = 0;
  if (!root->next) return PLL_FAILURE;

  *n_nodes = 0;

  /* we will traverse an unrooted tree in the following way

       3          1
        \        /
         * ---- *
        /        \
       4          2
   */

  utree_nodes_at_dist(root->back, outbuffer, n_nodes, distance, depth+1, fixed);
  utree_nodes_at_dist(root, outbuffer, n_nodes, distance, depth, fixed);

  return PLL_SUCCESS;
}
