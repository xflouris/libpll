#include "pll_tree.h"

/**
 * creates a new circular node
 */
static pll_utree_t * utree_create_node()
{
  pll_utree_t * new_node = (pll_utree_t *)calloc(1, sizeof(pll_utree_t));
  new_node->next         = (pll_utree_t *)calloc(1, sizeof(pll_utree_t));
  new_node->next->next   = (pll_utree_t *)calloc(1, sizeof(pll_utree_t));
  new_node->next->next->next = new_node;
  return new_node;
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
  parent->length = length;
  child->length = length;
  return PLL_SUCCESS;
}

PLL_EXPORT int pll_utree_bisect(pll_utree_t * edge,
                                pll_edge_t * parent_subtree,
                                pll_edge_t * child_subtree)
{
  assert(parent_subtree);
  assert(child_subtree);

  pll_utree_t c_edge = edge->back;

  /* connect parent subtree */
  parent_subtree->edge.utree.parent = edge->next->back;
  parent_subtree->edge.utree.child = edge->next->next->back;
  parent_subtree->length =
      parent_subtree->edge.utree.parent->length +
      parent_subtree->edge.utree.child->length;

  parent_subtree->edge.utree.parent->back   = parent_subtree->edge.utree.child;
  parent_subtree->edge.utree.child->back    = parent_subtree->edge.utree.parent;
  parent_subtree->edge.utree.parent->length = parent_subtree->length;
  parent_subtree->edge.utree.child->length  = parent_subtree->length;

  /* connect child subtree */
  child_subtree->edge.utree.parent = c_edge->next->back;
  child_subtree->edge.utree.child = c_edge->next->next->back;
  child_subtree->length =
      child_subtree->edge.utree.parent->length +
      child_subtree->edge.utree.child->length;

  child_subtree->edge.utree.parent->back   = child_subtree->edge.utree.child;
  child_subtree->edge.utree.child->back    = child_subtree->edge.utree.parent;
  child_subtree->edge.utree.parent->length = child_subtree->length;
  child_subtree->edge.utree.child->length  = child_subtree->length;

  /* remove edges */
  free(edge);
  free(c_edge);
  return c_edge;
}

PLL_EXPORT void pll_utree_reconnect(pll_edge_t * edge)
{
  /* create and connect 2 new nodes */
  pll_utree_t * parent_node = utree_create_node();
  pll_utree_t * child_node = utree_create_node();
  utree_connect_nodes(parent_node, child_node, edge->length);

  /* reconnect parent close to edge.parent */
  utree_connect_nodes(parent_node->next, edge->edge.utree.parent, 0);
  utree_connect_nodes(parent_node->next,
                      edge->edge.utree.parent->back,
                      edge->edge.utree.parent->back->length);

  /* reconnect child close to edge.child */
  utree_connect_nodes(child_node->next, edge->edge.utree.child, 0);
  utree_connect_nodes(child_node->next,
                      edge->edge.utree.child->back,
                      edge->edge.utree.child->back->length);
}
