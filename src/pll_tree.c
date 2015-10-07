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
  new_node->label = (char *) malloc(20);
  snprintf(new_node->label, 4, "NEW");
  new_node->next->label = new_node->next->next->label = new_node->label;
  new_node->next->data = new_node->next->next->data = new_node->data = NULL;
  new_node->next->length = new_node->next->next->length = new_node->length = 0.1;
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
  printf("\nConnect edges %s, %s\n", edge->next->back->label, edge->next->next->back->label);
  parent_subtree->edge.utree.parent = edge->next->back;
  parent_subtree->edge.utree.child = edge->next->next->back;
  parent_subtree->length =
      parent_subtree->edge.utree.parent->length +
      parent_subtree->edge.utree.child->length;

  utree_connect_nodes(parent_subtree->edge.utree.parent,
                      parent_subtree->edge.utree.child,
                      parent_subtree->length);

  /* connect child subtree */
  printf("\nConnect edges %s, %s\n", c_edge->next->back->label, c_edge->next->next->back->label);
  child_subtree->edge.utree.parent = c_edge->next->back;
  child_subtree->edge.utree.child = c_edge->next->next->back;
  child_subtree->length =
      child_subtree->edge.utree.parent->length +
      child_subtree->edge.utree.child->length;

  utree_connect_nodes(child_subtree->edge.utree.parent,
                      child_subtree->edge.utree.child,
                      child_subtree->length);

  /* remove edges */
  printf("\nFree edges %s, %s\n", edge->label, c_edge->label);
  utree_destroy_node(edge);
  utree_destroy_node(c_edge);
  return PLL_SUCCESS;
}

PLL_EXPORT void pll_utree_reconnect(pll_edge_t * edge)
{
  /* create and connect 2 new nodes */
  pll_utree_t * parent_node = utree_create_node();
  pll_utree_t * child_node = utree_create_node();

  printf("Reconnect edges %s, %s\n", parent_node->label, child_node->label);
  utree_connect_nodes(parent_node, child_node, edge->length);

  /* reconnect parent close to edge.parent */
  printf("Reconnect edges[p] %s, %s l: %f\n", parent_node->next->next->label, edge->edge.utree.parent->back->label,
         edge->edge.utree.parent->back->length);
  utree_connect_nodes(parent_node->next->next,
                      edge->edge.utree.parent->back,
                      edge->edge.utree.parent->back->length);
  printf("Reconnect edges[p] %s, %s\n", parent_node->next->label, edge->edge.utree.parent->label);
  utree_connect_nodes(parent_node->next, edge->edge.utree.parent, 0);

  /* reconnect child close to edge.child */
  printf("Reconnect edges[c] %s, %s l: %f\n", child_node->next->next->label, edge->edge.utree.child->back->label,
         edge->edge.utree.child->back->length);
  utree_connect_nodes(child_node->next->next,
                      edge->edge.utree.child->back,
                      edge->edge.utree.child->back->length);
  printf("Reconnect edges[c] %s, %s\n", child_node->next->label, edge->edge.utree.child->label);
  utree_connect_nodes(child_node->next, edge->edge.utree.child, 0);
}

PLL_EXPORT void pll_utree_TBR(pll_utree_t * b_edge, pll_edge_t * r_edge)
{
  pll_edge_t parent, child;
  unsigned int parent_clv_index, child_clv_index;
  unsigned int parent_scaler_index, child_scaler_index;

  /* keep removed CLVs */
  parent_clv_index = b_edge->clv_index;
  child_clv_index  = b_edge->back->clv_index;
  parent_scaler_index = b_edge->scaler_index;
  child_scaler_index = b_edge->back->scaler_index;

  /* bisect at b_edge */
  pll_utree_bisect(b_edge, &parent, &child);

  /* reconnect at r_edge */
  pll_utree_reconnect(r_edge);

  /* set CLV indices */
  pll_utree_t *p, *q;
  p = r_edge->edge.utree.parent->back;
  p->clv_index =
      p->next->clv_index =
      p->next->next->clv_index = parent_clv_index;
  p->scaler_index =
            p->next->scaler_index =
            p->next->next->scaler_index = parent_scaler_index;

  q = r_edge->edge.utree.child->back;
  q->clv_index =
        q->next->clv_index =
        q->next->next->clv_index = child_clv_index;
  q->scaler_index =
          q->next->scaler_index =
          q->next->next->scaler_index = child_scaler_index;
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

PLL_EXPORT int pll_utree_nodes_at_dist(pll_utree_t * root,
                                       pll_utree_t ** outbuffer,
                                       unsigned int * n_nodes,
                                       unsigned int distance,
                                       int fixed)
{
  unsigned int depth = 0;
  if (!root->next) return PLL_FAILURE;

  *n_nodes = 0;

  /* we will traverse an unrooted tree in the following way

              2
            /
      1  --*
            \
              3
   */

  utree_nodes_at_dist(root->back, outbuffer, n_nodes, distance, depth+1, fixed);
  utree_nodes_at_dist(root, outbuffer, n_nodes, distance, depth, fixed);

  return PLL_SUCCESS;
}
