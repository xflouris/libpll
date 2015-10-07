/*
 * pll_tree.h
 *
 *  Created on: Sep 28, 2015
 *      Author: diego
 */

#ifndef PLL_TREE_H_
#define PLL_TREE_H_

#include <pll.h>

typedef struct pll_edge
{
  union
  {
    struct
    {
      pll_utree_t * parent;
      pll_utree_t * child;
    } utree;
    struct
    {
      pll_rtree_t * parent;
      pll_rtree_t * child;
    } rtree;
  } edge;
  double length;
} pll_edge_t;

/**
 * Bisects the tree by removing one edge
 *
 * Removes the edge \p edge and frees the nodes defining that edge.
 * Reconnects the subtrees at the sides of the edge (figure below).
 * The branch lengths of the new edges are the sum of the removed ones.
 * Returns the new parent and child edges, where parent is the closest to \p edge.
 *
 *   A            C              A        C
 *    \___edge___/       ---->   |        |
 *    /          \               |        |
 *   B            D              B        D
 *   A,B,C,D are subtrees
 *
 * @param[in] edge edge to remove
 * @param[out] parent_subtree edge corresponding to the 'edge' subtree
 * @param[out] child_subtree  edge corresponding to the 'edge->back' subtree
 * @returns PLL_SUCCESS if OK
 */
PLL_EXPORT int pll_utree_bisect(pll_utree_t * edge,
                                pll_edge_t * parent_subtree,
                                pll_edge_t * child_subtree);

/**
 * Reconnects two subtrees by adding 2 new nodes and 1 edge.
 *
 * Addes 1 new edge connecting edges \p edge.parent and \p edge.child with
 * length \p edge.length.
 *
 *   A       C         A            C
 *   |       |  ---->   \___edge___/
 *   |       |          /          \
 *   B       D         B            D
 *   A,B,C,D are subtrees
 *
 * @param[in] edge new edge
 *
 * @returns PLL_SUCCESS if OK
 */
PLL_EXPORT void pll_utree_reconnect(pll_edge_t * edge);

/**
 * Returns the list of nodes at a certain distance from a specified root
 *
 * @param[in] root the root node
 * @param[out] outbuffer the list of nodes. Outbuffer should be allocated
 * @param[out] n_nodes the number of nodes returned in \p outbuffer
 * @paran[in] distance the maximum distance to check
 * @param[in] fixed if true, returns only the nodes at distance \p distance,
 *            otherwise, the nodes at distance <= \p distance.
 */
PLL_EXPORT int pll_utree_nodes_at_dist(pll_utree_t * root,
                                       pll_utree_t ** outbuffer,
                                       unsigned int * n_nodes,
                                       unsigned int distance,
                                       int fixed);

PLL_EXPORT void pll_utree_TBR(pll_utree_t * b_edge, pll_edge_t * r_edge);

#endif /* PLL_TREE_H_ */
