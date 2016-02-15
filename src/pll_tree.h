/*
 * pll_tree.h
 *
 *  Created on: Sep 28, 2015
 *      Author: diego
 */

#ifndef PLL_TREE_H_
#define PLL_TREE_H_

#ifndef PLL_H_
#define PLL_H_
#include "pll.h"
#endif

/* the bisection edge is a leaf */
#define PLL_ERROR_TBR_LEAF_BISECTION   5001
/* the bisection and reconnection points are overlapped */
#define PLL_ERROR_TBR_OVERLAPPED_NODES 5002
/* the reconnection branches belong to the same subtree */
#define PLL_ERROR_TBR_SAME_SUBTREE     5003
/* attempting to interchange a leaf */
#define PLL_ERROR_INTERCHANGE_LEAF     6001

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
  int additional_pmatrix_index;
  double length;
} pll_edge_t;

/**
 * Bisects the tree by removing one edge
 *
 * Removes the edge \p edge and frees the nodes defining that edge.
 * Reconnects the subtrees at the sides of the edge (figure below).
 * The branch lengths of the new edges are the sum of the removed ones.
 * The join branch contains the pmatrix index of the parent edges
 * The removed pmatrix indices are returned in the field
 *     'additional_pmatrix_index' of both output subtrees
 *
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
 * Adds 1 new edge connecting edges \p edge.parent and \p edge.child with
 * length \p edge.length.
 *
 *   A       C         A              C
 *   |       |  ---->   \            /
 *                       e1--edge--e2
 *   |       |          /            \
 *   B       D         B              D
 *   A,B,C,D are subtrees
 *
 * @param[in] edge                 new edge
 * @param[in] parent_pmatrix_index matrix index of e1-A
 * @param[in] parent_clv_index     clv index of e1
 * @param[in] parent_scaler_index  scaler index of e1
 * @param[in] child_pmatrix_index  matrix index of e2-C
 * @param[in] child_clv_index      clv index of e2
 * @param[in] child_scaler_index   scaler index of e2
 * @param[in] edge_pmatrix_index   matrix index of e1-e2
 *
 * @returns the new created edge
 */
PLL_EXPORT pll_edge_t pll_utree_reconnect(pll_edge_t * edge,
                                          unsigned int parent_pmatrix_index,
                                          unsigned int parent_clv_index,
                                          unsigned int parent_scaler_index,
                                          unsigned int child_pmatrix_index,
                                          unsigned int child_clv_index,
                                          unsigned int child_scaler_index,
                                          unsigned int edge_pmatrix_index);

/**
 * Returns the list of nodes at a certain distance from a specified edge
 *
 * @param[in] edge the root edge
 * @param[out] outbuffer the list of nodes. Outbuffer should be allocated
 * @param[out] n_nodes the number of nodes returned in \p outbuffer
 * @paran[in] distance the maximum distance to check
 * @param[in] fixed if true, returns only the nodes at distance \p distance,
 *            otherwise, the nodes at distance <= \p distance.
 */
PLL_EXPORT int pll_utree_nodes_at_edge_dist(pll_utree_t * edge,
                                            pll_utree_t ** outbuffer,
                                            unsigned int * n_nodes,
                                            unsigned int distance,
                                            int fixed);

/**
 * Returns the list of nodes at a certain distance from a specified node
 *
 * @param[in] node the root node
 * @param[out] outbuffer the list of nodes. Outbuffer should be allocated
 * @param[out] n_nodes the number of nodes returned in \p outbuffer
 * @paran[in] distance the maximum distance to check
 * @param[in] fixed if true, returns only the nodes at distance \p distance,
 *            otherwise, the nodes at distance <= \p distance.
 */
PLL_EXPORT int pll_utree_nodes_at_node_dist(pll_utree_t * node,
                                            pll_utree_t ** outbuffer,
                                            unsigned int * n_nodes,
                                            unsigned int distance,
                                            int fixed);

/**
 * Performs one TBR move by applying a bisection and a reconnection.
 * The CLV, scaler and pmatrix indices are updated.
 *
 * @param[in] b_edge bisection point
 * @param[in] r_edge reconnection point
 *
 * @returns true, if the move was applied correctly
 */
PLL_EXPORT int pll_utree_TBR(pll_utree_t * b_edge, pll_edge_t * r_edge);

/**
 * Interchanges 2 edges, represented by 2 internal nodes
 *
 * CLV and scaler indices, and labels are interchanged between nodes to match
 * the other 2 nodes in the triplet.
 *
 * @returns true, if the move was applied correctly
 */
PLL_EXPORT int pll_utree_interchange(pll_utree_t * edge1,
                                     pll_utree_t * edge2);
#endif /* PLL_TREE_H_ */
