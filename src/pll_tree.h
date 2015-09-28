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
 * @param[out] parent_subtree edge corresponding to the 'edge' subtree
 * @param[out] child_subtree  edge corresponding to the 'edge->back' subtree
 * @returns PLL_SUCCESS if OK
 */
PLL_EXPORT int pll_utree_bisect(pll_utree_t * edge,
                                        pll_edge_t * parent_subtree,
                                        pll_edge_t * child_subtree);
PLL_EXPORT void pll_utree_reconnect(pll_edge_t * edge);



#endif /* PLL_TREE_H_ */
