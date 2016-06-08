/*
    Copyright (C) 2015 Tomas Flouri, Diego Darriba

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

static int utree_find(pll_utree_t * start, pll_utree_t * target)
{
  /* checks whether the subtree rooted at 'start' (in the direction of 
     start->next and start->next->next) contains the node 'target' */

  if (!start) return 0;

  if (start == target) return 1;

  if (start->next)
  {
    if (start->next == target) return 1;
    if (utree_find(start->next->back, target)) return 1;
  }
  else
    return 0;

  if (start->next->next == target) return 1;
  if (utree_find(start->next->next->back, target)) return 1;

  return 0;
}

static void utree_link(pll_utree_t * a, pll_utree_t * b, double length)
{
  a->back = b;
  b->back = a;
  a->length = length;
  b->length = length;
}

/* drive fast and you will crash if not careful */
PLL_EXPORT int pll_utree_spr(pll_utree_t * p,
                             pll_utree_t * r,
                             pll_utree_rb_spr_t * rb)
{
  /* perform an SPR move in the following way:

      A           B          C             D           A          B
     ____        ____       ____          ____        ____       ____
     \  /        \  /       \  /          \  /        \  /       \  /
      \/          \/         \/            \/          \/         \/
       *          *          * p'           *          *          *
        \         |     q   /                \         |         /
         *'*_____.*._____*'* p     --->       *'*_____.*._____*'*
         '*'     *.*     '*'                  '*'     *.*     '*'
         / r       u    q' \                  /                 \
     r' *                   * v              *                   *
       /\                   /\              /\                   /\
      /__\                 /__\            /__\                 /__\
                                                                     
       D                    E               C                    E   
     
     node p must be part of an inner node (i.e. node with ->next set). The
     procedure prunes the subtree rooted at the opposite end-point of p
     (subtree C in our case) and regrafts it on the edge r'<->r. It is done
     in the following way:

     (a) prune the subtree rooted at the opposite end-point of p (p' on figure)
         by breaking the edges q<->u and q'<->v
     
     (b) connect node u with node v

     (c) break edge r<->r' by connecting node r with node q, and node r' with
         node q' 

     Node r must not be part of the subtree to be pruned (C in this case). Note
     that for speed reasons, the function *does not* check this property to save
     a tree traversal. A safer (albeit slower) function that checks this
     property is pll_utree_spr_safe
  */

  /* if p is a tip node then prompt an error */
  if (!p->next)
  {
    snprintf(pll_errmsg, 200, "Prune edge must be defined by an inner node");
    return PLL_FAILURE;
  }

  /* check whether the move will result in the same tree */
  if (r == p || r == p->back ||
      r == p->next || r == p->next->back ||
      r == p->next->next || r == p->next->next->back)
  {
    snprintf(pll_errmsg, 200, "Proposed move yields the same tree");
    return PLL_FAILURE;
  }

  /* check if rollback buffer is provided, and fill it up */
  if (rb)
  {
    rb->p = p;
    rb->r = r;
    rb->rb = r->back;
    rb->r_len = r->length;
    rb->pnb = p->next->back;
    rb->pnb_len = p->next->length;
    rb->pnnb = p->next->next->back;
    rb->pnnb_len = p->next->next->length;
  }

  /* (b) connect u and v */
  pll_utree_t * u = p->next->back;
  pll_utree_t * v = p->next->next->back;
  utree_link(u, v, u->length + v->length);

  /* (a) prune subtree C */
  p->next->back = p->next->next->back = NULL;

  /* (c) regraft C at r<->r' */
  double length = r->length / 2;
  utree_link(r->back, p->next->next, length);  /* r'<->q' */
  utree_link(r, p->next, length);              /* r<->q */
  
  return PLL_SUCCESS;
}

PLL_EXPORT void pll_utree_spr_rollback(pll_utree_rb_spr_t * rb)
{
  /* restore the tree topology from a previous SPR */
  utree_link(rb->pnb, rb->p->next, rb->pnb_len);
  utree_link(rb->pnnb, rb->p->next->next, rb->pnnb_len);
  utree_link(rb->r, rb->rb, rb->r_len);
}

/* drive slow and safe even with dope at hand */
PLL_EXPORT int pll_utree_spr_safe(pll_utree_t * p,
                                  pll_utree_t * r,
                                  pll_utree_rb_spr_t * rb)
{
  /* check all possible scenarios of failure */
  if (!p)
  {
    snprintf(pll_errmsg, 200, "Node p is set to NULL"); 
    return PLL_FAILURE;
  }

  if (!r)
  {
    snprintf(pll_errmsg, 200, "Node r is set to NULL"); 
    return PLL_FAILURE;
  }

  if (!p->next)
  {
    snprintf(pll_errmsg, 200, "Prune edge must be defined by an inner node");
    return PLL_FAILURE;
  }

  /* check whether the move results in the same tree */
  if (r == p || r == p->back ||
      r == p->next || r == p->next->back ||
      r == p->next->next || r == p->next->next->back)
  {
    snprintf(pll_errmsg, 200, "Proposed move yields the same tree");
    return PLL_FAILURE;
  }

  /* node r must not be in the same subtree as the one that is to be pruned */
  if (utree_find(p->back, r))
  {
    snprintf(pll_errmsg, 200, "Node r is part of the subtree to be pruned");
    return PLL_FAILURE;
  }

  return pll_utree_spr(p,r,rb);
}
