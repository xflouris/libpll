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
    Exelixis Lab, Heidelberg Instutute for Theoretical Studies
    Schloss-Wolfsbrunnenweg 35, D-69118 Heidelberg, Germany
*/

#include "pll.h"

static int dlist_insert(pll_dlist_t ** dlist, void * data, int insert_end)
{
  if (!*dlist)
  {
    *dlist = (pll_dlist_t *)malloc(sizeof(pll_dlist_t));
    if (!*dlist)
    {
      pll_errno = PLL_ERROR_MEM_ALLOC;
      snprintf(pll_errmsg, 200, "Unable to allocate enough memory.");
      return PLL_FAILURE;
    }

    (*dlist)->next = NULL;
    (*dlist)->prev = NULL;
    (*dlist)->data = data;
    return PLL_SUCCESS;
  }

  /* go to the last element if we chose to append */
  if (insert_end)
    for(; (*dlist)->next; dlist = &(*dlist)->next);

  (*dlist)->next = (pll_dlist_t *)malloc(sizeof(pll_dlist_t));
  if (!(*dlist)->next)
  {
    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf(pll_errmsg, 200, "Unable to allocate enough memory.");
    return PLL_FAILURE;
  }

  (*dlist)->next->next = NULL;
  (*dlist)->next->data = data;
  (*dlist)->next->prev = (*dlist);

  return PLL_SUCCESS;

}

PLL_EXPORT int pll_dlist_append(pll_dlist_t ** dlist, void * data)
{
  return dlist_insert(dlist, data, 1);
}

PLL_EXPORT int pll_dlist_prepend(pll_dlist_t ** dlist, void * data)
{
  return dlist_insert(dlist, data, 0);
}

PLL_EXPORT int pll_dlist_remove(pll_dlist_t ** dlist, void * data)
{
  for (; (*dlist) && (*dlist)->data != data; dlist = &((*dlist)->next));

  if (!*dlist) return PLL_FAILURE;

  if ((*dlist)->next)
    (*dlist)->next->prev = (*dlist)->prev;
  if ((*dlist)->prev)
    (*dlist)->prev->next = (*dlist)->next;

  free(*dlist);
  *dlist = NULL;

  return PLL_SUCCESS;
}
