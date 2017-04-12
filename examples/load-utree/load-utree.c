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
#include <stdarg.h>

static void fatal(const char * format, ...) __attribute__ ((noreturn));

static void fatal(const char * format, ...)
{
  va_list argptr;
  va_start(argptr, format);
  vfprintf(stderr, format, argptr);
  va_end(argptr);
  fprintf(stderr, "\n");
  exit(EXIT_FAILURE);
}

pll_utree_t * load_tree_unrooted(const char * filename)
{
  pll_utree_t * utree;
  pll_rtree_t * rtree;

  if (!(rtree = pll_rtree_parse_newick(filename)))
  {
    if (!(utree = pll_utree_parse_newick(filename)))
    {
      fprintf(stderr, "%s\n", pll_errmsg);
      return NULL;
    }
  }
  else
  {
    utree = pll_rtree_unroot(rtree);

    pll_unode_t * root = utree->nodes[utree->tip_count+utree->inner_count-1];

    /* optional step if using default PLL clv/pmatrix index assignments */
    pll_utree_reset_template_indices(root, utree->tip_count);

    pll_rtree_destroy(rtree,NULL);
  }

  return utree;
}

int main(int argc, char * argv[])
{
  if (argc != 2)
    fatal("syntax: %s [newick]", argv[0]);

  pll_utree_t * utree = load_tree_unrooted(argv[1]);
  if (!utree)
    fatal("Tree must be a rooted or unrooted binary.");

  /* select a random inner node */
  long int r = random() % utree->inner_count;
  pll_unode_t * root = utree->nodes[utree->tip_count + r];

  /* export the tree to newick format with the selected inner node as the root
     of the unrooted binary tree */
  char * newick = pll_utree_export_newick(root,NULL);

  printf("%s\n", newick);

  free(newick);

  pll_utree_destroy(utree,NULL);

  return 0;
}
