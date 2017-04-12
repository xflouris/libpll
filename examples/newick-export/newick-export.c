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
#include <time.h>

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

typedef struct nodedata_s
{
  double support;
  double rvalue;
} nodedata_t;

/* a callback function that prints the following information for the two types
 * of nodes, tips and inners:


   for tip nodes:
     label, branch length and support value
   
   for inner nodes:
     label, branch length, support value, and a random value 

   Note that this is not a standard format that can be read by all tree viewers

   Check the extended newick format for a more widespread notation:
   http://dmi.uib.es/~gcardona/BioInfo/enewick.html

*/
static char * cb_serialize(const pll_unode_t * node)
{
  char * s = NULL;
  nodedata_t * nd;


  /* inner node */
  if (node->next)
  {
    /* find which node of the round-about has the data element */
    if (node->data)
      nd = (nodedata_t *)(node->data);
    else if (node->next->data)
      nd = (nodedata_t *)(node->next->data);
    else
      nd = (nodedata_t *)(node->next->next->data);

    asprintf(&s,
             "%s[&support=%f,rvalue=%f]:%f",
             node->label ? node->label : "",
             nd->support,
             nd->rvalue,
             node->length);
  }
  else
  {
    nd = (nodedata_t *)(node->data);
    asprintf(&s,
             "%s[&support=%f]:%f",
             node->label ? node->label : "",
             nd->support,
             node->length);
  }

  return s;
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
  unsigned int i;

  if (argc != 2)
    fatal("syntax: %s [newick]", argv[0]);

  /* initialize pseudo-random number generator */
  srandom(time(NULL));

  /* load tree as unrooted */
  pll_utree_t * utree = load_tree_unrooted(argv[1]);
  if (!utree)
    fatal("Tree must be a rooted or unrooted binary.");

  /* set random support values for tip nodes */
  for (i = 0; i < utree->tip_count; ++i)
  {
    pll_unode_t * node = utree->nodes[i];
    nodedata_t * data;

    node->data = malloc(sizeof(nodedata_t));

    data = (nodedata_t *)(node->data);

    data->support = random() / (double)RAND_MAX;
  }

  /* set random support values and a random value to each inner node, but
     allocate the structure in only one of the three round-about nodes that make
     up an inner node */
  for (i = utree->tip_count; i < utree->tip_count + utree->inner_count; ++i)
  {
    pll_unode_t * node = utree->nodes[i];
    nodedata_t * data;

    node->data = malloc(sizeof(nodedata_t));

    data = (nodedata_t *)(node->data);

    data->support = random() / (double)RAND_MAX;
    data->rvalue = data->support * random();
  }

  /* select a random inner node */
  long int r = random() % utree->inner_count;
  pll_unode_t * root = utree->nodes[utree->tip_count + r];

  /* export the tree to newick format with the selected inner node as the root
     of the unrooted binary tree.
     
     If we pass NULL as the callback function, then the tree is printed in newick
     format with branch lengths only, i.e.
  
     char * newick = pll_utree_export_newick(root,NULL);
  */

  char * newick = pll_utree_export_newick(root,cb_serialize);

  printf("%s\n", newick);

  free(newick);

  pll_utree_destroy(utree,free);

  return 0;
}
