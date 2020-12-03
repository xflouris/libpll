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
#include <search.h>
#include <time.h>

#define RATE_CATS 4

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

int main(int argc, char * argv[])
{
  unsigned int i;
  unsigned int tip_nodes_count, inner_nodes_count, nodes_count, branch_count;
  pll_partition_t * partition;
  unsigned int arch = PLL_ATTRIB_ARCH_CPU;

  /* we accept only two arguments - the newick tree (unrooted binary) and the
     alignment in the form of FASTA reads */
  if (argc != 5)
    fatal(" syntax: %s [fasta] [seed] [attrib] [states]", argv[0]);

  /* open FASTA file */
  pll_fasta_t * fp = pll_fasta_open(argv[1], pll_map_fasta);
  if (!fp)
    fatal("Error opening file %s", argv[1]);

  if (!strcasecmp(argv[3], "avx"))
    arch = PLL_ATTRIB_ARCH_AVX;
  else if (!strcasecmp(argv[3], "avx2"))
    arch = PLL_ATTRIB_ARCH_AVX2;
  else if (!strcasecmp(argv[3], "sse"))
    arch = PLL_ATTRIB_ARCH_SSE;
  else if (!strcasecmp(argv[3], "tpavx"))
    arch = PLL_ATTRIB_ARCH_AVX | PLL_ATTRIB_PATTERN_TIP;
  else if (!strcasecmp(argv[3], "tpavx2"))
    arch = PLL_ATTRIB_ARCH_AVX2 | PLL_ATTRIB_PATTERN_TIP;
  else if (!strcasecmp(argv[3], "tpsse"))
    arch = PLL_ATTRIB_ARCH_SSE | PLL_ATTRIB_PATTERN_TIP;
  else if (!strcasecmp(argv[3], "tpcpu"))
    arch = PLL_ATTRIB_ARCH_CPU| PLL_ATTRIB_PATTERN_TIP;

  printf("Setting flags:\n");
  if (arch & PLL_ATTRIB_ARCH_CPU)
    printf("\tPLL_ATTRIB_ARCH_CPU\n");
  if (arch & PLL_ATTRIB_ARCH_SSE)
    printf("\tPLL_ATTRIB_ARCH_SSE\n");
  if (arch & PLL_ATTRIB_ARCH_AVX)
    printf("\tPLL_ATTRIB_ARCH_AVX\n");
  if (arch & PLL_ATTRIB_ARCH_AVX2)
    printf("\tPLL_ATTRIB_ARCH_AVX2\n");
  if (arch & PLL_ATTRIB_PATTERN_TIP)
    printf("\tPLL_ATTRIB_PATTERN_TIP\n");

  unsigned int states = atoi(argv[4]);
  unsigned const int * map = pll_map_nt;
  if (states == 20)
  {
    map = pll_map_aa;
  }

  char * seq = NULL;
  char * hdr = NULL;
  long seqlen;
  long hdrlen;
  long seqno;


  unsigned int max_alloc = 100;
  /* allocate arrays to store FASTA headers and sequences */
  char ** headers = (char **)calloc(max_alloc, sizeof(char *));
  char ** seqdata = (char **)calloc(max_alloc, sizeof(char *));

  /* read FASTA sequences and make sure they are all of the same length */
  int sites = -1;
  for (i = 0; pll_fasta_getnext(fp,&hdr,&hdrlen,&seq,&seqlen,&seqno); ++i)
  {
    if (i >= max_alloc)
    {
      max_alloc += 100;
      char ** tmpheaders = (char **)calloc(max_alloc, sizeof(char *));
      char ** tmpseqdata = (char **)calloc(max_alloc, sizeof(char *));
      memcpy(tmpheaders,headers,(max_alloc-100)*sizeof(char *));
      memcpy(tmpseqdata,seqdata,(max_alloc-100)*sizeof(char *));
      free(headers);
      free(seqdata);
      headers = tmpheaders;
      seqdata = tmpseqdata;
    }

    if (sites != -1 && sites != seqlen)
      fatal("FASTA file does not contain equal size sequences\n");

    if (sites == -1) sites = seqlen;

    headers[i] = hdr;
    seqdata[i] = seq;
  }

  /* did we stop reading the file because we reached EOF? */
  if (pll_errno != PLL_ERROR_FILE_EOF)
  {
    printf("errno %d %s\n", pll_errno, pll_errmsg);
    fatal("Error while reading file %s", argv[1]);
  }

  /* close FASTA file */
  pll_fasta_close(fp);

  if (sites == -1)
    fatal("Unable to read alignment");

  if (i <= 0)
    fatal("FASTA error");

  tip_nodes_count = i;
  inner_nodes_count = tip_nodes_count - 2;
  nodes_count = inner_nodes_count + tip_nodes_count;
  branch_count = nodes_count - 1;

  printf("Number of tip/leaf nodes in tree: %d\n", tip_nodes_count);
  printf("Number of inner nodes in tree: %d\n", inner_nodes_count);
  printf("Total number of nodes in tree: %d\n", nodes_count);
  printf("Number of branches in tree: %d\n", branch_count);
  printf("Number of sites: %d\n", sites);

  unsigned int * weight = pll_compress_site_patterns(seqdata,
                                                     map,
                                                     tip_nodes_count,
                                                     &sites);
  printf("Number of unique sites: %d\n", sites);


  /* create the PLL partition instance

  tip_nodes_count : the number of tip sequences we want to have
  inner_nodes_count : the number of CLV buffers to be allocated for inner nodes
  STATES : the number of states that our data have
  1 : number of different substitution models (or eigen decomposition)
      to use concurrently (i.e. 4 for LG4)
  branch_count: number of probability matrices to be allocated
  RATE_CATS : number of rate categories we will use
  inner_nodes_count : how many scale buffers to use
  PLL_ATTRIB_ARCH_SSE : list of flags for hardware acceleration (not yet implemented)

  */

  partition = pll_partition_create(tip_nodes_count,
                                   inner_nodes_count,
                                   states,
                                   (unsigned int)sites,
                                   1,
                                   branch_count,
                                   RATE_CATS,
                                   inner_nodes_count,
                                   arch);


  /* set pattern weights and free the weights array */
  pll_set_pattern_weights(partition,weight);
  free(weight);

  /* find sequences in hash table and link them with the corresponding taxa */
  for (i = 0; i < tip_nodes_count; ++i)
    pll_set_tip_states(partition, i, map, seqdata[i]);


  unsigned int score;

  pll_parsimony_t * parsimony = pll_fastparsimony_init(partition);

  pll_utree_t * tree = pll_fastparsimony_stepwise(&parsimony, headers, &score,1,atoi(argv[2]));
  printf("Score: %u\n", score);

  /* select a random inner node */
  long int r = random() % tree->inner_count;
  pll_unode_t * root = tree->nodes[tree->tip_count + r];

  /* export the tree to newick format with the selected inner node as the root
     of the unrooted binary tree */
  char * newick = pll_utree_export_newick(root,NULL);
  printf("%s\n", newick);
  free(newick);

  /* uncomment the following code to output an SVG file

  pll_svg_attrib_t * attr = pll_svg_attrib_create();
  attr->width = 1920;
  pll_utree_export_svg(tree,root,attr,"parsimony.svg");
  pll_svg_attrib_destroy(attr);

  */

  pll_parsimony_destroy(parsimony);

  /* destroy all structures allocated for the concrete PLL partition instance */
  pll_partition_destroy(partition);

  /* we will no longer need the tree structure */
  pll_utree_destroy(tree,NULL);
  /* ...neither the sequences and the headers as they are already
     present in the form of probabilities in the tip CLVs */
  for(i = 0; i < tip_nodes_count; ++i)
  {
    free(seqdata[i]);
    free(headers[i]);
  }
  free(seqdata);
  free(headers);


  return (EXIT_SUCCESS);
}
