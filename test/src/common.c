#include "common.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <search.h>

unsigned int get_attributes(int argc, char **argv)
{
  int i;
  unsigned int attributes = PLL_ATTRIB_ARCH_CPU;

  for (i=1; i<argc; ++i)
  {
    if (!strcmp (argv[i], "tv"))
    {
      /* tipvector */
      attributes |= PLL_ATTRIB_PATTERN_TIP;
    }
    else if (!strcmp (argv[i], "avx"))
    {
      /* avx vectorization */
      attributes |= PLL_ATTRIB_ARCH_AVX;
    }
    else if (!strcmp (argv[i], "sse"))
    {
      /* sse3 vectorization */
      attributes |= PLL_ATTRIB_ARCH_SSE;
    }
    else
    {
      printf("Unrecognised attribute: %s\n", argv[i]);
      exit(1);
    }
  }
    return attributes;
}

void skip_test ()
{
  printf ("Skip\n");
  exit (0);
}

/* note: There is no exhaustive error checking in this function.
         It is intended to use with the test datasets that were
         validated in advance. */
pll_partition_t * parse_msa(const char * filename,
                            unsigned int taxa_count,
                            unsigned int states,
                            unsigned int rate_cats,
                            unsigned int rate_matrices,
                            pll_utree_t * tree,
                            unsigned int attributes)
{
  unsigned int i;
  pll_partition_t * partition;
  long hdrlen, seqlen, seqno;
  char * seq = NULL,
       * hdr = NULL;

  /* open FASTA file */
  pll_fasta_t * fp = pll_fasta_open(filename, pll_map_fasta);
  if (!fp)
  {
    printf("Error opening file %s", filename);
    return NULL;
  }

  /* allocate arrays to store FASTA headers and sequences */
  char ** headers = (char **)calloc(taxa_count, sizeof(char *));
  char ** seqdata = (char **)calloc(taxa_count, sizeof(char *));

  /* read FASTA sequences and make sure they are all of the same length */
  int sites = -1;
  for (i = 0; pll_fasta_getnext(fp,&hdr,&hdrlen,&seq,&seqlen,&seqno); ++i)
  {
    if (sites == -1) sites = seqlen;

    headers[i] = hdr;
    seqdata[i] = seq;
  }

  if (pll_errno != PLL_ERROR_FILE_EOF)
  {
    printf("Error while reading file %s", filename);
    return NULL;
  }

  /* close FASTA file */
  pll_fasta_close(fp);

  if (sites == -1)
  {
    printf("Unable to read alignment");
    return NULL;
  }

  partition = pll_partition_create(taxa_count,           /* tip nodes */
                                   taxa_count - 2,       /* inner nodes */
                                   states,
                                   (unsigned int)sites,
                                   rate_matrices,        /* rate matrices */
                                   2*taxa_count - 3,     /* prob matrices */
                                   rate_cats,            /* rate categories */
                                   taxa_count - 2,       /* scale buffers */
                                   attributes);

  pll_utree_t ** tipnodes = (pll_utree_t  **)calloc(taxa_count,
                                                    sizeof(pll_utree_t *));
  pll_utree_query_tipnodes(tree, tipnodes);

  /* create a libc hash table of size tip_nodes_count */
  hcreate(taxa_count);

  /* populate a libc hash table with tree tip labels */
  unsigned int * data = (unsigned int *)malloc(taxa_count *
                                               sizeof(unsigned int));
  for (i = 0; i < taxa_count; ++i)
  {
    data[i] = i;
    ENTRY entry;
    entry.key = tipnodes[i]->label;
    entry.data = (void *)(data+i);
    hsearch(entry, ENTER);
  }

  free(tipnodes);

  for (i = 0; i < taxa_count; ++i)
  {
    ENTRY query;
    query.key = headers[i];
    ENTRY * found = NULL;

    found = hsearch(query,FIND);

    if (!found)
    {
      printf("Sequence with header %s does not appear in the tree", headers[i]);
      return NULL;
    }

    unsigned int tip_clv_index = *((unsigned int *)(found->data));

    const unsigned int * map = states == 4? pll_map_nt : pll_map_aa;
    pll_set_tip_states(partition, tip_clv_index, map, seqdata[i]);

    free(headers[i]);
    free(seqdata[i]);
  }

  /* clean */
  hdestroy();
  free(data);
  free(headers);
  free(seqdata);

  return partition;
}

int cb_full_traversal(pll_utree_t * node)
{
  return 1;
}

int cb_rfull_traversal(pll_rtree_t * node)
{
  return 1;
}
