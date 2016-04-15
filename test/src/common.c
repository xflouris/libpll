#include "common.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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
