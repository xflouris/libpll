/*
    Copyright (C) 2016 Tomas Flouri

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

static void vecswap(int i, int j, int n, char ** x)
{
  while (n--)
  {
    PLL_SWAP(x[i],x[j]);
    ++i; ++j;
  }
}

static void ssort1(char ** x, int n, int depth)
{
  int a,b,c,d,r,v;

  if (n <= 1) return;

  a = rand() % n;

  PLL_SWAP(x[0], x[a]);

  v = x[0][depth];

  a = b = 1;
  c = d = n-1;

  while (1)
  {
    while (b <= c && (r = x[b][depth]-v) <= 0)
    {
      if (r == 0)
      {
        PLL_SWAP(x[a], x[b]);
        ++a;
      }
      ++b;
    }
    while (b <= c && (r = x[c][depth]-v) >= 0)
    {
      if (r == 0)
      {
        PLL_SWAP(x[c], x[d]);
        --d;
      }
      --c;
    }
    if (b > c) break;
    PLL_SWAP (x[b], x[c]);
    ++b; --c;
  }

  r = PLL_MIN(a,b-a); vecswap(0,b-r,r,x);
  r = PLL_MIN(d-c,n-d-1); vecswap(b,n-r,r,x);
  r = b-a; ssort1(x,r,depth);

  if (x[r][depth] != 0)
    ssort1 (x + r, a + n - d - 1, depth + 1);

  r = d - c; ssort1(x+n-r,r,depth);
}

static void remap_range(const unsigned int * map,
                        unsigned char * charmap)
{
  unsigned int oldmap[PLL_ASCII_SIZE];
  unsigned int i,j;
  unsigned char k = 1;

  memcpy(oldmap, map, PLL_ASCII_SIZE * sizeof(unsigned int));
  memset(charmap, 0, PLL_ASCII_SIZE * sizeof(unsigned char));

  for (i = 0; i < PLL_ASCII_SIZE; ++i)
    if (oldmap[i])
    {
      charmap[i] = k;

      for (j = i+1; j < PLL_ASCII_SIZE; ++j)
        if (oldmap[i] == oldmap[j])
        {
          charmap[j] = k;
          oldmap[j] = 0;
        }

      ++k;
    }
}

static unsigned int findmax(const unsigned int * map)
{
  int i;
  unsigned int max = 0;

  for (i = 0; i < PLL_ASCII_SIZE; ++i)
    if (map[i] > max)
      max = map[i];

  return max;
}

static void encode(char ** sequence, const unsigned char * map, int count, int len)
{
  int i,j;
  char * p;

  for (i = 0; i < count; ++i)
  {
    p = sequence[i];
    j = len;
    while (j--)
    {
      *p = map[(int)(*p)];
      ++p;
    }
  }
}

PLL_EXPORT unsigned int * pll_compress_site_patterns(char ** sequence,
                                                     const unsigned int * map,
                                                     int count,
                                                     int * length)
{
  int i,j;
  char * memptr;
  char ** column;
  unsigned int * weight;

  unsigned char charmap[PLL_ASCII_SIZE];
  unsigned char inv_charmap[PLL_ASCII_SIZE];

  /* check that at least one sequence is given */
  if (!count) return NULL;

  /* a map must be given */
  if (!map) return NULL;

  /* a zero can never be used as a state */
  if (map[0]) return NULL;

  /* if map states are out of the BYTE range, remap */
  if (findmax(map) >= PLL_ASCII_SIZE)
  {
    remap_range(map,charmap);
  }
  else
  {
    for (i = 0; i < PLL_ASCII_SIZE; ++i)
      charmap[i] = (unsigned char)(map[i]);
  }

  /* create inverse charmap to decode states back to characters when
     compression is finished */
  for (i = 0; i < PLL_ASCII_SIZE; ++i)
    if (map[i])
      inv_charmap[charmap[i]] = (unsigned char)i;

  /* encode sequences using charmap */
  encode(sequence,charmap,count,*length);

  /* allocate memory for columns */
  column = (char **)malloc((size_t)(*length)*sizeof(char *));
  if (!column)
  {
    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf (pll_errmsg, 200,
              "Cannot allocate space for matrix columns.");
    return NULL;
  }

  /* allocate memory for the alignment */
  memptr = column[0] = (char *)malloc((size_t)(*length) *
                                      (size_t)(count+1) *
                                      sizeof(char));
  if (!memptr)
  {
    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf (pll_errmsg, 200,
              "Cannot allocate space for matrix data.");
    free(column);
    return NULL;
  }

  /* map memory to each column */
  for (i = 1; i < *length; ++i)
  {
    column[i] = column[i-1] + (count+1);
  }

  /* allocate space for weight vector */
  weight = (unsigned int *)malloc((size_t)(*length)*sizeof(unsigned int));
  if (!weight)
  {
    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf(pll_errmsg, 200,
             "Cannot allocate space for storing site weights.");
    free(column);
    free(memptr);
    return NULL;
  }

  /* split alignment into columns instead of rows */
  for (i = 0; i < (*length); ++i)
  {
    for (j = 0; j < count; ++j)
      column[i][j] = sequence[j][i];
    column[i][j] = 0;
  }

  /* sort the columns */
  ssort1(column, *length, 0);

  /* we have at least one unique site with weight 1 (the first site) */
  int compressed_length = 1;
  size_t ref = 0;
  weight[ref] = 1;

  /* find all unique columns and set their weights */
  for (i = 1; i < *length; ++i)
  {
    if (strcmp(column[i],column[i-1]))
    {
      column[ref+1] = column[i];
      ++ref;
      ++compressed_length;
      weight[ref] = 1;
    }
    else
      weight[ref]++;
  }

  /* copy the unique columns over the original sequences */
  for (i = 0; i < compressed_length; ++i)
    for (j = 0; j < count; ++j)
      sequence[j][i] = column[i][j];

  /* add terminating zero */
  for (j = 0; j < count; ++j)
    sequence[j][compressed_length] = 0;

  /* deallocate memory */
  free(memptr);
  free(column);

  /* adjust weight vector size to compressed length */
  unsigned int * mem = (unsigned int *)malloc((size_t)compressed_length *
                                              sizeof(unsigned int));
  if (mem)
  {
    /* copy weights */
    for (i = 0; i < compressed_length; ++i)
      mem[i] = weight[i];

    /* free and re-point */
    free(weight);
    weight = mem;
  }


  /* update length */
  *length = compressed_length;

  /* decode sequences using inv_charmap */
  encode(sequence,inv_charmap,count,compressed_length);

  return weight;
}
