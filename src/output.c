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

static void unscale(double * prob, unsigned int times);

PLL_EXPORT void pll_show_pmatrix(pll_partition_t * partition,
                                 unsigned int index,
                                 unsigned int float_precision)
{
  unsigned int i,j,k;
  double * pmatrix;
  unsigned int states = partition->states;
  unsigned int states_padded = partition->states_padded;

  for (k = 0; k < partition->rate_cats; ++k)
  {
    pmatrix = partition->pmatrix[index] + k*states*states_padded;
    for (i = 0; i < partition->states; ++i)
    {
      for (j = 0; j < states; ++j)
        printf("%+2.*f   ", float_precision, pmatrix[i*states_padded+j]);
      printf("\n");
    }
    printf("\n");
  }
}

static void unscale(double * prob, unsigned int times)
{
  unsigned int i;

  for (i = 0; i < times; ++i)
    *prob *= PLL_SCALE_THRESHOLD;
}

PLL_EXPORT void pll_show_clv(pll_partition_t * partition,
                             unsigned int clv_index,
                             int scaler_index,
                             unsigned int float_precision)
{
  unsigned int i,j,k;

  double * clv = partition->clv[clv_index];
  unsigned int * scaler = (scaler_index == PLL_SCALE_BUFFER_NONE) ?
                          NULL : partition->scale_buffer[scaler_index];
  unsigned int states = partition->states;
  unsigned int states_padded = partition->states_padded;
  unsigned int rates = partition->rate_cats;
  double prob;

  if ((clv_index < partition->tips) &&
      (partition->attributes & PLL_ATTRIB_PATTERN_TIP))
    return;

  printf ("[ ");
  for (i = 0; i < partition->sites; ++i)
  {
    printf("{");
    for (j = 0; j < rates; ++j)
    {
      printf("(");
      for (k = 0; k < states-1; ++k)
      {
        prob = clv[i*rates*states_padded + j*states_padded + k];
        if (scaler) unscale(&prob, scaler[i]);
        printf("%.*f,", float_precision, prob);
      }
      prob = clv[i*rates*states_padded + j*states_padded + k];
      if (scaler) unscale(&prob, scaler[i]);
      printf("%.*f)", float_precision, prob);
      if (j < rates - 1) printf(",");
    }
    printf("} ");
  }
  printf ("]\n");
}
