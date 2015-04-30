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

void pll_show_pmatrix(pll_partition_t * partition, int index)
{
  int i,j,k;
  double * pmatrix;
  int states = states; 

  for (k = 0; k < partition->rate_cats; ++k)
  {
    pmatrix = partition->pmatrix[index] + k*states*states;
    for (i = 0; i < partition->states; ++i)
    {
      for (j = 0; j < states; ++j)
        printf("%+2.5f   ", pmatrix[i*states+j]);
      printf("\n");
    }
    printf("\n");
  }
}

void pll_show_clv(pll_partition_t * partition, int index)
{
  int i,j,k;

  double * clv = partition->clv[index];
  int states = partition->states;
  int rates = partition->rate_cats;

  printf ("[ ");
  for (i = 0; i < partition->sites; ++i)
  {
    printf("{");
    for (j = 0; j < partition->rate_cats; ++j)
    {
      printf("(");
      for (k = 0; k < partition->states-1; ++k)
      {
        printf("%f,", clv[i*rates*states + j*states + k]);
      }
      printf("%f)", clv[i*rates*states + j*states + k]);
      if (j < rates - 1) printf(",");
    }
    printf("} ");
  }
  printf ("]\n");
}
