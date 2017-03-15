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
    Heidelberg Institute for Theoretical Studies,
    Schloss-Wolfsbrunnenweg 35, D-69118 Heidelberg, Germany
*/
%{
#include "pll.h"

extern int pll_phylip_lex();
extern FILE * pll_phylip_in;
extern void pll_phylip_lex_destroy();
extern int pll_phylip_lineno;

static int seq_count = 0;
static int g_error = 0;

struct sequence_s
{
  int len;
  char * label;
  char * sequence;
};

void pll_phylip_error(pll_msa_t ** msa, const char * s)
{
  fprintf(stderr, "%s.\n", s);
}

PLL_EXPORT void pll_msa_destroy(pll_msa_t * msa)
{
  if (!msa) return;

  int i;

  if (msa->label)
  {
    for (i = 0; i < msa->count; ++i)
      free(msa->label[i]);
    free(msa->label);
  }
  
  if (msa->sequence)
  {
    for (i = 0; i < msa->count; ++i)
      free(msa->sequence[i]);
    free(msa->sequence);
  }

  free(msa);
}

void msa_destroy(pll_msa_t * msa)
{
  if (msa)
  {
    msa->count = seq_count;
    pll_msa_destroy(msa);
  }
}
%}

%union
{
  char * s;
  pll_msa_t * msa;
  struct sequence_s ** sequence;
}

%error-verbose
%parse-param {pll_msa_t ** msa}
%destructor { msa_destroy($$); } msa

%token<s> STRING
%type<msa> msa
%type<sequence> sequence_list
%start input
%%

input: msa
{
  *msa = $1;
}

msa: STRING STRING sequence_list
{
  int i;

  if (!g_error)
    for (i = 0 ; $1[i]; ++i)
      if (!isdigit($1[i]))
      {
        pll_errno = PLL_ERROR_PHYLIP_SYNTAX;
        snprintf(pll_errmsg, 200, "Number of taxa must be an integer");
        g_error = 1;
        break;
      }

  if (!g_error)
    for (i = 0 ; $2[i]; ++i)
      if (!isdigit($2[i]))
      {
        pll_errno = PLL_ERROR_PHYLIP_SYNTAX;
        snprintf(pll_errmsg, 200, "Sequence length must be an integer");
        g_error = 1;
        break;
      }

  if (!g_error)
    if (atoi($1) != seq_count)
    {
      pll_errno = PLL_ERROR_PHYLIP_SYNTAX;
      snprintf(pll_errmsg, 200,
              "Number of sequences read is %d but expected %d\n",
              seq_count,
              atoi($1));
      g_error = 1;
    }

  if (!g_error)
    if (atoi($2) != $3[seq_count-1]->len)
    {
      pll_errno = PLL_ERROR_PHYLIP_SYNTAX;
      snprintf(pll_errmsg, 200,
              "Sequence length is %d but expected %d\n",
              $3[seq_count-1]->len,
              atoi($2));
      g_error = 1;
    }

  $$ = NULL;

  if (!g_error)
  {
    $$ = (pll_msa_t *)calloc(1,sizeof(pll_msa_t));

    $$->count = atoi($1);
    $$->length = atoi($2);
    $$->sequence = (char **)malloc(seq_count*sizeof(char *));
    $$->label = (char **)malloc(seq_count*sizeof(char *));

    for (i = 0; i < seq_count; ++i)
    {
      $$->sequence[i] = $3[seq_count-i-1]->sequence;
      $$->label[i]    = $3[seq_count-i-1]->label;
      free($3[seq_count-i-1]);
    }
  }
  else
  {
    for (i = 0; i < seq_count; ++i)
    {
      free($3[seq_count-i-1]->sequence);
      free($3[seq_count-i-1]->label);
      free($3[seq_count-i-1]);
    }
    seq_count = 0;
  }

  free($1);
  free($2);
  free($3);
}

sequence_list: STRING STRING sequence_list
{
  $$ = (struct sequence_s **)realloc((void *)$3,
                                     (seq_count+1)*sizeof(struct sequence_s *));
  $$[seq_count] = (struct sequence_s *)malloc(sizeof(struct sequence_s));
  $$[seq_count]->label = $1;
  $$[seq_count]->sequence = $2;
  $$[seq_count]->len = strlen($2);

  if ($$[seq_count]->len != $$[seq_count-1]->len)
  {
    pll_errno = PLL_ERROR_PHYLIP_SYNTAX;
    snprintf(pll_errmsg, 200, "Lengths of sequences %s and %s differ",
             $$[seq_count]->label, $$[seq_count-1]->label);
    g_error = 1;

  }
  seq_count++;

}
        | STRING STRING
{
  $$ = (struct sequence_s **)calloc(1,sizeof(struct sequence_s *));
  $$[seq_count] = (struct sequence_s *)malloc(sizeof(struct sequence_s));
  $$[seq_count]->label = $1;
  $$[seq_count]->sequence = $2;
  $$[seq_count]->len = strlen($2);
  seq_count++;
}

%%

PLL_EXPORT pll_msa_t * pll_phylip_parse_msa(const char * filename,
                                            unsigned int * msa_count)
{
  pll_msa_t * msa = NULL;

  pll_phylip_lineno = 1;
  seq_count= 0;
  g_error = 0;

  pll_phylip_in = fopen(filename, "r");
  if (!pll_phylip_in)
  {
    msa_destroy(msa);
    pll_errno = PLL_ERROR_FILE_OPEN;
    snprintf(pll_errmsg, 200, "Unable to open file (%s)", filename);
    return PLL_FAILURE;
  }
  else if (pll_phylip_parse(&msa))
  {
    msa_destroy(msa);
    msa = NULL;
    fclose(pll_phylip_in);
    pll_phylip_lex_destroy();
    return PLL_FAILURE;
  }

  if (pll_phylip_in) fclose(pll_phylip_in);

  pll_phylip_lex_destroy();

  *msa_count = seq_count;

  return msa;
}
