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

extern int pll_rtree_lex();
extern FILE * pll_rtree_in;
extern void pll_rtree_lex_destroy();

static int tip_cnt = 0;


int pll_rtree_wrap() 
{ 
  return 1;
}

void pll_destroy_rtree(pll_rtree_t * root)
{
  if (!root) return;

  pll_destroy_rtree(root->left);
  pll_destroy_rtree(root->right);

  free(root->label);
  free(root);
}


static void pll_rtree_error(pll_rtree_t * tree, const char * s) 
{
  fprintf(stderr, "%s.\n", s);
}

%}


%union
{
  char * s;
  char * d;
  struct pll_rtree * tree;
}

%define api.prefix {pll_rtree_}
%error-verbose
%parse-param {struct pll_rtree * tree}
%destructor { pll_destroy_rtree($$); } subtree

%token OPAR
%token CPAR
%token COMMA
%token COLON SEMICOLON 
%token<s> STRING
%token<d> NUMBER
%type<s> label optional_label
%type<d> number optional_length
%type<tree> subtree
%start input
%%

input: OPAR subtree COMMA subtree CPAR optional_label optional_length SEMICOLON
{
  tree->left   = $2;
  tree->right  = $4;
  tree->label  = $6;
  tree->length = $7 ? atof($7) : 0;
  free($7);
};

subtree: OPAR subtree COMMA subtree CPAR optional_label optional_length
{
  $$ = (pll_rtree_t *)calloc(1, sizeof(pll_rtree_t));
  $$->left   = $2;
  $$->right  = $4;
  $$->label  = $6;
  $$->length = $7 ? atof($7) : 0;
  free($7);
}
       | label optional_length
{
  $$ = (pll_rtree_t *)calloc(1, sizeof(pll_rtree_t));
  $$->label  = $1;
  $$->length = $2 ? atof($2) : 0;
  $$->left   = NULL;
  $$->right  = NULL;
  tip_cnt++;
  free($2);
};

 
optional_label:  { $$ = NULL;} | label  {$$ = $1;};
optional_length: { $$ = NULL;} | COLON number {$$ = $2;};
label: STRING    { $$=$1;};
number: NUMBER   { $$=$1;};

%%

pll_rtree_t * pll_parse_newick_rtree(const char * filename, int * tip_count)
{
  struct pll_rtree * tree;

  tree = (pll_rtree_t *)calloc(1, sizeof(pll_rtree_t));

  pll_rtree_in = fopen(filename, "r");
  if (!pll_rtree_in)
  {
    pll_destroy_rtree(tree);
    pll_errno = PLL_ERROR_FILE_OPEN;
    snprintf(pll_errmsg, 200, "Unable to open file (%s)", filename);
    return PLL_FAILURE;
  }
  else if (pll_rtree_parse(tree))
  {
    pll_destroy_rtree(tree);
    tree = NULL;
    fclose(pll_rtree_in);
    pll_rtree_lex_destroy();
    return PLL_FAILURE;
  }
  
  if (pll_rtree_in) fclose(pll_rtree_in);

  pll_rtree_lex_destroy();

  *tip_count = tip_cnt;

  return tree;
}
