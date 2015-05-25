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

extern int pll_utree_lex();
extern FILE * pll_utree_in;
extern void pll_utree_lex_destroy();

static int tip_cnt = 0;


int pll_utree_wrap()
{ 
  return 1;
}

static void dealloc_tree_recursive(pll_utree_t * node)
{
  if (!node->next)
  {
    free(node->label);
    free(node);
    return;
  }

  dealloc_tree_recursive(node->next->back);
  dealloc_tree_recursive(node->next->next->back);

  free(node->next->next);
  free(node->next);
  free(node->label);
  free(node);
}

void pll_destroy_utree(pll_utree_t * root)
{
  if (!root) return;
  if (!(root->next))
  {
    free(root->label);
    free(root);
    return;
  }

  if (root->next)
    dealloc_tree_recursive(root->next->back);
  if (root->next->next)
    dealloc_tree_recursive(root->next->next->back);
  if (root->back)
    dealloc_tree_recursive(root->back);

  free(root->label);
  free(root->next->next);
  free(root->next);
  free(root);
}

static void pll_utree_error(pll_utree_t * tree, const char * s) 
{
  fprintf(stderr, "%s.\n", s);
}

%}

%union
{
  char * s;
  char * d;
  struct pll_utree * tree;
}

%define api.prefix {pll_utree_}
%error-verbose
%parse-param {struct pll_utree * tree}
%destructor { pll_destroy_utree($$); } subtree

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

input: OPAR subtree COMMA subtree COMMA subtree CPAR optional_label optional_length SEMICOLON
{
  tree->next               = (pll_utree_t *)calloc(1, sizeof(pll_utree_t));

  tree->next->next         = (pll_utree_t *)calloc(1, sizeof(pll_utree_t));


  tree->back               = $2;
  tree->next->back         = $4;
  tree->next->next->back   = $6;

  tree->label              = $8;
  tree->next->label        = $8;
  tree->next->next->label  = $8;

  tree->length             = $2->length;
  tree->next->length       = $4->length;
  tree->next->next->length = $6->length;
  free($9);
};

subtree: OPAR subtree COMMA subtree CPAR optional_label optional_length
{
  $$                     = (pll_utree_t *)calloc(1, sizeof(pll_utree_t));

  $$->next               = (pll_utree_t *)calloc(1, sizeof(pll_utree_t));

  $$->next->next         = (pll_utree_t *)calloc(1, sizeof(pll_utree_t));


  $$->next->back         = $2;
  $$->next->next->back   = $4;

  $$->label              = $6;
  $$->next->label        = $6;
  $$->next->next->label  = $6;

  $$->length = $7 ? atof($7) : 0;
  free($7);

  $$->next->length       = $2->length;
  $$->next->next->length = $4->length;
}
       | label optional_length
{
  $$ = (pll_utree_t *)calloc(1, sizeof(pll_utree_t));

  $$->label  = $1;
  $$->length = $2 ? atof($2) : 0;
  $$->next   = NULL;
  tip_cnt++;
  free($2);
};

 
optional_label:  { $$ = NULL;} | label  {$$ = $1;};
optional_length: { $$ = NULL;} | COLON number {$$ = $2;};
label: STRING    { $$=$1;} | NUMBER {$$=$1;};
number: NUMBER   { $$=$1;};

%%

pll_utree_t * pll_parse_newick_utree(const char * filename, int * tip_count)
{
  struct pll_utree * tree;

  tree = (pll_utree_t *)calloc(1, sizeof(pll_utree_t));

  pll_utree_in = fopen(filename, "r");
  if (!pll_utree_in)
  {
    pll_destroy_utree(tree);
    pll_errno = PLL_ERROR_FILE_OPEN;
    snprintf(pll_errmsg, 200, "Unable to open file (%s)", filename);
    return PLL_FAILURE;
  }
  else if (pll_utree_parse(tree))
  {
    pll_destroy_utree(tree);
    tree = NULL;
    fclose(pll_utree_in);
    pll_utree_lex_destroy();
    return PLL_FAILURE;
  }
  
  if (pll_utree_in) fclose(pll_utree_in);

  pll_utree_lex_destroy();

  *tip_count = tip_cnt;

  return tree;
}
