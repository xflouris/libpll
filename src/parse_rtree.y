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
extern int pll_rtree_lineno;
extern int pll_rtree_colstart;
extern int pll_rtree_colend;

extern int pll_rtree_parse();
extern struct pll_rtree_buffer_state * pll_rtree__scan_string(char * str);
extern void pll_rtree__delete_buffer(struct pll_rtree_buffer_state * buffer);

static unsigned int tip_cnt = 0;

PLL_EXPORT void pll_rtree_destroy(pll_rtree_t * root,
                                  void (*cb_destroy)(void *))
{
  if (!root) return;

  pll_rtree_destroy(root->left, cb_destroy);
  pll_rtree_destroy(root->right, cb_destroy);

  if (root->data)
  {
    if (cb_destroy)
      cb_destroy(root->data);
  }

  free(root->label);
  free(root);
}


static void pll_rtree_error(pll_rtree_t * tree, const char * s)
{
  pll_errno = PLL_ERROR_NEWICK_SYNTAX;
  if (pll_rtree_colstart == pll_rtree_colend)
    snprintf(pll_errmsg, 200, "%s. (line %d column %d)\n",
             s, pll_rtree_lineno, pll_rtree_colstart);
  else
    snprintf(pll_errmsg, 200, "%s. (line %d column %d-%d)\n",
             s, pll_rtree_lineno, pll_rtree_colstart, pll_rtree_colend);
}

%}

%union
{
  char * s;
  char * d;
  struct pll_rtree * tree;
}

%error-verbose
%parse-param {struct pll_rtree * tree}
%destructor { pll_rtree_destroy($$,NULL); } subtree
%destructor { free($$); } STRING
%destructor { free($$); } NUMBER
%destructor { free($$); } label

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
  tree->parent = NULL;
  free($7);

  tree->left->parent  = tree;
  tree->right->parent = tree;

};

subtree: OPAR subtree COMMA subtree CPAR optional_label optional_length
{
  $$ = (pll_rtree_t *)calloc(1, sizeof(pll_rtree_t));
  $$->left   = $2;
  $$->right  = $4;
  $$->label  = $6;
  $$->length = $7 ? atof($7) : 0;
  free($7);

  $$->left->parent  = $$;
  $$->right->parent = $$;

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


optional_label:  {$$ = NULL;} | label  {$$ = $1;};
optional_length: {$$ = NULL;} | COLON number {$$ = $2;};
label: STRING    {$$=$1;} | NUMBER {$$=$1;};
number: NUMBER   {$$=$1;};

%%

static void recursive_assign_indices(pll_rtree_t * node,
                                     unsigned int * tip_clv_index,
                                     unsigned int * inner_clv_index,
                                     int * inner_scaler_index,
                                     unsigned int * inner_node_index)
{
  if (!node->left)
  {
    node->node_index = *tip_clv_index;
    node->clv_index = *tip_clv_index;
    node->pmatrix_index = *tip_clv_index;
    node->scaler_index = PLL_SCALE_BUFFER_NONE;
    *tip_clv_index = *tip_clv_index + 1;
    return;
  }

  recursive_assign_indices(node->left,
                           tip_clv_index,
                           inner_clv_index,
                           inner_scaler_index,
                           inner_node_index);

  recursive_assign_indices(node->right,
                           tip_clv_index,
                           inner_clv_index,
                           inner_scaler_index,
                           inner_node_index);

  node->node_index = *inner_node_index;
  node->clv_index = *inner_clv_index;
  node->scaler_index = *inner_scaler_index;
  node->pmatrix_index = *inner_clv_index;

  *inner_clv_index = *inner_clv_index + 1;
  *inner_scaler_index = *inner_scaler_index + 1;
  *inner_node_index = *inner_node_index + 1;
}

static void assign_indices(pll_rtree_t * root, unsigned int tip_count)
{
  unsigned int tip_clv_index = 0;
  unsigned int inner_clv_index = tip_count;
  unsigned int inner_node_index = tip_count;
  int inner_scaler_index = 0;

  recursive_assign_indices(root->left,
                           &tip_clv_index,
                           &inner_clv_index,
                           &inner_scaler_index,
                           &inner_node_index);

  recursive_assign_indices(root->right,
                           &tip_clv_index,
                           &inner_clv_index,
                           &inner_scaler_index,
                           &inner_node_index);

  root->node_index = inner_node_index;
  root->clv_index = inner_clv_index;
  root->scaler_index = inner_scaler_index;

  /* root gets any number for pmatrix since it will never be used */
  root->pmatrix_index = 0;
}

PLL_EXPORT pll_rtree_t * pll_rtree_parse_newick(const char * filename,
                                                unsigned int * tip_count)
{
  struct pll_rtree * tree;

  /* reset tip count */
  tip_cnt = 0;

  if (!(tree = (pll_rtree_t *)calloc(1, sizeof(pll_rtree_t))))
  {
    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf(pll_errmsg, 200, "Unable to allocate enough memory.");
    return PLL_FAILURE;
  }

  pll_rtree_in = fopen(filename, "r");
  if (!pll_rtree_in)
  {
    pll_rtree_destroy(tree,NULL);
    pll_errno = PLL_ERROR_FILE_OPEN;
    snprintf(pll_errmsg, 200, "Unable to open file (%s)", filename);
    return PLL_FAILURE;
  }
  else if (pll_rtree_parse(tree))
  {
    pll_rtree_destroy(tree,NULL);
    tree = NULL;
    fclose(pll_rtree_in);
    pll_rtree_lex_destroy();
    return PLL_FAILURE;
  }

  if (pll_rtree_in) fclose(pll_rtree_in);

  pll_rtree_lex_destroy();

  *tip_count = tip_cnt;

  /* initialize clv and scaler indices */
  assign_indices(tree, tip_cnt);

  return tree;
}

PLL_EXPORT pll_rtree_t * pll_rtree_parse_newick_string(char * s,
                                                       unsigned int * tip_count)
{
  int rc;
  struct pll_rtree * tree;

  /* reset tip count */
  tip_count = 0;

  if (!(tree = (pll_rtree_t *)calloc(1, sizeof(pll_rtree_t))))
  {
    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf(pll_errmsg, 200, "Unable to allocate enough memory.");
    return PLL_FAILURE;
  }

  struct pll_rtree_buffer_state * buffer = pll_rtree__scan_string(s);
  rc = pll_rtree_parse(tree);
  pll_rtree__delete_buffer(buffer);

  pll_rtree_lex_destroy();

  if (!rc)
  {
    *tip_count = tip_cnt;

    /* initialize clv and scaler indices */
    assign_indices(tree, tip_cnt);

    return tree;
  }

  free(tree);
  return NULL;
}
