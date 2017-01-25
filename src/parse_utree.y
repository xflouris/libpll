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
extern int pll_utree_lineno;
extern int pll_utree_colstart;
extern int pll_utree_colend;

extern int pll_utree_parse();
extern struct pll_utree_buffer_state * pll_utree__scan_string(char * str);
extern void pll_utree__delete_buffer(struct pll_utree_buffer_state * buffer);

static unsigned int tip_cnt = 0;

static void dealloc_data(pll_utree_t * node, void (*cb_destroy)(void *))
{
  if (node->data)
  {
    if (cb_destroy)
      cb_destroy(node->data);
  }
}

static void dealloc_tree_recursive(pll_utree_t * node,
                                   void (*cb_destroy)(void *))
{
  if (!node->next)
  {
    free(node->label);
    free(node);
    return;
  }

  dealloc_tree_recursive(node->next->back, cb_destroy);
  dealloc_tree_recursive(node->next->next->back, cb_destroy);

  /* deallocate any extra data */
  dealloc_data(node,cb_destroy);
  dealloc_data(node->next,cb_destroy);
  dealloc_data(node->next->next,cb_destroy);

  free(node->next->next);
  free(node->next);
  free(node->label);
  free(node);
}

PLL_EXPORT void pll_utree_destroy(pll_utree_t * root,
                                  void (*cb_destroy)(void *))
{
  if (!root) return;
  if (!(root->next))
  {
    free(root->label);
    free(root);
    return;
  }

  if (root->next)
    dealloc_tree_recursive(root->next->back,cb_destroy);
  if (root->next->next)
    dealloc_tree_recursive(root->next->next->back,cb_destroy);
  if (root->back)
    dealloc_tree_recursive(root->back,cb_destroy);

  /* deallocate any extra data */
  dealloc_data(root,cb_destroy);
  dealloc_data(root->next,cb_destroy);
  dealloc_data(root->next->next,cb_destroy);

  free(root->label);

  free(root->next->next);
  free(root->next);
  free(root);
}

static void pll_utree_error(pll_utree_t * tree, const char * s)
{
  pll_errno = PLL_ERROR_NEWICK_SYNTAX;
  if (pll_utree_colstart == pll_utree_colend)
    snprintf(pll_errmsg, 200, "%s. (line %d column %d)\n",
             s, pll_utree_lineno, pll_utree_colstart);
  else
    snprintf(pll_errmsg, 200, "%s. (line %d column %d-%d)\n",
             s, pll_utree_lineno, pll_utree_colstart, pll_utree_colend);
}

%}

%union
{
  char * s;
  char * d;
  struct pll_utree * tree;
}

%error-verbose
%parse-param {struct pll_utree * tree}
%destructor { pll_utree_destroy($$,NULL); } subtree
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

input: OPAR subtree COMMA subtree COMMA subtree CPAR optional_label optional_length SEMICOLON
{
  tree->next               = (pll_utree_t *)calloc(1, sizeof(pll_utree_t));

  tree->next->next         = (pll_utree_t *)calloc(1, sizeof(pll_utree_t));
  tree->next->next->next   = tree;


  tree->back               = $2;
  tree->next->back         = $4;
  tree->next->next->back   = $6;

  $2->back                 = tree;
  $4->back                 = tree->next;
  $6->back                 = tree->next->next;

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
  $$->next->next->next   = $$;


  $$->next->back         = $2;
  $$->next->next->back   = $4;

  $2->back               = $$->next;
  $4->back               = $$->next->next;

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

static void recursive_assign_indices(pll_utree_t * node,
                                    unsigned int * tip_clv_index,
                                    unsigned int * inner_clv_index,
                                    int * inner_scaler_index,
                                    unsigned int * inner_node_index)
{
  if (!node->next)
  {
    node->node_index = *tip_clv_index;
    node->clv_index = *tip_clv_index;
    node->pmatrix_index = *tip_clv_index;
    node->scaler_index = PLL_SCALE_BUFFER_NONE;
    *tip_clv_index = *tip_clv_index + 1;
    return;
  }

  recursive_assign_indices(node->next->back,
                           tip_clv_index,
                           inner_clv_index,
                           inner_scaler_index,
                           inner_node_index);

  recursive_assign_indices(node->next->next->back,
                           tip_clv_index,
                           inner_clv_index,
                           inner_scaler_index,
                           inner_node_index);

  node->node_index = *inner_node_index;
  node->next->node_index = *inner_node_index + 1;
  node->next->next->node_index = *inner_node_index + 2;

  node->clv_index = *inner_clv_index;
  node->next->clv_index = *inner_clv_index;
  node->next->next->clv_index = *inner_clv_index;

  node->pmatrix_index = *inner_clv_index;
  node->next->pmatrix_index = node->next->back->pmatrix_index;
  node->next->next->pmatrix_index = node->next->next->back->pmatrix_index;

  node->scaler_index = *inner_scaler_index;
  node->next->scaler_index = *inner_scaler_index;
  node->next->next->scaler_index = *inner_scaler_index;

  *inner_clv_index = *inner_clv_index + 1;
  *inner_scaler_index = *inner_scaler_index + 1;
  *inner_node_index = *inner_node_index + 3;
}

PLL_EXPORT void pll_utree_reset_template_indices(pll_utree_t * node,
                                                 unsigned int tip_count)
{
  unsigned int tip_clv_index = 0;
  unsigned int inner_clv_index = tip_count;
  unsigned int inner_node_index = tip_count;
  int inner_scaler_index = 0;

  recursive_assign_indices(node->back,
                           &tip_clv_index,
                           &inner_clv_index,
                           &inner_scaler_index,
                           &inner_node_index);

  recursive_assign_indices(node->next->back,
                           &tip_clv_index,
                           &inner_clv_index,
                           &inner_scaler_index,
                           &inner_node_index);

  recursive_assign_indices(node->next->next->back,
                           &tip_clv_index,
                           &inner_clv_index,
                           &inner_scaler_index,
                           &inner_node_index);

  node->node_index = inner_node_index;
  node->next->node_index = inner_node_index + 1;
  node->next->next->node_index = inner_node_index + 2;

  node->clv_index = inner_clv_index;
  node->next->clv_index = inner_clv_index;
  node->next->next->clv_index = inner_clv_index;

  node->scaler_index = inner_scaler_index;
  node->next->scaler_index = inner_scaler_index;
  node->next->next->scaler_index = inner_scaler_index;

  node->pmatrix_index = node->back->pmatrix_index;
  node->next->pmatrix_index = node->next->back->pmatrix_index;
  node->next->next->pmatrix_index = node->next->next->back->pmatrix_index;
}

PLL_EXPORT pll_utree_t * pll_utree_parse_newick(const char * filename,
                                                unsigned int * tip_count)
{
  struct pll_utree * tree;

  /* reset tip count */
  tip_cnt = 0;

  tree = (pll_utree_t *)calloc(1, sizeof(pll_utree_t));

  pll_utree_in = fopen(filename, "r");
  if (!pll_utree_in)
  {
    pll_utree_destroy(tree,NULL);
    pll_errno = PLL_ERROR_FILE_OPEN;
    snprintf(pll_errmsg, 200, "Unable to open file (%s)", filename);
    return PLL_FAILURE;
  }
  else if (pll_utree_parse(tree))
  {
    pll_utree_destroy(tree,NULL);
    tree = NULL;
    fclose(pll_utree_in);
    pll_utree_lex_destroy();
    return PLL_FAILURE;
  }

  if (pll_utree_in) fclose(pll_utree_in);

  pll_utree_lex_destroy();

  *tip_count = tip_cnt;

  /* initialize clv and scaler indices */
  pll_utree_reset_template_indices(tree, tip_cnt);

  return tree;
}

PLL_EXPORT pll_utree_t * pll_utree_parse_newick_string(char * s,
                                                       unsigned int * tip_count)
{
  int rc;
  struct pll_utree * tree;

  /* reset tip count */
  tip_count = 0;

  if (!(tree = (pll_utree_t *)calloc(1, sizeof(pll_utree_t))))
  {
    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf(pll_errmsg, 200, "Unable to allocate enough memory.");
    return PLL_FAILURE;
  }

  struct pll_utree_buffer_state * buffer = pll_utree__scan_string(s);
  rc = pll_utree_parse(tree);
  pll_utree__delete_buffer(buffer);

  pll_utree_lex_destroy();

  if (!rc)
  {
    *tip_count = tip_cnt;

    /* initialize clv and scaler indices */
    pll_utree_reset_template_indices(tree, tip_cnt);

    return tree;
  }

  free(tree);
  return NULL;
}
