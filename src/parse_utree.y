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
extern struct pll_utree_buffer_state * pll_utree__scan_string(const char * str);
extern void pll_utree__delete_buffer(struct pll_utree_buffer_state * buffer);

static unsigned int tip_cnt = 0;

static void dealloc_data(pll_unode_t * node, void (*cb_destroy)(void *))
{
  if (node->data)
  {
    if (cb_destroy)
      cb_destroy(node->data);
  }
}

static void dealloc_graph_recursive(pll_unode_t * node,
                                   void (*cb_destroy)(void *))
{
  if (!node->next)
  {
    dealloc_data(node,cb_destroy);
    free(node->label);
    free(node);
    return;
  }

  dealloc_graph_recursive(node->next->back, cb_destroy);
  dealloc_graph_recursive(node->next->next->back, cb_destroy);

  /* deallocate any extra data */
  dealloc_data(node,cb_destroy);
  dealloc_data(node->next,cb_destroy);
  dealloc_data(node->next->next,cb_destroy);

  free(node->next->next);
  free(node->next);
  free(node->label);
  free(node);
}

PLL_EXPORT void pll_utree_graph_destroy(pll_unode_t * root,
                                        void (*cb_destroy)(void *))
{
  if (!root) return;
  if (!(root->next))
  {
    dealloc_data(root,cb_destroy);
    free(root->label);
    free(root);
    return;
  }

  if (root->next)
    dealloc_graph_recursive(root->next->back,cb_destroy);
  if (root->next->next)
    dealloc_graph_recursive(root->next->next->back,cb_destroy);
  if (root->back)
    dealloc_graph_recursive(root->back,cb_destroy);

  /* deallocate any extra data */
  dealloc_data(root,cb_destroy);
  dealloc_data(root->next,cb_destroy);
  dealloc_data(root->next->next,cb_destroy);

  free(root->label);

  free(root->next->next);
  free(root->next);
  free(root);
}

PLL_EXPORT void pll_utree_destroy(pll_utree_t * tree,
                                  void (*cb_destroy)(void *))
{
  unsigned int i;
  pll_unode_t * node;

  /* deallocate tip nodes */
  for (i = 0; i < tree->tip_count; ++i)
  {
    dealloc_data(tree->nodes[i], cb_destroy);
    if (tree->nodes[i]->label)
      free(tree->nodes[i]->label);
    free(tree->nodes[i]);
  }

  /* deallocate inner nodes */
  for (i = tree->tip_count; i < tree->tip_count + tree->inner_count; ++i)
  {
    node = tree->nodes[i];

    dealloc_data(node, cb_destroy);
    dealloc_data(node->next, cb_destroy);
    dealloc_data(node->next->next, cb_destroy);

    if (node->label)
      free(node->label);

    free(node->next->next);
    free(node->next);
    free(node);
  }

  /* deallocate tree structure */
  free(tree->nodes);
  free(tree);
}

static void pll_utree_error(pll_unode_t * node, const char * s)
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
  struct pll_unode_s * tree;
}

%error-verbose
%parse-param {struct pll_unode_s * tree}
%destructor { pll_utree_graph_destroy($$,NULL); } subtree
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
  tree->next               = (pll_unode_t *)calloc(1, sizeof(pll_unode_t));

  tree->next->next         = (pll_unode_t *)calloc(1, sizeof(pll_unode_t));
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
  $$                     = (pll_unode_t *)calloc(1, sizeof(pll_unode_t));

  $$->next               = (pll_unode_t *)calloc(1, sizeof(pll_unode_t));

  $$->next->next         = (pll_unode_t *)calloc(1, sizeof(pll_unode_t));
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
  $$ = (pll_unode_t *)calloc(1, sizeof(pll_unode_t));

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

static void recursive_assign_indices(pll_unode_t * node,
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

PLL_EXPORT void pll_utree_reset_template_indices(pll_unode_t * node,
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

static void fill_nodes_recursive(pll_unode_t * node,
                                 pll_unode_t ** array,
                                 unsigned int * tip_index,
                                 unsigned int * inner_index)
{
  if (!node->next)
  {
    array[*tip_index] = node;
    *tip_index = *tip_index + 1;
    return;
  }

  fill_nodes_recursive(node->next->back,       array, tip_index, inner_index);
  fill_nodes_recursive(node->next->next->back, array, tip_index, inner_index);

  array[*inner_index] = node;
  *inner_index = *inner_index + 1;
}

static unsigned int utree_count_tips_recursive(pll_unode_t * node)
{
  unsigned int count = 0;

  if (!node->next)
    return 1;

  count += utree_count_tips_recursive(node->next->back);
  count += utree_count_tips_recursive(node->next->next->back);

  return count;
}

static unsigned int utree_count_tips(pll_unode_t * root)
{
  unsigned int count = 0;

  if (!root->next && !root->back->next)
    return 0;

  if (!root->next)
    root = root->back;

  count += utree_count_tips_recursive(root->back);
  count += utree_count_tips_recursive(root->next->back);
  count += utree_count_tips_recursive(root->next->next->back);

  return count;
}

/* wraps/encalupsates the unrooted tree graph into a tree structure
   that contains a list of nodes, number of tips and number of inner
   nodes. If 0 is passed as tip_count, then an additional recrursion
   of the tree structure is done to detect the number of tips */
PLL_EXPORT pll_utree_t * pll_utree_wraptree(pll_unode_t * root,
                                            unsigned int tip_count)
{
  pll_utree_t * tree = (pll_utree_t *)malloc(sizeof(pll_utree_t));
  if (!tree)
  {
    snprintf(pll_errmsg, 200, "Unable to allocate enough memory.");
    pll_errno = PLL_ERROR_MEM_ALLOC;
    return PLL_FAILURE;
  }

  if (tip_count < 3 && tip_count != 0)
  {
    snprintf(pll_errmsg, 200, "Invalid tip_count value (%u).", tip_count);
    pll_errno = PLL_ERROR_PARAM_INVALID;
    return PLL_FAILURE;
  }

  if (tip_count == 0)
  {
    tip_count = utree_count_tips(root);
    if (!tip_count)
    {
      snprintf(pll_errmsg, 200, "Input tree contains no inner nodes.");
      pll_errno = PLL_ERROR_PARAM_INVALID;
      return PLL_FAILURE;
    }
  }

  tree->nodes = (pll_unode_t **)malloc((2*tip_count-2)*sizeof(pll_unode_t *));
  if (!tree->nodes)
  {
    snprintf(pll_errmsg, 200, "Unable to allocate enough memory.");
    pll_errno = PLL_ERROR_MEM_ALLOC;
    return PLL_FAILURE;
  }
  
  unsigned int tip_index = 0;
  unsigned int inner_index = tip_count;

  fill_nodes_recursive(root->back,             tree->nodes, &tip_index, &inner_index);
  fill_nodes_recursive(root->next->back,       tree->nodes, &tip_index, &inner_index);
  fill_nodes_recursive(root->next->next->back, tree->nodes, &tip_index, &inner_index);

  tree->nodes[inner_index] = root;
  tree->tip_count = tip_count;
  tree->edge_count = 2*tip_count-3;
  tree->inner_count = tip_count-2;

  return tree;
}

PLL_EXPORT pll_utree_t * pll_utree_parse_newick(const char * filename)
{
  pll_utree_t * tree;

  struct pll_unode_s * root;

  /* reset tip count */
  tip_cnt = 0;

  pll_utree_in = fopen(filename, "r");
  if (!pll_utree_in)
  {
    pll_errno = PLL_ERROR_FILE_OPEN;
    snprintf(pll_errmsg, 200, "Unable to open file (%s)", filename);
    return PLL_FAILURE;
  }

  if (!(root = (pll_unode_t *)calloc(1, sizeof(pll_unode_t))))
  {
    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf(pll_errmsg, 200, "Unable to allocate enough memory.");
    return PLL_FAILURE;
  }

  if (pll_utree_parse(root))
  {
    pll_utree_graph_destroy(root,NULL);
    root = NULL;
    fclose(pll_utree_in);
    pll_utree_lex_destroy();
    return PLL_FAILURE;
  }

  if (pll_utree_in) fclose(pll_utree_in);

  pll_utree_lex_destroy();

  /* initialize clv and scaler indices to the default template */
  pll_utree_reset_template_indices(root, tip_cnt);

  /* wrap tree */
  tree = pll_utree_wraptree(root,tip_cnt);

  return tree;
}

PLL_EXPORT pll_utree_t * pll_utree_parse_newick_string(const char * s)
{
  int rc;
  struct pll_unode_s * root;
  pll_utree_t * tree = NULL;

  /* reset tip count */
  tip_cnt = 0;

  if (!(root = (pll_unode_t *)calloc(1, sizeof(pll_unode_t))))
  {
    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf(pll_errmsg, 200, "Unable to allocate enough memory.");
    return PLL_FAILURE;
  }

  struct pll_utree_buffer_state * buffer = pll_utree__scan_string(s);
  rc = pll_utree_parse(root);
  pll_utree__delete_buffer(buffer);

  pll_utree_lex_destroy();

  if (!rc)
  {
    /* initialize clv and scaler indices */
    pll_utree_reset_template_indices(root, tip_cnt);

    tree = pll_utree_wraptree(root,tip_cnt);
  }
  else
    free(root);

  return tree;
}
