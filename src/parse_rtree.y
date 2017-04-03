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
extern struct pll_rtree_buffer_state * pll_rtree__scan_string(const char * str);
extern void pll_rtree__delete_buffer(struct pll_rtree_buffer_state * buffer);

static unsigned int tip_cnt = 0;

static void dealloc_data(pll_rnode_t * node, void (*cb_destroy)(void *))
{
  if (node->data)
  {
    if (cb_destroy)
      cb_destroy(node->data);
  }
}

PLL_EXPORT void pll_rtree_graph_destroy(pll_rnode_t * root,
                                        void (*cb_destroy)(void *))
{
  if (!root) return;

  pll_rtree_graph_destroy(root->left, cb_destroy);
  pll_rtree_graph_destroy(root->right, cb_destroy);

  dealloc_data(root, cb_destroy);
  free(root->label);
  free(root);
}

PLL_EXPORT void pll_rtree_destroy(pll_rtree_t * tree,
                                  void (*cb_destroy)(void *))
{
  unsigned int i;
  pll_rnode_t * node;

  /* deallocate all nodes */
  for (i = 0; i < tree->tip_count + tree->inner_count; ++i)
  {
    node = tree->nodes[i];
    dealloc_data(node, cb_destroy);

    if (node->label)
      free(node->label);

    free(node);
  }

  /* deallocate tree structure */
  free(tree->nodes);
  free(tree);
}

static void pll_rtree_error(pll_rnode_t * node, const char * s)
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
  struct pll_rnode_s * tree;
}

%error-verbose
%parse-param {struct pll_rnode_s * tree}
%destructor { pll_rtree_graph_destroy($$,NULL); } subtree
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
  $$ = (pll_rnode_t *)calloc(1, sizeof(pll_rnode_t));
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
  $$ = (pll_rnode_t *)calloc(1, sizeof(pll_rnode_t));
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

static void recursive_assign_indices(pll_rnode_t * node,
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

PLL_EXPORT void pll_rtree_reset_template_indices(pll_rnode_t * root,
                                                 unsigned int tip_count)
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

static void fill_nodes_recursive(pll_rnode_t * node,
                                 pll_rnode_t ** array,
                                 unsigned int * tip_index,
                                 unsigned int * inner_index)
{
  if (!node->left)
  {
    array[*tip_index] = node;
    *tip_index = *tip_index + 1;
    return;
  }

  fill_nodes_recursive(node->left,  array, tip_index, inner_index);
  fill_nodes_recursive(node->right, array, tip_index, inner_index);

  array[*inner_index] = node;
  *inner_index = *inner_index + 1;
}

static unsigned int rtree_count_tips(pll_rnode_t * root)
{
  unsigned int count = 0;

  if (root->left)
    count += rtree_count_tips(root->left);
  if (root->right)
    count += rtree_count_tips(root->right);

  if (!root->left && !root->right)
    return 1;

  return count;
}

PLL_EXPORT pll_rtree_t * pll_rtree_wraptree(pll_rnode_t * root,
                                            unsigned int tip_count)
{
  pll_rtree_t * tree = (pll_rtree_t *)malloc(sizeof(pll_rtree_t));
  if (!tree)
  {
    snprintf(pll_errmsg, 200, "Unable to allocate enough memory.");
    pll_errno = PLL_ERROR_MEM_ALLOC;
    return PLL_FAILURE;
  }

  if (tip_count < 2 && tip_count != 0)
  {
    snprintf(pll_errmsg, 200, "Invalid tip_count value (%u).", tip_count);
    pll_errno = PLL_ERROR_PARAM_INVALID;
    return PLL_FAILURE;
  }

  if (tip_count == 0)
  {
    tip_count = rtree_count_tips(root);
    if (tip_count < 2)
    {
      snprintf(pll_errmsg, 200, "Input tree contains no inner nodes.");
      pll_errno = PLL_ERROR_PARAM_INVALID;
      return PLL_FAILURE;
    }
  }

  tree->nodes = (pll_rnode_t **)malloc((2*tip_count-1)*sizeof(pll_rnode_t *));
  if (!tree->nodes)
  {
    snprintf(pll_errmsg, 200, "Unable to allocate enough memory.");
    pll_errno = PLL_ERROR_MEM_ALLOC;
    return PLL_FAILURE;
  }
  
  unsigned int tip_index = 0;
  unsigned int inner_index = tip_count;

  fill_nodes_recursive(root->left, tree->nodes, &tip_index, &inner_index);
  fill_nodes_recursive(root->right,tree->nodes, &tip_index, &inner_index);
  tree->nodes[inner_index] = root;

  tree->tip_count = tip_count;
  tree->edge_count = 2*tip_count-2;
  tree->inner_count = tip_count-1;
  tree->root = root;

  return tree;
}

PLL_EXPORT pll_rtree_t * pll_rtree_parse_newick(const char * filename)
{
  pll_rtree_t * tree;

  struct pll_rnode_s * root;

  /* reset tip count */
  tip_cnt = 0;

  /* open input file */
  pll_rtree_in = fopen(filename, "r");
  if (!pll_rtree_in)
  {
    pll_errno = PLL_ERROR_FILE_OPEN;
    snprintf(pll_errmsg, 200, "Unable to open file (%s)", filename);
    return PLL_FAILURE;
  }

  /* create root node */
  if (!(root = (pll_rnode_t *)calloc(1, sizeof(pll_rnode_t))))
  {
    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf(pll_errmsg, 200, "Unable to allocate enough memory.");
    return PLL_FAILURE;
  }

  if (pll_rtree_parse(root))
  {
    pll_rtree_graph_destroy(root,NULL);
    root = NULL;
    fclose(pll_rtree_in);
    pll_rtree_lex_destroy();
    return PLL_FAILURE;
  }

  if (pll_rtree_in) fclose(pll_rtree_in);

  pll_rtree_lex_destroy();

  /* initialize clv and scaler indices */
  pll_rtree_reset_template_indices(root, tip_cnt);

  /* wrap tree */
  tree = pll_rtree_wraptree(root,tip_cnt);

  return tree;
}

PLL_EXPORT pll_rtree_t * pll_rtree_parse_newick_string(const char * s)
{
  int rc;
  struct pll_rnode_s * root;
  pll_rtree_t * tree = NULL;

  /* reset tip count */
  tip_cnt = 0;

  if (!(root = (pll_rnode_t *)calloc(1, sizeof(pll_rnode_t))))
  {
    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf(pll_errmsg, 200, "Unable to allocate enough memory.");
    return PLL_FAILURE;
  }

  struct pll_rtree_buffer_state * buffer = pll_rtree__scan_string(s);
  rc = pll_rtree_parse(root);
  pll_rtree__delete_buffer(buffer);

  pll_rtree_lex_destroy();

  if (!rc)
  {
    /* initialize clv and scaler indices */
    pll_rtree_reset_template_indices(root, tip_cnt);

    tree = pll_rtree_wraptree(root,tip_cnt);
  }
  else
    free(root);

  return tree;
}
