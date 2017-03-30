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

#include "pll.h"


typedef struct pll_svg_data_s 
{
  int height;
  double x;
  double y;
} pll_svg_data_t; 

typedef struct pll_svg_aux_s
{
  int tip_occ;
  double scaler;
  double canvas_width;
  double max_font_len;
  double max_tree_len;
} pll_svg_aux_t;

static pll_svg_data_t * create_data(int height, double x, double y)
{
  pll_svg_data_t * data = (pll_svg_data_t *)malloc(sizeof(pll_svg_data_t));
  if (!data) return NULL;

  data->height = height;
  data->x = x;
  data->y = y;

  return data;
}

static int utree_height_recursive(pll_unode_t * node)
{
  if (!node->next)
  {
    node->data = (void *)create_data(0,0,0);
    if (!node->data) return 0;
    return 1;
  }

  if (!utree_height_recursive(node->next->back)) return 0;
  if (!utree_height_recursive(node->next->next->back)) return 0;

  pll_svg_data_t * d1 = (pll_svg_data_t *)(node->next->back->data);
  pll_svg_data_t * d2 = (pll_svg_data_t *)(node->next->next->back->data);
  pll_svg_data_t * d  = create_data(0,0,0);
  if (!d) return 0;

  if (d1->height > d2->height)
    d->height = d1->height+1;
  else
    d->height = d2->height+1;

  node->data = node->next->data = node->next->next->data = d;

  return 1;
}

static int utree_set_height(pll_unode_t * root)
{
  if (!root->next) return PLL_FAILURE;

  if (!utree_height_recursive(root->back)) return PLL_FAILURE;
  if (!utree_height_recursive(root)) return PLL_FAILURE;

  pll_svg_data_t * db = (pll_svg_data_t *)(root->back->data);
  pll_svg_data_t * d = (pll_svg_data_t *)(root->data);
  
  if (db->height >= d->height)
    d->height = db->height+1;
    

  return PLL_SUCCESS;
}

static void draw_line(FILE * fp,
                      double x1, double y1,
                      double x2, double y2,
                      double stroke_width)
{
  fprintf(fp,
          "<line x1=\"%f\" y1=\"%f\" x2=\"%f\" y2=\"%f\" "
          "stroke=\"#31a354\" stroke-width=\"%f\" />\n",
          x1, y1, x2, y2, stroke_width);
}

static void draw_circle(FILE * fp, double cx, double cy, double r)
{
  fprintf(fp,
          "<circle cx=\"%f\" cy=\"%f\" r=\"%f\" fill=\"#31a354\" "
          "stroke=\"#31a354\" />\n",
          cx, cy, r);
}

static void utree_set_offset(pll_unode_t * node,
                             const pll_svg_attrib_t * attr,
                             const pll_svg_aux_t * aux)
{
  pll_unode_t * parent = NULL;

  /* scale node's branch length (edge towards parent) */
  pll_svg_data_t * data = (pll_svg_data_t *)(node->data);
  data->x = node->length * aux->scaler;
  
  /* did we reach the root ? */
  pll_svg_data_t * parent_data = (pll_svg_data_t *)(node->back->data);
  if (parent_data->height > data->height)
    parent = node->back;

  /* if node has a parent then add up the parent's x coord such that
     the branch is shifted towards right. Otherwise, if node is the root,
     align it with the left margin */
  if (parent)
    data->x += parent_data->x;
  else
    data->x= attr->margin_left;

  /* if it's a tip then stop here */
  if (!node->next) 
    return;
  
  /* otherwise recursively set coordinates for the other nodes in a
     pre-order fashion */
  utree_set_offset(node->next->back, attr, aux);
  utree_set_offset(node->next->next->back, attr, aux);
  if (!parent)
    utree_set_offset(node->back,attr, aux);
}

static void utree_plot(FILE * fp,
                       pll_unode_t * node,
                       const pll_svg_attrib_t * attr,
                       pll_svg_aux_t * aux)
{
  double y;
//  static int tip_occ = 0;
  pll_unode_t * parent = NULL;

  pll_svg_data_t * data = (pll_svg_data_t *)(node->data);
  pll_svg_data_t * parent_data = (pll_svg_data_t *)(node->back->data);

  if (parent_data->height > data->height)
    parent = node->back;

  if (node->next)
  {
    utree_plot(fp, node->next->back, attr, aux);
    utree_plot(fp, node->next->next->back, attr, aux);
    if (!parent)
      utree_plot(fp, node->back, attr, aux);
  }

  if (parent)
  {
    double x,px;

    x = data->x;
    px = parent_data->x;

    if (!node->next)
    {
      y = aux->tip_occ * attr->tip_spacing + 
          attr->margin_top +
          attr->legend_spacing;
      aux->tip_occ = aux->tip_occ + 1;
    }
    else
    {
      double ly,ry;

      pll_svg_data_t * nb_data  = node->next->back->data;
      pll_svg_data_t * nnb_data = node->next->next->back->data;

      ly = nb_data->y;
      ry = nnb_data->y;
      y = (ly + ry) / 2.0;

      draw_line(fp,x,ly,x,ry,attr->stroke_width);
      draw_circle(fp,x,y,attr->node_radius);
    }
    /* horizontal line */
    draw_line(fp,
              px,
              y,
              x,
              y,
              attr->stroke_width);
    data->y = y;

    if (!node->next)
    {
      fprintf(fp, "<text x=\"%f\" y=\"%f\" "
                  "font-size=\"%ld\" font-family=\"Arial;\">%s</text>\n",
              x+5,
              y+attr->font_size/3.0,
              attr->font_size,
              node->label);
    }
    else
      fprintf(fp, "\n");
  }
  else
  {
    double ly,ry,x;
    pll_svg_data_t * nb_data  = (pll_svg_data_t *)(node->next->back->data);

    ly = nb_data->y;
    ry = parent_data->y;
    y = (ly + ry) / 2.0;
    x = attr->margin_left;

    draw_line(fp,x,ly,x,ry,attr->stroke_width);
    draw_circle(fp,x,y,attr->node_radius);
  }
}

static void utree_scaler_init(const pll_svg_attrib_t * attr,
                              pll_svg_aux_t * aux,
                              pll_utree_t * tree)
{
  unsigned int i;
  double len = 0;
  double label_len;

  /* compute the length of all tip-to-root paths and store the longest one in
     max_tree_len */
  for (i = 0; i < tree->tip_count; ++i)
  {
    pll_unode_t * node = tree->nodes[i];

    len = node->length;
    node = node->back;
    while(1)
    {
      pll_svg_data_t * data = (pll_svg_data_t *)(node->data);
      pll_svg_data_t * nb_data = (pll_svg_data_t *)(node->next->back->data);
      pll_svg_data_t * nnb_data = (pll_svg_data_t *)(node->next->next->back->data);

      if (nb_data->height > data->height)
        node = node->next->back;
      else if (nnb_data->height > data->height)
        node = node->next->next->back;
      else
        break;

      len += node->length;
    }
    
    if (len > aux->max_tree_len) 
      aux->max_tree_len = len;

    label_len = (attr->font_size / 1.5) * 
                (tree->nodes[i]->label ? strlen(tree->nodes[i]->label) : 0);

    len = (aux->canvas_width - label_len) / len;
    if (i == 0)
    {
      aux->scaler = len;
      aux->max_font_len = label_len;
    }
    else
      if (len < aux->scaler)
      {
        aux->scaler = len;
        aux->max_font_len = label_len;
      }
  }
}

static void print_header(FILE * fp,
                         pll_utree_t * tree,
                         const pll_svg_attrib_t * attr,
                         pll_svg_aux_t * aux)
{
  long svg_height;

  aux->canvas_width = attr->width - attr->margin_left - attr->margin_right;

  /* initialize pixel scaler (scaler) and compute max tree 
     length (max_tree_len) */
  utree_scaler_init(attr, aux, tree);

  svg_height = attr->margin_top + attr->legend_spacing + attr->margin_bottom + 
               attr->tip_spacing * tree->tip_count;

  /* print svg header tag with dimensions and grey border */
  fprintf(fp, "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"%ld\" "
          "height=\"%ld\" style=\"border: 1px solid #cccccc;\">\n",
          attr->width,
          svg_height);

  /* draw legend */
  if (attr->legend_show)
  {
    draw_line(fp,
              attr->margin_left,
              10,
              (aux->canvas_width - aux->max_font_len)*attr->legend_ratio + 
                attr->margin_left,
              10,
              3);

    fprintf(fp, "<text x=\"%f\" y=\"%f\" font-size=\"%ld\" "
                "font-family=\"Arial;\">%.*f</text>\n",
            (aux->canvas_width - aux->max_font_len)*attr->legend_ratio +
             attr->margin_left + 5,
            20-attr->font_size/3.0,
            attr->font_size,
            attr->precision,
            aux->max_tree_len * attr->legend_ratio);
  }

  /* uncomment to print a dashed border to indicate margins */
  
  /*
  fprintf(svg_fp, "<rect x=\"%ld\" y=\"%ld\" width=\"%ld\" fill=\"none\" "
          "height=\"%ld\" stroke=\"#999999\" stroke-dasharray=\"5,5\" "
          "stroke-width=\"1\" />\n",
          opt_svg_marginleft, 
          opt_svg_margintop + legend_spacing, 
          svg_width - opt_svg_marginleft - opt_svg_marginright,
          svg_height - opt_svg_margintop - legend_spacing - opt_svg_marginbottom);
  */
  
}

static void svg_make(FILE * fp,
                     pll_utree_t * tree,
                     pll_unode_t * root,
                     const pll_svg_attrib_t * attr)
{

  /* initialize auxiliary variables */
  pll_svg_aux_t aux;
  aux.max_font_len = 0;
  aux.max_tree_len = 0;
  aux.canvas_width = 0;
  aux.tip_occ = 0;

  /* print SVG header */
  print_header(fp,tree,attr,&aux);

  /* compute position for each node */
  utree_set_offset(root,attr,&aux);

  /* plot tree */
  utree_plot(fp, root, attr, &aux);

  /* closing svg tag */
  fprintf(fp, "</svg>\n");
}

PLL_EXPORT pll_svg_attrib_t * pll_svg_attrib_create()
{
  pll_svg_attrib_t * x;

  x = (pll_svg_attrib_t *)malloc(sizeof(pll_svg_attrib_t));
  if (!x) return NULL;

  /* set some defaults */
  x->width = 1920;
  x->font_size = 12;
  x->tip_spacing = 20;
  x->stroke_width = 3;
  x->legend_show = 1;
  x->legend_spacing = 10;
  x->legend_ratio = 0.1;
  x->margin_left = 20;
  x->margin_right = 20;
  x->margin_bottom = 20;
  x->margin_top = 20;
  x->node_radius = 0;
  x->precision = 7;

  return x;
}

PLL_EXPORT void pll_svg_attrib_destroy(pll_svg_attrib_t * attrib)
{
  free(attrib);
}

PLL_EXPORT int pll_utree_export_svg(pll_utree_t * tree,
                                    pll_unode_t * root,
                                    const pll_svg_attrib_t * attribs,
                                    const char * filename)
{
  unsigned int i;

  /* clone the tree */
  int rc = PLL_SUCCESS;

  if (!root || !(root->next))
    return PLL_FAILURE;

  /* open output file for writing */
  FILE * fp = fopen(filename, "w");
  if (!fp)
  {
    return PLL_FAILURE;
  }


  /* backup data */
  void ** data_old = (void **)malloc((tree->tip_count+tree->inner_count) *
                                      sizeof(void *));
  if (!data_old)
  {
    return PLL_FAILURE;
  }

  /* copy old data */
  for (i = 0; i < tree->tip_count+tree->inner_count; ++i)
  {
    data_old[i] = tree->nodes[i]->data;
    tree->nodes[i]->data = NULL;
  }


  /* treat unrooted tree as rooted binary with a ternary root
     and compute the height of each node */
  //if (!utree_set_height(cloned))
  if (!utree_set_height(root))
    rc = PLL_FAILURE;
  else
    svg_make(fp, tree, root, attribs);

  fclose(fp);

  /* restore old data */
  for (i = 0; i < tree->tip_count+tree->inner_count; ++i)
  {
    if (tree->nodes[i]->data)
      free(tree->nodes[i]->data);
    tree->nodes[i]->data = data_old[i];
  }
  free(data_old);


  return rc;
}
