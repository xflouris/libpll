/*
    Copyright (C) 2015 Diego Darriba, Tomas Flouri

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

    Contact: Diego Darriba <Diego.Darriba@h-its.org>,
    Exelixis Lab, Heidelberg Instutute for Theoretical Studies
    Schloss-Wolfsbrunnenweg 35, D-69118 Heidelberg, Germany
*/
#ifndef COMMON_H_
#define COMMON_H_

#ifndef PLL_H_
#define PLL_H_
#include "pll.h"
#endif

extern const unsigned int odd5_map[256];

unsigned int get_attributes(int argc, char **argv);
void skip_test();

pll_partition_t * parse_msa(const char * filename,
                            unsigned int states,
                            unsigned int rate_cats,
                            unsigned int rate_matrices,
                            pll_utree_t * tree,
                            unsigned int attributes);

pll_partition_t * parse_msa_reduced(const char * filename,
                            unsigned int states,
                            unsigned int rate_cats,
                            unsigned int rate_matrices,
                            pll_utree_t * tree,
                            unsigned int attributes,
                            unsigned int max_sites);
int cb_full_traversal(pll_unode_t * node);
int cb_rfull_traversal(pll_rnode_t * node);

/* print error and exit */
void fatal(const char * format, ...) __attribute__ ((noreturn));
char * xstrdup(const char * s);
void * xmalloc(size_t size);

#endif /* COMMON_H_ */
