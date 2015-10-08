/*
 * utils.h
 *
 *  Created on: Sep 28, 2015
 *      Author: diego
 */

#ifndef UTILS_H_
#define UTILS_H_

#include "pll_optimize.h"

#include <assert.h>
#include <stdarg.h>
#include <search.h>
#include <time.h>

#define STATES    4
#define SUBST_PARAMS (STATES*(STATES-1)/2)

void fatal (const char * format, ...) __attribute__ ((noreturn));

pll_partition_t * partition_fasta_create (const char *file,
                                          unsigned int states,
                                          unsigned int n_rate_matrices,
                                          unsigned int n_rate_cats,
                                          int attributes,
                                          int rooted,
                                          unsigned int tip_count,
                                          const char **tipnames);

void set_missing_branch_length (pll_utree_t * tree, double length);

int * build_model_symmetries (const char * modelmatrix);

/* a callback function for performing a full traversal */
int cb_full_traversal(pll_utree_t * node);

#endif /* UTILS_H_ */
