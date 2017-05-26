/*
    Copyright (C) 2017 Alexey Kozlov, Diego Darriba, Tomas Flouri

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

    Contact: Alexey Kozlov <Alexey.Kozlov@h-its.org>,
    Exelixis Lab, Heidelberg Instutute for Theoretical Studies
    Schloss-Wolfsbrunnenweg 35, D-69118 Heidelberg, Germany
*/
#include "common.h"

#define N_STATES_NT 4
#define N_STATES_AA 20
#define N_STATES_ODD 5
#define N_CAT_GAMMA 4
#define N_BASE_FREQS 3
#define N_SUBST_MAT 3
#define N_SITES 5
#define N_BRANCHES 5
#define FLOAT_PRECISION 4

#define MIN(a,b) (a<b?a:b)

#define DATATYPE_NT 0
#define DATATYPE_AA  1
#define DATATYPE_ODD 2

static double cat_rates[N_CAT_GAMMA] = {1e-31, 1e-6, 1.0, 100.0 };
static unsigned int params_indices[16] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

static double branch_lengths[N_BRANCHES] = { 1e-6, 1e-2, 0.2, 1.0, 100. };
static unsigned int matrix_indices[N_BRANCHES] = { 0, 1, 2, 3, 4 };

static pll_partition_t *part_nt, *part_aa, *part_odd;

void check_matrix(pll_partition_t * p, unsigned int matrix_index)
{
  const double * mat = p->pmatrix[matrix_index];
  unsigned int i;
  for (i = 0; i < p->rate_cats * p->states_padded * p->states; ++i)
  {
    if (isnan(mat[i]) || !isfinite(mat[i]) || mat[i] < 0.)
      fatal("ERROR invalid value in p-matrix %u: %lf\n", matrix_index, mat[i]);
  }
}

pll_partition_t * init_partition(unsigned int attrs, int datatype)
{
  unsigned int i;

  unsigned int states = 0;

  switch(datatype)
  {
    case DATATYPE_NT:
      states = N_STATES_NT;
      break;
    case DATATYPE_AA:
      states = N_STATES_AA;
      break;
    case DATATYPE_ODD:
      states = N_STATES_ODD;
      break;
    default:
      assert(0);
  }

  pll_partition_t * p = pll_partition_create(N_BRANCHES-1,  /* tips */
                                             0,             /* clv vectors */
                                             states,
                                             N_SITES,
                                             1,             /* rate matrices */
                                             N_BRANCHES,    /* prob matrices */
                                             N_CAT_GAMMA,
                                             0,             /* scalers */
                                             attrs);

  if (!p)
    fatal("ERROR creating partition: %s\n", pll_errmsg);

  pll_set_category_rates(p, cat_rates);

  printf("category rates(%d): [", N_CAT_GAMMA);
  for (i = 0; i < N_CAT_GAMMA; ++i)
    printf("%lf ", cat_rates[i]);
  printf("]\n");

  return p;
}

void init(unsigned int attrs)
{
  part_nt = init_partition(attrs, DATATYPE_NT);
  part_aa = init_partition(attrs, DATATYPE_AA);
  part_odd = init_partition(attrs, DATATYPE_ODD);
}

double * init_freqs(unsigned int n_states)
{
  double * base_freqs =
      (double *) malloc(N_BASE_FREQS * n_states * sizeof(double));

  unsigned int k;
  unsigned int offset = 0;

  for (k = 0; k < n_states; ++k)   /* equal freqs */
    base_freqs[offset + k] = 1.0 / (double) n_states;

  offset += n_states;
  for (k = 0; k < n_states; ++k)   /* skewed freqs */
  {
    base_freqs[offset + k] = 1.0 / (double) n_states;
    const double skew = 1.0 / (3. * n_states);
    if (k % 2 == 0)
      base_freqs[offset + k] += skew;
    else if (k != n_states-1)
      base_freqs[offset + k] -= skew;
  }

  const double minfreq = 1e-3;
  const double maxfreq = (1. - 0.5 * n_states * minfreq) / (0.5 * n_states);

  offset += n_states;
  for (k = 0; k < n_states; ++k)   /* extreme freqs */
    base_freqs[offset + k] = (k % 2 == 0) ? minfreq : maxfreq;

  return base_freqs;
}

double * init_rates(unsigned int n_rates)
{
  double * subst_rates =
      (double *) malloc(N_SUBST_MAT * n_rates * sizeof(double));

  unsigned int k;
  unsigned int offset = 0;
  for (k = 0; k < n_rates; ++k)   /* equal rates */
    subst_rates[offset + k] = 1.0;

  offset += n_rates;
  for (k = 0; k < n_rates; ++k)   /* skewed rates */
  {
    subst_rates[offset + k] = 1.0;
    const double skew = 5.;
    if (k % 2 == 0)
      subst_rates[offset + k] *= skew;
    else if (k != n_rates-1)
      subst_rates[offset + k] /= skew;
  }

  offset += n_rates;
  for (k = 0; k < n_rates-1; ++k)   /* extreme rates */
    subst_rates[offset + k] = (k % 2 == 0) ? 1e-3 : 1e3;
  subst_rates[offset + n_rates-1] = 1.0;

  return subst_rates;
}

int eval(pll_partition_t * partition, double * base_freqs, double * subst_rates)
{
  unsigned int i;

  pll_set_frequencies(partition, 0, base_freqs);
  pll_set_subst_params(partition, 0, subst_rates);

  printf("datatype = ");
  if (partition->states == 4)
    printf("DNA");
  else if (partition->states == 20)
    printf("PROT");
  else
    printf("ODD");
  printf("\n");

  printf("base freqs: [");
  for (i = 0; i < partition->states; ++i)
    printf("%lf ", partition->frequencies[0][i]);
  printf("]\n");

  printf("subst rates: [");
  for (i = 0; i < partition->states * (partition->states -1) / 2; ++i)
    printf("%lf ", partition->subst_params[0][i]);
  printf("]\n");

  pll_update_prob_matrices(partition,
                           params_indices,
                           matrix_indices,
                           branch_lengths,
                           N_BRANCHES);

  for (i = 0; i < N_BRANCHES; ++i)
  {
    printf("P-matrix: %d, brlen = %lf\n", i, branch_lengths[i]);
    pll_show_pmatrix(partition, i, 9);
    check_matrix(partition, i);
  }

  return 1;
}

void cleanup()
{
  pll_partition_destroy(part_nt);
  pll_partition_destroy(part_aa);
  pll_partition_destroy(part_odd);
}

int main(int argc, char * argv[])
{
  unsigned int i, j, k;

  /* check attributes */
  unsigned int attributes = get_attributes(argc, argv);

  /* pattern tip is not relevant for pmatrix computation */
  if (attributes & PLL_ATTRIB_PATTERN_TIP)
    skip_test();

  init(attributes);

  pll_partition_t * parts[] = {part_nt, part_aa, part_odd};

  for (i = 0; i < 3; ++i)
  {
    pll_partition_t * p = parts[i];

    const unsigned int n_states = p->states;
    const unsigned int n_rates = n_states * (n_states-1) / 2;

    double * base_freqs = init_freqs(n_states);
    double * subst_rates = init_rates(n_rates);

    for (j = 0; j < N_BASE_FREQS; ++j)
    {
      for (k = 0; k < N_SUBST_MAT; ++k)
      {
        eval(p, base_freqs + j * n_states, subst_rates + k * n_rates);
      }
    }

    free(base_freqs);
    free(subst_rates);
  }

  cleanup();

  return (0);
}
