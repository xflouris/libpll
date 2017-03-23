lib_LTLIBRARIES=libpll.la

AM_CFLAGS=-Wall -Wsign-compare -D_GNU_SOURCE -std=c99 -O3 -fPIC

LIBPLLHEADERS=\
pll.h

AM_YFLAGS = -d -p `${SED} -n 's/.*_\(.*\)/pll_\1_/p' <<<"$*"`
AM_LFLAGS = -P `${SED} -n 's/.*_\(.*\)/pll_\1_/p' <<<"$*"` -o lex.yy.c

libpll_la_SOURCES=\
fasta.c \
gamma.c \
likelihood.c \
list.c \
maps.c \
models.c \
pll.c \
output.c \
utree.c \
rtree.c \
derivatives.c \
partials.c \
compress.c \
utree_moves.c \
utree_svg.c \
parsimony.c \
core_derivatives.c \
core_partials.c \
core_pmatrix.c \
core_likelihood.c \
parse_utree.y \
parse_rtree.y \
parse_phylip.y \
lex_utree.l \
lex_rtree.l \
lex_phylip.l \
fast_parsimony.c \
stepwise.c \
random.c

libpll_la_CFLAGS = $(AM_CFLAGS)

# To allow cross-compilation, those SIMD flags will be used for the respective source files only  
AVX2FLAGS=-mfma -mavx2
AVXFLAGS=-mavx
SSEFLAGS=-msse3

SIMD_KERNELS=

if HAVE_AVX2
 SIMD_KERNELS+=libsimd_avx2.la
 libsimd_avx2_la_CFLAGS=$(AM_CFLAGS) $(AVX2FLAGS)
 libsimd_avx2_la_SOURCES=\
 core_partials_avx2.c \
 core_derivatives_avx2.c \
 core_pmatrix_avx2.c \
 core_likelihood_avx2.c \
 fast_parsimony_avx2.c
endif

if HAVE_AVX
SIMD_KERNELS+=libsimd_avx.la
libsimd_avx_la_CFLAGS=$(AM_CFLAGS) $(AVXFLAGS)
libsimd_avx_la_SOURCES=\
core_partials_avx.c \
core_derivatives_avx.c \
core_pmatrix_avx.c \
core_likelihood_avx.c \
fast_parsimony_avx.c
endif

if HAVE_SSE3
SIMD_KERNELS+=libsimd_sse.la
libsimd_sse_la_CFLAGS=$(AM_CFLAGS) $(SSEFLAGS)
libsimd_sse_la_SOURCES=\
core_partials_sse.c \
core_derivatives_sse.c \
core_likelihood_sse.c \
core_pmatrix_sse.c \
fast_parsimony_sse.c
endif

libpll_la_LIBADD = $(SIMD_KERNELS)
noinst_LTLIBRARIES = $(SIMD_KERNELS) 

# This is the version information for libtool and not the libpll version. For
#  more information see:
# https://www.gnu.org/software/libtool/manual/html_node/Updating-version-info.html
libpll_la_LDFLAGS = -version-info 0:0:0

pkgincludedir=$(includedir)/libpll
pkginclude_HEADERS = pll.h
