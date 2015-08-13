/*
 * L-BFGS-B is released under the "New BSD License" (aka "Modified BSD License"
 * or "3-clause license")
 * Please read attached file License.txt
 *
 * ===========   L-BFGS-B (version 3.0.  April 25, 2011  ===================
 *
 *    This is a modified version of L-BFGS-B. Minor changes in the updated
 *    code appear preceded by a line comment as follows
 *
 *    c-jlm-jn
 *
 *    Major changes are described in the accompanying paper:
 *
 *        Jorge Nocedal and Jose Luis Morales, Remark on "Algorithm 778:
 *        L-BFGS-B: Fortran Subroutines for Large-Scale Bound Constrained
 *        Optimization"  (2011). To appear in  ACM Transactions on
 *        Mathematical Software,
 */
#ifndef LBFGSB_H
#define LBFGSB_H

#define PLL_LBFGSB_ERROR             1.0e-4

typedef int ftnlen;
typedef int logical;
#define TRUE_  (1)
#define FALSE_ (0)

#ifdef DEBUG
    #define DBG(fmt, ...) printf(fmt, ##__VA_ARGS__);
#else
    #define DBG(fmt, ...)
#endif

#include <math.h>
#include <stdio.h>

#ifdef MATLAB_MEX_FILE
#include "mex.h"
#endif

#ifndef abs
#define abs(x) ((x) >= 0 ? (x) : -(x))
#endif
#ifndef min
#define min(a,b) ((a) <= (b) ? (a) : (b))
#endif
#ifndef max
#define max(a,b) ((a) >= (b) ? (a) : (b))
#endif

/* Modified L-BFGS-B to use ints instead
 * of strings, for testing the "task". Way
 * simpler this way. For ease-of-use, use
 * these aliases. For each class of task,
 * make sure numbering is contiguous
 * so that I have easy tests for class
 * membership */

#define START    1
#define NEW_X    2
#define ABNORMAL 3 /* message: ABNORMAL_TERMINATION_IN_LNSRCH. */
#define RESTART  4 /* message: RESTART_FROM_LNSRCH. */

#define FG        10
#define FG_END    15
#define IS_FG(x) ( ((x)>=FG) ?  ( ((x)<=FG_END) ? 1 : 0 ) : 0 )
#define FG_LN     11
#define FG_LNSRCH FG_LN
#define FG_ST     12
#define FG_START  FG_ST

#define CONVERGENCE      20
#define CONVERGENCE_END  25
#define IS_CONVERGED(x) ( ((x)>=CONVERGENCE) ?  ( ((x)<=CONVERGENCE_END) ? 1 : 0 ) : 0 )
#define CONV_GRAD        21 /* message: CONVERGENCE: NORM_OF_PROJECTED_GRADIENT_<=_PGTOL. */
#define CONV_F           22 /* message: CONVERGENCE: REL_REDUCTION_OF_F_<=_FACTR*EPSMCH. */

#define STOP      30
#define STOP_END  40
#define IS_STOP(x) ( ((x)>=STOP) ?  ( ((x)<=STOP_END) ? 1 : 0 ) : 0 )
#define STOP_CPU  31 /* message: STOP: CPU EXCEEDING THE TIME LIMIT. */
#define STOP_ITER 32 /* message: STOP: TOTAL NO. of f AND g EVALUATIONS EXCEEDS LIM.  */
#define STOP_GRAD 33 /* message: STOP: THE PROJECTED GRADIENT IS SUFFICIENTLY SMALL. */

#define WARNING        100
#define WARNING_END    110
#define IS_WARNING(x) ( ((x)>=WARNING) ?  ( ((x)<=WARNING_END) ? 1 : 0 ) : 0 )
#define WARNING_ROUND  101  /* WARNING: ROUNDING ERRORS PREVENT PROGRESS */
#define WARNING_XTOL   102  /* WARNING: XTOL TEST SATISIED */
#define WARNING_STPMAX 103 /* WARNING: STP = STPMAX */
#define WARNING_STPMIN 104 /* WARNING: STP = STPMIN */

#define ERROR          200
#define ERROR_END      240
#define IS_ERROR(x) ( ((x)>=ERROR) ?  ( ((x)<=ERROR_END) ? 1 : 0 ) : 0 )
/* More specific conditions below */
#define ERROR_SMALLSTP 201 /* message: ERROR: STP .LT. STPMIN  */
#define ERROR_LARGESTP 202 /* message: ERROR: STP .GT. STPMAX  */
#define ERROR_INITIAL  203 /* message: ERROR: INITIAL G .GE. ZERO */
#define ERROR_FTOL     204 /* message: ERROR: FTOL .LT. ZERO   */
#define ERROR_GTOL     205 /* message: ERROR: GTOL .LT. ZERO   */
#define ERROR_XTOL     206 /* message: ERROR: XTOL .LT. ZERO   */
#define ERROR_STP0     207 /* message: ERROR: STPMIN .LT. ZERO */
#define ERROR_STP1     208 /* message: ERROR: STPMAX .LT. STPMIN */
#define ERROR_N0       209 /* ERROR: N .LE. 0 */
#define ERROR_M0       210 /* ERROR: M .LE. 0 */
#define ERROR_FACTR    211 /* ERROR: FACTR .LT. 0 */
#define ERROR_NBD      212 /* ERROR: INVALID NBD */
#define ERROR_FEAS     213 /* ERROR: NO FEASIBLE SOLUTION */

#if !defined( _USE_OPTIMIZED_BLAS )
/* make alias to our reference implementation */
#define daxpy daxpyRef
#define dcopy dcopyRef
#define ddot  ddotRef
#define dscal dscalRef
#else
#if defined( _USE_32_BIT_BLAS )
#define daxpy daxpy32
#define dcopy dcopy32
#define ddot  ddot32
#define dscal dscal32
#else
#define daxpy daxpy
#define dcopy dcopy
#define ddot  ddot
#define dscal dscal
#endif
#endif

/* and "word" was a char that was one fo these: */
#define WORD_DEFAULT 0 /* aka "---".  */
#define WORD_CON 1 /*  the subspace minimization converged. */
#define WORD_BND 2 /* the subspace minimization stopped at a bound. */
#define WORD_TNT 3 /*  the truncated Newton step has been used. */

/* First decision: use the included BLAS code
 * (un-optimized version taken from Netlib;
 * this is the "reference" implementation, or
 * "Ref" for short; this is your best option
 * since these routines should not be a bottleneck
 * in your computation);
 * or, use your own BLAS library, such 
 * as Intel MKL, ATLAS, GotoBLAS, etc.
 * (For example, you could use libmwblas
 * or libmwblas_compat32, which are included
 * with Mathworks and usually based on Intel MKL).
 *
 * The reason not to use your own BLAS library
 * is that (1) it may be faster but it won't be
 * a bottleneck, and (2) some BLAS libraries
 * use 32-bit/4byte ints, and others use 64-bit/
 * 8byte ints, and you pass by reference,
 * so if you get it wrong, you crash.
 * In Matlab, useing -lmwblascompat32
 *   uses names like ddot32
 *
 * In short, use the Ref version unless you feel
 * lucky.
 * */

#define PLL_LBFGSB_DEFAULT_FTOL .001
#define PLL_LBFGSB_DEFAULT_GTOL .9
#define PLL_LBFGSB_DEFAULT_XTOL .1
#define PLL_LBFGSB_DEFAULT_STEPMIN 0.
#define PLL_LBFGSB_DBL_EPSILON 2.2e-16

#ifdef __cplusplus
    extern "C" {
#endif

extern double ddot(int *, double *, int *, double *, 
        int *);

extern  int daxpy(int *, double *, double *, 
        int *, double *, int *);

extern  int dscal(int *, double *, double *, 
        int *);

extern  int dcopy(int *, double *, int *, 
	    double *, int *);

extern int setulb(int *n, int *m, double *x, 
	double *l, double *u, int *nbd, double *f, double 
	*g, double *factr, double *pgtol, double *wa, int *
	iwa, int *task, int *iprint, int *csave, logical *lsave, 
	int *isave, double *dsave); /* ftnlen task_len, ftnlen csave_len); */

extern int mainlb(int *n, int *m, double *x, 
        double *l, double *u, int *nbd, double *f, double 
        *g, double *factr, double *pgtol, double *ws, double *
        wy, double *sy, double *ss, double *wt, double *wn, 
        double *snd, double *z__, double *r__, double *d__, 
        double *t, double *xp, double *wa, int *index, 
        int *iwhere, int *indx2, int *task, int *iprint, 
        int *csave, logical *lsave, int *isave, double *dsave);


extern  int freev(int *, int *, int *, 
        int *, int *, int *, int *, logical *, logical *, 
        logical *, int *, int *);

extern int formk(int *, 
        int *, int *, int *, int *, int *, int *, 
        logical *, double *, double *, int *, double *, 
        double *, double *, double *, int *, int *, 
        int *);

extern  int formt(int *, double *, double *, 
        double *, int *, double *, int *);

extern int subsm(int *, int *, int *, int *, double *, double *, 
        int *, double *, double *, double *, double *,
        double *, double *, double *, double *, int *
        , int *, int *, double *, double *, int *, 
        int *);

extern  int active(int *, double *, double *,
        int *, double *, int *, int *, logical *, 
        logical *, logical *);

extern int cauchy(int *, double *, 
        double *, double *, int *, double *, int *, 
        int *, double *, double *, double *, int *, 
        double *, double *, double *, double *, 
        double *, int *, int *, double *, double *, 
        double *, double *, int *, int *, double *, 
        int *, double *);

extern  int cmprlb(int *, int *, double *, 
        double *, double *, double *, double *, 
        double *, double *, double *, double *, int *,
        double *, int *, int *, int *, logical *, 
        int *);

extern  int matupd(int *, int *, double *, 
        double *, double *, double *, double *, 
        double *, int *, int *, int *, int *, 
        double *, double *, double *, double *, 
        double *);

extern int lnsrlb(int *n, double *l, double *u, 
        int *nbd, double *x, double *f, double *fold, 
        double *gd, double *gdold, double *g, double *d__, 
        double *r__, double *t, double *z__, double *stp, 
        double *dnorm, double *dtd, double *xstep, double *
        stpmx, int *iter, int *ifun, int *iback, int *nfgv, 
        int *info, int *task, logical *boxed, logical *cnstnd, int *
        csave, int *isave, double *dsave); /* ftnlen task_len, ftnlen 
        csave_len); */

extern  int projgr(int *, double *, double *,
        int *, double *, double *, double *);

/* in linesearch.c */
extern int dcsrch(double *f, double *g, double *stp, 
        double *ftol, double *gtol, double *xtol, double *
        stpmin, double *stpmax, int *task, int *isave, double *
        dsave);/* ftnlen task_len);*/

extern  int dcstep(double *, double *,
        double *, double *, double *, double *,
        double *, double *, double *, logical *, double *,
        double *);

extern int dpofa(double *, int *, int *, 
		int *);

extern int dtrsl(double *, int *, int *,
		double *, int *, int *);

#ifdef __cplusplus
    }   /* extern "C" */
#endif

#endif /* LBFGSB_H */
