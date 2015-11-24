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
#include "lbfgsb.h"

static int c__1 = 1;

int dpofa(double *a, int *lda, int *n, int *
	info)
{
    /* System generated locals */
    int a_dim1, a_offset, i__1, i__2, i__3;

    /* Local variables */
    int j, k;
    double s, t;
    int jm1;

/*
    dpofa factors a double precision symmetric positive definite 
    matrix. 

    dpofa is usually called by dpoco, but it can be called 
    directly with a saving in time if  rcond  is not needed. 
    (time for dpoco) = (1 + 18/n)*(time for dpofa) . 

    on entry 

       a       double precision(lda, n) 
               the symmetric matrix to be factored.  only the 
               diagonal and upper triangle are used. 

       lda     int 
               the leading dimension of the array  a . 

       n       int 
               the order of the matrix  a . 

    on return 

       a       an upper triangular matrix  r  so that  a = trans(r)*r 
               where  trans(r)  is the transpose. 
               the strict lower triangle is unaltered. 
               if  info .ne. 0 , the factorization is not complete. 

       info    int 
               = 0  for normal return. 
               = k  signals an error condition.  the leading minor 
                    of order  k  is not positive definite. 

    linpack.  this version dated 08/14/78 . 
    cleve moler, university of new mexico, argonne national lab. 

    subroutines and functions 

    blas ddot 
    fortran sqrt 
    */

/*     internal variables */

/*     begin block with ...exits to 40 */


    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;

    /* Function Body */
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	*info = j;
	s = 0.;
	jm1 = j - 1;
	if (jm1 < 1) {
	    goto L20;
	}
	i__2 = jm1;
	for (k = 1; k <= i__2; ++k) {
	    i__3 = k - 1;
	    t = a[k + j * a_dim1] - ddot(&i__3, &a[k * a_dim1 + 1], &c__1, &
		    a[j * a_dim1 + 1], &c__1);
	    t /= a[k + k * a_dim1];
	    a[k + j * a_dim1] = t;
	    s += t * t;
/* L10: */
	}
L20:
	s = a[j + j * a_dim1] - s;
/*     ......exit */
	if (s <= 0.) {
	    goto L40;
	}
	a[j + j * a_dim1] = sqrt(s);
/* L30: */
    }
    *info = 0;
L40:
    return 0;
} /* dpofa */

/* ====================== The end of dpofa =============================== */

int dtrsl(double *t, int *ldt, int *n, 
	double *b, int *job, int *info)
{
    /* System generated locals */
    int t_dim1, t_offset, i__1, i__2;

    /* Local variables */
    int j, jj, case__;
    double temp;
    /*
    extern double ddot(int *, double *, int *, double *, 
	    int *);
    extern  int daxpy(int *, double *, double *, 
	    int *, double *, int *);
    */


/*     dtrsl solves systems of the form */

/*                   t * x = b */
/*     or */
/*                   trans(t) * x = b */

/*     where t is a triangular matrix of order n. here trans(t) */
/*     denotes the transpose of the matrix t. */

/*     on entry */

/*         t         double precision(ldt,n) */
/*                   t contains the matrix of the system. the zero */
/*                   elements of the matrix are not referenced, and */
/*                   the corresponding elements of the array can be */
/*                   used to store other information. */

/*         ldt       int */
/*                   ldt is the leading dimension of the array t. */

/*         n         int */
/*                   n is the order of the system. */

/*         b         double precision(n). */
/*                   b contains the right hand side of the system. */

/*         job       int */
/*                   job specifies what kind of system is to be solved. */
/*                   if job is */

/*                        00   solve t*x=b, t lower triangular, */
/*                        01   solve t*x=b, t upper triangular, */
/*                        10   solve trans(t)*x=b, t lower triangular, */
/*                        11   solve trans(t)*x=b, t upper triangular. */

/*     on return */

/*         b         b contains the solution, if info .eq. 0. */
/*                   otherwise b is unaltered. */

/*         info      int */
/*                   info contains zero if the system is nonsingular. */
/*                   otherwise info contains the index of */
/*                   the first zero diagonal element of t. */

/*     linpack. this version dated 08/14/78 . */
/*     g. w. stewart, university of maryland, argonne national lab. */

/*     subroutines and functions */

/*     blas daxpy,ddot */
/*     fortran mod */

/*     internal variables */


/*     begin block permitting ...exits to 150 */

    /*        check for zero diagonal elements. */

    /* Parameter adjustments */
    t_dim1 = *ldt;
    t_offset = 1 + t_dim1;
    t -= t_offset;
    --b;

    /* Function Body */
    i__1 = *n;

    for (*info = 1; *info <= i__1; ++(*info)) {
        /*     ......exit */
        if (t[*info + *info * t_dim1] == 0.) {
            goto L150;
        }
        /* L10: */
    }
    *info = 0;
    /* SRB: there is a problem with this. Not sure */
    /* It keeps returning with info==2.  Feb 17 2015 */

    /*        determine the task and go to it. */

    case__ = 1;
    if (*job % 10 != 0) {
        case__ = 2;
    }
    if (*job % 100 / 10 != 0) {
        case__ += 2;
    }
    switch (case__) {
        case 1:  goto L20;
        case 2:  goto L50;
        case 3:  goto L80;
        case 4:  goto L110;
    }

    /*        solve t*x=b for t lower triangular */

L20:
    b[1] /= t[t_dim1 + 1];
    if (*n < 2) {
        goto L40;
    }
    i__1 = *n;
    for (j = 2; j <= i__1; ++j) {
        temp = -b[j - 1];
        i__2 = *n - j + 1;
        daxpy(&i__2, &temp, &t[j + (j - 1) * t_dim1], &c__1, &b[j], &c__1);
        b[j] /= t[j + j * t_dim1];
        /* L30: */
    }
L40:
    goto L140;

/*        solve t*x=b for t upper triangular. */

L50:
    b[*n] /= t[*n + *n * t_dim1];
    if (*n < 2) {
        goto L70;
    }
    i__1 = *n;
    for (jj = 2; jj <= i__1; ++jj) {
        j = *n - jj + 1;
        temp = -b[j + 1];
        daxpy(&j, &temp, &t[(j + 1) * t_dim1 + 1], &c__1, &b[1], &c__1);
        b[j] /= t[j + j * t_dim1];
        /* L60: */
    }
L70:
    goto L140;

/*        solve trans(t)*x=b for t lower triangular. */

L80:
    b[*n] /= t[*n + *n * t_dim1];
    if (*n < 2) {
        goto L100;
    }
    i__1 = *n;
    for (jj = 2; jj <= i__1; ++jj) {
        j = *n - jj + 1;
        i__2 = jj - 1;
        b[j] -= ddot(&i__2, &t[j + 1 + j * t_dim1], &c__1, &b[j + 1], &c__1);
        b[j] /= t[j + j * t_dim1];
        /* L90: */
    }
L100:
    goto L140;

/*        solve trans(t)*x=b for t upper triangular. */

L110:
    b[1] /= t[t_dim1 + 1];
    if (*n < 2) {
        goto L130;
    }
    i__1 = *n;
    for (j = 2; j <= i__1; ++j) {
        i__2 = j - 1;
        b[j] -= ddot(&i__2, &t[j * t_dim1 + 1], &c__1, &b[1], &c__1);
        b[j] /= t[j + j * t_dim1];
        /* L120: */
    }
L130:
L140:
L150:
    return 0;
} /* dtrsl */

