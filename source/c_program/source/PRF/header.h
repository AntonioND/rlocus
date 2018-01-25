/******************************************************************/
/*                                                                */
/* file:          header.h                                        */
/*                                                                */
/* version:       1.2                                             */
/*                                                                */
/* author:        B. Frenzel                                      */
/*                                                                */
/* date:          Jan 7 1993                                      */
/*                                                                */
/* description:                                                   */
/* This headerfile contains all definitions and declarations      */
/* for the files:                                                 */
/* test.c    test enivironment for testing null()                 */
/* null.c    main rootfinding routine                             */
/* muller.c  rootfinding routine using muller method              */
/* newton.c  rootfinding routine using Newton method              */
/*                                                                */
/* further necessary files:                                       */
/* complex.c contains complex operations                          */
/* tools.c   contains tools:                                      */
/*           fdvalue()   computes P(x) and optional P'(x)         */
/*           poldef()    deflates polynomial                      */
/*                                                                */
/* Copyright:                                                     */
/* Lehrstuhl fuer Nachrichtentechnik Erlangen                     */
/* Cauerstr. 7, 8520 Erlangen, FRG, 1993                          */
/* e-mail: int@nt.e-technik.uni-erlangen.de                       */
/*                                                                */
/******************************************************************/










/***** common declarations and defines for all files *****/
#include <hpmath.h>           /* mathematical functions like log10()..*/
/* #include <float.h>           defines for DBL_MIN ...              */

                            /* type definition for a complex number */
typedef  struct DCOMPLEX {
                double r,
                       i;
                } dcomplex;

#define DBL_EPSILON     2.2204460492503131E-16
#define DBL_MAX         1.7976931348623157E+308
#ifndef TRUE
#define  TRUE     1
#define  FALSE    0
#endif
#define  PI  3.14159265358979323846 /* circular transcendental nb.  */
#define  MAXCOEFF 5001     /* max. number of coefficients          */
#define  NLL     0L         /* pointer to 0                         */

                            /* macros for complex operations        */
                            /* for fast operations                  */
#define CADD(x,a,b)     {(x).r = (a).r+(b).r; (x).i = (a).i+(b).i;}
#define CMUL(x,a,b)     {(x).r = (a).r*(b).r - (a).i*(b).i; \
                         (x).i = (a).i*(b).r + (a).r*(b).i;}
#define COMPLEXM(x,a,b) {(x).r = (a); (x).i = (b);}

#ifndef TOOLS
                            /* deflates polynomial with Horner     */
extern int    poldef();
                            /* delivers P(x) and optional P'(x)    */
extern void   fdvalue();
#endif

#ifndef COMPL    /* Access to extern operations with complex numbers */
                 /* stored in file complex.c                         */
                           /* add two complex numbers                */
extern dcomplex Cadd();
                           /* subtract two complex numbers           */
extern dcomplex Csub();
                           /* multiply two complex numbers           */
extern dcomplex Cmul();
                           /* divide two complex numbers             */
extern dcomplex Cdiv();
                           /* compute a complex number out of two    */
                           /* real numbers                           */
extern dcomplex Complex();
                           /* absolute value of a complex number     */
extern double   Cabs();
                           /* argument of a complex number           */
extern double   Carg();
                           /* conjugate complex number               */
extern dcomplex Conjg();
                           /* root of a complex number               */
extern dcomplex Csqrt();
                           /* multiply a double and a complex number */
extern dcomplex RCmul();
                           /* add a real and a complex number        */
extern dcomplex RCadd();
#endif










#ifdef NUL
                     /* rootfinding routine using Muller's method */
extern dcomplex muller();
                     /* rootfinding routine using Newton's method */
extern dcomplex newton_();

unsigned char poly_check();
unsigned char lin_or_quad();
void monic();
#endif










#ifdef MULLER
#define ITERMAX   150   /* max. number of iteration steps                 */
#define CONVERGENCE 100/* halve q2, when |P(x2)/P(x1)|^2 > CONVERGENCE   */
#define MAXDIST   1e3  /* max. relative change of distance between       */
                       /* x-values allowed in one step                   */
#define FACTOR    1e5  /* if |f2|<FACTOR*macc and (x2-x1)/x2<FACTOR*macc */
                       /* then root is determined; end routine           */
#define KITERMAX  1e3  /* halve distance between old and new x2 max.     */
                       /* KITERMAX times in case of possible overflow    */
#define FVALUE    1e36 /* initialisation of |P(x)|^2                     */

#define BOUND1    1.01 /* improve convergence in case of small changes   */
#define BOUND2    0.99 /* of |P(x)|^2                                    */
#define BOUND3    0.01

#define BOUND4    sqrt(DBL_MAX)/1e4 /* if |P(x2).r|+|P(x2).i|>BOUND4 =>  */
                       /* suppress overflow of |P(x2)|^2                 */
#define BOUND6    log10(BOUND4)-4   /* if |x2|^nred>10^BOUND6 =>         */
                       /* suppress overflow of P(x2)                     */
#define BOUND7    1e-5 /* relative distance between determined root and  */
                       /* real root bigger than BOUND7 => 2. iteration   */
#define NOISESTART DBL_EPSILON*1e2 /* when noise starts counting         */
#define NOISEMAX  5     /* if noise>NOISEMAX: terminate iteration        */

dcomplex x0,x1,x2,    /* common points [x0,f(x0)=P(x0)], ... [x2,f(x2)]  */
         f0,f1,f2,    /* of parabola and polynomial                      */
         h1,h2,       /* distance between x2 and x1                      */
         q2;          /* smaller root of parabola                        */
int      iter;        /* iteration counter                               */

void initialize();
void root_of_parabola();
void iteration_equation();
void suppress_overflow();
void too_big_functionvalues();
void convergence_check();
void compute_function();
void check_x_value();
void root_check();
#endif










#ifdef NEWTON
#define  ITERMAX  20   /* max. number of iterations                      */
#define  FACTOR   5    /* calculate new dx, when change of x0 is smaller */
                       /* than FACTOR*(old change of x0)                 */
#define  FVALUE   1E36 /* initialisation of |P(xmin)|                    */
#define  BOUND    sqrt(DBL_EPSILON)
                       /* if the imaginary part of the root is smaller   */
                       /* than BOUND5 => real root                       */
#define  NOISEMAX 5    /* max. number of iterations with no better value */
#endif


