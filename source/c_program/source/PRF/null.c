/******************************************************************/
/*                                                                */
/* file:          null.c                                          */
/*                                                                */
/* main function: null()                                          */
/*                                                                */
/* version:       1.3                                             */
/*                                                                */
/* author:        B. Frenzel                                      */
/*                                                                */
/* date:          Jan 19 1993                                     */
/*                                                                */
/* input:         p[]     coefficient vector of the original      */
/*                        polynomial                              */
/*                pred[]  coefficient vector of the deflated      */
/*                        polynomial                              */
/*                *n      the highest exponent of the original    */
/*                        polynomial                              */
/*                flag    flag = TRUE  => complex coefficients    */
/*                        flag = FALSE => real    coefficients    */
/*                                                                */
/* output:        root[]  vector of determined roots              */
/*                *maxerr estimation of max. error of all         */
/*                        determined roots                        */
/*                                                                */
/* subroutines:   poly_check(),quadratic(),lin_or_quad(),monic(), */
/*                muller(),newton_()                               */
/*                                                                */
/* description:                                                   */
/* main rootfinding routine for polynomials with complex or real  */
/* coefficients using Muller's method combined with Newton's m.   */
/*                                                                */
/* Copyright:                                                     */
/* Lehrstuhl fuer Nachrichtentechnik Erlangen                     */
/* Cauerstr. 7, 8520 Erlangen, FRG, 1993                          */
/* e-mail: int@nt.e-technik.uni-erlangen.de                       */
/*                                                                */
/******************************************************************/
#define  NUL
#include "header.h"










/***** main function null() *****/
unsigned char null(dcomplex *p,dcomplex *pred,int *n,dcomplex *root,
                   double *maxerr,unsigned char flag)
/*dcomplex *p,     coefficient vector of the original polynomial   */
/*         *pred,  coefficient vector of the deflated polynomial   */
/*         *root;  determined roots                                */
/*int      *n;     the highest exponent of the original polynomial */
/*double   *maxerr; max. error of all determined roots             */
/*unsigned char flag;  flag = TRUE  => complex coefficients        */
/*                     flag = FALSE => real    coefficients        */
{
     double   newerr; /* error of actual root                      */
     dcomplex ns;     /* root determined by Muller's method        */
     int      nred,   /* highest exponent of deflated polynomial   */
              i;      /* counter                                   */
     unsigned char error; /* indicates an error in poly_check      */
     int      red,
              diff;   /* number of roots at 0                      */

     *maxerr = 0.;    /* initialize max. error of determined roots */
     nred = *n;       /* At the beginning: degree defl. polyn. =   */
                      /* degree of original polyn.                 */

                      /* check input of the polynomial             */
     error = poly_check(p,&nred,n,root);
     diff  = (*n-nred); /* reduce polynomial, if roots at 0        */
     p    += diff;
     *n   =  nred;

     if (error)
          return error; /* error in poly_check(); return error     */

                        /* polynomial is linear or quadratic       */
     if (lin_or_quad(p,nred,root)==0) {
          *n += diff;     /* remember roots at 0                   */
          *maxerr = DBL_EPSILON;
          return 0;       /* return no error                       */
     }

     monic(p,n);          /* get monic polynom                     */

     for (i=0;i<=*n;i++)  pred[i]=p[i];  /* original polynomial    */
                           /* = deflated polynomial at beginning   */
                           /* of Muller                            */

     do {                  /* main loop of null()                  */
                           /* Muller method                        */
          ns = muller(pred,nred);



                           /* Newton method                        */
          root[nred-1] = newton_(p,*n,ns,&newerr,flag);

                           /* stores max. error of all roots       */
          if (newerr>*maxerr)
               *maxerr=newerr;
                           /* deflate polynomial                   */
          red = poldef(pred,nred,root,flag);
          pred += red;        /* forget lowest coefficients        */
          nred -= red;        /* reduce degree of polynomial       */
     } while (nred>2);
                              /* last one or two roots             */
     (void) lin_or_quad(pred,nred,root);
     if (nred==2) {
          root[1] = newton_(p,*n,root[1],&newerr,flag);
          if (newerr>*maxerr)
               *maxerr=newerr;
     }
     root[0] = newton_(p,*n,root[0],&newerr,flag);
     if (newerr>*maxerr)
          *maxerr=newerr;

    *n += diff;              /* remember roots at 0               */
    return 0;                /* return no error                   */
}










/***** poly_check() check the formal correctness of input *****/
unsigned char poly_check(dcomplex *pred,int *nred,int *n,dcomplex *root)
/*dcomplex *pred,  coefficient vector of the original polynomial   */
/*         *root;  determined roots                                */
/*int      *nred,  highest exponent of the deflated polynomial     */
/*         *n;     highest exponent of the original polynomial     */
{
     int  i = -1, /* i stores the (reduced) real degree            */
          j;      /* counter variable                              */
     unsigned char
           notfound=TRUE; /* indicates, whether a coefficient      */
                          /* unequal zero was found                */

     if (*n<0) return 1;  /* degree of polynomial less than zero   */
                          /* return error                          */

     for (j=0;j<=*n;j++) {      /* determines the "real" degree of       */
          if(Cabs(pred[j])!=0.) /* polynomial, cancel leading roots      */
               i=j;
     }
     if (i==-1) return 2;   /* polynomial is a null vector; return error */
     if (i==0) return 3;    /* polynomial is constant unequal null;      */
                            /* return error                              */

     *n=i;                  /* set new exponent of polynomial            */
     i=0;                   /* reset variable for exponent               */
     do {                   /* find roots at 0                           */
          if (Cabs(pred[i])==0.)
               i++;
          else
               notfound=FALSE;
     } while (i<=*n && notfound);

          if (i==0) {            /* no root determined at 0              */
               *nred = *n;       /* original degree=deflated degree and  */
               return 0;         /* return no error                      */
          } else {               /* roots determined at 0:               */
               for (j=0;j<=i-1;j++) /* store roots at 0                  */
                    root[*n-j-1] = Complex(0.,0.);
               *nred = *n-i;  /* reduce degree of deflated polynomial    */
               return 0;      /* and return no error                     */
          }
}










/***** quadratic() calculates the roots of a quadratic polynomial *****/
void quadratic(dcomplex *pred,dcomplex *root)
/*dcomplex *pred,  coefficient vector of the deflated polynomial   */
/*         *root;  determined roots                                */

{
     dcomplex discr,       /* discriminate                         */
              Z1,Z2,       /* numerators of the quadratic formula  */
              N;           /* denominator of the quadratic formula */

                           /* discr = p1^2-4*p2*p0                 */
     discr   = Csub(Cmul(pred[1],pred[1]),
               RCmul(4.,Cmul(pred[2],pred[0])));
                           /* Z1 = -p1+sqrt(discr)                 */
     Z1      = Cadd(RCmul(-1.,pred[1]),Csqrt(discr));
                           /* Z2 = -p1-sqrt(discr)                 */
     Z2      = Csub(RCmul(-1.,pred[1]),Csqrt(discr));
                           /* N  = 2*p2                            */
     N       = RCmul(2.,pred[2]);
     root[0] = Cdiv(Z1,N); /* first root  = Z1/N                   */
     root[1] = Cdiv(Z2,N); /* second root = Z2/N                   */
}










/***** lin_or_quad() calculates roots of lin. or quadratic equation *****/
unsigned char lin_or_quad(dcomplex *pred,int nred,dcomplex *root)
/*dcomplex *pred,  coefficient vector of the deflated polynomial   */
/*         *root;  determined roots                                */
/*int      nred;   highest exponent of the deflated polynomial     */
{
     if (nred==1) {     /* root = -p0/p1                           */
          root[0] = Cdiv(RCmul(-1.,pred[0]),pred[1]);
          return 0;     /* and return no error                     */
     } else if (nred==2) { /* quadratic polynomial                 */
          quadratic(pred,root);
          return 0;        /* return no error                      */
     }

     return 1; /* nred>2 => no roots were calculated               */
}










/***** monic() computes monic polynomial for original polynomial *****/
void monic(dcomplex *p,int *n)
/*dcomplex *p;     coefficient vector of the original polynomial   */
/*int      *n;     the highest exponent of the original polynomial */
{
     double factor;   /* stores absolute value of the coefficient  */
                      /* with highest exponent                     */
     int    i;        /* counter variable                          */

     factor=1./Cabs(p[*n]);     /* factor = |1/pn|                 */
     if ( factor!=1.)           /* get monic pol., when |pn| != 1  */
         for (i=0;i<=*n;i++)
               p[i]=RCmul(factor,p[i]);
}
