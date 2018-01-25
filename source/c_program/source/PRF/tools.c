/******************************************************************/
/*                                                                */
/* file:          tools.c                                         */
/*                                                                */
/* version:       1.2                                             */
/*                                                                */
/* author:        B. Frenzel                                      */
/*                                                                */
/* date:          Jan 7 1993                                      */
/*                                                                */
/* function:      description:                                    */ 
/*                                                                */
/* hornc()        hornc() deflates the polynomial with coeff.     */
/*                stored in pred[0],...,pred[n] by Horner's       */
/*                method (one root)                               */
/*                                                                */
/* horncd()       horncd() deflates the polynomial with coeff.    */
/*                stored in pred[0],...,pred[n] by Horner's       */
/*                method (two roots)                              */
/*                                                                */
/* poldef()       decides whether to call hornc() or horncd()     */
/*                                                                */
/* fdvalue()      fdvalue() computes f=P(x0) by Horner's method;  */
/*                if flag=TRUE, it additionally computes the      */
/*                derivative df=P'(x0)                            */
/*                                                                */
/* Copyright:                                                     */
/* Lehrstuhl fuer Nachrichtentechnik Erlangen                     */
/* Cauerstr. 7, 8520 Erlangen, FRG, 1993                          */
/* e-mail: int@nt.e-technik.uni-erlangen.de                       */
/*                                                                */
/******************************************************************/
#define  TOOLS
#include "header.h"










/***** Horner method to deflate one root *****/ 
void hornc(dcomplex *pred,int nred,dcomplex x0,unsigned char flag) 
/*dcomplex      *pred,  coefficient vector of the polynomial           */
/*              x0;     root to be deflated                            */
/*int           nred;   the highest exponent of (deflated) polynomial  */
/*unsigned char flag;   indicates how to reduce polynomial             */

{
     int      i;          /* counter                        */
     dcomplex help1;      /* help variable                  */

     if ((flag&1)==0)     /* real coefficients              */
          for(i=nred-1; i>0; i--) 
               pred[i].r += (x0.r*pred[i+1].r);
     else                 /* complex coefficients           */
          for (i=nred-1; i>0; i--) {
               CMUL(help1,pred[i+1],x0);
               CADD(pred[i],help1,pred[i]);
          }
}










/***** Horner method to deflate two roots *****/
void horncd(dcomplex *pred,int nred,double a,double b)
/*dcomplex *pred;     coefficient vector of the polynomial           */
/*double   a,         coefficients of the quadratic polynomial       */
/*         b;         x^2+ax+b                                       */
/*int      nred;      the highest exponent of (deflated) polynomial  */
{
     int i;        /* counter */

     pred[nred-1].r += pred[nred].r*a; 
     for (i=nred-2; i>1; i--)
          pred[i].r += (a*pred[i+1].r+b*pred[i+2].r);
}










/***** main routine to deflate polynomial *****/
int poldef(dcomplex *pred,int nred,dcomplex *root,unsigned char flag)
/*dcomplex *pred,     coefficient vector of the polynomial           */
/*         *root;     vector of determined roots                     */
/*int       nred;     the highest exponent of (deflated) polynomial  */
/*unsigned char flag;  indicates how to reduce polynomial            */
{
     double   a,   /* coefficients of the quadratic polynomial       */
              b;   /* x^2+ax+b                                       */
     dcomplex x0;  /* root to be deflated                            */


     x0 = root[nred-1];
     if (x0.i!=0.)             /* x0 is complex                      */
          flag |=2; 

     if (flag==2) {            /* real coefficients and complex root */
          a = 2*x0.r;          /* => deflate x0 and Conjg(x0)        */
          b = -(x0.r*x0.r+x0.i*x0.i); 
          root[nred-2]=Conjg(x0); /* store second root = Conjg(x0)   */
          horncd(pred,nred,a,b);
          return 2;            /* two roots deflated                 */
     } else {
          hornc(pred,nred,x0,flag); /* deflate only one root         */
          return 1;            /* one root deflated                  */
     }
}










/***** fdvalue computes P(x0) and optional P'(x0) *****/
void fdvalue(dcomplex *p,int n,dcomplex *f,dcomplex *df,dcomplex x0,
             unsigned char flag) 
/*dcomplex      *p,     coefficient vector of the polynomial P(x)   */
/*              *f,     the result f=P(x0)                          */
/*              *df,    the result df=P'(x0), if flag=TRUE          */
/*              x0;     polynomial will be computed at x0           */
/*int           n;      the highest exponent of p                   */
/*unsigned char flag;   flag==TRUE => compute P'(x0)                */
{
     int      i;     /* counter                                     */
     dcomplex help1; /* help variable                               */ 

     *f  = p[n];
     if (flag==TRUE) {              /* if flag=TRUE, compute P(x0)  */
          COMPLEXM(*df,0.,0.);      /* and P'(x0)                   */
          for (i=n-1; i>=0; i--) {
               CMUL(help1,*df,x0);  /* *df = *f   + *df * x0        */
               CADD(*df,help1,*f);
               CMUL(help1,*f,x0);   /* *f  = p[i] + *f * x0         */
               CADD(*f,help1,p[i]);
          }
     } else                       /* otherwise: compute only P(x0)  */
          for (i=n-1; i>=0; i--) {
               CMUL(help1,*f,x0);   /* *f = p[i] + *f * x0          */
               CADD(*f,help1,p[i]);
          }
}
