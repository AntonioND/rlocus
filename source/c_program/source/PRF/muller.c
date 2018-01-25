/******************************************************************/
/*                                                                */
/* file:          muller.c                                        */
/*                                                                */
/* main function: muller()                                        */
/*                                                                */
/* version:       1.2                                             */
/*                                                                */
/* author:        B. Frenzel                                      */
/*                                                                */
/* date:          Jan 7 1993                                      */ 
/*                                                                */
/* input:         pred[]  coefficient vector of the deflated      */
/*                        polynomial                              */
/*                nred    the highest exponent of the deflated    */
/*                        polynomial                              */  
/*                                                                */
/* return:        xb      determined root                         */
/*                                                                */
/* subroutines:   initialize(),root_of_parabola(),                */
/*                iteration_equation(), compute_function(),       */
/*                check_x_value(), root_check()                   */
/*                                                                */
/* description:                                                   */
/* muller() determines the roots of a polynomial with complex     */
/* coefficients with the Muller method; these roots are the       */
/* initial estimations for the following Newton method            */
/*                                                                */
/* Copyright:                                                     */
/* Lehrstuhl fuer Nachrichtentechnik Erlangen                     */
/* Cauerstr. 7, 91054 Erlangen, FRG, 1993                         */
/* e-mail: int@nt.e-technik.uni-erlangen.de                       */
/*                                                                */
/******************************************************************/
#define  MULLER
#include "header.h"











/***** main routine of Mueller's method *****/
dcomplex muller(dcomplex *pred,int nred)
/*dcomplex *pred;        coefficient vector of the deflated polynomial      */
/*int      nred;         the highest exponent of the deflated polynomial    */
{
     double   f1absq=FVALUE,  /* f1absq=|f1|^2                              */
              f2absq=FVALUE,  /* f2absq=|f2|^2                              */
              f2absqb=FVALUE, /* f2absqb=|P(xb)|^2                          */
              h2abs,          /* h2abs=|h2|                                 */
              epsilon;        /* bound for |q2|                             */
     int      seconditer=0,   /* second iteration, when root is too bad     */
              noise=0,        /* noise counter                              */
              rootd=FALSE;    /* rootd = TRUE  => root determined           */
                              /* rootd = FALSE => no root determined        */
     dcomplex xb;             /* best x-value                               */

                              /* initializing routine                       */
     initialize(pred,&xb,&epsilon); 

     fdvalue(pred,nred,&f0,&f0,x0,FALSE);   /* compute exact function value */
     fdvalue(pred,nred,&f1,&f1,x1,FALSE);   /* oct-29-1993 ml               */
     fdvalue(pred,nred,&f2,&f2,x2,FALSE);


do {                          /* loop for possible second iteration         */
     do {                     /* main iteration loop                        */
                              /* calculate the roots of the parabola        */
          root_of_parabola();
  
                              /* store values for the next iteration        */
          x0 = x1;  
          x1 = x2; 
          h2abs = Cabs(h2);   /* distance between x2 and x1                 */
          
                              /* main iteration-equation                    */
          iteration_equation(&h2abs);

                              /* store values for the next iteration        */ 
          f0 = f1;
          f1 = f2;
          f1absq = f2absq;

                              /* compute P(x2) and make some checks         */
          compute_function(pred,nred,f1absq,&f2absq,epsilon);

          /* printf("betrag %10.5e  %4.2d  %4.2d\n",f2absq,iter,seconditer);  */

                              /* is the new x-value the best approximation? */
          check_x_value(&xb,&f2absqb,&rootd,f1absq,f2absq,
                        epsilon,&noise);

                              /* increase noise counter                     */
         if (fabs((Cabs(xb)-Cabs(x2))/Cabs(xb))<NOISESTART)
              noise++;
     } while ((iter<=ITERMAX) && (!rootd) && (noise<=NOISEMAX));

     seconditer++;            /* increase seconditer                        */

                              /* check, if determined root is good enough   */
     root_check(pred,nred,f2absqb,&seconditer,&rootd,&noise,xb); 
} while (seconditer==2);

     return xb;               /* return best x value                        */
}










/***** initializing routine *****/
void initialize(dcomplex *pred,dcomplex *xb,double *epsilon)
/*dcomplex *pred,     coefficient vector of the deflated polynomial */    
/*         *xb;       best x-value                                  */
/*double   *epsilon;  bound for |q2|                                */
{
     /* initial estimations for x0,...,x2 and its values            */
     /* ml, 12-21-94 changed                                        */

     x0 = Complex(0.,0.);                 /* x0 = 0 + j*1           */ 
     x1 = Complex(-1./sqrt(2),-1./sqrt(2));                /* x1 = 0 - j*1           */
     x2 = Complex(1./sqrt(2),1./sqrt(2)); /* x2 = (1 + j*1)/sqrt(2) */

     h1 = Csub(x1,x0);                         /* h1 = x1 - x0      */
     h2 = Csub(x2,x1);                         /* h2 = x2 - x1      */
     q2 = Cdiv(h2,h1);                         /* q2 = h2/h1        */

     *xb      = x2;    /* best initial x-value = zero   */
     *epsilon = FACTOR*DBL_EPSILON;/* accuracy for determined root  */ 
     iter     = 0;                 /* reset iteration counter       */
}










/***** root_of_parabola() calculate smaller root of Muller's parabola *****/
void root_of_parabola(void)
{
     dcomplex A2,B2,C2,  /* variables to get q2                */
              discr,     /* discriminante                      */
              N1,N2;     /* denominators of q2                 */

                     /* A2 = q2(f2 - (1+q2)f1 + f0q2)          */
                     /* B2 = q2[q2(f0-f1) + 2(f2-f1)] + (f2-f1)*/
                     /* C2 = (1+q2)f[2]                        */
     A2   = Cmul(q2,Csub(Cadd(f2,Cmul(q2,f0)),
                         Cmul(f1,RCadd(1.,q2))));
     B2   = Cadd(Csub(f2,f1),Cmul(q2,Cadd(Cmul(q2,
                         Csub(f0,f1)),RCmul(2.,Csub(f2,f1)))));
     C2   = Cmul(f2,RCadd(1.,q2));
                     /* discr = B2^2 - 4A2C2                   */
     discr = Csub(Cmul(B2,B2),RCmul(4.,Cmul(A2,C2)));
                     /* denominators of q2                     */
     N1 = Csub(B2,Csqrt(discr));  
     N2 = Cadd(B2,Csqrt(discr));  
                 /* choose denominater with largest modulus    */
     if (Cabs(N1)>Cabs(N2) && Cabs(N1)>DBL_EPSILON)
          q2 = Cdiv(RCmul(-2.,C2),N1);  
     else if (Cabs(N2)>DBL_EPSILON)
          q2 = Cdiv(RCmul(-2.,C2),N2);  
     else 
          q2 = Complex(cos(iter),sin(iter));  
}










/***** main iteration equation: x2 = h2*q2 + x2 *****/
void iteration_equation(double *h2abs)
/*double *h2abs;                  Absolute value of the old distance        */
{
     double h2absnew,          /* Absolute value of the new h2              */
            help;              /* help variable                             */

     h2 = Cmul(h2,q2);         
     h2absnew = Cabs(h2);      /* distance between old and new x2           */

     if (h2absnew > (*h2abs*MAXDIST)) { /* maximum relative change          */ 
          help = MAXDIST/h2absnew;
          h2 = RCmul(help,h2);
          q2 = RCmul(help,q2); 
     } 

     *h2abs = h2absnew; /* actualize old distance for next iteration        */

     x2 = Cadd(x2,h2);
}










/**** suppress overflow *****/
void suppress_overflow(int nred)
/*int nred;        the highest exponent of the deflated polynomial        */
{
     int           kiter;            /* internal iteration counter        */
     unsigned char loop;             /* loop = FALSE => terminate loop    */
     double        help;             /* help variable                     */

     kiter = 0;                      /* reset iteration counter           */
     do { 
          loop=FALSE;                /* initial estimation: no overflow   */
          help = Cabs(x2);           /* help = |x2|                       */
          if (help>1. && fabs(nred*log10(help))>BOUND6) {
               kiter++;              /* if |x2|>1 and |x2|^nred>10^BOUND6 */ 
               if (kiter<KITERMAX) { /* then halve the distance between   */
                    h2=RCmul(.5,h2); /* new and old x2                    */
                    q2=RCmul(.5,q2); 
                    x2=Csub(x2,h2);
                    loop=TRUE;
               } else 
                    kiter=0;  
          }
     } while(loop);
}










/***** check of too big function values *****/
void too_big_functionvalues(double *f2absq)
/*double *f2absq;                                   f2absq=|f2|^2          */
{
     if ((fabs(f2.r)+fabs(f2.i))>BOUND4)         /* limit |f2|^2, when     */
          *f2absq = fabs(f2.r)+fabs(f2.i);       /* |f2.r|+|f2.i|>BOUND4   */
     else                                  
          *f2absq = (f2.r)*(f2.r)+(f2.i)*(f2.i); /* |f2|^2 = f2.r^2+f2.i^2 */
}










/***** Muller's modification to improve convergence *****/
void convergence_check(int *overflow,double f1absq,double f2absq,
                       double epsilon)
/*double f1absq,       f1absq = |f1|^2                          */
/*       f2absq,       f2absq = |f2|^2                          */
/*       epsilon;      bound for |q2|                           */
/*int    *overflow;    *overflow = TRUE  => overflow occures    */ 
/*                     *overflow = FALSE => no overflow occures */
{
     if ((f2absq>(CONVERGENCE*f1absq)) && (Cabs(q2)>epsilon) && 
     (iter<ITERMAX)) {
          q2 = RCmul(.5,q2); /* in case of overflow:            */
          h2 = RCmul(.5,h2); /* halve q2 and h2; compute new x2 */
          x2 = Csub(x2,h2);
          *overflow = TRUE;
     }
}










/***** compute P(x2) and make some checks *****/
void compute_function(dcomplex *pred,int nred,double f1absq,
                      double *f2absq,double epsilon)
/*dcomplex *pred;      coefficient vector of the deflated polynomial   */    
/*int      nred;       the highest exponent of the deflated polynomial */
/*double   f1absq,     f1absq = |f1|^2                                 */
/*         *f2absq,    f2absq = |f2|^2                                 */
/*         epsilon;    bound for |q2|                                  */
{
     int    overflow;    /* overflow = TRUE  => overflow occures       */ 
                         /* overflow = FALSE => no overflow occures    */

     do {
         overflow = FALSE; /* initial estimation: no overflow          */ 

                   /* suppress overflow                                */
         suppress_overflow(nred); 

                   /* calculate new value => result in f2              */
         fdvalue(pred,nred,&f2,&f2,x2,FALSE);

                   /* check of too big function values                 */
         too_big_functionvalues(f2absq);

                   /* increase iterationcounter                        */
         iter++;

                   /* Muller's modification to improve convergence     */
         convergence_check(&overflow,f1absq,*f2absq,epsilon);
     } while (overflow);
}










/***** is the new x2 the best approximation? *****/
void check_x_value(dcomplex *xb,double *f2absqb,int *rootd,
                   double f1absq,double f2absq,double epsilon,
                   int *noise)
/*dcomplex *xb;        best x-value                                    */
/*double   *f2absqb,   f2absqb |P(xb)|^2                               */
/*         f1absq,     f1absq = |f1|^2                                 */
/*         f2absq,     f2absq = |f2|^2                                 */
/*         epsilon;    bound for |q2|                                  */
/*int      *rootd,     *rootd = TRUE  => root determined               */
/*                     *rootd = FALSE => no root determined            */
/*         *noise;     noisecounter                                    */
{
     if ((f2absq<=(BOUND1*f1absq)) && (f2absq>=(BOUND2*f1absq))) {
                                  /* function-value changes slowly     */
          if (Cabs(h2)<BOUND3) {  /* if |h[2]| is small enough =>      */
              q2 = RCmul(2.,q2);  /* double q2 and h[2]                */
              h2 = RCmul(2.,h2);     
          } else {                /* otherwise: |q2| = 1 and           */
                                  /*            h[2] = h[2]*q2         */
              q2 = Complex(cos(iter),sin(iter));
              h2 = Cmul(h2,q2);
          }
     } else if (f2absq<*f2absqb) {
          *f2absqb = f2absq;      /* the new function value is the     */
          *xb      = x2;          /* best approximation                */
          *noise   = 0;           /* reset noise counter               */
          if ((sqrt(f2absq)<epsilon) && 
          (Cabs(Cdiv(Csub(x2,x1),x2))<epsilon))
               *rootd = TRUE;     /* root determined                   */
     }
}










/***** check, if determined root is good enough. *****/
void root_check(dcomplex *pred,int nred,double f2absqb,int *seconditer,
                int *rootd,int *noise,dcomplex xb)
/*dcomplex *pred,        coefficient vector of the deflated polynomial      */
/*         xb;           best x-value                                       */
/*int      nred,         the highest exponent of the deflated polynomial    */
/*         *noise,       noisecounter                                       */
/*         *rootd,       *rootd = TRUE  => root determined                  */
/*                       *rootd = FALSE => no root determined               */
/*         *seconditer;  *seconditer = TRUE  => start second iteration with */
/*                                              new initial estimations     */
/*                       *seconditer = FALSE => end routine                 */
/*double   f2absqb;      f2absqb |P(xb)|^2                                  */
{

     dcomplex df;     /* df=P'(x0)                                          */

     if ((*seconditer==1) && (f2absqb>0)) { 
          fdvalue(pred,nred,&f2,&df,xb,TRUE); /* f2=P(x0), df=P'(x0)        */
         if (Cabs(f2)/(Cabs(df)*Cabs(xb))>BOUND7) {
              /* start second iteration with new initial estimations        */
/*              x0 = Complex(-1./sqrt(2),1./sqrt(2)); 
              x1 = Complex(1./sqrt(2),-1./sqrt(2)); 
              x2 = Complex(-1./sqrt(2),-1./sqrt(2)); */
/*ml, 12-21-94: former initial values: */
              x0 = Complex(1.,0.);                 
              x1 = Complex(-1.,0.);                
              x2 = Complex(0.,0.);       /*   */
              fdvalue(pred,nred,&f0,&df,x0,FALSE); /* f0 =  P(x0)           */
              fdvalue(pred,nred,&f1,&df,x1,FALSE); /* f1 =  P(x1)           */
              fdvalue(pred,nred,&f2,&df,x2,FALSE); /* f2 =  P(x2)           */
              iter = 0;                /* reset iteration counter           */
              (*seconditer)++;         /* increase seconditer               */
              *rootd = FALSE;          /* no root determined                */
              *noise = 0;              /* reset noise counter               */
          }
     }
}
