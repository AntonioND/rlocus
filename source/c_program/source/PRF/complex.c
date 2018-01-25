/******************************************************************/
/*                                                                */
/* file:          complex.c                                       */
/*                                                                */
/* source:        Numerical Recipies in C/The Art of Scientific   */
/*                Computing; William H. Press et al.;             */
/*                Cambridge 1990; Cambridge University Press      */
/*                                                                */
/* function:      description:                                    */
/*                                                                */
/* Cadd(a,b)      add complex numbers a and b and return result   */
/* Csub(a,b)      subtract complex a and b and return result      */
/* Cmul(a,b)      multiply complex a and b and return result      */
/* Cdiv(a,b)      divide complex a by b and return result         */
/* Complex(a,b)   create complex number with real part a and      */
/*                imaginary part b and return result              */
/* Cabs(z)        compute absolute value of complex z and return  */
/*                result                                          */
/* Conjg(z)       compute conjugate complex value of z and return */
/*                result                                          */
/* Csqrt(z)       compute square root of complex z and return res.*/
/* RCmul(x,a)     multiply real x with complex a and return res.  */
/* RCadd(x,a)     add real x and complex a and return result      */
/*                                                                */
/* further functions:                                             */
/* author:        B. Frenzel                                      */
/* date:          Jan 15 1993                                     */
/* function:      description:                                    */
/*                                                                */
/* Carg(z)        compute argument of complex z and return result */
/*                                                                */
/* Copyright:                                                     */
/* Lehrstuhl fuer Nachrichtentechnik Erlangen                     */
/* Cauerstr. 7, 8520 Erlangen, FRG, 1993                          */
/* e-mail: int@nt.e-technik.uni-erlangen.de                       */
/*                                                                */
/******************************************************************/
#define COMPL
#include "header.h"









/***** Cadd(a,b)add complex numbers a and b and return result *****/
dcomplex Cadd(dcomplex a,dcomplex b)
{
     dcomplex c;

     c.r = a.r + b.r;
     c.i = a.i + b.i;
     return c;
}





/***** Csub(a,b) subtract complex a and b and return result *****/
dcomplex Csub(dcomplex a,dcomplex b)
{
     dcomplex c;

     c.r = a.r - b.r;
     c.i = a.i - b.i;
     return c;
}





/***** Cmul(a,b) multiply complex a and b and return result *****/
dcomplex Cmul(dcomplex a,dcomplex b)
{
     dcomplex c;

     c.r = a.r * b.r - a.i * b.i;
     c.i = a.i * b.r + a.r * b.i;
     return c;
}





/***** Cdiv(a,b) divide complex a by b and return result *****/
dcomplex Cdiv(dcomplex a,dcomplex b)
{
      dcomplex c;
      double r,den;

      if (fabs(b.r) >= fabs(b.i)) {
           r   = b.i/b.r;
           den = b.r+r*b.i;
           c.r = (a.r+r*a.i)/den;
           c.i = (a.i-r*a.r)/den;
      } else {
           r   = b.r/b.i;
           den = b.i + r*b.r;
           c.r = (a.r*r+a.i)/den;
           c.i = (a.i*r-a.r)/den;
      }
      return c;
}





/***** Complex(a,b) create complex number with real part a and *****/
/***** imaginary part b and return result                      *****/
dcomplex Complex(double a,double b)
{
     dcomplex c;

     c.r = a;
     c.i = b;
     return c;
}





/***** Cabs(z) compute absolute value of complex z *****/
double Cabs(dcomplex z)
{
     double x,y,ans,temp;

     x = fabs(z.r);
     y = fabs(z.i);
     if (x == 0.0)
          ans=y;
     else if (y == 0.0)
          ans=x;
     else if (x>y) {
          temp=y/x;
          ans=x*sqrt(1.0+temp*temp);
     } else {
          temp=x/y;
          ans=y*sqrt(1.0+temp*temp);
     }
     return ans;
}





/***** Carg(z) compute argument of complex z *****/
double Carg(dcomplex z)
{
     if (z.r==0.)
     {
          if (z.i>0)
               return PI/2.;
          else if (z.i<0)
               return (-PI/2.);
          else
               return 0;
     }

     return atan2(z.i,z.r);
}





/***** Conjg(z) compute conjugate complex value of z *****/
dcomplex Conjg(dcomplex z)
{
     dcomplex c;

     c.r = z.r;
     c.i = -z.i;
     return c;
}





/***** Csqrt(z) compute square root of complex z *****/
dcomplex Csqrt(dcomplex z)
{
     dcomplex c;
     double x,y,w,r;

     if ((z.r == 0.0) && (z.i == 0.0)) {
          c.r = c.i = 0.0;
          return c;
     } else {
          x = fabs(z.r);
          y = fabs(z.i);
          if (x >= y) {
               r = y/x;
               w = sqrt(x)*sqrt(0.5*(1.0+sqrt(1.0+r*r)));
          } else {
               r = x/y;
               w = sqrt(y)*sqrt(0.5*(r+sqrt(1.0+r*r)));
          }
          if (z.r >= 0.0) {
               c.r = w;
               c.i = z.i/(2.0*w);
          } else {
               c.i = (z.i >= 0) ? w : -w;
               c.r = z.i/(2.0*c.i);
          }
          return c;
     }
}





/***** RCmul(x,a) multiply real x with complex a *****/
dcomplex RCmul(double x,dcomplex a)
{
     dcomplex c;

     c.r = x*a.r;
     c.i = x*a.i;
     return c;
}





/***** RCadd(x,a) add real x and complex a and return result *****/
dcomplex RCadd(double x,dcomplex a)
{
     dcomplex c;

     c.r = x+a.r;
     c.i = a.i;
     return c;
}
