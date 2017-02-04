#include "fourier.h"
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <complex>

using namespace std;

bool is_pot2(int n);
bool FFT(int dir,long n,complex<double> *c);
bool DFT(int dir,int m,complex<double> *c);
bool DFT2D(int dir,int n,int m,complex<double> **c);
bool FFT2D(int dir,int n,int m,complex<double> **c);

/***
* APlicar Transformada de Fourier luego DEP
*/
double* Fourier::get_coef_fourier(double* v, int n,int num_coef)
{
    if(v==NULL) return NULL;

    complex<double> *c=new complex<double>[n];
    for(int i=0;i<n;i++)
        c[i]=complex<double>(v[i],0);

    if(is_pot2(n))
        FFT(1,n,c);
    else
        DFT(1,n,c);

    num_coef =  std::min(num_coef, n);
    double* transf = new double[num_coef];

    double fact = 2.0/n; // normalize 
    
    for(int i=0; i< num_coef; i++)
    {
      // Y = magnitud * sin( Y' + phase)
      transf[i] = std::abs(c[i])*fact; // magnitud
      //transf[num_coef/2 + i] = std::arg(c[i]); // phase
    }
 
    delete[] c;

    return transf;
}

double** Fourier::get_coef_fourier2D(double** ma, int n,int m,int num_coef)
{
    if(ma==NULL) return NULL;

    complex<double> **c=new complex<double>*[n];
    for(int i=0;i<n;i++)
    {
        c[i]=new complex<double>[m];
        for(int j=0;j<m;j++)
            c[i][j]=complex<double>(ma[i][j],0);
    }

    if(is_pot2(n) && is_pot2(m))
    {
        FFT2D(1,n,m,c);
    }
    else
    {
        DFT2D(1,n,m,c);
    }


    double** transf=new double*[std::min(num_coef,n)];


    for(int i=0; i<std::min(num_coef,n); i++)
    {
        transf[i]=new double[m];
        for(int j=0;j<m;j++)
            transf[i][j]=abs(c[i][j]);

       delete[] c[i];
    }

    delete[] c;

    return transf;
}

/***************************************************************************************/
/***
* valida n si es multiplo de 2
*/
bool is_pot2(int n)
{
    int i=0;
    double num=0.0;
    bool pot2=false;

    while(num < n )
    {
        num=pow(2,++i);
        pot2= num == n;
    }

    return pot2;
}

/*
    Direct fourier transform
    dir =  1 gives forward transform
    dir = -1 gives reverse transform

     Formula: reverse
                  N-1
                  ---
              1   \           j k 2 pi n / N
      X(n) = ---   >   x(k) e                    = reverse transform
              N   /                                n=0..N-1
                  ---
                  k=0

      Formula: forward
                  N-1
                  ---
                  \          - j k 2 pi n / N
      X(n) =       >   x(k) e                    = forward transform
                  /                                n=0..N-1
                  ---
                  k=0
*/
bool DFT(int dir,int m,complex<double> *c)
{
   long i,k;
   double arg;
   double cosarg,sinarg;
   double *x2=NULL,*y2=NULL;

   x2 = new double[m];
   y2 = new double[m];

   for (i=0;i<m;i++)
   {
      x2[i] = 0;
      y2[i] = 0;
      arg = - dir * 2.0 * 3.141592654 * (double)i / (double)m;
      for (k=0;k<m;k++) {
         cosarg = cos(k * arg);
         sinarg = sin(k * arg);
         x2[i] += (real(c[k]) * cosarg - imag(c[k]) * sinarg);
         y2[i] += (real(c[k]) * sinarg + imag(c[k]) * cosarg);
      }
   }

   /* Copy the data back */
   if (dir == -1) { // reverse (IDFT)
      for (i=0;i<m;i++)
      {
         c[i]=complex<double>( x2[i]/(double)m, y2[i]/(double)m );
      }
   } else {
      for (i=0;i<m;i++)
      {
         c[i]=complex<double>(x2[i], y2[i]);
      }
   }

   delete[] x2;
   delete[] y2;

   return true;
}

/*
   This computes an in-place complex-to-complex FFT
   x and y are the real and imaginary arrays of n=2^m points.
*/
bool FFT(int dir,long n,complex<double> *c)
{
   long m,i,i1,j,k,i2,l,l1,l2;
   double c1,c2,t1,t2,u1,u2,z;
   complex<double> t;

   /* Calculate m: pow(2,m)=n */
   m=1;
   while(pow(2,m)!=n) m++;

   /* Do the bit reversal */
   i2 = n >> 1;
   j = 0;
   for (i=0;i<n-1;i++) {
      if (i < j) {
         t=c[i];
         c[i]=c[j];
         c[j]=t;
      }
      k = i2;
      while (k <= j) {
         j -= k;
         k >>= 1;
      }
      j += k;
   }

   /* Compute the FFT */
   c1 = -1.0;
   c2 = 0.0;
   l2 = 1;
   for (l=0;l<m;l++) {
      l1 = l2;
      l2 <<= 1;
      u1 = 1.0;
      u2 = 0.0;
      for (j=0;j<l1;j++) {
         for (i=j;i<n;i+=l2) {
            i1 = i + l1;
            t1 = u1 * real(c[i1]) - u2 * imag(c[i1]);
            t2 = u1 * imag(c[i1]) + u2 * real(c[i1]);
            t=complex<double>(t1,t2);
            c[i1]=c[i]-t;
            c[i]+=t;
         }
         z =  u1 * c1 - u2 * c2;
         u2 = u1 * c2 + u2 * c1;
         u1 = z;
      }
      c2 = sqrt((1.0 - c1) / 2.0);
      if (dir == 1)
         c2 = -c2;
      c1 = sqrt((1.0 + c1) / 2.0);
   }

   /* Scaling for forward transform */
   t=complex<double>(n,0);
   if (dir == 1) {
      for (i=0;i<n;i++) {
         c[i]/=t;
      }
   }

   return true;
}

/*-------------------------------------------------------------------------
   Perform a 2D DFT inplace given a complex 2D array
   The direction dir, 1 for forward, -1 for reverse
   The size of the array (nx,ny)
   Return false if there are memory problems
*/
bool DFT2D(int dir,int n,int m,complex<double> **c)
{
   long i,j,u,v;
   double arg1,arg2;
   double cosarg,sinarg;
   double **x2=NULL,**y2=NULL;

   x2 = new double*[n];
   y2 = new double*[n];
   if (x2 == NULL || y2 == NULL)
      return false;

   for (i=0;i<n;i++)
   {
	   x2[i]=new double[m];
	   y2[i]=new double[m];
   }

   for(u=0;u<n;u++)
	for(v=0;v<m;v++)
	{
		x2[u][v] = 0;
		y2[u][v] = 0;
		arg1 = - dir * 2.0 * 3.141592654;
		for(i=0;i<n;i++)
			for(j=0;j<m;j++)
			{
				arg2=arg1*(u*i*1.0/(double)n + v*j*1.0/(double)m);
				cosarg = cos(arg2);
				sinarg = sin(arg2);
				x2[u][v] += (real(c[i][j]) * cosarg - imag(c[i][j]) * sinarg);
				y2[u][v] += (real(c[i][j]) * sinarg + imag(c[i][j]) * cosarg);
			}
	}

   /* Copy the data back */
   if (dir==1) {
      for (i=0;i<n;i++)
		for (j=0;j<m;j++)
		  {
			 c[i][j]=complex<double>( x2[i][j]/(double)(n*m), y2[i][j]/(double)(n*m) );
		  }
   } else {
      for (i=0;i<n;i++)
		for (j=0;j<m;j++)
	    {
			c[i][j]=complex<double>(x2[i][j], y2[i][j]);
		}
   }


   for (i=0;i<n;i++)
   {
	   delete[] x2[i];
	   delete[] y2[i];
   }

   delete[] x2;
   delete[] y2;

   return true;
}

/*-------------------------------------------------------------------------
   Perform a 2D FFT inplace given a complex 2D array
   The direction dir, 1 for forward, -1 for reverse
   The size of the array (nx,ny)
   Return false if there are memory problems or
      the dimensions are not powers of 2
*/
bool FFT2D(int dir,int n,int m,complex<double> **c)
{
   int i,j;
   complex<double> *row,*col;

   if (!is_pot2(n) || !is_pot2(m))
      return false;

   /* Transform the rows */
   row = new complex<double>[n];
   col = new complex<double>[m];
   if (row == NULL || col == NULL)
      return false;


   for (j=0;j<m;j++) {
      for (i=0;i<n;i++) {
         row[i] = c[i][j];
      }
      FFT(dir,n,row);
      for (i=0;i<n;i++) {
         c[i][j] = row[i];
      }
   }

   /* Transform the columns */
   for (i=0;i<n;i++) {
      for (j=0;j<m;j++) {
         col[j] = c[i][j];
      }
      FFT(dir,m,col);
      for (j=0;j<m;j++) {
         c[i][j] = col[j];
      }
   }

   delete[] row;
   delete[] col;

   return true;
}
