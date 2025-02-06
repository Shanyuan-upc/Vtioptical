/* Copyright (c) China University of Petroleum, 2015.*/
/* All rights reserved.                       */

/* siig_yjd.h - include file for general purpose siig research group stuff */
#ifndef SIIG_H
#define SIIG_H

/* INCLUDES */
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <limits.h>
#include <float.h>
#include <stdarg.h>
#include <errno.h>
#include<time.h>

#include <fcntl.h>      /* non-ANSI */
#include <unistd.h>     /* non-ANSI */
#include <sys/types.h>  /* non-ANSI */

#define NINT(x) ((int)((x)>0.0?(x)+0.5:(x)-0.5))
#define	MAX(x,y) ((x) > (y) ? (x) : (y))
#define	MIN(x,y) ((x) < (y) ? (x) : (y))
#define PI (3.141592653589793)

#ifndef ABS
	#define ABS(x) ((x) < 0 ? -(x) : (x))
#endif


/* complex struct*/
typedef struct _complexStruct {
	float r;
	float i;
} complex;


#ifndef __cplusplus/* if not C++, declare  C complex function*/
	
	/* complex a plus complex b and return a complex structure c=a+b*/
	complex cadd (complex a, complex b);

	/* complex a sub complex b and return a complex structure c=a-b*/
	complex csub (complex a, complex b);

	/* complex a multiply complex b and return a complex structure c=a*b*/
	complex cmul (complex a, complex b);

	/* complex a divide complex b and return a complex structure c=a/b*/
	complex cdiv (complex a, complex b);

	/* cmplx	make a complex number from two real numbers */
	complex cmplx (float re, float im);

	/* conjg	complex conjugate of a complex number */
	complex conjg (complex z);

	/* cneg	negate a complex number */
	complex cneg (complex z);

	/* cinv	invert a complex number*/
	complex cinv (complex z);

	/* csqrt	complex square root of a complex number*/
	//complex csqrt (complex z);

	/* cexp	complex exponential of a complex number*/
	//complex cexp (complex z);

	/* crmul	multiply a complex number by a real number */
	complex crmul (complex a, float x);

	/* rcabs	real magnitude of a complex number*/
	float rcabs (complex z);
#endif /*endif C++ */

#ifndef ALLOC
	#define ALLOC
	/* allocate and free multi-dimensional arrays */
	void *alloc1(int n1,int size);
	void *realloc1(void *v,int n1,int size);
	void **alloc2(int n1,int n2,int size);
	void ***alloc3(int n1,int n2,int n3,int size);
	void free1 (void *p);
	void free2 (void **p);
	void free3 (void ***p);

	int *alloc1int(int n1);
	void free1int(int *p);
	int *realloc1int(int *v, int n1);
	int **alloc2int(int n1, int n2);
	void free2int(int **p);
	int ***alloc3int(int n1, int n2, int n3);
	void free3int(int ***p);

	float *alloc1float(int n1);
	float *realloc1float(float *v, int n1);
	void free1float(float *p);
	float **alloc2float(int n1, int n2);
	void free2float(float **p);
	float ***alloc3float(int n1, int n2, int n3);
	void free3float(float ***p);

	double *alloc1double(int n1);
	double *realloc1double(double *v, int n1);
	void free1double(double *p);
	double **alloc2double(int n1, int n2);
	void free2double(double **p);
	double ***alloc3double(int n1, int n2, int n3);
	void free3double(double ***p);

	complex *alloc1complex(int n1);
	void free1complex(complex *p);
	complex **alloc2complex(int n1, int n2);
	void free2complex(complex **p);
	complex ***alloc3complex(int n1, int n2, int n3);
	void free3complex(complex ***p);

#endif /*ALLOC*/

#ifndef RANDOM
	#define RANDOM
	/* pseudo-random numbers */
	float franuni (void);
	void sranuni (int seed);
#endif/* RANDOM*/

#ifndef SIN_FUN
	#define SIN_FUN

	/* other linear system solvers */
	void stoepd (int n, double r[], double g[], double f[], double a[]);
	void stoepf (int n, float r[], float g[], float f[], float a[]);

	/* interpolation */
	float fsinc (float x);
	double dsinc (double x);
	void mksinc (float d, int lsinc, float sinc[]);
	void ints8r (int nxin, float dxin, float fxin, float yin[], 
		float yinl, float yinr, int nxout, float xout[], float yout[]);
	void ints8c (int nxin, float dxin, float fxin, complex yin[], 
		complex yinl, complex yinr, int nxout, float xout[], complex yout[]);
	void intt8r (int ntable, float table[][8],
		int nxin, float dxin, float fxin, float yin[], 
		float yinl, float yinr, int nxout, float xout[], float yout[]);
	void intt8c (int ntable, float table[][8],
		int nxin, float dxin, float fxin, complex yin[], 
		complex yinl, complex yinr, int nxout, float xout[], complex yout[]);

#endif/* SIN_FUN*/

#ifndef FFT
	#define FFT
	/* Prime Factor FFTs */
	int npfa (int nmin);
	int npfao (int nmin, int nmax);
	int npfar (int nmin);
	int npfaro (int nmin, int nmax);
	void pfacc (int isign, int n, complex z[]);
	void pfarc (int isign, int n, float rz[], complex cz[]);
	void pfacr (int isign, int n, complex cz[], float rz[]);
	void pfa2cc (int isign, int idim, int n1, int n2, complex z[]);
	void pfa2rc (int isign, int idim, int n1, int n2, float rz[], complex cz[]);
	void pfa2cr (int isign, int idim, int n1, int n2, complex cz[], float rz[]);
	void pfamcc (int isign, int n, int nt, int k, int kt, complex z[]);
#endif /* FFT*/

/* exit function */
int SIIG_Exit();



#endif /* SIIG_H */
