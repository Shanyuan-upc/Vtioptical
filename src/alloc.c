/* Copyright (c) China University of Petroelum, 2015.*/
/* All rights reserved.                       */

/*********************** self documentation **********************/
/*****************************************************************************
ALLOC - Allocate and free multi-dimensional arrays
void *alloc1(int n1,int size)
void *realloc1(void *v,int n1,int size)
void **alloc2(int n1,int n2,int size)
void ***alloc3(int n1,int n2,int n3,int size)
void free1 (void *p)
void free2 (void **p)
void free3 (void ***p)
int *alloc1int(int n1)
void free1int(int *p)
int *realloc1int(int *v, int n1)
int **alloc2int(int n1, int n2)
void free2int(int **p)
int ***alloc3int(int n1, int n2, int n3)
void free3int(int ***p)
float *alloc1float(int n1)
float *realloc1float(float *v, int n1)
void free1float(float *p)
float **alloc2float(int n1, int n2)
void free2float(float **p)
float ***alloc3float(int n1, int n2, int n3)
void free3float(float ***p)
double *alloc1double(int n1)
double *realloc1double(double *v, int n1)
void free1double(double *p)
double **alloc2double(int n1, int n2)
void free2double(double **p)
double ***alloc3double(int n1, int n2, int n3)
void free3double(double ***p)
complex *alloc1complex(int n1)
void free1complex(complex *p)
complex **alloc2complex(int n1, int n2)
void free2complex(complex **p)
complex ***alloc3complex(int n1, int n2, int n3)
void free3complex(complex ***p)

******************************************************************************
Author:    	Jidong Yang, SWPI research Group, /25/11/2015     
*****************************************************************************/
/**************** end self doc ********************************/

#include "siig_yjd.h"

/* allocate a 1-d array */
void *alloc1 (int n1, int size)
{
	void *p;

	if ((p=malloc(n1*size))==NULL)
		return NULL;
	return p;
}

/* re-allocate a 1-d array */
void *realloc1(void *v, int n1, int size)
{
	void *p;

	if ((p=realloc(v,n1*size))==NULL)
		return NULL;
	return p;
}

/* free a 1-d array */
void free1 (void *p)
{
	free(p);
}

/* allocate a 2-d array */
void **alloc2 (int n1, int n2, int size)
{
	int i2;
	void **p;

	if ((p=(void**)malloc(n2*sizeof(void*)))==NULL) 
		return NULL;
	if ((p[0]=(void*)malloc(n2*n1*size))==NULL) {
		free(p);
		return NULL;
	}
	for (i2=0; i2<n2; i2++)
		p[i2] = (char*)p[0]+size*n1*i2;
	return p;
}

/* free a 2-d array */
void free2 (void **p)
{
	free(p[0]);
	free(p);
}

/* allocate a 3-d array */
void ***alloc3 (int n1, int n2, int n3, int size)
{
// int mat3d1[3][2][2];
//  memset(mat3d1, 0, Rows * Cols * Dims * sizeof(int));
//
// double ***mat3d2 = NULL;
// mat3d2 = (double***)malloc(Dims * sizeof(double**));
// for (int d = 0; d < Dims; d++)
// 	mat3d2[d] = (double**)malloc(Rows * sizeof(double*));
// 	for (int d = 0; d < Dims; d++)
// 	for (int r = 0; r < Rows; r++){
// 		mat3d2[d][r] = (double*)malloc(Cols * sizeof(double));	
// 			memset(mat3d2[d][r], 0, Cols * sizeof(double));
// 			}
// 			// 查看略
//
//#if 0
	int i3,i2;
	void ***p;

	if ((p=(void***)malloc(n3*sizeof(void**)))==NULL)
		return NULL;
	if ((p[0]=(void**)malloc(n3*n2*sizeof(void*)))==NULL) {
		free(p);
		return NULL;
	}
	if ((p[0][0]=(void*)malloc(n3*n2*n1*size))==NULL) {
		free(p[0]);
		free(p);
		return NULL;
	}

	for (i3=0; i3<n3; i3++) {
		p[i3] = p[0]+n2*i3;
		for (i2=0; i2<n2; i2++)
			p[i3][i2] = (char*)p[0][0]+size*n1*(i2+n2*i3);
	}
	return p;
//#endif
}

/* free a 3-d array */
void free3 (void ***p)
{
	free(p[0][0]);
	free(p[0]);
	free(p);
}


/* allocate a 1-d array of ints */
int *alloc1int(int n1)
{
	return (int*)alloc1(n1,sizeof(int));
}

/* re-allocate a 1-d array of ints */
int *realloc1int(int *v, int n1)
{
	return (int*)realloc1(v,n1,sizeof(int));
}

/* free a 1-d array of ints */
void free1int(int *p)
{
	free1(p);
}

/* allocate a 2-d array of ints */
int **alloc2int(int n1, int n2)
{
	return (int**)alloc2(n1,n2,sizeof(int));
}

/* free a 2-d array of ints */
void free2int(int **p)
{
	free2((void**)p);
}

/* allocate a 3-d array of ints */
int ***alloc3int(int n1, int n2, int n3)
{
	return (int***)alloc3(n1,n2,n3,sizeof(int));
}

/* free a 3-d array of ints */
void free3int(int ***p)
{
	free3((void***)p);
}

/* allocate a 1-d array of floats */
float *alloc1float(int n1)
{
	return (float*)alloc1(n1,sizeof(float));
}

/* re-allocate a 1-d array of floats */
float *realloc1float(float *v, int n1)
{
	return (float*)realloc1(v,n1,sizeof(float));
}

/* free a 1-d array of floats */
void free1float(float *p)
{
	free1(p);
}

/* allocate a 2-d array of floats */
float **alloc2float(int n1, int n2)
{
	return (float**)alloc2(n1,n2,sizeof(float));
}

/* free a 2-d array of floats */
void free2float(float **p)
{
	free2((void**)p);
}

/* allocate a 3-d array of floats */
float ***alloc3float(int n1, int n2, int n3)
{
	return (float***)alloc3(n1,n2,n3,sizeof(float));
}

/* free a 3-d array of floats */
void free3float(float ***p)
{
	free3((void***)p);
}

/* allocate a 1-d array of doubles */
double *alloc1double(int n1)
{
	return (double*)alloc1(n1,sizeof(double));
}

/* re-allocate a 1-d array of doubles */
double *realloc1double(double *v, int n1)
{
	return (double*)realloc1(v,n1,sizeof(double));
}


/* free a 1-d array of doubles */
void free1double(double *p)
{
	free1(p);
}

/* allocate a 2-d array of doubles */
double **alloc2double(int n1, int n2)
{
	return (double**)alloc2(n1,n2,sizeof(double));
}

/* free a 2-d array of doubles */
void free2double(double **p)
{
	free2((void**)p);
}

/* allocate a 3-d array of doubles */
double ***alloc3double(int n1, int n2, int n3)
{
	return (double***)alloc3(n1,n2,n3,sizeof(double));
}

/* free a 3-d array of doubles */
void free3double(double ***p)
{
	free3((void***)p);
}

/* allocate a 1-d array of complexs */
complex *alloc1complex(int n1)
{
	return (complex*)alloc1(n1,sizeof(complex));
}

/* re-allocate a 1-d array of complexs */
complex *realloc1complex(complex *v, int n1)
{
	return (complex*)realloc1(v,n1,sizeof(complex));
}

/* free a 1-d array of complexs */
void free1complex(complex *p)
{
	free1(p);
}

/* allocate a 2-d array of complexs */
complex **alloc2complex(int n1, int n2)
{
	return (complex**)alloc2(n1,n2,sizeof(complex));
}

/* free a 2-d array of complexs */
void free2complex(complex **p)
{
	free2((void**)p);
}

/* allocate a 3-d array of complexs */
complex ***alloc3complex(int n1, int n2, int n3)
{
	return (complex***)alloc3(n1,n2,n3,sizeof(complex));
}

/* free a 3-d array of complexs */
void free3complex(complex ***p)
{
	free3((void***)p);
}
