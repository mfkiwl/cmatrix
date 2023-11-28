#include<stdio.h>
#include<stdlib.h>
#include"cmatrix.h"
/* new matrix ------------------------------------------------------------------
* allocate memory of matrix
* args   : int    row,column       I   number of rows and columns of matrix
* return : matrix pointer (if row<=0 or column<=0, return NULL)
*-----------------------------------------------------------------------------*/
extern double* CMatrixd(int row, int column) {
	double* p;
	if (row <= 0 || column <= 0) return;
	if (!(p = (double*)malloc(sizeof(double) * row * column))) {
		printf("matrix memory allocation error: n=%d,m=%d\n", row, column);
		return NULL;
	}
	return p;
}
extern float* CMatrixf(int row, int column) {
	float* p;
	if (row <= 0 || column <= 0) return;
	if (!(p = (float*)malloc(sizeof(float) * row * column))) {
		printf("matrix memory allocation error: n=%d,m=%d\n", row, column);
		return NULL;
	}
	return p;
}
extern int* CMatrixi(int row, int column) {
	int* p;
	if (row <= 0 || column <= 0) return;
	if (!(p = (int*)malloc(sizeof(int) * row * column))) {
		printf("matrix memory allocation error: n=%d,m=%d\n", row, column);
		return NULL;
	}
	return p;
}

/* zero matrix -----------------------------------------------------------------
* generate new zero matrix
* args   : int    n,m       I   number of rows and columns of matrix
* return : matrix pointer (if n<=0 or m<=0, return NULL)
*-----------------------------------------------------------------------------*/
extern double* CMatrixd_Zero(int row, int column)
{
	double* p;

#if NOCALLOC
	if ((p = CMatrixd(row, column))) for (row = row * column - 1; row >= 0; row--) p[row] = 0.0;
#else
	if (row <= 0 || column <= 0) return NULL;
	if (!(p = (double*)calloc(sizeof(double), row * column))) {
		printf("matrix memory allocation error: n=%d,m=%d\n", row, column);
		return NULL;
	}
#endif
	return p;
}
extern float* CMatrixf_Zero(int row, int column)
{
	float* p;

#if NOCALLOC
	if ((p = CMatrixf(row, column))) for (row = row * column - 1; row >= 0; row--) p[row] = 0.0;
#else
	if (row <= 0 || column <= 0) return NULL;
	if (!(p = (float*)calloc(sizeof(float), row * column))) {
		printf("matrix memory allocation error: n=%d,m=%d\n", row, column);
		return NULL;
	}
#endif
	return p;
}
extern int* CMatrixi_Zero(int row, int column)
{
	int* p;

#if NOCALLOC
	if ((p = CMatrixi(row, column))) for (row = row * column - 1; row >= 0; row--) p[row] = 0.0;
#else
	if (row <= 0 || column <= 0) return NULL;
	if (!(p = (int*)calloc(sizeof(int), row * column))) {
		printf("matrix memory allocation error: n=%d,m=%d\n", row, column);
		return NULL;
	}
#endif
	return p;
}

/* identity matrix -------------------------------------------------------------
* generate new identity matrix
* args   : int    n         I   number of rows and columns of matrix
* return : matrix pointer (if n<=0, return NULL)
*-----------------------------------------------------------------------------*/
extern double* CMatrixd_Eye(int order)
{
	double* p;
	int i;

	if ((p = CMatrixd_Zero(order, order))) {
		for (i = 0; i < order; i++) p[i + i * order] = 1.0;
		return p;
	}
	else {
		printf("matrix memory allocation error: n=%d,m=%d\n", order, order);
		return NULL;
	}
}
extern float* CMatrixf_Eye(int order)
{
	float* p;
	int i;

	if ((p = CMatrixf_Zero(order, order))) {
		for (i = 0; i < order; i++) p[i + i * order] = 1.0;
		return p;
	}
	else {
		printf("matrix memory allocation error: n=%d,m=%d\n", order, order);
		return NULL;
	}
}
extern int* CMatrixi_Eye(int order)
{
	int* p;
	int i;

	if ((p = CMatrixi_Zero(order, order))) {
		for (i = 0; i < order; i++) p[i + i * order] = 1.0;
		return p;
	}
	else {
		printf("matrix memory allocation error: n=%d,m=%d\n", order, order);
		return NULL;
	}
}

/* inner product ---------------------------------------------------------------
* inner product of vectors
* args   : double *a,*b     I   vector a,b (n x 1)
*          int    n         I   size of vector a,b
* return : a'*b
*-----------------------------------------------------------------------------*/
extern double Dotd(const double* vector1, const double* vector2, int order)
{
	double result = 0.0;
	while (--order >= 0) result += vector1[order] * vector2[order];
	return result;
}
extern float Dotf(const float* vector1, const float* vector2, int order)
{
	float result = 0.0;
	while (--order >= 0) result += vector1[order] * vector2[order];
	return result;
}
extern int Doti(const int* vector1, const int* vector2, int order)
{
	int result = 0.0;
	while (--order >= 0) result += vector1[order] * vector2[order];
	return result;
}

/* euclid norm -----------------------------------------------------------------
* euclid norm of vector
* args   : double *a        I   vector a (n x 1)
*          int    n         I   size of vector a
* return : || a ||
*-----------------------------------------------------------------------------*/
extern double Normd(const double* vector, int order)
{
	return sqrt(Dotd(vector, vector, order));
}
extern float Normf(const float* vector, int order)
{
	return sqrt(Dotd(vector, vector, order));
}
extern int Normi(const int* vector, int order)
{
	return sqrt(Dotd(vector, vector, order));
}

/* outer product of 3d vectors -------------------------------------------------
* outer product of 3d vectors
* args   : double *a,*b     I   vector a,b (3 x 1)
*          double *c        O   outer product (a x b) (3 x 1)
* return : none
*-----------------------------------------------------------------------------*/
extern void Cross3d(const double* vector1, const double* vector2, double* vector_result)
{
	vector_result[0] = vector1[1] * vector2[2] - vector1[2] * vector2[1];
	vector_result[1] = vector1[2] * vector2[0] - vector1[0] * vector2[2];
	vector_result[2] = vector1[0] * vector2[1] - vector1[1] * vector2[0];
}
extern void Cross3f(const float* vector1, const float* vector2, float* vector_result)
{
	vector_result[0] = vector1[1] * vector2[2] - vector1[2] * vector2[1];
	vector_result[1] = vector1[2] * vector2[0] - vector1[0] * vector2[2];
	vector_result[2] = vector1[0] * vector2[1] - vector1[1] * vector2[0];
}
extern void Cross3i(const int* vector1, const int* vector2, int* vector_result)
{
	vector_result[0] = vector1[1] * vector2[2] - vector1[2] * vector2[1];
	vector_result[1] = vector1[2] * vector2[0] - vector1[0] * vector2[2];
	vector_result[2] = vector1[0] * vector2[1] - vector1[1] * vector2[0];
}