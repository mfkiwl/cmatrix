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

extern int CMatrixd_Print(double* matrix, int row, int column) {
	printf(">>Matrix_%p:\n", matrix);
	int i, j;
	for (i = 0; i < row; i++) {
		for (j = 0; j < column; j++) {
			printf("%.3lf\t", matrix[i * column + j]);
		}
		printf("\n");
	}
	return 0;
}
extern int CMatrixf_Print(float* matrix, int row, int column) {
	printf(">>Matrix_%p:\n", matrix);
	int i, j;
	for (i = 0; i < row; i++) {
		for (j = 0; j < column; j++) {
			printf("%.3f\t", matrix[i * column + j]);
		}
		printf("\n");
	}
	return 0;
}
extern int CMatrixi_Print(int* matrix, int row, int column) {
	printf(">>Matrix_%p:\n", matrix);
	int i, j;
	for (i = 0; i < row; i++) {
		for (j = 0; j < column; j++) {
			printf("%d\t", matrix[i * column + j]);
		}
		printf("\n");
	}
	return 0;
}

extern void CMatrixd_copy(double* mat_dest, const double* mat_src, int row, int col) {
	memcpy(mat_dest, mat_src, sizeof(double) * row * col);
}
extern void CMatrixf_copy(float* mat_dest, const float* mat_src, int row, int col) {
	memcpy(mat_dest, mat_src, sizeof(float) * row * col);
}
extern void CMatrixi_copy(int* mat_dest, const int* mat_src, int row, int col) {
	memcpy(mat_dest, mat_src, sizeof(int) * row * col);
}

/* multiply matrix (wrapper of blas dgemm) -------------------------------------
* multiply matrix by matrix (C=alpha*A*B+beta*C)
* args   : char   *tr       I  transpose flags ("N":normal,"T":transpose)
*          int    n,k,m     I  size of (transposed) matrix A,B
*          double alpha     I  alpha
*          double *A,*B     I  (transposed) matrix A (n x m), B (m x k)
*          double beta      I  beta
*          double *C        IO matrix C (n x k)
* return : none
*-----------------------------------------------------------------------------*/
extern void matmul(const char* transpose, int left_mat_row, int l_row_r_col, int right_mat_col, double alpha,
	const double* left_matrix, const double* right_matrix, double beta, double* result_matrix)
{
	double d;
	int i, j, x, f = transpose[0] == 'N' ? (transpose[1] == 'N' ? 1 : 2) : (transpose[1] == 'N' ? 3 : 4);

	for (i = 0; i < left_mat_row; i++) for (j = 0; j < right_mat_col; j++) {
		d = 0.0;
		switch (f) {
		case 1: for (x = 0; x < l_row_r_col; x++) d += left_matrix[i + x * left_mat_row] * right_matrix[x + j * l_row_r_col]; break;
		case 2: for (x = 0; x < l_row_r_col; x++) d += left_matrix[i + x * left_mat_row] * right_matrix[j + x * right_mat_col]; break;
		case 3: for (x = 0; x < l_row_r_col; x++) d += left_matrix[x + i * l_row_r_col] * right_matrix[x + j * l_row_r_col]; break;
		case 4: for (x = 0; x < l_row_r_col; x++) d += left_matrix[x + i * l_row_r_col] * right_matrix[j + x * right_mat_col]; break;
		}
		if (beta == 0.0) result_matrix[i + j * left_mat_row] = alpha * d; 
		else result_matrix[i + j * left_mat_row] = alpha * d + beta * result_matrix[i + j * left_mat_row];
	}
}
