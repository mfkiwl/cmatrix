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

extern void CMatrixd_Copy(double* mat_dest, const double* mat_src, int row, int col) {
	memcpy(mat_dest, mat_src, sizeof(double) * row * col);
}
extern void CMatrixf_Copy(float* mat_dest, const float* mat_src, int row, int col) {
	memcpy(mat_dest, mat_src, sizeof(float) * row * col);
}
extern void CMatrixi_Copy(int* mat_dest, const int* mat_src, int row, int col) {
	memcpy(mat_dest, mat_src, sizeof(int) * row * col);
}

/**
 * @brief ¾ØÕó³Ë·¨(C=alpha4lr*left_matrix(transpose?)*right_matrix(transpose?)+beta4re*C).
 * 
 * @param transpose		[in] ³ËÖ®Ç°ÊÇ·ñ×ªÖÃ
 * @param left_matrix		[in] ×ó¾ØÕó
 * @param right_matrix		[in] ÓÒ¾ØÕó
 * @param left_mat_row		[in] ×ó¾ØÕóÐÐ
 * @param l_row_r_col		[in] ×ó¾ØÕóÁÐ£¬ÓÒ¾ØÕóÐÐ
 * @param right_mat_col	[in] ÓÒ¾ØÕóÁÐ
 * @param alpha4lr			[in] ´ý³Ë¾ØÕóÏµÊý
 * @param beta4re			[in] ½á¹û¾ØÕóÏµÊý
 * @param result_matrix	[in/out] ½á¹û¾ØÕó
 */
extern void CMatrixd_Mul(const char* transpose, const double* left_matrix, const double* right_matrix, 
	int left_mat_row, int l_row_r_col, int right_mat_col, double alpha4lr, double beta4re, double* result_matrix)
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
		if (beta4re == 0.0|| beta4re == 1.0) result_matrix[i + j * left_mat_row] = alpha4lr * d;
		else result_matrix[i + j * left_mat_row] = alpha4lr * d + beta4re * result_matrix[i + j * left_mat_row];
	}
}
extern void CMatrixf_Mul(const char* transpose, const float* left_matrix, const float* right_matrix,
	int left_mat_row, int l_row_r_col, int right_mat_col, float alpha4lr, float beta4re, float* result_matrix)
{
	float d;
	int i, j, x, f = transpose[0] == 'N' ? (transpose[1] == 'N' ? 1 : 2) : (transpose[1] == 'N' ? 3 : 4);

	for (i = 0; i < left_mat_row; i++) for (j = 0; j < right_mat_col; j++) {
		d = 0.0;
		switch (f) {
		case 1: for (x = 0; x < l_row_r_col; x++) d += left_matrix[i + x * left_mat_row] * right_matrix[x + j * l_row_r_col]; break;
		case 2: for (x = 0; x < l_row_r_col; x++) d += left_matrix[i + x * left_mat_row] * right_matrix[j + x * right_mat_col]; break;
		case 3: for (x = 0; x < l_row_r_col; x++) d += left_matrix[x + i * l_row_r_col] * right_matrix[x + j * l_row_r_col]; break;
		case 4: for (x = 0; x < l_row_r_col; x++) d += left_matrix[x + i * l_row_r_col] * right_matrix[j + x * right_mat_col]; break;
		}
		if (beta4re == 0.0 || beta4re == 1.0) result_matrix[i + j * left_mat_row] = alpha4lr * d;
		else result_matrix[i + j * left_mat_row] = alpha4lr * d + beta4re * result_matrix[i + j * left_mat_row];
	}
}
extern void CMatrixi_Mul(const char* transpose, const int* left_matrix, const int* right_matrix,
	int left_mat_row, int l_row_r_col, int right_mat_col, int alpha4lr, int beta4re, int* result_matrix)
{
	int d;
	int i, j, x, f = transpose[0] == 'N' ? (transpose[1] == 'N' ? 1 : 2) : (transpose[1] == 'N' ? 3 : 4);

	for (i = 0; i < left_mat_row; i++) for (j = 0; j < right_mat_col; j++) {
		d = 0.0;
		switch (f) {
		case 1: for (x = 0; x < l_row_r_col; x++) d += left_matrix[i + x * left_mat_row] * right_matrix[x + j * l_row_r_col]; break;
		case 2: for (x = 0; x < l_row_r_col; x++) d += left_matrix[i + x * left_mat_row] * right_matrix[j + x * right_mat_col]; break;
		case 3: for (x = 0; x < l_row_r_col; x++) d += left_matrix[x + i * l_row_r_col] * right_matrix[x + j * l_row_r_col]; break;
		case 4: for (x = 0; x < l_row_r_col; x++) d += left_matrix[x + i * l_row_r_col] * right_matrix[j + x * right_mat_col]; break;
		}
		if (beta4re == 0 || beta4re == 1) result_matrix[i + j * left_mat_row] = alpha4lr * d;
		else result_matrix[i + j * left_mat_row] = alpha4lr * d + beta4re * result_matrix[i + j * left_mat_row];
	}
}



/* LU decomposition ----------------------------------------------------------*/
static int ludcmp(double* A, int n, int* indx, double* d)
{
	double big, s, tmp, * vv = mat(n, 1);
	int i, imax = 0, j, k;

	*d = 1.0;
	for (i = 0; i < n; i++) {
		big = 0.0; for (j = 0; j < n; j++) if ((tmp = fabs(A[i + j * n])) > big) big = tmp;
		if (big > 0.0) vv[i] = 1.0 / big; else { free(vv); return -1; }
	}
	for (j = 0; j < n; j++) {
		for (i = 0; i < j; i++) {
			s = A[i + j * n]; for (k = 0; k < i; k++) s -= A[i + k * n] * A[k + j * n]; A[i + j * n] = s;
		}
		big = 0.0;
		for (i = j; i < n; i++) {
			s = A[i + j * n]; for (k = 0; k < j; k++) s -= A[i + k * n] * A[k + j * n]; A[i + j * n] = s;
			if ((tmp = vv[i] * fabs(s)) >= big) { big = tmp; imax = i; }
		}
		if (j != imax) {
			for (k = 0; k < n; k++) {
				tmp = A[imax + k * n]; A[imax + k * n] = A[j + k * n]; A[j + k * n] = tmp;
			}
			*d = -(*d); vv[imax] = vv[j];
		}
		indx[j] = imax;
		if (A[j + j * n] == 0.0) { free(vv); return -1; }
		if (j != n - 1) {
			tmp = 1.0 / A[j + j * n]; for (i = j + 1; i < n; i++) A[i + j * n] *= tmp;
		}
	}
	free(vv);
	return 0;
}
/* LU back-substitution ------------------------------------------------------*/
static void lubksb(const double* A, int n, const int* indx, double* b)
{
	double s;
	int i, ii = -1, ip, j;

	for (i = 0; i < n; i++) {
		ip = indx[i]; s = b[ip]; b[ip] = b[i];
		if (ii >= 0) for (j = ii; j < i; j++) s -= A[i + j * n] * b[j]; else if (s) ii = i;
		b[i] = s;
	}
	for (i = n - 1; i >= 0; i--) {
		s = b[i]; for (j = i + 1; j < n; j++) s -= A[i + j * n] * b[j]; b[i] = s / A[i + i * n];
	}
}
/**
 * @brief ¾ØÕóÇóÄæ matrix = matrix^-1.
 * 
 * @param A
 * @param n
 * @return 
 */
extern int CMatrix_Inv(double* matrix, int order)
{
	double d, * B;
	int i, j, * indx;

	indx = CMatrixi(order, 1); B = CMatrixd(order, order); CMatrix_Cpy(B, matrix, order, order);
	if (ludcmp(B, order, indx, &d)) { free(indx); free(B); return -1; }
	for (j = 0; j < order; j++) {
		for (i = 0; i < order; i++) matrix[i + j * order] = 0.0;
		matrix[j + j * order] = 1.0;
		lubksb(B, order, indx, matrix + j * order);
	}
	free(indx); free(B);
	return 0;
}

