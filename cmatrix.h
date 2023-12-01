#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#ifndef CMATRIX_H
#define CMATRIX_H
extern double* CMatrixd(int row, int column);
extern float* CMatrixf(int row, int column);
extern int* CMatrixi(int row, int column);

extern double* CMatrixd_Zero(int row, int column);
extern float* CMatrixf_Zero(int row, int column);
extern int* CMatrixi_Zero(int row, int column);

extern double* CMatrixd_Eye(int order);
extern float* CMatrixf_Eye(int order);
extern int* CMatrixi_Eye(int order);

extern double Dotd(const double* vector1, const double* vector2, int order);
extern float Dotf(const float* vector1, const float* vector2, int order);
extern int Doti(const int* vector1, const int* vector2, int order);

extern double Normd(const double* vector, int order);
extern float Normf(const float* vector, int order);
extern int Normi(const int* vector, int order);

extern void Cross3d(const double* vector1, const double* vector2, double* vector_result);
extern void Cross3f(const float* vector1, const float* vector2, float* vector_result);
extern void Cross3i(const int* vector1, const int* vector2, int* vector_result);


extern int CMatrixd_Print(double* matrix, int row, int column);
extern int CMatrixf_Print(float* matrix, int row, int column);
extern int CMatrixi_Print(int* matrix, int row, int column);

extern void CMatrixd_Copy(double* mat_dest, const double* mat_src, int row, int col);
extern void CMatrixf_Copy(float* mat_dest, const float* mat_src, int row, int col);
extern void CMatrixi_Copy(int* mat_dest, const int* mat_src, int row, int col);

extern void CMatrixd_Mul(const char* transpose, const double* left_matrix, const double* right_matrix,
	int left_mat_row, int l_row_r_col, int right_mat_col, double alpha, double beta, double* result_matrix);
extern void CMatrixf_Mul(const char* transpose, const float* left_matrix, const float* right_matrix,
	int left_mat_row, int l_row_r_col, int right_mat_col, float alpha, float beta, float* result_matrix);
extern void CMatrixi_Mul(const char* transpose, const int* left_matrix, const int* right_matrix,
	int left_mat_row, int l_row_r_col, int right_mat_col, int alpha, int beta, int* result_matrix);

extern int CMatrix_Inv(double* matrix, int order);

#endif // !CMATRIX_H

