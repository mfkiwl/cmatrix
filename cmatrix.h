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

#endif // !CMATRIX_H

