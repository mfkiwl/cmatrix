#include"cmatrix.h"
int main() {

	double* mat1 = CMatrixd(3, 3);
	double* mat2 = CMatrixd(3, 3);
	double* mat_r = CMatrixd(3, 3);
	double m_temp1[9] = { 1,2,3,4,5,6,7,8,9 };
	double m_temp2[9] = { 1,2,3,4,5,6,7,8,9 };
	
	CMatrixd_copy(mat1, m_temp1, 3, 3);
	CMatrixd_copy(mat2, m_temp2, 3, 3);
	//memcpy(mat1, m_temp, sizeof(double) * 9);
	
	CMatrixd_Print(mat1, 3, 3);
	CMatrixd_Print(mat2, 3, 3);

	CMatrixd_mul("NN", mat1, mat2, 3, 3, 3, 1.0, 0.0, mat_r);


	CMatrixd_Print(mat_r, 3, 3);

	double matt[9] = {0,1,2,1,0,3,4,-3,8};
	CMatrixd_Copy(mat_r, matt, 3, 3);

	CMatrix_Inv(mat_r, 3);

	CMatrixd_Print(mat_r,3,3);

	free(mat1);
	free(mat2);
	free(mat_r);

	

	return 10086;
}
