#include"cmatrix.h"
int main() {

	double* mat1 = CMatrixd(3, 3);
	double m_temp[9] = { 1,2,3,4,5,6,7,8,9 };

	CMatrixd_copy(mat1, m_temp, 3, 3);
	//memcpy(mat1, m_temp, sizeof(double) * 9);
	
	CMatrixd_Print(mat1, 3, 3);

	free(mat1);

	

	return 10086;
}
