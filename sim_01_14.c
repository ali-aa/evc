#include <math.h>
#define addr(i, j) ((i)*n+(j))

void comp_phi(double *A, double *cos_phi, double *sin_phi, int n, int i, int j)
{
	double xi, xj, t;
	xi = A[addr(i, i-1)];
	xj = A[addr(j, i-1)];
	t = sqrt(xi*xi+xj*xj);
	*cos_phi = xi / t;
	*sin_phi = -xj / t;
}

int is_symmetric(int n, double *A)
{
	int i, j;
	for(i=0; i<n; i++)
		for(j=0; j<i; j++)
			if(A[addr(i, j)]!=A[addr(j, i)]) return -1;
	return 0;
}

int sim_memsize_01_14(int n)
{
	return 0;
}

void left_mult_T_v(double *A, int k, double cos_phi, double sin_phi, int i, int j, int n)
{
	double ti, tj;
	ti = A[addr(i, k)]; tj = A[addr(j, k)];
	A[addr(i, k)] = ti*cos_phi - tj*sin_phi;
	A[addr(j, k)] = ti*sin_phi + tj*cos_phi;
}

void left_mult_T_A(double *A, int n, double cos_phi, double sin_phi, int i, int j)
{
	int k;
	
	for(k=i-1; k<n; k++)
		left_mult_T_v(A, k, cos_phi, sin_phi, i, j, n);
}

void right_mult_T_v(double *A, int k, double cos_phi, double sin_phi, int i, int j, int n)
{
	double ti, tj;
	ti = A[addr(k, i)]; tj = A[addr(k, j)];
	A[addr(k, i)] = ti*cos_phi - tj*sin_phi;
	A[addr(k, j)] = ti*sin_phi + tj*cos_phi;
}

void right_mult_T_A(double *A, int n, double cos_phi, double sin_phi, int i, int j)
{
	int k;
	
	for(k=i-1; k<n; k++)
		right_mult_T_v(A, k, cos_phi, sin_phi, i, j, n);
}

// returns
// 0: success
// -1: fail
int sim_01_14(int n, double* A, double* tmp, double precision)
{
	int i, j;
	double cos_phi, sin_phi;

	if(is_symmetric(n, A)!=0) return -1;

	for(i=1; i<n-1; i++)
		for(j=i+1; j<n; j++)
		{
			if(A[addr(i, i-1)]*A[addr(i, i-1)] + A[addr(j, i-1)]* A[addr(j, i-1)]<precision) continue;
			comp_phi(A, &cos_phi, &sin_phi, n, i, j);
			// умножаем слева
			left_mult_T_A(A, n, cos_phi, sin_phi, i, j);
			// умножаем справа
			right_mult_T_A(A, n, cos_phi, sin_phi, i, j);
		}
	return 0;
}


