#include <math.h>
#define addr(i, j) ((i)*n+(j))
#define SWAP(A, B) { double t = A; A = B; B = t; }

void sort(double *A, int n)
{
  int i, j;
 
  for (i=n-1; i>0; i--)
  {
    for (j=0; j<i; j++)
    {
      if (A[j]>A[j+1]) 
        SWAP(A[j], A[j+1]);
    }
  }
}

//Норма A бесконечность
double norm(double *A, int n)
{
	int i, j;
	double t, max=0;
	for(i=0; i<n; i++)
	{
		t=0;
		for(j=0; j<n; j++) t += fabs(A[addr(i, j)]);
		if(t>max) max = t;
	}
	return max;
}

void add_shift(double *A, double s_k, int n, int cur_n)
{
	int i;
	for(i=0; i<cur_n; i++) A[addr(i, i)] += s_k;
}

void sub_shift(double *A, double s_k, int n, int cur_n)
{
	int i;
	for(i=0; i<cur_n; i++) A[addr(i, i)] -= s_k;
}

void comp(double *A, double *cos_phi, double *sin_phi, int n, int i, int j)
{
	double xi, xj, t;
	xi = A[addr(i, i)]; 
	xj = A[addr(j, i)];
	t = sqrt(xi*xi+xj*xj);
	*cos_phi = xi / t;
	*sin_phi = -xj / t;
}

int evc_memsize_01_14(int n)
{
	return 2*n*sizeof(double);
}

void left_mul_T_v(double *A, int k, double cos_phi, double sin_phi, int i, int j, int n)
{
	double ti, tj;
	ti = A[addr(i, k)]; tj = A[addr(j, k)];
	A[addr(i, k)] = ti*cos_phi - tj*sin_phi;
	A[addr(j, k)] = ti*sin_phi + tj*cos_phi;
}

void left_mul_T_A(double *A, int n, int cur_n, double cos_phi, double sin_phi, int i, int j)
{
	int k;
	
	for(k=i; k<cur_n; k++)
		left_mult_T_v(A, k, cos_phi, sin_phi, i, j, n);
}

void right_mul_T_v(double *A, int k, double cos_phi, double sin_phi, int i, int j, int n)
{
	double ti, tj;
	ti = A[addr(k, i)]; tj = A[addr(k, j)];
	A[addr(k, i)] = ti*cos_phi - tj*sin_phi;
	A[addr(k, j)] = ti*sin_phi + tj*cos_phi;
}

void right_mul_T_A(double *A, int n, int cur_n, double cos_phi, double sin_phi, int i, int j)
{
	int k;
	
	for(k=i; k<cur_n; k++)
		right_mul_T_v(A, k, cos_phi, sin_phi, i, j, n);
	for(k=i; k<n; k++) A[addr(i, k)] = A[addr(k, i)];
}

int evc_01_14(int n, int max_iterations, double epsilon, double* A, double* E, double* tmp, double precision)
{
	int i, iter, cur_n;
	double *Q_cos, *Q_sin, cos_phi, sin_phi, s_k, D, a00, a01, a10, a11, norm_a;
	
	// Норма преобразованной трехдиаганальной матрицы, до работы QR алгоритма
	norm_a = norm(A, n);
	Q_cos = &tmp[0];
	Q_sin = &tmp[n];
	cur_n = n;

	for(iter=0; (iter<max_iterations || max_iterations==0) && cur_n>2; iter++)
	{
		// Отнимаем сдвиг
		if(iter%101==0) 
		{
			s_k = A[addr(cur_n-1, cur_n-1)];
			sub_shift(A, s_k, n, cur_n);
		}

		// QR разложение
		for(i=0; i<cur_n-1; i++)
		{
if(A[addr(i, i)]*A[addr(i, i)] + A[addr(i+1, i)]* A[addr(i+1, i)]<precision) continue;

			comp(A, &cos_phi, &sin_phi, n, i, i+1);
			Q_cos[i] = cos_phi;
			Q_sin[i] = sin_phi;
			// умножаем слева
			left_mul_T_A(A, n, cur_n, cos_phi, sin_phi, i, i+1);
		}

		// теперь A есть матрица R, умножаем ее на матрицу Q
		for(i=0; i<cur_n-1; i++)
			right_mul_T_A(A, n, cur_n, Q_cos[i], Q_sin[i], i, i+1);

		//Добавляем сдвиг
		if(iter%101==0) add_shift(A, s_k, n, cur_n);

		if(fabs(A[addr(cur_n-1, cur_n-2)])<epsilon*norm_a)
		{
			E[cur_n-1] = A[addr(cur_n-1, cur_n-1)];
			cur_n--;
		}
	}
	// Для матрицы 2x2 собственные значения находим из решения квадратного уравнения
	if(cur_n==2)
	{
		a00 = A[addr(0, 0)]; a01 = A[addr(0, 1)];
		a10 = A[addr(1, 0)]; a11 = A[addr(1, 1)];
		D = (a00 + a11)*(a00 + a11) - 4*(a00*a11 - a01*a10);
		E[0] = (a00 + a11 + sqrt(D))/2;
		E[1] = (a00 + a11 - sqrt(D))/2;
	}
	else if(cur_n>2) return 1;

	// Сортируем собсвенные значения
	sort(E, n);
	
	return 0;
	
}

