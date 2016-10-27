#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "task_01_14.h"

#define EXIT -1
#define ERR_MEM -2
#define ERR_FILE -3
#define ERR_FILE_FORMAT -4
#define OK 0
#define addr(i, j) ((i)*n+(j))
#define MAX(a, b) ((a)>(b) ? (a) : (b) )

extern char dbg_mode = 0, err_mode = 0;

void _(char *msg, int line)
{
	if(dbg_mode) printf("%s:%d: %s\n", "file", line, msg);
}

void err(char *msg)
{
	if(err_mode) fprintf(stderr, "Error: %s\n", msg);
}

void free_res(double *A, double *E, double *tmp, FILE *fin, FILE *fout)
{
	_("Deallocating resources", __LINE__);
	if(A!=NULL) free(A);
	if(E!=NULL) free(E);
	if(tmp!=NULL) free(tmp);
	if(fin!=NULL) free(fin);
	if(fout!=NULL) free(fout);
}

int get_param(char *param, int *max_iter, double *eps, double *prec)
{
	int i = 0, iter_n = 9, eps_n = 4, prec_n = 5;
	char iter_s[] = "max_iter=";
	char eps_s[] = "eps=";
	char prec_s[] = "prec=";

	while(param[i]!=0 && iter_s[i]!=0)
		if(param[i]!=iter_s[i]) break; else i++;
	if(i==iter_n) { *max_iter = atoi(param+iter_n); return OK; }

	i=0;
	while(param[i]!=0 && eps_s[i]!=0)
		if(param[i]!=eps_s[i]) break; else i++;
	if(i==eps_n) { *eps = atof(param+eps_n); return OK; }

	i=0;
	while(param[i]!=0 && prec_s[i]!=0)
		if(param[i]!=prec_s[i]) break; else i++;
	if(i==prec_n) { *prec = atof(param+prec_n); return OK; }

	return EXIT;
}

void print_matrix(int n, double *A)
{
	int i, j;
	for(i=0; i<n; i++)
	{
		for(j=0; j<n; j++) printf ("%lf ", A[addr(i, j)]);
		printf("\n");
	}
	printf("\n");
}

void print_eig_val(int n, double *E)
{
	int i;
	for(i=0; i<n; i++) printf("%lf ", E[i]);
	printf("\n");
}

int parse_cmdline(int argc, char **argv, int *in_file, int *out_file, char *out_mode, char *stat_mode, int *max_iter, double *eps, double *prec)
{
	char help[] = "Usage: evc [input_file_name] [output_file_name] [options]\n \
Where options include:\n \
  -d    print debug messages [default OFF]\n \
  -e    print errors [default OFF]\n \
  -p    print matrix [default OFF]\n \
  -t    print execution time [default OFF]\n \
  -prec=<num>	precision [default - 1e-14]\n \
  -eps=<num>	epsilon [default - 1e-10]\n \
  -max_iter=<num>	limit number of iterations \
			[default - 0, i.e. no limit]\n \
  -h, -?    print this and exit\n \
Default input_file_name value is 01_14_in.txt, default output_file_name value is 01_14_out.txt.\n";

	int i, count=0;

	for(i=1; i<argc; i++)
	{
		if(argv[i][0]=='-' && argv[i][2]==0) {
			char t = argv[i][1];
			if(t=='d') dbg_mode=1;
			else if(t=='e') err_mode=1;
			else if(t=='p') *out_mode=1;
			else if(t=='t') *stat_mode=1;
			else if(t=='h' || t=='?') { printf(help); return OK; }
			else { printf("Error: Wrong usage!\n\n%s", help); return EXIT; }
		}
		else if(argv[i][0]=='-' && argv[i][2]!=0) {
			if(get_param(argv[i]+1, max_iter, eps, prec)==-1) {
				printf("Error: Wrong usage!\n\n%s", help); return EXIT;
			}
		}
		else {
			if(count==0) { *in_file=i; count++; }
			else if(count==1) { *out_file=i; count++; }
			else { printf("Error: Wrong usage!\n\n%s", help); return EXIT; }
		}
	}
	return OK;
}

int main(int argc, char **argv)
{
	double *A=NULL, *E=NULL, *tmp=NULL;
	FILE *fin, *fout;
	int n, i, j, sz_tmp, max_iter=0, res=0, in_file=-1, out_file=-1;
	double prec=1e-14, eps=1e-10;
	char out_mode = 0, stat_mode = 0;
	clock_t t;	

	if(parse_cmdline(argc, argv, &in_file, &out_file, &out_mode, &stat_mode, &max_iter, &eps, &prec)==EXIT) return EXIT;

	_("Opening files", __LINE__);
	fin = fopen(in_file==-1 ? "01_14_in.txt" : argv[in_file], "r");
	fout = fopen(out_file==-1 ? "01_14_out.txt" : argv[out_file], "w");
	
	if(fin==NULL || fout==NULL) {
		_("Can't open input or/and output file!", __LINE__);
		err("Can't open input or/and output file!");
		free_res(A, E, tmp, fin, fout);
		return ERR_FILE;
	}
	
	_("Reading input file", __LINE__);
	if(fscanf(fin, "%d", &n)==-1) {
		err("Wrong file format!");
		free_res(A, E, tmp, fin, fout);
		return ERR_FILE_FORMAT;
	}

	sz_tmp = MAX(evc_memsize_01_14(n), sim_memsize_01_14(n));

	_("Allocating memory", __LINE__);
	A=(double *)malloc(sizeof(double)*n*n);
	tmp=(double *)malloc(sz_tmp);
	E=(double *)malloc(sizeof(double)*n);

	if(A==NULL || E==NULL || (tmp==NULL && sz_tmp!=0))
	{
		_("Memory allocation failed", __LINE__);
		err("Can't allocate memory!");
		free_res(A, E, tmp, fin, fout);
		return ERR_MEM;
	}

	_("Reading rest of input file", __LINE__);
	for(i=0; i<n; i++)
		for(j=0; j<n; j++)
			if(fscanf(fin, "%lf", &A[i*n+j])==-1) 
			{ 
				err("Wrong file format!"); 
				free_res(A, E, tmp, fin, fout);
				return ERR_FILE_FORMAT;
			}

	if(out_mode==1) print_matrix(n, A);
	
	_("Simplifying matrix", __LINE__);
	t = clock();
	res = sim_01_14(n, A, tmp, prec);
	t = clock() - t;
	if(stat_mode) printf("Execution time of simplification part %f seconds\n", ((float)t)/CLOCKS_PER_SEC);

	if(res==-1) { fprintf(fout, "0\n"); free_res(A, E, tmp, fin, fout); return OK; }
	if(out_mode==1) print_matrix(n, A);
	
	_("Computing eigen values", __LINE__);
	t = clock();
	res = evc_01_14(n, max_iter, prec, A, E, tmp, eps);
	t = clock() - t;
	if(stat_mode) printf("Execution time of computin eigen values %f seconds\n", ((float)t)/CLOCKS_PER_SEC);

	if(res==1) { fprintf(fout, "0\n"); free_res(A, E, tmp, fin, fout); return OK; }
	if(out_mode==1) print_eig_val(n, E);
	
	_("Writing to output file", __LINE__);
	fprintf(fout, "%d\n", n);
	for(i=0; i<n; i++) fprintf(fout, "%1.9lf\n", E[i]);

	free_res(A, E, tmp, fin, fout);

	return OK;
}


