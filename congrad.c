#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<time.h>
#define MAX_ITERS 10000
#define tol 1e-6

double error = 1000;
int iter = 0;

void check_tol(double *r1, double *error, int n){
    int i;
    double e = 0;
    for(i=0;i<n;i++){
	e +=pow(*(r1+i),2);
    }
    *error = sqrt(e);
}

void init_A_b(double **A, double *b, int n){
    int i, j;
srand(time(NULL));
    for(i=0;i<n;i++){
	*(b+i) = rand()%10;
	for(j=0;j<n;j++){
	    *(*(A+i)+j) = rand()%10 ;
	}
    }

}

double run_congrad(double **A, double *b, double *x0, double *x1, double *p0, double *p1, double *r0, double*r1, int n, double * error){
    double e = 0;
    int i, j;
    double a = 0;
    double a_den = 0;
    double a_num = 0;
    double beta = 0;
    double beta_num = 0;
    double beta_den = 0;
    double * temp = malloc(sizeof(double)*n);
    /*clear temp*/
    for(i=0;i<n;i++){
	*(temp+i) = 0;
    }
    /* get p0*A  */
    for(i=0;i<n;i++){
	for(j=0;j<n;j++){
	    *(temp+i) += *(p0+j) * *(*(A+i)+j);
	}		
    }
    /* get (p0A) * p0 */
    for(i = 0;i<n;i++){
	a_den += (*(temp+i)) * (*(p0+i));
    }
    /*get numerator of a  = r' * r*/
    for(i = 0;i<n;i++){
	a_num += (*(r0+i)) * (*(r0+i));
    }
    a = a_num / a_den;
    //  printf("a: %f\n",a);
    /*clear temp*/
    for(i=0;i<n;i++){
	*(temp+i) = 0;
    }

    /* get a * p0 */
    for(i = 0;i<n;i++){
	*(temp+i)  = a * (*(p0+i));
    }
    /* get x1 = x0 + ap0 */
    for(i = 0;i<n;i++){
	*(x1+i) = (*(x0+i)) + *(temp+i);
	//	printf("x1[%d]: %f\n",i, *(x1+i));
    }
    /* get a*A  */
    for(i=0;i<n;i++){
	for(j=0;j<n;j++){
	    *(*(A+i)+j) = a * (*(*(A+i)+j));
	}		
    }
    /*clear temp*/
    for(i=0;i<n;i++){
	*(temp+i) = 0;
    }
    /* get A*p0  */
    for(i=0;i<n;i++){
	for(j=0;j<n;j++){
	    *(temp+i) += *(*(A+i)+j) * (*(p0 + j));
	}		
    }
    /*finish r1 =  r0 - aAp0*/
    for(i = 0;i<n;i++){
	(*(r1+i))  =  (*(r0+i)) - (*(temp+i));
	//	printf("r1[%d]: %f\n",i, *(r1+i));
    }
    /* clear  a*A  */
    for(i=0;i<n;i++){
	for(j=0;j<n;j++){
	    *(*(A+i)+j) =  (*(*(A+i)+j)) / a;
	}		
    }

    /*clear temp*/
    for(i=0;i<n;i++){
	*(temp+i) = 0;
    }
    /*now beta = r1*r1/r0*r0*/
    for(i = 0;i<n;i++){
	beta_num +=  (*(r1+i)) * (*(r1+i));
    }
    /*cont beta = r1*r1/r0*r0*/
    for(i = 0;i<n;i++){
	beta_den +=  (*(r0+i)) * (*(r0+i));
    }

    beta = beta_num / beta_den;
    /*clear temp*/
    for(i=0;i<n;i++){
	*(temp+i) = 0;
    }

    /*p1 time*/
    for(i = 0;i<n;i++){
	*(temp+i) =  beta  * *(p0+i);
    }
    /*p1 time*/
    for(i = 0;i<n;i++){
	*(p1+i) =  *(r1+i)  + *(temp+i);
    }
    check_tol(r1, error, n);	
    e = *error;
    /*clear temp*/
    for(i=0;i<n;i++){
	*(temp+i) = 0;
    }


    return e;
}


int main(int argc, char *argv[]){
    int i, n = 2;
    if(argc == 2){
	n = atoi(argv[1]);
    }
    double **A = malloc(sizeof(double)*n);
    double *data = malloc(sizeof(double)*n*n);
    for(i=0;i<n;i++)
    {
	*(A+i) = data + n*i;
    }

    double *b = malloc(sizeof(double)*n);
    init_A_b(A, b, n);
    double *x0 = malloc(sizeof(double)*n);
    double *x = malloc(sizeof(double)*n);
    double *p0 = malloc(sizeof(double)*n);
    double *r0 = malloc(sizeof(double)*n);
    double *x1 = malloc(sizeof(double)*n);
    double *p1 = malloc(sizeof(double)*n);
    double *r1 = malloc(sizeof(double)*n);
    /*init x0, r0, p0*/
    for(i=0;i<n;i++){
	*(x0+i) = 0;
	*(r0+i) = *(b+i);
	*(p0+i) = *(r0+i);
    }
    clock_t start = clock();
    while(error > tol  && iter < MAX_ITERS){
	iter++;
	error = run_congrad(A, b, x0, x1, p0, p1, r0, r1, n, &error);
	p0 = p1;
	x0 = x1;
	r0 = r1;	
    }
    clock_t end = clock();
    x = x0;
    
/*    printf("x = \n",*(x+i));*/
    for(i=0;i<n;i++){
//	printf("    |  %.5f  |\n",*(x+i));
    }
    if(*(x+0) != *(x+0)){
//	printf("*************************************************************************************\n");
//	printf("Try running the program again, random number generator produced a degenerate A matrix\n");
//	printf("*************************************************************************************\n");
    }
	printf("Running time: %f\n",(double)(end-start)/(double)CLOCKS_PER_SEC);

    return 0;
}
