
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "qwsim.h"


/*Complex numbers and martrices*/
complex c_sum(complex z1, complex z2){
	complex r;
	r.re = z1.re + z2.re;
	r.im = z1.im + z2.im;

	return r;
}

complex c_mult(complex z1, complex z2){
/*Receivestwo complex numbers z1 and z2 and returns r = z_1 + z_2*/
	double a,b,c,d;
	complex r;
	a = z1.re;
	b = z1.im;
	c = z2.re;
	d = z2.im;
	r.re = a*c-b*d;
	r.im = b*c+a*d;
	return r;
}

double c_angle(complex z){
/*Receives complex number z and returns arctan(Re(z)/Im(z))*/
	double r = 0.0;
	if(z.re > 0.0){
		r = atan(z.im/z.re);
	}
	if(z.re < 0.0 && z.im >= 0.0){
		r = atan(z.im/z.re) + M_PI;
	}
	if(z.re < 0.0 && z.im < 0.0){
		r = atan(z.im/z.re) - M_PI;
	}
	if(z.re == 0.0 && z.im > 0.0){
		r = M_PI/2.0;
	}
	if(z.re == 0.0 && z.im < 0.0){
		r = -M_PI/2.0;
	}
	return r;
}

complex c_div(complex z1, complex z2){
/*Receives complex numbers z1 and z2 and returns z1/z2*/
	double a,b,c,d;
	complex r;
	a = z1.re;
	b = z1.im;
	c = z2.re;
	d = z2.im;

	r.re = (a*c+b*d)/(c*c + d*d);
	r.im = (b*c-a*d)/(c*c + d*d);

	return r;
}

complex c_conj(complex c){
/*Receives complex number z and retunrs its complex conjugate*/
	complex z;

	z.re = c.re;
	z.im = -c.im;

	return z;
}


double c_abs_sq(complex z){
/*Receives complex number z and outputs it squared modulus*/
double a = z.re;
double b = z.im;

return a*a+b*b;

}

/*complex matrix functions*/

double vector_norm_sq(c_matrix a){
/*Receives vector a and returns squared norm of a*/
	int i;
	double r = 0.0;
	for(i=0; i < a.rows; i++){
		r = r + c_abs_sq(M(a,i,0));
	}
	return r;
}

void new_matrix(c_matrix *r,int rows,int cols)
{
/*Initializes matrix r with "rowns" rows and "cols" columns with all entries equal to 0.0*/
	int i;
	r->rows=rows;
	r->cols=cols;
	r->t=(complex *)malloc(sizeof(complex)*cols*rows);
	for(i=0;i< rows*cols; i++){
		r->t[i].re = 0.0;
		r->t[i].im = 0.0;
	}
	N_malloc = N_malloc + 1;

}

void reset_matrix(c_matrix *r){
/*Receives a matrix r and resets its entries to 0.0*/
	int i;
	for(i=0;i< (r->rows)*(r->cols); i++){
		r->t[i].re = 0.0;
		r->t[i].im = 0.0;
	}
	
}

void erase_matrix(c_matrix *a){
/*Erases matrix a*/
	free(a->t);
	a->t = NULL;
	a->cols = 0;
	a->rows = 0;
	N_free = N_free + 1;
}

void print_matrix(c_matrix t) {
/*Prints matrix t on screen*/
	int i,j;
	for(i=0;i<t.rows;i++) {
		printf("| ");
		for(j=0;j<t.cols;j++)
			printf("%f +i%f ",M(t,i,j).re,M(t,i,j).im);
		printf("|\n");
	}
	printf("\n");
}

void matrix_mult(c_matrix r,c_matrix a,c_matrix b){
/*multiplies matrix a with b and storees result on r*/
	int m,n,k;
	for(m=0;m<a.rows;m++){
		for(k=0;k<b.cols;k++) {
			M(r,m,k).re = 0.0;
			M(r,m,k).im = 0.0;
			for(n=0;n<a.cols;n++){
				M(r,m,k)=c_sum(M(r,m,k),c_mult(M(a,m,n),M(b,n,k)));
			}
		}
	}
}

void matrix_sum(c_matrix r, c_matrix a, c_matrix b) {
/*Sums matrix a with b and storees result on r*/
	int i,j;
	for(i=0;i<a.rows;i++){
		for(j=0;j<b.cols;j++) {
			M(r,i,j)=c_sum(M(a,i,j),M(b,i,j));
		}
	}
}

void c_times_matrix(c_matrix r, complex c, c_matrix a){
/*multiplies constant c with matrix a and the result is stored on r*/
	int i,j;
	for(i=0;i<r.rows;i++){
		for(j=0;j<r.cols;j++) {
			M(r,i,j)=c_mult(c,M(a,i,j));
		}
	}
}


void matrix_conj(c_matrix r,c_matrix a){
/*stores in matrix r the complex conjugate of the entries of matrix a */
	int i,j;

	for(i = 0; i < a.cols; i++){
		for(j=0; j < a.rows; j ++){
			M(r,i,j) = c_conj(M(a,i,j));
		}
	}
}

void matrix_transpose(c_matrix r, c_matrix a){
/* stores in matrix r the transpose of matrix a */
	int i,j;

	for(i=0;i<a.rows;i++){
		for(j=0;j<a.cols;j++){
			M(r,j,i) = M(a,i,j);
		}
	}
}

void matrix_adjoint(c_matrix *r,c_matrix a){
/*stores in matrix r the adjoint of matrix a*/
	int i,j;
	for(i=0;i<r->rows;i++){
		for(j=0;j<r->cols;j++){
			M((*r),i,j) = c_conj(M(a,j,i));
		}
	}
}

void position_trace_2D(state *s,int t){
/*Given state |s> this function determines the partial trace over the coin space*/
	int c1,c2,i_1,j_1, i_2, j_2;


	for(c1=0; c1 < Edges; c1++)
		for(c2=0; c2 < Edges; c2++){
	M(rho_c,c1,c2).re = 0.0;
	M(rho_c,c1,c2).im = 0.0;
	}

	for(c1=0; c1 < Edges; c1++)
		for(c2=0; c2 < Edges; c2++)
			for(i_1=-t; i_1 <= t; i_1++)
				for(j_1=-t; j_1 <= t; j_1++){

	M(rho_c,c1,c2) = c_sum(M(rho_c,c1,c2), c_mult(M(s[biject(i_1,j_1)].direction,c1,0), c_conj(M(s[biject(i_1, j_1)].direction,c2,0))));
	}
}



void identity(c_matrix I){
/*Initializes matrix I with entries of the identity matrix */
	int i;

	for(i=0; i < I.cols; i++){
		M(I,i,i).re = 1.0;
		M(I,i,i).im = 0.0;
	}
}

void copy_matrix(c_matrix r, c_matrix a){
/*Copies entries of matrix a into matrix r*/
	int i,j;

	for(i = 0; i < a.rows; i++){
		for(j=0; j < a.cols; j ++){
			M(r,i,j) = M(a,i,j);
		}
	}
}

c_matrix I_diag_B(c_matrix a, int dim){
/*Returns a block diagonal matrix r of dimension "dim" where the first block is the diagonal matrix and the second is matrix a */
	c_matrix r;
	int i,j,ir,jr;
	i = -1;

	new_matrix(&r,dim,dim);
	identity(r);
	for(ir= dim-a.rows; ir < dim; ir++){
		i++;
		j=-1;
		for(jr=dim-a.cols; jr < dim; jr++){
			j++;
			M(r,ir,jr) = M(a,i,j);
		}
	}
	return r;
}

c_matrix copy_sub_matrix(c_matrix a,int *lines, int *cols){
/*Copies the entries of matrix a from lines[0] to lines[1] and cols[0] to cols[1] into matrix r*/
	c_matrix r;
	int i,j,ir,jr;
	ir = -1;

	new_matrix(&r,1+lines[1]-lines[0],1+cols[1]-cols[0]);
	for(i=lines[0]; i <= lines[1]; i++){
		ir++;
		jr=-1;
		for(j=cols[0]; j <= cols[1]; j++){
			jr++;
			M(r,ir,jr) = M(a,i,j);
		}
	}
	return r;
}
 
c_matrix QR_decomposition(c_matrix A){
/*Returns matrix Q from QR decompositionn of A, A = QR*/
	int i,lins[2],cols[2];
	int T_malloc = 0;
	int T_free = 0;
	c_matrix x, P_i,B_i, Q_aux,PA, PA_aux,Q;

	lins[1] = A.cols-1;
	new_matrix(&Q_aux,A.rows, A.cols);
	T_malloc++;
	new_matrix(&PA, A.rows, A.cols);
	T_malloc++;
	new_matrix(&PA_aux,A.rows, A.cols);
	T_malloc++;
	new_matrix(&Q,A.rows, A.cols);
	//T_malloc++;
	copy_matrix(PA,A);
	identity(Q);

	for(i=0; i < A.cols; i++){
		lins[0] = i;
		cols[0] = i;
		cols[1] = i;

		x = copy_sub_matrix(PA,lins,cols);
		T_malloc++;
		new_matrix(&B_i, 1+lins[1]-lins[0],1+ lins[1]-lins[0]);
		T_malloc++;
		Householder_matrix(B_i,x);
		P_i = I_diag_B(B_i,A.cols);
		T_malloc++;
		matrix_mult(PA_aux,P_i,PA);
		copy_matrix(PA, PA_aux);

		matrix_mult(Q_aux,P_i,Q);
		copy_matrix(Q,Q_aux);

		erase_matrix(&B_i);
		T_free++;
		erase_matrix(&P_i);
		T_free++;
		erase_matrix(&x);
		T_free++;
	}
	erase_matrix(&PA);
	T_free++;
	erase_matrix(&PA_aux);
	T_free++;
	erase_matrix(&Q_aux);
	T_free++;

	return Q;
}



c_matrix I_diag_B2(c_matrix a, int *lins, int *cols, int dim){
/*Returns block matrix dim x dim r initialized with identity matrix and with
 entries of matrix ''a'' on lines ''lins[0]'' to ''lins[1]'' and columns ''cols[0]'' to ''cols[1]''*/
	c_matrix r;
	int i,j,ir,jr;
	i = -1;

	new_matrix(&r,dim,dim);
	identity(r);
	for(ir=lins[0]; ir < lins[1]; ir++){
		i++;
		j=-1;
		for(jr=cols[0]; jr < cols[1]; jr++){
			j++;
			M(r,ir,jr) = M(a,i,j);
		}
	}
	return r;
}



void Householder_matrix(c_matrix P, c_matrix x){
/*Computes the Householder matrix for QR decomposition*/
	c_matrix v, v_ad;
	int i;
	double gamma = c_angle(M(x,0,0));
	double v_norm_sq;
	double x_norm = sqrt(vector_norm_sq(x));
	complex e_gamma;
	complex aux, aux_2;
	int tm = 0;
	int tf = 0;

	e_gamma.re = cos(gamma);

	e_gamma.im = sin(gamma);


	new_matrix(&v,x.rows,x.cols);
	tm++;
	new_matrix(&v_ad, x.cols, x.rows);
	tm++;

	copy_matrix(v,x);
	aux_2.re = x_norm;
	aux_2.im = 0.0;
	aux = c_mult(e_gamma, aux_2);
	M(v,0,0) = c_sum(M(v,0,0) , aux);

	v_norm_sq = vector_norm_sq(v);
	matrix_adjoint(&v_ad, v);

	matrix_mult(P,v,v_ad);

	aux.re = 0.0;
	aux.im = 0.0;
	if(v_norm_sq != 0.0)
		aux.re = -2.0/v_norm_sq;

	c_times_matrix(P,aux, P);
	aux.re = 1.0;
	for(i=0; i < P.rows; i++){
		M(P,i,i) = c_sum(M(P,i,i), aux);
	}

	erase_matrix(&v_ad);
	tf++;
	erase_matrix(&v);
	tf++;

}

double cost_function(c_matrix A){
/*returns the sum of squared absolute values of matrix A's diagonal matrix*/
	int i;
	double r = 0.0;

	for(i=0;i < A.cols;i++){
		r = r + c_abs_sq(M(A,i,i));
	}
	return r;
}

void Householder_matrix2(c_matrix P, c_matrix x){
/*Computes the Householder matrix for QR decomposition*/
	c_matrix v, v_ad;
	int i;
	double gamma = c_angle(M(x,0,0));
	double v_norm_sq;
	double x_norm = sqrt(vector_norm_sq(x));
	complex e_gamma;
	complex aux, aux_2;

	e_gamma.re = cos(gamma);

	e_gamma.im = sin(gamma);

	new_matrix(&v,x.rows,x.cols);
	new_matrix(&v_ad, x.cols, x.rows);

	copy_matrix(v,x);

	aux_2.re = x_norm;
	aux_2.im = 0.0;
	aux = c_mult(e_gamma, aux_2);
	
	M(v,1,0) = c_sum(M(v,1,0) , aux);

	v_norm_sq = vector_norm_sq(v);
	matrix_adjoint(&v_ad, v);

	matrix_mult(P,v,v_ad);

	aux.re = 0.0;
	aux.im = 0.0;
	if(v_norm_sq != 0.0)
		aux.re = -2.0/v_norm_sq;

	c_times_matrix(P,aux, P);
	aux.re = 1.0;
	for(i=0; i < P.rows; i++){
		M(P,i,i) = c_sum(M(P,i,i), aux);
	}
	
	erase_matrix(&v_ad);
	erase_matrix(&v);
}


c_matrix QR_algorithm(c_matrix A){
/*Returns matrix R of QR factorization A = QR*/
	int k;
	int T_malloc = 0;
	int T_free = 0;
	complex c;
	c_matrix Q,Qt,R,Ak, Ak_1,Aux, Aux2;
	c.re = -1.0;
	c.im = 0.0;
	new_matrix(&R,A.cols,A.cols);
	//T_malloc++;
	new_matrix(&Ak,A.cols,A.cols);
	T_malloc++;
	new_matrix(&Ak_1,A.cols,A.cols);
	T_malloc++;
	new_matrix(&Aux,A.cols,A.cols);
	T_malloc++;
	new_matrix(&Aux2, A.cols, A.cols);
	T_malloc++;
	new_matrix(&Q, A.cols, A.cols);
	T_malloc++;
	copy_matrix(Ak_1,A);

	for(k=1;k<=50;k++){

		Qt = QR_decomposition(Ak_1);
		T_malloc++;
		matrix_adjoint(&Q,Qt);
		matrix_mult(R,Qt,Ak_1);
		matrix_mult(Ak,R,Q);



		c_times_matrix(Aux, c, Ak_1);
		copy_matrix(Ak_1,Ak);

		erase_matrix(&Qt);
		T_free++;
		matrix_sum(Aux2, Aux, Ak);
		if(cost_function(Aux2) <= pow(10.0,-5.0)){
			break;
		}
		
		
	}

	erase_matrix(&Ak_1);
	T_free++;
	erase_matrix(&Aux);
	T_free++;
	erase_matrix(&Aux2);
	T_free++;
	erase_matrix(&Ak);
	T_free++;
	erase_matrix(&Q);
	T_free++;
	
	for(k=0;k<R.cols; k++){
		M(R,k,k).re = - M(R,k,k).re;
		M(R,k,k).im = - M(R,k,k).im;
	}

	return R;

}

void QR_algorithm2(c_matrix A,c_matrix R ,c_matrix Qn){
/*Same as QR, but changes the values of R and Qn such that A = QRn */
	int k;
	int tm = 0;
	int tf = 0;
	complex c;
	c_matrix Qn_aux,Q,Qt,Ak, Ak_1,Aux, Aux2;
	c.re = -1.0;
	c.im = 0.0;
	new_matrix(&Ak,A.cols,A.cols);
	tm++;
	new_matrix(&Ak_1,A.cols,A.cols);
	tm++;
	new_matrix(&Aux,A.cols,A.cols);
	tm++;
	new_matrix(&Aux2, A.cols, A.cols);
	tm++;
	new_matrix(&Q, A.cols, A.cols);
	tm++;
	new_matrix(&Qn_aux, A.cols, A.rows);
	tm++;
	identity(Qn_aux);
	copy_matrix(Qn,Qn_aux);
	copy_matrix(Ak_1,A);

	for(k=1;k<=500;k++){

		Qt = QR_decomposition(Ak_1);
	tm++;
		matrix_adjoint(&Q,Qt);
		matrix_mult(Qn_aux, Qt,Qn);
		copy_matrix(Qn,Qn_aux);
		matrix_mult(R,Qt,Ak_1);
		matrix_mult(Ak,R,Q);

		c_times_matrix(Aux, c, Ak_1);
		copy_matrix(Ak_1,Ak);

		matrix_sum(Aux2, Aux, Ak);

		erase_matrix(&Qt);
		tf++;
	}

	erase_matrix(&Ak_1);
	tf++;
	erase_matrix(&Aux);
	tf++;
	erase_matrix(&Aux2);
	tf++;
	erase_matrix(&Ak);
	tf++;
	erase_matrix(&Qn_aux);
	tf++;
	erase_matrix(&Q);
	tf++;
	for(k=0;k<R.cols; k++){
		M(R,k,k).re = - M(R,k,k).re;
		M(R,k,k).im = - M(R,k,k).im;
	}

}






complex entropy(c_matrix D){
/*Determines the von Neumann entropy/Shannon entropy of diagonal matrix D*/
	int i;
	complex S;
	double aux;
	S.re = 0.0;
	S.im = 0.0;

	for(i=0; i < D.cols; i++){
		if(c_abs_sq(M(D,i,i)) <= pow(10.0,-10.0)){
			S.re = S.re + 0.0;
			S.im = S.im + 0.0;
		}
		else
		{
		aux = sqrt(c_abs_sq(M(D,i,i)));

		S.re = S.re - aux*log(aux);
		S.im = S.im - c_angle(M(D,i,i));
		}
	}

return S;
}
/*Quantum Walk functions*/

int biject(int a, int b){

	a = a+N;
	b = b+N;
	return b + (Nodes1D)*a;
}

void bijinv(int *r,int ab){
	int a,b;
	b = ab%(Nodes1D);
	a = (ab-b)/(Nodes1D);
	r[0] = a-N;
	r[1] = b-N;
}


state * initialize_state(){
/*Allocates memory for the overall quantum state*/
	state *s = (state *)malloc(sizeof(state) * Nodes);
	int i,j;
	for(i=0;i< Nodes ; i++)
	{
		s[i].direction.cols = 1;
		s[i].direction.rows = Edges;
		s[i].direction.t=(complex *)malloc(sizeof(complex)*Edges);
		for(j=0;j<Edges; j++){
			M(s[i].direction,j,0).re = 0.0;
			M(s[i].direction,j,0).im = 0.0;
		}
	}
	return s;
}

void erase_state(state *s){
/*Erases quantum state from memory*/
	int i;

	for(i=0; i < Nodes; i++)
		{
			s[i].direction.cols = 0;
			s[i].direction.rows = 0;
			free(s[i].direction.t);
		}

}


void copy_state(state *s1, state*s2){
/*Copies all the entries from quantum state s1 to quantum state s2*/
int i,e,j;

	for(i=0;i<Nodes ;i++){
		for(e=0; e < Edges; e++){
			M(s1[i].direction,e,0)= M(s2[i].direction,e,0);
		}
	}


}

void shift2D(state *s1, state *s2, int t){

	int i,j;

	for(i=-N;i <= N; i++){
		for(j=-N;j<= N;j++){
			shift_cond[biject(i,j)][0] = biject(i+1,j+1);
			shift_cond[biject(i,j)][1] = 0;

			shift_cond[biject(i,j)][2] = biject(i+1,j-1);
			shift_cond[biject(i,j)][3] = 1;

			
			shift_cond[biject(i,j)][4] = biject(i-1,j+1);
			shift_cond[biject(i,j)][5] = 2;


			shift_cond[biject(i,j)][6] = biject(i-1,j-1);
			shift_cond[biject(i,j)][7] = 3;
			
		}
	}

}

void circ(){
/*Sets the circle topology for both walkers */
	int i,j;
	for(j=-N; j< N; j++){
		shift_cond[biject(N,j)][0] = biject(-N,j+1);
		shift_cond[biject(N,j)][1] = 0;

		shift_cond[biject(-N,j)][4] = biject(N,j+1);
		shift_cond[biject(-N,j)][5] = 2;
		
	}


	for(j=-N+1; j<= N; j++){


		shift_cond[biject(N,j)][2] = biject(-N,j-1);
		shift_cond[biject(N,j)][3] = 1;

		shift_cond[biject(-N,j)][6] = biject(N,j-1);
		shift_cond[biject(-N,j)][7] = 3;
		
	}

	for(i=-N; i<N; i++){		
		shift_cond[biject(i,N)][0] = biject(i+1,-N);
		shift_cond[biject(i,N)][1] = 0;
        

		shift_cond[biject(i,-N)][2] = biject(i+1,N);
		shift_cond[biject(i,-N)][3] = 1;
	}

	for(i=-N+1; i<= N; i++){		

		shift_cond[biject(i,N)][4] = biject(i-1,-N);
		shift_cond[biject(i,N)][5] = 2;


		shift_cond[biject(i,-N)][6] = biject(i-1,N);
		shift_cond[biject(i,-N)][7] = 3;
	}

	shift_cond[biject(N,N)][0] = biject(-N,-N);
	shift_cond[biject(N,N)][1] = 0;


	shift_cond[biject(N,N)][2] = biject(-N,N-1);
	shift_cond[biject(N,N)][3] = 1;

	shift_cond[biject(N,N)][4] = biject(N-1,-N);
	shift_cond[biject(N,N)][5] = 2;

	shift_cond[biject(N,N)][6] = biject(N-1,N-1);
	shift_cond[biject(N,N)][7] = 3;

//__________________________________________________________________________________


	shift_cond[biject(N,-N)][0] = biject(-N,-N+1);
	shift_cond[biject(N,-N)][1] = 0;


	shift_cond[biject(N,-N)][2] = biject(-N,N);
	shift_cond[biject(N,-N)][3] = 1;

	shift_cond[biject(N,-N)][4] = biject(N-1,-N+1);
	shift_cond[biject(N,-N)][5] = 2;

	shift_cond[biject(N,-N)][6] = biject(N-1,N);
	shift_cond[biject(N,-N)][7] = 3;

//__________________________________________________________________________________

	shift_cond[biject(-N,N)][0] = biject(-N+1,-N);
	shift_cond[biject(-N,N)][1] = 0;


	shift_cond[biject(-N,N)][2] = biject(-N+1,N-1);
	shift_cond[biject(-N,N)][3] = 1;

	shift_cond[biject(-N,N)][4] = biject(N,-N);
	shift_cond[biject(-N,N)][5] = 2;

	shift_cond[biject(-N,N)][6] = biject(N,N-1);
	shift_cond[biject(-N,N)][7] = 3;

//_____________________________________________________________________________

	shift_cond[biject(-N,-N)][0] = biject(-N+1,-N+1);
	shift_cond[biject(-N,-N)][1] = 0;


	shift_cond[biject(-N,-N)][2] = biject(-N+1,N);
	shift_cond[biject(-N,-N)][3] = 1;

	shift_cond[biject(-N,-N)][4] = biject(N,-N+1);
	shift_cond[biject(-N,-N)][5] = 2;

	shift_cond[biject(-N,-N)][6] = biject(N,N);
	shift_cond[biject(-N,-N)][7] = 3;
}

void reflect_i(){
/*Sets to the reflecting coin the boundary positions where ''i = N'' and ''i = -N''*/
	int i,j;
	
		pos_coin[N+N][0] = -1;
		pos_coin[-N+N][0] = -1;
	
}
void reflect_j(){
/*Sets to the reflecting coin the boundary positions where ''j = N'' and ''j = -N''*/
	int i,j;
			
		pos_coin[N+N][1] = -1;
		pos_coin[-N+N][1] = -1;
	
}


void broken_link_1st_walker(int i, int j, int direction){
/*Sets broken link between points ''(i,j)'' and ''(i+1,j)'' where "i" is the position of walker 1 and "j" the position of walker 2 */
	int c, v[2];

	for(j=-N; j <= N; j++){
		
		for(c=0; c < Edges; c++){
			switch(shift_cond_altered[biject(i, j)][2*c+1]){
				case 0:
					shift_cond_altered[biject(i, j)][2*c+1] = 2;
					bijinv(v,shift_cond_altered[biject(i, j)][2*c]);
					shift_cond_altered[biject(i, j)][2*c] = biject(i,v[1]);
					break;
				case 1:
					shift_cond_altered[biject(i, j)][2*c+1] = 3;
					bijinv(v,shift_cond_altered[biject(i, j)][2*c]);
					shift_cond_altered[biject(i, j)][2*c] = biject(i,v[1]);
					break;
				case 2:
					shift_cond_altered[biject(i, j)][2*c+1] = 2;
					bijinv(v,shift_cond_altered[biject(i, j)][2*c]);
					break;
				case 3:
					shift_cond_altered[biject(i, j)][2*c+1] = 3;
					bijinv(v,shift_cond_altered[biject(i, j)][2*c]);
					break;
			}		
		}
		for(c=0; c < Edges; c++){
			
			switch(shift_cond_altered[biject(i+1, j)][2*c+1]){
				case 0:
					shift_cond_altered[biject(i+1, j)][2*c+1] = 0;
					bijinv(v,shift_cond_altered[biject(i+1, j)][2*c]);
					break;
				case 1:
					shift_cond_altered[biject(i+1, j)][2*c+1] = 1;
					bijinv(v,shift_cond_altered[biject(i+1, j)][2*c]);
					break;
				case 2:
					shift_cond_altered[biject(i+1, j)][2*c+1] = 0;
					bijinv(v,shift_cond_altered[biject(i+1, j)][2*c]);
					shift_cond_altered[biject(i+1, j)][2*c] = biject(i+1,v[1]);
					break;
				case 3:
					shift_cond_altered[biject(i+1, j)][2*c+1] = 1;
					bijinv(v,shift_cond_altered[biject(i+1, j)][2*c]);
					shift_cond_altered[biject(i+1, j)][2*c] = biject(i+1,v[1]);
					break;
			}		
		}
	}

	
}

void permanent_broken_link_1st_walker(int i, int j, int direction){
/*Sets broken link between points ''(i,j)'' and ''(i+1,j)'' where "i" is the position of walker 1 and "j" the position of walker 2 */
	int c, v[2];

	for(j=-N; j <= N; j++){
		
		for(c=0; c < Edges; c++){
			switch(shift_cond[biject(i, j)][2*c+1]){
				case 0:
					shift_cond[biject(i, j)][2*c+1] = 2;
					bijinv(v,shift_cond[biject(i, j)][2*c]);
					shift_cond[biject(i, j)][2*c] = biject(i,v[1]);
					break;
				case 1:
					shift_cond[biject(i, j)][2*c+1] = 3;
					bijinv(v,shift_cond[biject(i, j)][2*c]);
					shift_cond[biject(i, j)][2*c] = biject(i,v[1]);
					break;
				case 2:
					shift_cond[biject(i, j)][2*c+1] = 2;
					bijinv(v,shift_cond[biject(i, j)][2*c]);
					break;
				case 3:
					shift_cond[biject(i, j)][2*c+1] = 3;
					bijinv(v,shift_cond[biject(i, j)][2*c]);
					break;
			}		
		}
		for(c=0; c < Edges; c++){
			
			switch(shift_cond[biject(i+1, j)][2*c+1]){
				case 0:
					shift_cond[biject(i+1, j)][2*c+1] = 0;
					bijinv(v,shift_cond[biject(i+1, j)][2*c]);
					break;
				case 1:
					shift_cond[biject(i+1, j)][2*c+1] = 1;
					bijinv(v,shift_cond[biject(i+1, j)][2*c]);
					break;
				case 2:
					shift_cond[biject(i+1, j)][2*c+1] = 0;
					bijinv(v,shift_cond[biject(i+1, j)][2*c]);
					shift_cond[biject(i+1, j)][2*c] = biject(i+1,v[1]);
					break;
				case 3:
					shift_cond[biject(i+1, j)][2*c+1] = 1;
					bijinv(v,shift_cond[biject(i+1, j)][2*c]);
					shift_cond[biject(i+1, j)][2*c] = biject(i+1,v[1]);
					break;
			}		
		}
	}

	
}


void broken_link_2nd_walker(int l, int j, int direction){
/*Sets broken link between points ''(i,j)'' and ''(i,j+1)'' where "i" is the position of walker 1 and "j" the position of walker 2 */
	int v[2],c, i;

	
	for(i=-N; i <= N; i++){
		
		for(c=0; c < Edges; c++){
			switch(shift_cond_altered[biject(i, j)][2*c+1]){
				case 0:
					shift_cond_altered[biject(i, j)][2*c+1] = 1;
					bijinv(v,shift_cond_altered[biject(i, j)][2*c]);
					shift_cond_altered[biject(i, j)][2*c] = biject(v[0],j);
					break;
				case 1:
					shift_cond_altered[biject(i, j)][2*c+1] = 1;
					bijinv(v,shift_cond_altered[biject(i, j)][2*c]);
					break;
				case 2:
					shift_cond_altered[biject(i, j)][2*c+1] = 3;
					bijinv(v,shift_cond_altered[biject(i, j)][2*c]);
					shift_cond_altered[biject(i, j)][2*c] = biject(v[0],j);
					break;
				case 3:
					shift_cond_altered[biject(i, j)][2*c+1] = 3;
					bijinv(v,shift_cond_altered[biject(i, j)][2*c]);
					break;
			}		
		}
		for(c=0; c < Edges; c++){
			
			switch(shift_cond_altered[biject(i, j+1)][2*c+1]){
				case 0:
					shift_cond_altered[biject(i, j+1)][2*c+1] = 0;
					bijinv(v,shift_cond_altered[biject(i, j+1)][2*c]);
					break;
				case 1:
					shift_cond_altered[biject(i, j+1)][2*c+1] = 0;
					bijinv(v,shift_cond_altered[biject(i, j+1)][2*c]);
					shift_cond_altered[biject(i, j+1)][2*c] = biject(v[0],j+1);
					break;
				case 2:
					shift_cond_altered[biject(i, j+1)][2*c+1] = 2;
					bijinv(v,shift_cond_altered[biject(i, j+1)][2*c]);
					break;
				case 3:
					shift_cond_altered[biject(i, j+1)][2*c+1] = 2;
					bijinv(v,shift_cond_altered[biject(i, j+1)][2*c]);
					shift_cond_altered[biject(i, j+1)][2*c] = biject(v[0],j+1);
					break;
			}		
		}

	}
}

void permanent_broken_link_2nd_walker(int l, int j, int direction){
/*Sets broken link between points ''(i,j)'' and ''(i,j+1)'' where "i" is the position of walker 1 and "j" the position of walker 2 */
	int v[2],c, i;

	
	for(i=-N; i <= N; i++){
		
		for(c=0; c < Edges; c++){
			switch(shift_cond[biject(i, j)][2*c+1]){
				case 0:
					shift_cond[biject(i, j)][2*c+1] = 1;
					bijinv(v,shift_cond[biject(i, j)][2*c]);
					shift_cond[biject(i, j)][2*c] = biject(v[0],j);
					break;
				case 1:
					shift_cond[biject(i, j)][2*c+1] = 1;
					bijinv(v,shift_cond[biject(i, j)][2*c]);
					break;
				case 2:
					shift_cond[biject(i, j)][2*c+1] = 3;
					bijinv(v,shift_cond[biject(i, j)][2*c]);
					shift_cond[biject(i, j)][2*c] = biject(v[0],j);
					break;
				case 3:
					shift_cond[biject(i, j)][2*c+1] = 3;
					bijinv(v,shift_cond[biject(i, j)][2*c]);
					break;
			}		
		}
		for(c=0; c < Edges; c++){
			
			switch(shift_cond[biject(i, j+1)][2*c+1]){
				case 0:
					shift_cond[biject(i, j+1)][2*c+1] = 0;
					bijinv(v,shift_cond[biject(i, j+1)][2*c]);
					break;
				case 1:
					shift_cond[biject(i, j+1)][2*c+1] = 0;
					bijinv(v,shift_cond[biject(i, j+1)][2*c]);
					shift_cond[biject(i, j+1)][2*c] = biject(v[0],j+1);
					break;
				case 2:
					shift_cond[biject(i, j+1)][2*c+1] = 2;
					bijinv(v,shift_cond[biject(i, j+1)][2*c]);
					break;
				case 3:
					shift_cond[biject(i, j+1)][2*c+1] = 2;
					bijinv(v,shift_cond[biject(i, j+1)][2*c]);
					shift_cond[biject(i, j+1)][2*c] = biject(v[0],j+1);
					break;
			}		
		}

	}
}



void step(state *s1, state *s2, int t){
/*Performs one step of the quantum walk*/
int k,i,j;
	copy_shift_cond();
	link_impurities();

	if(varying_coin_each_step_1 == 1){
		reset_matrix(&random_coin);
		random_matrix_1D_31(0,randomt_1);
		//copy_matrix(randomt_1, random_coin);
	}

	if(varying_coin_each_step_2 == 1){
		reset_matrix(&random_coin);
		random_matrix_1D_32(0,randomt_2);
		//copy_matrix(randomt_2, random_coin);
	}

	coin2D(s1,s2, N);
	
	
	for(i=-N;i <= N; i++){
		for(j = -N; j <= N; j++){
			k = biject(i,j);
				M(   s2[  shift_cond_altered[k][0]].direction  ,  shift_cond_altered[k][1]   ,    0) = M(s1[k].direction,0,0);	
				M(   s2[  shift_cond_altered[k][2]].direction  ,  shift_cond_altered[k][3]   ,    0) = M(s1[k].direction,1,0);
				M(   s2[  shift_cond_altered[k][4]].direction  ,  shift_cond_altered[k][5]   ,    0) = M(s1[k].direction,2,0);
				M(   s2[  shift_cond_altered[k][6]].direction  ,  shift_cond_altered[k][7]   ,    0) = M(s1[k].direction,3,0);
		}
	}


copy_state(s1,s2);
}


void initialize_measure_ij(){
/*Initializes the structure ''initialize_measure_ij'' with entries equal to 0*/
	int i;
		if(dim > 0){
		i_0 = (int *)malloc(sizeof(int) * dim);
		
		j_0 = (int *)malloc(sizeof(int) * dim);
	
		
		for(i=0; i<dim; i++){
			i_0[i] = 0;
			j_0[i] = 0;
		}
	}
}

void initialize_shift_cond(){
/*Initializes the structure ''initialize_shift_cond'' with entries equal to 0*/
	int i,j;

	shift_cond = calloc(Nodes,sizeof(int *));

	for(i=0; i<Nodes; i++){
		shift_cond[i] = calloc(2*Edges,sizeof(int));
		for(j=0;j<2*Edges;j++)
			shift_cond[i][j] = 0;
	}

}

void initialize_shift_cond_altered(){
/*Initializes the structure ''initialize_shift_cond_altered'' with entries equal to 0*/
	int i,j;

	shift_cond_altered = calloc(Nodes,sizeof(int *));

	for(i=0; i<Nodes; i++){
		shift_cond_altered[i] = calloc(2*Edges,sizeof(int));
		for(j=0;j<2*Edges;j++)
			shift_cond_altered[i][j] = 0;
	}

}

void choose_coin_walker_1(k){
/*Resets matrix ''A_C'' according to the entry $k$ of ''pos_coin_altered''*/
	reset_matrix(&A_C);

	switch(pos_coin_altered[k][0]){
	case -2:
		identity(A_C);
		break;
	case -1:
		reflect_coin_1D(A_C);
		break;
	case 0: 

		Hadamard1D(A_C);

		break;
	case 1: 

		Grover1D(A_C);
		break;
	case 3: 

		random_matrix_1D(A_C);
		break;
	case 4:
		random_matrix_1D_21(k, A_C);
		break;
	case 5:
		copy_matrix(A_C, randomt_1);
		break;
	}
}

void choose_coin_walker_2(k){
/*Resets matrix ''B_C'' according to the entry $k$ of ''pos_coin_altered''*/
	reset_matrix(&B_C);
	switch(pos_coin_altered[k][1]){
	case -2:
		identity(B_C);
		break;
	case -1:
		reflect_coin_1D(B_C);
		break;
	case 0: 

		Hadamard1D(B_C);

		break;
	case 1: 

		Grover1D(B_C);
		break;
	case 3: 

		random_matrix_1D(B_C);
		break;
	case 4:
		random_matrix_1D_22(k, B_C);
		break;
	case 5:
		copy_matrix(B_C, randomt_2);
		break;
	}

}



void initialize_coin_w1(){
/*Allocate memory to pointer coins_w1. coins_w1 will store the coin operators of walker 1.*/
	int i;
	coins_w1 = (c_matrix *)malloc(sizeof(c_matrix) * Nodes1D);
		
	for(i=-N;i <= N; i++)
		new_matrix(&coins_w1[i+N],2,2);

}

void initialize_coin_w2(){
/*Allocate memory to pointer coins_w2. coins_w2 will store the coin operators of walker 2.*/
	int i;
	coins_w2 = (c_matrix *)malloc(sizeof(c_matrix) * Nodes1D);
	
		
	for(i=-N;i <= N; i++)
		new_matrix(&coins_w2[i+N],2,2);

}

void copy_position_coin(){
/*Copies all the entries of ''pos_coin'' to ''pos_coin_altered''*/
	int i;

	for(i=0; i< Nodes1D; i++){
		pos_coin_altered[i][0] = pos_coin[i][0];
		pos_coin_altered[i][1] = pos_coin[i][1];
	}
}

void coin2D(state *s1, state *s2, int t){
/*Applies the coin operator in each position of the line*/
	int i,j,k;
	t = N-1;
	copy_position_coin();
	coin_impurities();
	
	if(different_line == 0){
		for(i=-N; i <= N; i++){
			choose_coin_walker_1(i+N);
			copy_matrix(coins_w1[i+N], A_C);
			copy_matrix(coins_w2[i+N], A_C);
		}
	}
	else{
		for(i=-N; i <= N; i++){
			choose_coin_walker_1(i+N);
			copy_matrix(coins_w1[i+N], A_C);
			choose_coin_walker_2(i+N);
			copy_matrix(coins_w2[i+N], B_C);
		}
	}

	for(i=-N;i <= N; i++){
		for(j=-N;j <= N; j++){
			k = biject(i,j);

			reset_matrix(&random_coin);
			kron_A_B(coins_w1[i+N], coins_w2[j+N]);

			matrix_mult(s2[k].direction, random_coin, s1[k].direction);
		}
	}
	copy_state(s1,s2);
	
}



void initialize_position_coin(){
/*Initializes the structure ''initialize_position_coin'' with entries equal to 0*/
	int i;
	pos_coin = calloc(Nodes1D,sizeof(int *));

	for(i=0; i< Nodes1D; i++){
		pos_coin[i] = calloc(2,sizeof(int));
		pos_coin[i][0] = 0;//coin1;
		pos_coin[i][1] = 0;//coin2;

	}
}

void initialize_position_coin_altered(){
/*Initializes the structure ''initialize_position_coin_altered'' with entries equal to 0*/
	int i;
	pos_coin_altered = calloc(Nodes1D,sizeof(int *));;

	for(i=0; i< Nodes1D; i++){
		pos_coin_altered[i] = calloc(2,sizeof(int));
		pos_coin_altered[i][0] = 0; //coin1;
		pos_coin_altered[i][1] = 0; //coin2;
	}

}


void plot2D(char *filename,double *xvals, double *yvals, int NUM_POINTS){
/*Writes the data contained in xvals and yvals in file named in string "filename", at folder Output.*/
        char *aux = (char *)malloc(sizeof(char) * (strlen(filename) + 8));
        strcpy(aux,"Output/");
        strcat(aux, filename);
        FILE * temp = fopen(aux, "w");
        
        
        int i;
        
    	for (i=0; i < NUM_POINTS; i++)
    	{
    		fprintf(temp, "%lf %lf \n", xvals[i], yvals[i]); //Write the data to a temporary file
    	}
        
    	fclose(temp);
        free(aux);
        
    }
    
void plot3D(char *filename, double *zvals, int NUM_POINTS){
/*Writes the data contained in zvals (values at positions (i,j)) in file named in string "filename", at folder Output.*/
        
        char *aux = (char *)malloc(sizeof(char) * (strlen(filename) + 8));
        strcpy(aux,"Output/");
        strcat(aux, filename);
        FILE * temp = fopen(aux, "w");
        
        
        if(temp == NULL)
            printf("\n não abriu ficheiro nenhum \n");
        
        int i;
    	for (i=0; i < NUM_POINTS; i++)
    		fprintf(temp, "%lf %lf %lf \n", x[i], y[i], zvals[i]);
        
    	fclose(temp);
        
        free(aux);
    }

void Hadamard1D(c_matrix H){
/*Sets matrix "H" to Hadamard matrix*/
	M(H,0,0).re = 1.0/sqrt(2.0);
	M(H,0,1).re = 1.0/sqrt(2.0);
	M(H,1,0).re = 1.0/sqrt(2.0);
	M(H,1,1).re = -1.0/sqrt(2.0);

}

void Grover1D(c_matrix H){
/*Sets matrix "H" to Grover matrix*/
	M(H,0,0).re = 0.0;
	M(H,0,1).re = 1.0;
	M(H,1,0).re = 1.0;
	M(H,1,1).re = 0.0;
}


void Hadamard2D(){
/*Sets ''random_coin'' to Hadamard matrix*/
		M(random_coin,0,0).re = 1.0/2;
		M(random_coin,0,1).re = 1.0/2;
		M(random_coin,1,0).re = 1.0/2;
		M(random_coin,1,1).re = -1.0/2;

		M(random_coin,2,0).re = 1.0/2;
		M(random_coin,2,1).re = 1.0/2;
		M(random_coin,3,0).re = 1.0/2;
		M(random_coin,3,1).re = -1.0/2;

		M(random_coin,0,2).re = 1.0/2;
		M(random_coin,0,3).re = 1.0/2;
		M(random_coin,1,2).re = 1.0/2;
		M(random_coin,1,3).re = -1.0/2;

		M(random_coin,2,2).re = -1.0/2;
		M(random_coin,2,3).re = -1.0/2;
		M(random_coin,3,2).re = -1.0/2;
		M(random_coin,3,3).re = 1.0/2;


}


void Grover2D(){
/*Sets ''random_coin'' to Grover matrix*/
	M(random_coin,0,0).re = -1.0/2;
	M(random_coin,0,1).re =  1.0/2;
	M(random_coin,1,0).re =  1.0/2;
	M(random_coin,1,1).re = -1.0/2;

	M(random_coin,2,0).re =  1.0/2;
	M(random_coin,2,1).re =  1.0/2;
	M(random_coin,3,0).re =  1.0/2;
	M(random_coin,3,1).re =  1.0/2;

	M(random_coin,0,2).re =  1.0/2;
	M(random_coin,0,3).re =  1.0/2;
	M(random_coin,1,2).re =  1.0/2;
	M(random_coin,1,3).re =  1.0/2;

	M(random_coin,2,2).re = -1.0/2;
	M(random_coin,2,3).re =  1.0/2;
	M(random_coin,3,2).re =  1.0/2;
	M(random_coin,3,3).re = -1.0/2;

}

double random_0_1(){
/*Returns a random number from 0 to 1*/
	double a,b,c;

	a = rand();
	b = RAND_MAX;
	c = a/b;
	return c;

}

double random_0_Nodes(){
/*Returns a random integer number between 0 and the ''Nodes''*/
	double a,b,c;

	a = rand();
	b = RAND_MAX;
	c = a/b*Nodes;
	return c;

}

double random_01_1(){
/*Returns a random number from -1 to 1*/
	double a,b,c,d;

	a = rand();
	b = RAND_MAX;
	c = a/b;
	a = rand();
	d = a/b;
	if(d < 0.5)
		c = -c;

	return c;
}

void random_matrix_1D(c_matrix R){
/*Returns a random 2x2 unitary matrix*/
	double theta, zeta, xi;
	complex e_xi, e_zeta, cos_theta, sin_theta,c;

	xi = random_0_1()*M_PI/2.0;
	zeta = random_0_1()*M_PI/2.0;
	theta = random_0_1()*M_PI/2.0;

	e_xi.re = cos(xi);
	e_xi.im = sin(xi);
	e_zeta.re = cos(zeta);
	e_zeta.im = sin(zeta);
	cos_theta.re = cos(theta);
	cos_theta.im = 0.0;
	sin_theta.re = sin(theta);
	sin_theta.im = 0.0;

	c.re = -1.0;
	c.im = 0.0;

	M(R,0,0) = c_mult(e_xi, cos_theta);
	M(R,0,1) = c_mult(e_zeta, sin_theta);
	e_zeta.im = sin(-zeta);
	M(R,1,0) = c_mult(e_zeta, sin_theta);
	e_xi.im = sin(-xi);
	M(R,1,1) = c_mult(e_xi, c_mult(c,cos_theta));

}

void random_matrix_1D_21(int pos, c_matrix R){
/*Returns a random 2x2 unitary matrix*/
	double theta, zeta, xi;
	complex e_xi, e_zeta, cos_theta, sin_theta,c;

	theta = range1[pos][0] + random_0_1()*(range1[pos][1] - range1[pos][0]);
	zeta = range1[pos][2] + random_0_1()*(range1[pos][3] - range1[pos][2]);
	xi = range1[pos][4] + random_0_1()*(range2[pos][5] - range1[pos][4]);

	e_xi.re = cos(xi);
	e_xi.im = sin(xi);
	e_zeta.re = cos(zeta);
	e_zeta.im = sin(zeta);
	cos_theta.re = cos(theta);
	cos_theta.im = 0.0;
	sin_theta.re = sin(theta);
	sin_theta.im = 0.0;

	c.re = -1.0;
	c.im = 0.0;

	M(R,0,0) = c_mult(e_xi, cos_theta);
	M(R,0,1) = c_mult(e_zeta, sin_theta);
	e_zeta.im = sin(-zeta);
	M(R,1,0) = c_mult(e_zeta, sin_theta);
	e_xi.im = sin(-xi);
	M(R,1,1) = c_mult(e_xi, c_mult(c,cos_theta));

}

void random_matrix_1D_22(int pos, c_matrix R){
/*Returns a random 2x2 unitary matrix*/
	double theta, zeta, xi;
	complex e_xi, e_zeta, cos_theta, sin_theta,c;

	theta = range2[pos][0] + random_0_1()*(range2[pos][1] - range2[pos][0]);
	zeta = range2[pos][2] + random_0_1()*(range2[pos][3] - range2[pos][2]);
	xi = range2[pos][4] + random_0_1()*(range2[pos][5] - range2[pos][4]);

	e_xi.re = cos(xi);
	e_xi.im = sin(xi);
	e_zeta.re = cos(zeta);
	e_zeta.im = sin(zeta);
	cos_theta.re = cos(theta);
	cos_theta.im = 0.0;
	sin_theta.re = sin(theta);
	sin_theta.im = 0.0;

	c.re = -1.0;
	c.im = 0.0;

	M(R,0,0) = c_mult(e_xi, cos_theta);
	M(R,0,1) = c_mult(e_zeta, sin_theta);
	e_zeta.im = sin(-zeta);
	M(R,1,0) = c_mult(e_zeta, sin_theta);
	e_xi.im = sin(-xi);
	M(R,1,1) = c_mult(e_xi, c_mult(c,cos_theta));

}

void random_matrix_1D_31(int pos, c_matrix R){
    /*Returns a random 2x2 unitary matrix*/
	double theta, zeta, xi;
	complex e_xi, e_zeta, cos_theta, sin_theta,c;
    
	theta = range1_fixed[pos][0] + random_0_1()*(range1_fixed[pos][1] - range1_fixed[pos][0]);
	zeta = range1_fixed[pos][2] + random_0_1()*(range1_fixed[pos][3] - range1_fixed[pos][2]);
	xi = range1_fixed[pos][4] + random_0_1()*(range1_fixed[pos][5] - range1_fixed[pos][4]);
    
	e_xi.re = cos(xi);
	e_xi.im = sin(xi);
	e_zeta.re = cos(zeta);
	e_zeta.im = sin(zeta);
	cos_theta.re = cos(theta);
	cos_theta.im = 0.0;
	sin_theta.re = sin(theta);
	sin_theta.im = 0.0;
    
	c.re = -1.0;
	c.im = 0.0;
    
	M(R,0,0) = c_mult(e_xi, cos_theta);
	M(R,0,1) = c_mult(e_zeta, sin_theta);
	e_zeta.im = sin(-zeta);
	M(R,1,0) = c_mult(e_zeta, sin_theta);
	e_xi.im = sin(-xi);
	M(R,1,1) = c_mult(e_xi, c_mult(c,cos_theta));
    
}

void random_matrix_1D_32(int pos, c_matrix R){
    /*Returns a random 2x2 unitary matrix*/
	double theta, zeta, xi;
	complex e_xi, e_zeta, cos_theta, sin_theta,c;
    
	theta = range2_fixed[pos][0] + random_0_1()*(range2_fixed[pos][1] - range2_fixed[pos][0]);
	zeta = range2_fixed[pos][2] + random_0_1()*(range2_fixed[pos][3] - range2_fixed[pos][2]);
	xi = range2_fixed[pos][4] + random_0_1()*(range2_fixed[pos][5] - range2_fixed[pos][4]);
    
	e_xi.re = cos(xi);
	e_xi.im = sin(xi);
	e_zeta.re = cos(zeta);
	e_zeta.im = sin(zeta);
	cos_theta.re = cos(theta);
	cos_theta.im = 0.0;
	sin_theta.re = sin(theta);
	sin_theta.im = 0.0;
    
	c.re = -1.0;
	c.im = 0.0;
    
	M(R,0,0) = c_mult(e_xi, cos_theta);
	M(R,0,1) = c_mult(e_zeta, sin_theta);
	e_zeta.im = sin(-zeta);
	M(R,1,0) = c_mult(e_zeta, sin_theta);
	e_xi.im = sin(-xi);
	M(R,1,1) = c_mult(e_xi, c_mult(c,cos_theta));
    
}


void random_matrix_2D(){
/*Returns a random 4x4 unitary matrix*/
	c_matrix A,B;

	new_matrix(&A,2,2);
	new_matrix(&B,2,2);
	random_matrix_1D(A);
	random_matrix_1D(B);

	M(random_coin,0,0) = c_mult(M(A,0,0),M(B,0,0));
	M(random_coin,0,1) = c_mult(M(A,0,0),M(B,0,1));
	M(random_coin,1,0) = c_mult(M(A,0,0),M(B,1,0));
	M(random_coin,1,1) = c_mult(M(A,0,0),M(B,1,1));

	M(random_coin,0,2) = c_mult(M(A,0,1),M(B,0,0));
	M(random_coin,0,3) = c_mult(M(A,0,1),M(B,0,1));
	M(random_coin,1,2) = c_mult(M(A,0,1),M(B,1,0));
	M(random_coin,1,3) = c_mult(M(A,0,1),M(B,1,1));

	M(random_coin,2,0) = c_mult(M(A,1,0),M(B,0,0));
	M(random_coin,2,1) = c_mult(M(A,1,0),M(B,0,1));
	M(random_coin,3,0) = c_mult(M(A,1,0),M(B,1,0));
	M(random_coin,3,1) = c_mult(M(A,1,0),M(B,1,1));

	M(random_coin,2,2) = c_mult(M(A,1,1),M(B,0,0));
	M(random_coin,2,3) = c_mult(M(A,1,1),M(B,0,1));
	M(random_coin,3,2) = c_mult(M(A,1,1),M(B,1,0));
	M(random_coin,3,3) = c_mult(M(A,1,1),M(B,1,1));
	
	erase_matrix(&A);
	erase_matrix(&B);
}






void kron_A_B(c_matrix A,c_matrix B){
/*Performs kronecker product of matrices "A" and "B" and stores it in matrix "random_coin"*/
	M(random_coin,0,0) = c_mult(M(A,0,0),M(B,0,0));
	M(random_coin,0,1) = c_mult(M(A,0,0),M(B,0,1));
	M(random_coin,1,0) = c_mult(M(A,0,0),M(B,1,0));
	M(random_coin,1,1) = c_mult(M(A,0,0),M(B,1,1));

	M(random_coin,0,2) = c_mult(M(A,0,1),M(B,0,0));
	M(random_coin,0,3) = c_mult(M(A,0,1),M(B,0,1));
	M(random_coin,1,2) = c_mult(M(A,0,1),M(B,1,0));
	M(random_coin,1,3) = c_mult(M(A,0,1),M(B,1,1));

	M(random_coin,2,0) = c_mult(M(A,1,0),M(B,0,0));
	M(random_coin,2,1) = c_mult(M(A,1,0),M(B,0,1));
	M(random_coin,3,0) = c_mult(M(A,1,0),M(B,1,0));
	M(random_coin,3,1) = c_mult(M(A,1,0),M(B,1,1));

	M(random_coin,2,2) = c_mult(M(A,1,1),M(B,0,0));
	M(random_coin,2,3) = c_mult(M(A,1,1),M(B,0,1));
	M(random_coin,3,2) = c_mult(M(A,1,1),M(B,1,0));
	M(random_coin,3,3) = c_mult(M(A,1,1),M(B,1,1));
	
}


void rho_x_coin(c_matrix phi_k_1, c_matrix v){
/*partial trace over y of the density matrix of partial scalar product of the overall quantum state and an eigen vector ''v'' of the coin state */
	int i,j,k,c1, c2;
	complex c, ijkc, kc, kc_;
	double test;
	c.re=0.0;
	c.im =0.0;
	reset_matrix(&phi_k_1);

	for(i=-N; i <= N; i++){
		for(k=-N; k<=N; k++){

			ijkc.re = 0.0;
			ijkc.im = 0.0;
			for(j=-N; j<=N; j++){
				kc.re = 0.0;
				kc.im = 0.0;
				kc_.re = 0.0;
				kc_.im = 0.0;
				for(c1=0; c1 < Edges; c1++){
					kc = c_sum(kc, c_mult(c_conj(M(v,c1,0)) , M(s1[biject(i,j)].direction,c1,0 )));
					kc_ = c_sum(kc_, c_mult(M(v,c1,0) , c_conj(M(s1[biject(k,j)].direction,c1,0 ))));	
				}
				ijkc = c_sum(ijkc , c_mult(kc,kc_));			
 			}
			M(phi_k_1,(i+N),(k+N)) = c_sum(M(phi_k_1,(i+N),(k+N)) , ijkc);
		}
	}
	

	for(i=-N; i<=N; i++){
		c = c_sum (c ,M(phi_k_1,(i+N),(i+N)));
	}

}




double Entanglement_of_formation(){
/*Upper bound of the overall state entanglement of formation*/
	int i, j,l,h;
	c_matrix lambda, V,vi, pi_pit, R, p1_tracet;
	complex c[Edges],a;
	double b = 0.0;
	double aux;
	int tm = 0; 
	int tf = 0;


	new_matrix(&p1_tracet, 1,2*N+1);
	tm++;
	new_matrix(&vi,Edges,1);
	tm++;
	new_matrix(&lambda,Edges,Edges);
	tm++;	
	new_matrix(&V, Edges, Edges);
	tm++;	
	new_matrix(&pi_pit, 2*N+1, 2*N+1);
	tm++;
	QR_algorithm2(rho_c,lambda, V);
	for(j=0; j< Edges; j++){
		aux = sqrt(sqrt(c_abs_sq(M(lambda,j,j))));

		if(c_abs_sq(M(lambda,j,j)) <= pow(10.0, -10.0))
			aux = 1.0;

		c[j].re = 0.0;
		c[j].im = 0.0;
		for(i=0; i < Edges; i++){
			M(vi,i,0) = c_conj(M(V,j,i));
			M(vi,i,0).re = M(vi,i,0).re/aux;
			M(vi,i,0).im = M(vi,i,0).im/aux;
		}
	
		rho_x_coin(pi_pit,vi);
		R = QR_algorithm(pi_pit);

		tm++;
		c[j] = entropy(R);
		erase_matrix(&R);
		tf++;
	}

	
	for(j=0; j< Edges; j++){
		b = b + c[j].re*M(lambda,j,j).re;
	}
	erase_matrix(&p1_tracet);
	tf++;
	erase_matrix(&lambda);
	tf++;
	erase_matrix(&V);
	tf++;
	erase_matrix(&vi);
	tf++;
	erase_matrix(&pi_pit);
	tf++;
	return b;
	
}



void mean_value_xy(state *s){
/*Calculates the mean value of particles 1 and 2 at each step*/
	int i,c,x_y[2];
	double aux;

	mean_x[t_counter] = 0.0;
	mean_y[t_counter] = 0.0;
	for(i=0;i< Nodes; i++){
		bijinv(x_y,i);
		for(c=0; c < Edges; c++){
			aux = (double) x_y[0];
			mean_x[t_counter] = mean_x[t_counter] + aux * c_abs_sq(M(s[i].direction,c,0));
			aux = (double) x_y[1];
			mean_y[t_counter] = mean_y[t_counter] + aux * c_abs_sq(M(s[i].direction,c,0));
		}
		
	}
}

void copy_shift_cond(){
/*Copies ''shift_cond'' into ''shift_cond_altered''*/
	int i,j;

	for(i=0;i< Nodes; i++)
		for(j=0;j< 2*Edges; j++){
			shift_cond_altered[i][j] = shift_cond[i][j];
		}
}


void link_impurities(){
/*Randomly choose broken link according to "brokenlink_imp_1", "brokenlink_imp_2" and "different_line"*/
	int k=0;
	int i,j,c,c2, x_,y_, coin_state1, coin_state2;// i_0[2*N+1];
	x_ = N;
	y_ = N;


	if(different_line == 0){
		for(i=-N+1;i<N-1; i++){
			//	i_0[i+N] = 0;
				if((random_0_1() < index_brokenlink_imp_1) || (ib[i] == 1) ){

			//		i_0[i+N]=1;
					broken_link_1st_walker(i,k, 0);
                    			broken_link_2nd_walker(k,i,0);

				}
		}
		for(i=-N+1;i<N-1; i++){
				//if(i_0[i+N] == 1){
				//	broken_link_2nd_walker(k,i,0);
				//}
		}
	}
	else{
		for(i=-N+1;i<N-1; i++){
				if((random_0_1() < index_brokenlink_imp_1) ||(ib[i] == 1) )
					broken_link_1st_walker(i,k, 0);		
		}
		for(i=-N+1;i<N-1; i++){
				if(random_0_1() < index_brokenlink_imp_2 || (jb[i] == 1))
					broken_link_2nd_walker(k,i,0);
		
		}
	}

}


void coin_impurities(){
/*Randomly chooses positions for random coins with respect to ''index_coin_imp1'', ''index_coin_imp2'' and "different_line"*/
	int i,k, x_,y_;
	x_ = N+1;
	y_ = N+1;
	if(different_line == 1){
		for(i=-N+1;i<N; i++){
			if(random_0_1() < index_coin_imp_1)
				pos_coin_altered[i+N][0] = 3;

			if(random_0_1() < index_coin_imp_2)
				pos_coin_altered[i+N][1] = 3;
		}
	}
	else{
		for(i=-N+1;i<N; i++){
				if(random_0_1() < index_coin_imp_1){
					pos_coin_altered[i+N][0] = 3;
					pos_coin_altered[i+N][1] = 3;
				}
		}
	}

}


void covariance_xy(state *s){
/*Determines the covariance between particles 1 and 2 (<ij>-<i> <j>) in each step*/
	int i,c,x_y[2];
	double aux1, aux2;


	for(i=0;i< Nodes; i++){
		bijinv(x_y,i);
		aux1 = (double) x_y[0];
		aux2 = (double) x_y[1];
		xy_cov[t_counter] = xy_cov[t_counter] + aux1*aux2*P[i];
	}
	xy_cov[t_counter] = xy_cov[t_counter] - mean_x[t_counter]*mean_y[t_counter];

}

void mean_distance_xy(state *s){
/*Determines <x-y> in each step*/
	int i,j;
	for(i=-N; i <= N; i++)
		for(j=-N; j <= N; j++)
			mean_dist[t_counter ] = mean_dist[t_counter ] + ((double) abs( i-j))*P[biject(i,j)];
}

void mutual_info_Shannon(){
/*Shannon mutual information of variables i and j*/
	I_xy[t_counter] = H_y[t_counter] + H_x[t_counter] - H_xy[t_counter];

}

void mutual_info_von_Neumann(){
/*von Neuman mutual information of the positions of particles 1 and 2*/
	Iv_xy[t_counter] = S_y[t_counter] + S_x[t_counter] - S_c[t_counter];

}

void entropy_xy(){
/*Shannon entropy for position probability distribution*/
	int i;

	for(i=0; i< Nodes; i++)
		if(P[i] != 0.0)
			H_xy[t_counter] = H_xy[t_counter] - (P[i])*log(P[i])/log(2);
}

void entropy_x(){
/*Shannon entropy for particle 1*/
	int i,j;
	double px;
	for(i=-N; i<= N; i++){
		px = 0.0;
		for(j = -N; j <= N; j++ ){
			px = px + P[biject(i,j)];
		}
		if(px != 0.0)
			H_x[t_counter] = H_x[t_counter] - px*log(px)/log(2);
	}
}


void entropy_y(){
/*Shannon entropy for particle 2*/
	int i,j;
	double py;
	for(j=-N; j<= N; j++){
		py = 0.0;
		for(i = -N; i <= N; i++ ){
			py = py + P[biject(i,j)];
		}
		if(py != 0.0)
			H_y[t_counter] = H_y[t_counter] - py*log(py)/log(2);
	}
}


void initialize_global_variables(){
/*Initializes all the global variables and structures*/
	int i;

	Nodes1D = 2*N + 1;
	Nodes2D = Nodes1D * Nodes1D;
	Nodes = Nodes2D;

	new_matrix(&Y_pos,2*N+1,2*N+1);
	new_matrix(&Y_lambda,2*N+1,2*N+1);
	new_matrix(&random_coin,4,4);
	new_matrix(&A_C,2,2);
	new_matrix(&B_C,2,2);
	new_matrix(&randomt_1,2,2);
	new_matrix(&randomt_2,2,2);
	new_matrix(&rho_c,Edges, Edges);
	new_matrix(&X_p, Nodes1D, Nodes1D);
	new_matrix(&Y_p, Nodes1D, Nodes1D);
    
	P_av = (double *)malloc(sizeof(double) * Nodes);
	P = (double *)malloc(sizeof(double) * Nodes);
	S_c = (double *)malloc(sizeof(double) * steps);
	H_c = (double *)malloc(sizeof(double) * steps);
	S_x = (double *)malloc(sizeof(double) * steps);
	S_y = (double *)malloc(sizeof(double) * steps);
	S_y_pi = (double *)malloc(sizeof(double) * steps);
	E_f = (double *)malloc(sizeof(double) * steps);
	H_x = (double *)malloc(sizeof(double) * steps);
	H_y = (double *)malloc(sizeof(double) * steps);
	H_xy = (double *)malloc(sizeof(double) * steps);
	I_xy = (double *)malloc(sizeof(double) * steps);
	quantum_discord = (double *)malloc(sizeof(double) * steps);
	Iv_xy = (double *)malloc(sizeof(double) * steps);
	mean_x = (double *)malloc(sizeof(double) * steps);
	mean_y = (double *)malloc(sizeof(double) * steps);
	xy_cov = (double *)malloc(sizeof(double ) * steps);
	x = (double *)malloc(sizeof(double ) * Nodes2D);
	y = (double *)malloc(sizeof(double ) * Nodes2D);
	xy_cov = (double *)malloc(sizeof(double ) * steps);
	mean_dist = (double *)malloc(sizeof(double )* steps);
	measure_prob = (double *)malloc(sizeof(double )* steps);
	norm_factor = (double *)malloc(sizeof(double )* steps);
	one_shot_prob_hit = (double *)malloc(sizeof(double )* steps);
	av_hitting_time = (double *)malloc(sizeof(double )* steps);
	first_time_hit = (double *)malloc(sizeof(double )* steps);
	concurrent_hit_time = (double *)malloc(sizeof(double )* steps);

	initialize_coin_w1();
	initialize_coin_w2();
	initialize_shift_cond();
	initialize_shift_cond_altered();
	shift2D(s1,s2,t_counter);
	initialize_position_coin();
	initialize_position_coin_altered();

	for(i=0;i< Nodes; i++){
		P_av[i] = 0.0;
		P[i] = 0.0;
	}
	for(i=0;i< steps; i++){
		S_c[i] = 0.0;
		H_c[i] = 0.0;
		S_x[i] = 0.0;
		S_y[i] = 0.0;
		H_x[i] = 0.0;
		H_y[i] = 0.0;
		H_xy[i] = 0.0;
		I_xy[i] = 0.0;
		S_y_pi[i] = 0.0;
		mean_x[i] = 0.0;
		mean_y[i] = 0.0;
		xy_cov[i] = 0.0;
		mean_dist[i] = 0.0;
		measure_prob[i] = 0.0;
		norm_factor[i] = 1.0;
		one_shot_prob_hit[i] = 0.0;
		av_hitting_time[i] = 0.0;
		first_time_hit[i] = 0.0;
		concurrent_hit_time[i] = 0.0;
		quantum_discord[i] = 0.0;

	}

	for(i=0;i< Nodes2D; i++){
		x[i]=0;
		y[i]=0;
	}

	s1 = initialize_state();
	s2 = initialize_state();


}

void plots(){
/*Plots in GNUplot the quantities determined in each step of the quantum walk*/
	int i,j, r[2];

	for(i=0; i < Nodes2D; i++)
		P_av[i] = P_av[i]/((double) steps);

	for(i=0;i < Nodes2D; i++){
		bijinv(r,i);
		x[i] = r[0];
		y[i] = r[1];
	}
	plot3D("Average_Probability_Distribution\0",P_av, Nodes);
	plot3D("Probability_Distribution\0",P, Nodes);
	double *k = (double *)malloc(sizeof(double ) * steps);

	for(i=0 ; i<steps ; i++)
		k[i] = (double) i;

	plot2D("Von_Newman_entropy of coin state\0",k,S_c,steps);

	plot2D("Shannon_entropy of coin state\0",k,H_c,steps);
	

	plot2D("mean_x\0",k,mean_x,steps);
	plot2D("mean_y\0",k,mean_y,steps);
	plot2D("xy_cov\0",k, xy_cov, steps);
	plot2D("mean_dist\0",k, mean_dist, steps);
	plot2D("S_x\0",k, S_x, steps);
	plot2D("S_y\0",k, S_y, steps);
	plot2D("H_x\0",k, H_x, steps);
	plot2D("H_y\0",k, H_y, steps);
	plot2D("E_f\0",k, E_f, steps);
	plot2D("I_xy\0",k, I_xy, steps);
	plot2D("Iv_xy\0",k, Iv_xy, steps);
	plot2D("H_xy\0",k, H_xy, steps);

	plot2D("one Shot probability to hit\0",k,one_shot_prob_hit, steps);
	plot2D("Average hitting time \0",k, av_hitting_time, steps);
	plot2D("First time to hit\0",k, first_time_hit, steps);
	plot2D("Concurrence hitting time\0",k, concurrent_hit_time, steps);
	plot2D("entropy of truncated measurements for particle 2 given particle 1\0",k, S_y_pi, steps);
	plot2D("Quantum Discord\0",k, quantum_discord, steps);

	
}


void sum_prop_pos(state *s){
/*Determines the entries of ''P'', probability of the particle beeing at each position*/
	int i,j;
	for(i=0;i < Nodes; i++){
		P[i] = 0.0;
		for(j=0;j < Edges;j++)
			P[i] = P[i] + c_abs_sq(M(s[i].direction,j,0));
	}
}

double prob_test(state *s){
/*Sums all the entries of ''P''*/
	int i,j;
	double p=0.0;
	for(i=0;i < Nodes; i++)
		for(j=0;j < Edges;j++)
			p = p + c_abs_sq(M(s[i].direction,j,0));

	return p;
}

void P_av_actualization(){
/*Calculates at each step the average probability distribution of the particles positions*/
	int i;
	for(i=0;i< Nodes; i++)
		P_av[i] = P_av[i] + P[i];
}


void qw_on_plane(state *s1, state *s2){
/*Performs the quantum walk and determines all the quantities*/
	complex c;
	int t;

	c_matrix diag_rho_c, diag_x, diag_y;

	for(t_counter=0;t_counter < steps; t_counter++){
        	printf("step = %d \n ", t_counter);
			if(t_counter > N)
				t = N;
			else
				t = t_counter;


			if(t_counter > 0){

				step(s1, s2, t);
				measure_options(m_options, t_counter, 2);
				copy_state_normalized(s1,s2, t_counter);
			}

            
			
			mean_distance_xy(s1);
			mean_value_xy(s1);
			covariance_xy(s1);

			
			sum_prop_pos( s1);
			P_av_actualization();

			position_trace_2D(s1,N);
			diag_rho_c = QR_algorithm(rho_c);
			H_c[t_counter] = entropy(rho_c).re;
			S_c[t_counter] = entropy(diag_rho_c).re;

			erase_matrix(&diag_rho_c);

			entropy_x();
			entropy_y();
			entropy_xy();
			mutual_info_Shannon();



		if(vn_x == 1){
			rho_x();
			diag_x = QR_algorithm(X_p);
			S_x[t_counter] = entropy(diag_x).re;
			erase_matrix(&diag_x);
		}

		if(vn_y == 1){
			rho_y();
			diag_y = QR_algorithm(Y_p);
			S_y[t_counter] = entropy(diag_y).re;
			erase_matrix(&diag_y);
		}
		
		if(vn_m_info == 1 && vn_x == 1 && vn_y == 1){
			mutual_info_von_Neumann();
		}
		if(q_discord == 1 && vn_m_info == 1 && vn_x == 1 && vn_y == 1){
			entropy_truncated_measurements(t_counter);
			quantum_discord_calc(t_counter);
		}
		if(ent_form == 1)
			E_f[t_counter] = Entanglement_of_formation();

	}
}



void kron(c_matrix R, c_matrix A, c_matrix B){
/*Performs kronecker product of matrices "A" and "B" and stores it in matrix "R"*/
int i1, j1, i2, j2, i,j;

	for(i1=0; i1< A.rows; i1++)
		for(j1=0; j1< A.cols; j1++)
			for(i2=0; i2< B.rows; i2++)
				for(j2=0; j2< B.rows; j2++){
					i = i1*(B.rows)+i2;
					j = j1*(B.cols)+j2;
					M(R,i,j) = c_mult(M(A,i1,j1), M(B,i2,j2));
				}

}

void reflect_coin_1D(c_matrix R){
/*The reflecting coin is set on matrix "R" of dimensions 2x2.*/
	M(R,0,0).re = 0.0;
	M(R,0,1).re = 1.0;
	M(R,1,0).re = 1.0;
	M(R,1,1).re = 0.0;
	M(R,0,0).im = 0.0;
	M(R,0,1).im = 0.0;
	M(R,1,0).im = 0.0;
	M(R,1,1).im = 0.0;

}

void reflect_coin_2D(c_matrix R){
/*The reflecting coin is set on matrix "R" of dimensions 4x4.*/
	c_matrix R1D, id, R1D_id, id_R1D;
	
	new_matrix(&R1D,2,2);
	new_matrix(&id,2,2);
	new_matrix(&R1D_id,4,4);
	new_matrix(&id_R1D,4,4);

	identity(id);
	reflect_coin_1D(R1D);
	kron(R1D_id, R1D, id);
	//print_matrix(R1D_id);
	kron(id_R1D, id, R1D);
	//print_matrix(R1D_id);
	matrix_mult(R, id_R1D, R1D_id);
	////////printf("R\n");
	//print_matrix(R);
	erase_matrix(&R1D);
	erase_matrix(&id);
	erase_matrix(&R1D_id);
	erase_matrix(&id_R1D); 
}




void rho_x(){
/*Determines the density matrix of particle 1.*/
	int i,j,k,c1;
	complex c;
	double test;
	c.re=0.0;
	c.im =0.0;
	reset_matrix(&X_p);
	for(i=-N; i <= N; i++){
		for(k=-N; k<=N; k++){
			for(j=-N; j<=N; j++){
				for(c1=0; c1 < Edges; c1++){
					M(X_p,(i+N),(k+N)) = c_sum(M(X_p,(i+N),(k+N)) , c_mult( M( s1[biject(i,j)].direction,c1, 0) , c_conj(M( s1[biject(k,j)].direction,c1, 0 ))));
					
				}

 			}
		}
	}
	

	for(i=-N; i<=N; i++){
		c = c_sum (c ,M(X_p,(i+N),(i+N)));

	}

}

void rho_y(){
/*Determines the density matrix of particle 2.*/
	int i,j,k,c1;
	complex c;
	double test;
	c.re=0.0;
	c.im =0.0;
	reset_matrix(&Y_p);

	for(j=-N; j <= N; j++){
		for(k=-N; k<=N; k++){
			for(i=-N; i<=N; i++){
				for(c1=0; c1 < Edges; c1++){
					M(Y_p,(j+N),(k+N)) = c_sum(M(Y_p,(j+N),(k+N)) , c_mult( M( s1[biject(i,j)].direction,c1, 0) , c_conj(M( s1[biject(i,k)].direction,c1, 0 ))));
					
				}

 			}
		}
	}
	

	for(i=-N; i<=N; i++){
		c = c_sum (c ,M(Y_p,(i+N),(i+N)));
	}

}

void copy_state_normalized(state *s1, state*s2,int t){
/*Copies to ''s1'' the normalized state of ''s2''*/
int i,e,j;
double sr_norm_factor = sqrt(norm_factor[t]);
	for(i=0;i<Nodes ;i++){
		for(e=0; e < Edges; e++){
			M(s1[i].direction,e,0).re = M(s2[i].direction,e,0).re/sr_norm_factor;
			M(s1[i].direction,e,0).im = M(s2[i].direction,e,0).im/sr_norm_factor;
		}
	}
}

void P_0(int t){
/*Updates the probability of the particle being in position specified by structures ''i_0'' and ''j_0'' and stores it in ''measure_prob'' */
	int c;
	int i,j;
	measure_prob[t] = 0.0;
	for(i=0; i < dim; i++)
		
			for(c=0; c< Edges; c++){
				measure_prob[t] = measure_prob[t] + c_abs_sq(M(s2[biject(i_0[i], j_0[i])].direction,c,0));
			}
}

void norm_factor_act(int t ){
/*Updates norm_factor*/
	int c;
	
	norm_factor[t] = 1.0 - measure_prob[t];

}

void P_1(int t){
/*Performs quantum measurement on positions specified in ''i_0'' and ''j_0''*/
	int c;
	int i,j;
	measure_prob[t] = 0.0;
	for(i=0; i < dim; i++)
		
		for(c=0; c< Edges; c++){
			measure_prob[t] = measure_prob[t] + c_abs_sq(M(s2[biject(i_0[i], j_0[i])].direction,c,0));
			M(s2[biject(i_0[i], j_0[i])].direction,c,0).re = 0.0;
			M(s2[biject(i_0[i], j_0[i])].direction,c,0).im = 0.0;
		}
}

void one_shot_probability_hit(int t){
/*Determines one shot probability to hit.*/
	P_0(t);
	one_shot_prob_hit[t] = measure_prob[t];
}

void first_time_prob_hit(int t_0, int t){
/*Calculates first time probability to hit.*/
	if( t >= t_0){
		P_1( t);
		norm_factor_act( t);
		//P_1( t);
	}
	first_time_hit[t] = measure_prob[t];
}

void average_hitting_time(int t){
/*Determines the average hitting time*/
	if(t > 0)
		av_hitting_time[t] = av_hitting_time[t-1] + first_time_hit[t]*((double) t);
	else
		av_hitting_time[t] = first_time_hit[t]*((double) t);

}


void concurrent_hitting_time(int t){
/*Determines the concurrent hitting time.*/
	if(t > 0)
		concurrent_hit_time[t] = concurrent_hit_time[t-1] + first_time_hit[t];
	else
		concurrent_hit_time[t] = first_time_hit[t];
}

void measure_options(int m, int t, int t_0){
/*According to the value of "m", especific hitting time is set.*/
int z;
	if(m == 1){
			one_shot_probability_hit(t);
	}

	if (m == 2){
			first_time_prob_hit( t_0, t);
	}
	
	if(m == 3){
			first_time_prob_hit( t_0, t);
			average_hitting_time(t);
			concurrent_hitting_time(t);
	}

}

void Truncated_position_measurement(int i, c_matrix Y_pos){
/*Determines the density matrix ("Y_pos") of variable y given measurements in variable x*/
	int j1,j2, c1;

	for(j1=-N; j1 <= N ; j1++){
		for(j2=-N; j2 <= N ; j2++){
			M(Y_pos,(j1+N),(j2+N)).re = 0.0;
			M(Y_pos,(j1+N),(j2+N)).im = 0.0;
			for(c1=0; c1 < Edges; c1++){
				M(Y_pos,(j1+N),(j2+N)) = c_sum(M(Y_pos,(j1+N),(j2+N)) , c_mult( M( s1[biject(i,j1)].direction,c1, 0) , c_conj(M( s1[biject(i,j2)].direction,c1, 0 ))));
			}
		}
 	}
}

void normalize_truncated_position_measurement(c_matrix Y_pos, double prob){
/*Normalizes the truncated measurement*/
	int j1,j2, c1;

	
	for(j1=-N; j1<=N; j1++){
		for(j2=-N; j2 <= N; j2++){
			M(Y_pos,(j1+N),(j2+N)).re = M(Y_pos,(j1+N),(j2+N)).re/prob;
			M(Y_pos,(j1+N),(j2+N)).im = M(Y_pos,(j1+N),(j2+N)).im/prob;
		}
	}
}

void entropy_truncated_measurements(int t){
/*Calculates the entropy of the density matrix after truncated measurement*/
	int i,k,c;

	complex aux;
	double prob;

	for(i=-N; i<= N; i++){
		Truncated_position_measurement(i,Y_pos);
		
		prob = 0.0;

	for(k=-N; k < N+1; k++)
			for(c=0; c< 4; c++)
			prob = prob + c_abs_sq(M(s1[biject(i,k)].direction,c,0));//M(Y_pos,k,k).re;

		if(prob != 0.0)
			normalize_truncated_position_measurement(Y_pos, prob);


		Y_lambda = QR_algorithm(Y_pos);
		aux = entropy(Y_lambda);
		aux.re = prob*(aux.re);
		aux.im = prob*(aux.im);
		S_y_pi[t] = S_y_pi[t] + aux.re;
		erase_matrix(&Y_lambda);
	}

}


void quantum_discord_calc(int t){
/*Calculates the quantum discord*/
	quantum_discord[t] = S_x[t]-S_c[t] + S_y_pi[t];
}


void gnu_plot2D(char *filename,double *xvals, double *yvals, int NUM_POINTS){
/*Plots the data contained in file "filename", using GNUplot.*/
    FILE * gnuplotPipe = popen ("gnuplot -persistent", "w");
    
    int i;
    
	fprintf(gnuplotPipe,"set title '%s'\n", filename);
    fprintf(gnuplotPipe, "plot '%s' with l\n", filename);
	pclose(gnuplotPipe);
}

void gnu_plot3D(char *filename, double *zvals, int NUM_POINTS){
/*Plots the data contained in file "filename", using GNUplot*/
    FILE * gnuplotPipe = popen ("gnuplot -persistent", "w");
    
    
    fprintf(gnuplotPipe, "set dgrid3d %d %d %d \n",2*Nodes, 2*Nodes, 1);
	fprintf(gnuplotPipe, "set hidden3d \n",filename);
    fprintf(gnuplotPipe, "splot '%s'  u 1:2:3 with lines \n",filename);
	fprintf(gnuplotPipe, "set output \"color.eps\"\n",filename);
    pclose(gnuplotPipe);
    
}


void gnu_plots(){
/*Plots every data file contained in folder Output*/
	int i,j, r[2];
    
	for(i=0; i < Nodes2D; i++)
		P_av[i] = P_av[i]/((double) steps);
    
	for(i=0;i < Nodes2D; i++){
		bijinv(r,i);
		x[i] = r[0];
		y[i] = r[1];
	}
	gnu_plot3D("Output/Average_Probability_Distribution\0",P_av, Nodes);
	gnu_plot3D("Output/Probability_Distribution\0",P, Nodes);
	double k[steps];
    
	for(i=0 ; i<steps ; i++)
		k[i] = (double) i;
    
	gnu_plot2D("Output/Von_Newman_entropy of coin state\0",k,S_c,steps);
    
	gnu_plot2D("Output/Shannon_entropy of coin state\0",k,H_c,steps);
	
    
	gnu_plot2D("Output/mean_x\0",k,mean_x,steps);
	gnu_plot2D("Output/mean_y\0",k,mean_y,steps);
	gnu_plot2D("Output/xy_cov\0",k, xy_cov, steps);
	gnu_plot2D("Output/mean_dist\0",k, mean_dist, steps);
	gnu_plot2D("Output/S_x\0",k, S_x, steps);
	gnu_plot2D("Output/S_y\0",k, S_y, steps);
	gnu_plot2D("Output/H_x\0",k, H_x, steps);
	gnu_plot2D("Output/H_y\0",k, H_y, steps);
	gnu_plot2D("Output/E_f\0",k, E_f, steps);
	gnu_plot2D("Output/I_xy\0",k, I_xy, steps);
	gnu_plot2D("Output/Iv_xy\0",k, Iv_xy, steps);
	gnu_plot2D("Output/H_xy\0",k, H_xy, steps);
    
	gnu_plot2D("Output/one Shot probability to hit\0",k,one_shot_prob_hit, steps);
	gnu_plot2D("Output/Average hitting time \0",k, av_hitting_time, steps);
	gnu_plot2D("Output/First time to hit\0",k, first_time_hit, steps);
	gnu_plot2D("Output/Concurrence hitting time\0",k, concurrent_hit_time, steps);
	gnu_plot2D("Output/entropy of truncated measurements for particle 2 given particle 1\0",k, S_y_pi, steps);
	gnu_plot2D("Output/Quantum Discord\0",k, quantum_discord, steps);
    
	
	printf("\n test = %lf\n",prob_test(s1));
    
    
}







