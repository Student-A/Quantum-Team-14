#ifndef _sim_H_
#define	_sim_H_

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

/*********This file contains all the global variables and functions/subroutines prototypes of file qwsim.c**********/


/******************Structures*****************************/
typedef struct complex{
double re;
double im;
}complex;

typedef struct c_matrix{
	int rows;
	int cols;
	complex *t;
}c_matrix;

/*This macro is to be used with the structure c_matrix
 *m(A,i,j) accesses the line i and column j entry of matrix A
 * */
#define M(m,x,y) m.t[x+y*m.rows]


typedef struct state{
c_matrix direction;
}state;



/*******************************Global Variables/structures********************************/

//s1 and s2 are the overall state. s1 is the main structure and the operations on s1 are stored in s2. at the end of a step, s2 is copied to s1. This will avoid erasing values when performing shift operation.
state *s1, *s2;
//Number of steps for Quantum Walk
int steps;
int fixed_coin_1, fixed_coin_2;
//Dimension of the lattice
int N;
//Number of nodes on the lattice
int Nodes;
//Maximum number of positions of each particle
int Nodes1D;
//Total number of nodes on the lattice
int Nodes2D;
//Number of edges per node (Only 4 is used)
int Edges;
//Counter for the number of steps performed
int t_counter;
//Pointer for matrix that stores the code for the coin operator at each position
int **pos_coin;
//Pointer for matrix that stores the code for the coin operator at each position and per step alterations
int **pos_coin_altered;
//Memory test variables
int N_malloc;
int N_free;
//Choice for different line or same line
int different_line;
//Stores the options for the mixing times.
int m_options;
//Choice for calculating von Neumann Entropy of first particle's position
int vn_x;
//Choice for calculating von Neumann Entropy of second particle's position
int vn_y;
//Choice for calculating von Neumann mutual information of the particles positions
int vn_m_info;
//Choice for calculating quantum discord
int q_discord;
//Choice for calculating upper bound for entanglement of formation
int ent_form;
//Choice of coin to apply on all positions of the line for the first particle
int coin1;
//Choice of coin to apply on all positions of the line for the second particle
int coin2;
//Choice for reflecting boundary of first particle
int rf1;
//Choice for reflecting boundary of second particle
int rf2;
//Choice for ploting output data with gnuplot
int gnu;
//number of fixed broken links
int fixed_broken_link_i, fixed_broken_link_j;
//Average position probability over time
double *P_av;
//Position probability
double *P;
double **range1, **range2, **range1_fixed, **range2_fixed;
//Pointer for vector that von Neumann entropy of coin state at each step
double *S_c;
//Pointer for vector that von Neumann entropy of first particle's position at each step
double *S_x;
//Pointer for vector that von Neumann entropy of second particle's position at each step
double *S_y;
//Pointer for vector that von Neumann mutual information of variables x and y at each step
double *Iv_xy;
//Pointer for vector that von Neumann entropy of truncated measured density matrix at each step
double *S_y_pi;
//Pointer for vector that quantum discord at each step
double *quantum_discord;
//Pointer for vector storing Shannon entropy of coin state at each step
double *H_c;
//Pointer for vector storing Shannon entropy of first particle's position at each step
double *H_x;
//Pointer for vector storing Shannon entropy of second particle's position at each step
double *H_y;
//Pointer for vector storing Shannon entropy of particles positions at each step
double *H_xy;
//Pointer for vector storing Shannon mutual information of particles positions at each step
double *I_xy;
//Pointer for vector storing mean value of first particle's position at each step
double *mean_x;
//Pointer for vector storing mean value of second particle's position at each step
double *mean_y;
//x and y are pointers to vectors with values from -N to N. They are used in plot3D
double *x, *y;
//Pointer for vector storing covariance of particles positions at each step
double *xy_cov;
//Pointer for vector storing mean distance between particles positions at each step
double *mean_dist;
//Fraction of coin impurities for first particle
double index_coin_imp_1;
//Fraction of coin impurities for second particle
double index_coin_imp_2;
//Fraction of link impurities of first particle
double index_brokenlink_imp_1;
//Fraction of link impurities of second particle
double index_brokenlink_imp_2;
//Pointer for vector storing upper bound of entanglement of formation at each step 
double *E_f;
int *ib, *jb;
//Number of measuring points
int dim;
//Pointers for vectors storing measuring coordinates
int *i_0, *j_0;

c_matrix randomt_1, randomt_2;

int varying_coin_each_step_1, varying_coin_each_step_2;

//Stores the coin operator to be applied to first particle at each position and each step
c_matrix A_C;
//Stores the coin operator to be applied to second particle at each position and each step
c_matrix B_C;
//Coin state of the system (trace over position)
c_matrix rho_c; 
//Matrix used to perform coin operation on the overall state
c_matrix random_coin;
//Density matrix of first particle's position
c_matrix X_p;
//Density matrix of second particle's position
c_matrix Y_p;
//2 by 2 matrices used to generate other 4 by 4 matrices on Grover2D, Hadamard2D and in kron_A_B
c_matrix AA,BB;
//Truncated measured density matrix
c_matrix Y_pos;
//Stores eigenvalues of Y_pos
c_matrix Y_lambda;
// Stores all the coin operations to be applied at every position at each step for the first particle
c_matrix *coins_w1;
// Stores all the coin operations to be applied at every position at each step for the second particle
c_matrix *coins_w2;
//Pointer to a matrix that stores the topology of the lattice
int **shift_cond;
//Pointer for matrix that stores topology of the lattice and per step alterations
int **shift_cond_altered;
//Stores probability of particle in certain positions, determined by i_0 and j_0, at each step.
double *measure_prob;
//Normalization factor for state after measurement
double *norm_factor;
//Pointer for vector that stores one shot probability to hit at each step
double *one_shot_prob_hit;
//Pointer for vector that stores the average hitting time at each step
double *av_hitting_time;
//Pointer for vector that stores the first time probability to hit at each step
double *first_time_hit;
//Pointer for vector that stores the concurrent probability to hit at each step
double *concurrent_hit_time;


/*Global Variables initialization functions*/
void initialize_global_variables();


/*Complex number operations*/
complex c_sum(complex, complex);
complex c_mult(complex, complex);
complex c_div(complex, complex);
complex c_conj(complex);
double c_abs_sq(complex);
double c_angle(complex);


/*complex matrix operations*/
double random_0_1();
double vector_norm_sq(c_matrix);
void new_matrix(c_matrix * ,int ,int );
void print_matrix(c_matrix );
void erase_matrix(c_matrix *);
void matrix_mult( c_matrix, c_matrix, c_matrix);
void matrix_sum(c_matrix ,c_matrix ,c_matrix );
void c_times_matrix(c_matrix, complex , c_matrix);
void matrix_conj(c_matrix ,c_matrix );
void matrix_transpose(c_matrix, c_matrix);
void matrix_adjoint(c_matrix *, c_matrix );
void identity(c_matrix);
void copy_matrix(c_matrix, c_matrix);
c_matrix copy_sub_matrix(c_matrix ,int *, int *);
void random_matrix_1D(c_matrix );
void random_matrix_2D();
void kron_A_B(c_matrix ,c_matrix);
c_matrix random_rotation_2D();
void Grover1D(c_matrix);
void kron(c_matrix, c_matrix, c_matrix);
void reflect_coin_1D(c_matrix );
void reflect_coin_2D(c_matrix );
void random_matrix_2D_same_line();
void random_matrix_1D_21(int, c_matrix);
void random_matrix_1D_22(int, c_matrix);
void random_matrix_1D_31(int, c_matrix);
void random_matrix_1D_32(int, c_matrix);

/*Eigenvalues*/
void Householder_matrix(c_matrix, c_matrix);
c_matrix matrix_P(c_matrix , c_matrix );
c_matrix I_diag_B(c_matrix , int );
c_matrix QR_decomposition(c_matrix);
c_matrix QR_algorithm(c_matrix);
void QR_algorithm2(c_matrix, c_matrix, c_matrix);
void Householder_matrix2(c_matrix, c_matrix);
int test_eigen_values(c_matrix, c_matrix);

/****************************************Quantum Walk***************************************/
int biject(int, int);
int biject3D(int, int, int);
void bijinv(int* ,int);
void bijinv3D(int*, int);
state * initialize_state();
void shift(state *,state*, int );
void shift2D(state *, state *, int );
void step(state *, state *, int );
void erase_state(state *);
void coin1D(state *, state *, c_matrix , int);
void coin2D(state *, state *, int);
void copy_state_normalized(state *,state *, int);
c_matrix position_trace_1D(state *);
void position_trace_2D(state *, int );
void coin_trace_2D(state *,int );
void Hadamard1D();
void Hadamard2D();
void Fourier2D();
void Grover2D();


/*Broken links, and graph related functions*/
void initialize_shift_cond();
void initialize_shift_cond_altered();
void reflect_i();
void reflect_j();
void circ();
void moebius2D_i();
void moebius2D_j();
void broken_link2D(int , int, int , int);
void copy_shift_cond();
void link_impurities();
void initialize_position_coin();
void initialize_position_coin_altered();
void coin_impurities();
void permanent_broken_link2D(int , int, int , int);
void permanent_broken_link_1st_walker(int  , int , int );
void permanent_broken_link_2nd_walker(int  , int , int );

/*Relevant Quantities*/
void mean_value_xy(state *);
void covariance_xy(state *);
complex entropy(c_matrix);
void rho_x();
void rho_y();
void entropy_xy();
void entropy_x();
void entropy_y();
void mutual_info_Shannon();
void mutual_info_von_neumann();
void sum_prop_pos(state *);
double prob_test(state *);
void P_av_actualization();

/*Quantum Walk Function*/
void qw_on_plane(state *s1, state *s2);


/*Measurement subroutines*/
void initialize_measure_ij();
void measure_options(int , int , int );
void entropy_truncated_measurements(int );
void quantum_discord_calc(int );


/*subroutines related to file generation and storage*/
void plot2D(char *,double *, double *, int);
void plot3D(char *, double *, int);
void plots();


/*subroutines related to GNUplot*/
void gnu_plot2D(char *, double *, double *, int);
void gnu_plot3D(char *, double *, int);
void gnu_plots();

#endif
