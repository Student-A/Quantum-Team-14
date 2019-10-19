#include <stdio.h>
#include <stdlib.h>
#include "qwsim.h"
#include <math.h>
#include <time.h>


void parse_aux(FILE *fp,double * a){
/*Processes the information regarding a field of the parse file*/
	int ch;
	double test;
	int k=1;
	
	ch = fgetc(fp);
	while( ch != EOF && ch != ' ' && ch != '\n'){
			ch = fgetc(fp);		
	}

	k = 0;
	while( ch != EOF && ch != '\n'){
		fscanf(fp,"%lf",&a[k]);
		k++;
		ch = fgetc(fp);
	}

	while(ch == ' '){
		fscanf(fp,"%lf",&a[k]);

		ch = fgetc(fp);
	}
	ch = fgetc(fp);
	while( ch != EOF && ch != ' ' && ch != '\n'){
			ch = fgetc(fp);		
	}

}

void parser(){
/*Extracts all the information from the parse file*/
	FILE *fp;
	int points;
	double test, pp, *aux, normalizer = 0.0;
	int *i, *j, *c,l;
	complex *v;
	double *aux_i,*aux_j;
	double *aux_di, *aux_dj;

	fp = fopen("parse_file.txt","r");
	
	parse_aux(fp,&pp);
	N = (int )pp;
	printf(" M = %d \n", N);

	Edges = 4;

	parse_aux(fp,&pp);
	steps = (int )pp;printf(" steps = %d \n", steps);
	steps = steps +1;

	parse_aux(fp,&pp);
	points = (int )pp;printf(" points = %d \n", points);

	aux = (double *)malloc(sizeof(double) * points);
	i = (int *)malloc(sizeof(int) * points);
	j = (int *)malloc(sizeof(int) * points);
	c = (int *)malloc(sizeof(int) * points);
	v = (complex *)malloc(sizeof(complex) * points);

	parse_aux(fp,aux); for(l=0; l < points; l++){ i[l] =(int ) aux[l]; printf("i[%d] = %d\n", l, i[l]);}
	parse_aux(fp,aux); for(l=0; l < points; l++){ j[l] =(int ) aux[l]; printf("j[%d] = %d\n", l, j[l]);}
	parse_aux(fp,aux); for(l=0; l < points; l++){ c[l] =(int ) aux[l]; printf("c[%d] = %d\n", l, c[l]);}
	parse_aux(fp,aux); for(l=0; l < points; l++){ v[l].re = aux[l]; normalizer = normalizer + aux[l]*aux[l]; printf("v.re[%d] = %lf\n", l, v[l].re);}
	parse_aux(fp,aux); for(l=0; l < points; l++){ v[l].im = aux[l]; normalizer = normalizer + aux[l]*aux[l]; printf("v.im[%d] = %lf\n", l, v[l].im);}
	free(aux);

	parse_aux(fp,&pp);
	different_line = (int ) pp;printf("different_line = %d \n", different_line);
	
	initialize_global_variables();
   	circ();
	parse_aux(fp,&pp);printf("index_brokenlink_imp_1 = %lf \n", pp);
	index_brokenlink_imp_1 = pp;

	int ai=0;
	int aj=0;
	int ad=0;


	double *aux_2;
	parse_aux(fp,&pp);printf("number of fixed broken links 1st particle = %lf \n", pp);
	fixed_broken_link_i = (int ) pp;
	aux_i = (double *)malloc(sizeof(double) * fixed_broken_link_i);
   ib = (int *)malloc(sizeof(int) * 2*N+1);
   jb = (int *)malloc(sizeof(int) * 2*N+1);
    
    
    for(l=0; l <= 2*N; l++){
        ib[l] = 0;
        jb[l] = 0;
    }
    
	parse_aux(fp,aux_i);
	for(l=0; l < fixed_broken_link_i; l++)printf("fixed broken link i = %lf\n",aux_i[l]);
	for(l=0; l < fixed_broken_link_i; l++) ib[(int) aux_i[l] + N] = 1; //permanent_broken_link_1st_walker((int ) aux_i[l] , aj , ad );
	free(aux_i);


	parse_aux(fp,&pp);
	coin1 = (int) pp;printf("coin1 = %lf \n", pp);

		
	for(l=0; l< Nodes1D; l++){
		pos_coin[l][0] = coin1;
		pos_coin_altered[l][0] = coin1;		
	}

	parse_aux(fp,&pp);printf("varying_coin_each_step_1 = %lf \n", pp);
	varying_coin_each_step_1 =(int ) pp;
	
	
	if(varying_coin_each_step_1 == 1){
		
	//	for(l=0; l < Nodes; l++) free(range1[l]);
	//	free(range1);
				
		for(l=0; l < Nodes1D; l++) pos_coin[l][0] = 5;


		range1_fixed = (double **)malloc(sizeof(double * )*1);
		range1_fixed[0] = (double *)malloc(sizeof(double ) * 6);
		
		parse_aux(fp,&pp);
		range1_fixed[0][0] = pp;

		parse_aux(fp,&pp);
		range1_fixed[0][1] = pp;

		parse_aux(fp,&pp);
		range1_fixed[0][2] = pp;

		parse_aux(fp,&pp);
		range1_fixed[0][3] = pp;

		parse_aux(fp,&pp);
		range1_fixed[0][4] = pp;

		parse_aux(fp,&pp);
		range1_fixed[0][5] = pp;

	}
	else{
		parse_aux(fp,&pp);
		parse_aux(fp,&pp);
		parse_aux(fp,&pp);
		parse_aux(fp,&pp);
		parse_aux(fp,&pp);
		parse_aux(fp,&pp);

	}
	

	



	parse_aux(fp,&pp);printf("index_coin_imp_1 = %lf \n", pp);
	index_coin_imp_1 =  pp;

	parse_aux(fp,&pp);printf("number of fixed coins for first walker = %lf \n", pp);
	fixed_coin_1 =(int ) pp;

	range1 = (double **)malloc(sizeof(double * ) * ((2*N+1)));
	for(l=0; l < (2*N+1); l++) range1[l] = (double *)malloc(sizeof(double ) * 6);

	aux_i = (double *)malloc(sizeof(double) * fixed_coin_1);
	parse_aux(fp,aux_i);

	aux_2 = (double *)malloc(sizeof(double) * fixed_coin_1);

	for(l=0; l < fixed_coin_1; l++) {pos_coin[((int ) aux_i[l])+N][0] = 4;
	
}

	parse_aux(fp,aux_2);
	
	for(l=0; l < fixed_coin_1; l++)range1[((int ) aux_i[l])+N][0] = aux_2[l];
	
	parse_aux(fp,aux_2);
	for(l=0; l < fixed_coin_1; l++)range1[((int ) aux_i[l])+N][1] = aux_2[l];
	
	parse_aux(fp,aux_2);
	for(l=0; l < fixed_coin_1; l++)range1[((int ) aux_i[l])+N][2] = aux_2[l];
	
	parse_aux(fp,aux_2);
	for(l=0; l < fixed_coin_1; l++)range1[((int ) aux_i[l])+N][3] = aux_2[l];
	
	parse_aux(fp,aux_2);
	for(l=0; l < fixed_coin_1; l++)range1[((int ) aux_i[l])+N][4] = aux_2[l];
	
	parse_aux(fp,aux_2);
	for(l=0; l < fixed_coin_1; l++)range1[((int ) aux_i[l])+N][5] = aux_2[l];
	
	free(aux_2);
	free(aux_i);

	

	
	parse_aux(fp,&pp);
	rf1 = (int) pp;


	parse_aux(fp,&pp);printf("index_brokenlink_imp_2 = %lf \n", pp);
	index_brokenlink_imp_2 = pp;


	parse_aux(fp,&pp);printf("number of fixed broken links 2nd particle = %lf \n", pp);
	fixed_broken_link_j = (int ) pp;
	aux_i = (double *)malloc(sizeof(double) * fixed_broken_link_j);
	parse_aux(fp,aux_i);
	for(l=0; l < fixed_broken_link_j; l++)printf("fixed broken link j = %lf\n",aux_i[l]);
	for(l=0; l < fixed_broken_link_j; l++) jb[(int) aux_i[l]+N] = 1;
	free(aux_i);


	parse_aux(fp,&pp);
	coin2 = (int) pp;printf("coin2 = %lf \n", pp);

	
	for(l=0; l< Nodes1D; l++){
		pos_coin[l][1] = coin2;
		pos_coin_altered[l][1] = coin2;
	}

	parse_aux(fp,&pp);printf("varying_coin_each_step_2 = %lf \n", pp);
	varying_coin_each_step_2 =(int ) pp;
	printf("What's up 1?\n");

	if(varying_coin_each_step_2 == 1){
		
	//	for(l=0; l < Nodes; l++) free(range1[l]);
	//	free(range1);
		printf("What's up 2?\n");		
		for(l=0; l < Nodes1D; l++) pos_coin[l][1] = 5;

		printf("pos_coin[0] = %d\n", pos_coin[0][1]);

		range2_fixed = (double **)malloc(sizeof(double * )*1);
		range2_fixed[0] = (double *)malloc(sizeof(double ) * 6);
		
		parse_aux(fp,&pp);
		range2_fixed[0][0] = pp; printf("range2_fixed[0][0] = %lf \n",range2_fixed[0][0]);

		parse_aux(fp,&pp);
		range2_fixed[0][1] = pp; printf("range2_fixed[0][1] = %lf \n",range2_fixed[0][1]);

		parse_aux(fp,&pp);
		range2_fixed[0][2] = pp;printf("range2_fixed[0][2] = %lf \n",range2_fixed[0][2]);

		parse_aux(fp,&pp);
		range2_fixed[0][3] = pp;printf("range2_fixed[0][3] = %lf \n",range2_fixed[0][3]);

		parse_aux(fp,&pp);
		range2_fixed[0][4] = pp;printf("range2_fixed[0][4] = %lf \n",range2_fixed[0][4]);

		parse_aux(fp,&pp);
		range2_fixed[0][5] = pp;printf("range2_fixed[0][5] = %lf \n",range2_fixed[0][5]);

	}
	else{
		parse_aux(fp,&pp);
		parse_aux(fp,&pp);
		parse_aux(fp,&pp);
		parse_aux(fp,&pp);
		parse_aux(fp,&pp);
		parse_aux(fp,&pp);

	}
	



	parse_aux(fp,&pp);printf("index_coin_imp_2 = %lf \n", pp);
	index_coin_imp_1 =  pp;

	parse_aux(fp,&pp);printf("number of fixed coins for second walker = %lf \n", pp);
	fixed_coin_2 =(int ) pp;

	range2 = (double **)malloc(sizeof(double * ) * ((2*N+1)));
	for(l=0; l < (2*N+1); l++) range2[l] = (double *)malloc(sizeof(double ) * 6);

	aux_i = (double *)malloc(sizeof(double) * fixed_coin_2);
	parse_aux(fp,aux_i);

	aux_2 = (double *)malloc(sizeof(double) * fixed_coin_2);

	for(l=0; l < fixed_coin_2; l++) {pos_coin[((int ) aux_i[l])+N][1] = 4;
}

	parse_aux(fp,aux_2);
	
	for(l=0; l < fixed_coin_2; l++)range2[((int ) aux_i[l])+N][0] = aux_2[l];
	
	parse_aux(fp,aux_2);
	for(l=0; l < fixed_coin_2; l++)range2[((int ) aux_i[l])+N][1] = aux_2[l];
	
	parse_aux(fp,aux_2);
	for(l=0; l < fixed_coin_2; l++)range2[((int ) aux_i[l])+N][2] = aux_2[l];
	
	parse_aux(fp,aux_2);
	for(l=0; l < fixed_coin_2; l++)range2[((int ) aux_i[l])+N][3] = aux_2[l];
	
	parse_aux(fp,aux_2);
	for(l=0; l < fixed_coin_2; l++)range2[((int ) aux_i[l])+N][4] = aux_2[l];
	
	parse_aux(fp,aux_2);
	for(l=0; l < fixed_coin_2; l++)range2[((int ) aux_i[l])+N][5] = aux_2[l];
	
	free(aux_2);
	free(aux_i);

	
	parse_aux(fp,&pp);
	rf2 = (int) pp;
		
	normalizer = sqrt(normalizer);
	for(l=0; l < points; l++){
		printf("normalizer = %lf\n", normalizer);
		M(s1[biject(i[l],j[l])].direction,c[l],0).re = v[l].re/normalizer;
		M(s1[biject(i[l],j[l])].direction,c[l],0).im = v[l].im/normalizer;
		printf("prob test = %lf\n", prob_test(s1));
	}
	
	parse_aux(fp,&pp);
	dim = (int) pp; printf("dim = %d \n", dim);
	
	initialize_measure_ij();
	aux = (double *)malloc(sizeof(double) * dim);
	parse_aux(fp,aux); for(l=0; l < dim; l++){ i_0[l] =(int ) aux[l]; printf("i_0[%d] = %d\n", l, i_0[l]);}
	parse_aux(fp,aux); for(l=0; l < dim; l++){ j_0[l] =(int ) aux[l]; printf("j_0[%d] = %d\n", l, j_0[l]);}
	free(aux);

	parse_aux(fp,&pp);
	m_options = (int) pp;
    printf("m_options = %d\n", m_options);


	parse_aux(fp,&pp);
	vn_x = (int) pp;printf("vn_x = %lf \n", pp);

	parse_aux(fp,&pp);
	vn_y = (int) pp;printf("vn_y = %lf \n", pp);

	parse_aux(fp,&pp);
	vn_m_info = (int) pp;printf("vn_m_info = %lf \n", pp);

	parse_aux(fp,&pp);
	q_discord = (int) pp;printf("q_discord = %lf \n", pp);

	parse_aux(fp,&pp);
	ent_form = (int) pp;printf("ent_form = %lf \n", pp);

	

	
	if(rf1 == 1)reflect_i();
	if(rf2 == 1)reflect_j();

	parse_aux(fp,&pp);
	gnu = (int) pp;
	
	fclose(fp);


}


int main(void){
	double p, *pp, xy;
	N_malloc = 0;
	N_free = 0;
	
	// setting random number generator seed 
	time_t t;
	t = time( NULL );
	srand(t);

	t_counter = 1;
	parser();

	qw_on_plane(s1, s2);

	p = prob_test(s1);


	printf("\n %f \n",p);
	
	plots();
	if(gnu == 1) gnu_plots();

	return 0;
}

