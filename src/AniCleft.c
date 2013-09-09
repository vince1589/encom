#include "STeM.h"

int main(int argc, char *argv[]) {
	int all; /*Nombre d'atomes dans pdb*/

	int all_t; /*Nombre d'atomes dans pdb*/

 	int help_flag = 1;
 	char file_name[500];
 	char check_name[500];
 	char out_name[500];
 	int verbose = 0;
	int i;
	long seed;
	seed = time_seed();
	int nconn;
	int print_flag = 0;
	int it = 1000;
	int dens_flag = 1;
	char wholeinit[500] = "UNDEF";
	char wholetarg[500] = "UNDEF";
	
 	for (i = 1;i < argc;i++) {
 		if (strcmp("-i",argv[i]) == 0) {strcpy(file_name,argv[i+1]);--help_flag;}
 		if (strcmp("-h",argv[i]) == 0) {help_flag = 1;}
 		if (strcmp("-dist",argv[i]) == 0) {dens_flag = 0;}
 		if (strcmp("-v",argv[i]) == 0) {verbose = 1;}
		if (strcmp("-pi",argv[i]) == 0) {strcpy(wholeinit,argv[i+1]);}
		if (strcmp("-pt",argv[i]) == 0) {strcpy(wholetarg,argv[i+1]);}
 		if (strcmp("-t",argv[i]) == 0) {strcpy(check_name,argv[i+1]);help_flag = 0;}
 		if (strcmp("-o",argv[i]) == 0) {strcpy(out_name,argv[i+1]);++print_flag;}
 		if (strcmp("-it",argv[i]) == 0) {int temp;sscanf(argv[i+1],"%d",&temp);it = temp;}
 	
 	}
	 	
 	if (help_flag == 1) {
 		printf("****************************\nHelp Section\n-i\tFile Input (PDB)\n-v\tVerbose\n-w\tWeight Vector\n-t\tInitial Value of template (negative value for random)\n\tIf Load Template, multiply the template\n-lt\tLoad input template\n-sp\tSuper Node Mode (CA, N, C)\n-kt\tPoid de l'angle entre les nodes (1)\n-kr\tPoid de la distance entre les nodes (1)\n-f\tFile to fit\n****************************\n");
 		return(0); 
 	} 

 	//***************************************************
 	//*													*
 	//*Build a structure contaning information on the pdb
 	//*													*
 	//***************************************************
 	
 	// Première strucutre
 	
 	all = count_atom(file_name);
 	nconn = count_connect(file_name);
 	
 	if (verbose == 1) {printf("Connect:%d\n",nconn);}
 	
	if (verbose == 1) {printf("Assigning Structure\n\tAll Atom\n");}
	
	// Array qui comprend tous les connects
	
	int **connect_h=(int **)malloc(nconn*sizeof(int *)); 
    for(i=0;i<nconn;i++) { connect_h[i]=(int *)malloc(6*sizeof(int));}
    
    assign_connect(file_name,connect_h);

	// Assign tous les atoms
	
	struct pdb_atom strc_all[all];
	build_all_strc(file_name,strc_all); // Retourne le nombre de Node
	
	if (verbose == 1) {printf("	Atom:%d\n",all);}
	check_lig(strc_all,connect_h,nconn,all);
	
	//Construit la structure a comparer
 	nconn = 0;
 	all_t = count_atom(check_name);
 	nconn = count_connect(check_name);
 	
 	if (verbose == 1) {printf("Connect:%d\n",nconn);}
 	
	if (verbose == 1) {printf("Assigning Structure\n\tAll Atom\n");}
	
	// Array qui comprend tous les connects
	
	int **connect_t=(int **)malloc(nconn*sizeof(int *)); 
    for(i=0;i<nconn;i++) { connect_t[i]=(int *)malloc(6*sizeof(int));}

  assign_connect(check_name,connect_t);
	
	// Assign tous les atoms
	
	struct pdb_atom strc_all_t[all_t];
	build_all_strc(check_name,strc_all_t); // Retourne le nombre de Node
	

 	
 	assign_atom_type(strc_all,all);
 	assign_atom_type(strc_all_t,all_t);
 	
 	printf("I found %d Anisou with %d atom !\n",load_anisou(strc_all,file_name,all),all);
  printf("I found %d Anisou with %d atom !\n",load_anisou(strc_all_t,check_name,all_t),all_t);
	

	
	
	center_strc(all_t,strc_all_t);
	center_strc(all,strc_all);
	
	struct pdb_atom tstrc[all_t];
	
	
	// Rotation axe !
	gsl_vector *xaxis = gsl_vector_alloc(3);
	gsl_vector_set_all(xaxis,0);
	gsl_vector_set(xaxis,0,1);
	
	gsl_vector *yaxis = gsl_vector_alloc(3);
	gsl_vector_set_all(yaxis,0);
	gsl_vector_set(yaxis,1,1);
	
	gsl_vector *zaxis = gsl_vector_alloc(3);
	gsl_vector_set_all(zaxis,0);
	gsl_vector_set(zaxis,2,1);
	
	// Find the simlar atom type
	
	int two[8][all];
	int twom[8];
	for(i=0;i<8;++i){twom[i] = 0;}
	for(i=0;i<all;++i) {
			
			if (strc_all[i].type > 8) {printf("Cannot assign type\n");continue;}
			two[strc_all[i].type-1][twom[strc_all[i].type-1]] = i;
			++twom[strc_all[i].type-1];
	}

	gsl_matrix *rota = gsl_matrix_alloc(3,3);
	gsl_vector *trans = gsl_vector_alloc(3);
	// Optimize the best super
	
	float range[6];
	for(i=0;i<3;++i) {range[i] = 3.1416*2;range[i+3] = 10;}
	
	
	float gene[6];
	float bgen[6];
	float qgen[6];
	int j;
	float old = -1000;
	float best = -1000;
	double conj_dens12;
	gsl_vector *delr = gsl_vector_alloc(3);
	gsl_matrix *incov12 = gsl_matrix_alloc(3,3);
	gsl_vector *pos = gsl_vector_alloc(3);
		gsl_vector_set_all(pos,0);
	for(i=0;i<6;++i) {gene[i] = 0;bgen[i] = 0;qgen[i] = 0;}
	
	
	// The Monte Carlo
	printf("Starting MC !\n");
	for(i =0;i<it;++i) {
		copy_strc(tstrc,strc_all_t,all_t);
		
		// Apply rotation
		int change = 0;
		if (i != 0) {
			while(change == 0) {
				for(j=0;j<6;++j) {
					gene[j] = qgen[j];
					if(ran2(&seed) > 0.0) {
						gene[j] = ran2(&seed)*range[j];
						if(j>2) {gene[j]-=5;}
						++change;
					}
				}
			}
		}
		rotate_matrix(rota,gene[0],xaxis); // In radian !
		rotate_all(rota,tstrc,all_t);
		rotate_matrix(rota,gene[1],yaxis); // In radian !
		rotate_all(rota,tstrc,all_t);
		rotate_matrix(rota,gene[2],zaxis); // In radian !
		rotate_all(rota,tstrc,all_t);
		
		for(j=0;j<3;++j) {gsl_vector_set(trans,j,gene[j+3]);}
		translate_strc(tstrc,all_t,trans);
	
		// Compare minimum distance
		
		// Comapre prob density
		float dens_sum = 0;
		float avg = 0;
		for(j=0;j<all_t;++j) {
			int type = tstrc[j].type-1;
			if (type < 0 || type > 7) {printf("I didn't have a type:\n");continue;}
			int k;
			float max = 20*20;
			for(k=0;k<twom[type];++k) {
				
				float dist = pow(tstrc[j].x_cord-strc_all[two[type][k]].x_cord,2)+pow(tstrc[j].y_cord-strc_all[two[type][k]].y_cord,2)+pow(tstrc[j].z_cord-strc_all[two[type][k]].z_cord,2);
			
				if (dist < max) {max = dist;}

				if (dist < 10*10) {
					
					if (conj_prob_init(&tstrc[j], &strc_all[two[type][k]], incov12,delr,&conj_dens12) != -1 && dens_flag == 1) {
						float dens = density_prob(incov12, delr, conj_dens12,pos);
					//	float prob =proxim_prob(incov12, delr, conj_dens12, 0.00, 4, 50);
					//	printf("Dist:%f Dens:%.12f Prob:%.12f\n",sqrt(dist),dens,prob);
						dens_sum += dens;
					}
				}
			
			}
			avg += (max);
		//	printf("Type:%d J:%d Dist:%f Dens:%f\n",type,j,sqrt(max),dens_sum);
			
		}
		

		
		if (i % 1000 == 0) {printf("IT:%d Best:%.6f Dens sum:%.6f Dist:%.6f Gene:",i,best,dens_sum,avg/all);
		for(j=0;j<6;++j) {printf("%.6f ",gene[j]);}printf("\n");
		}
		if (dens_flag == 1) {
			if (dens_sum > old) {
				for(j=0;j<6;++j) {qgen[j] = gene[j];}
				old = dens_sum;		
			}
			if (dens_sum > best) {
				for(j=0;j<6;++j) {bgen[j] = gene[j];}
				best = dens_sum;		
			}
		} else {
			if (-avg/all > old) {
				for(j=0;j<6;++j) {qgen[j] = gene[j];}
				old = -avg/all;		
			}
			if (-avg/all > best) {
				for(j=0;j<6;++j) {bgen[j] = gene[j];}
				best = -avg/all;		
			}
		}
	}
	
	// Apply rotation
	
	
		
	rotate_matrix(rota,bgen[0],xaxis); // In radian !
	rotate_all(rota,tstrc,all_t);
	rotate_matrix(rota,bgen[1],yaxis); // In radian !
	rotate_all(rota,tstrc,all_t);
	rotate_matrix(rota,bgen[2],zaxis); // In radian !
	rotate_all(rota,strc_all_t,all_t);
	
	for(j=0;j<3;++j) {gsl_vector_set(trans,j,bgen[j+3]);}
	translate_strc(strc_all_t,all_t,trans);
	
	printf("My best = %f\n",best);
	for(j=0;j<6;++j) {printf("%.5f ",bgen[j]);}printf("\n");
	write_strc("match.pdb",strc_all_t,all_t,1.00);
	write_strc("init.pdb",strc_all,all,1.00);
	
	// Load the whole structure
	
	if (strcmp(wholeinit,"UNDEF") != 0) {
		printf("WholeInit:%s\n",wholeinit);	// Première strucutre
 		float tocenter[3];
 		build_all_strc(file_name,strc_all);
 		for(i=0;i<3;++i) {tocenter[i] = 0.0;}
 		for(i=0;i<all;++i) {
 			tocenter[0] -= strc_all[i].x_cord/all;
 			tocenter[1] -= strc_all[i].y_cord/all;
 			tocenter[2] -= strc_all[i].z_cord/all;
 		}
 		for(j=0;j<3;++j) {gsl_vector_set(trans,j,tocenter[j]);}
 		
	 	all = count_atom(wholeinit);
	 	nconn = count_connect(wholeinit);
	 	
	 	if (verbose == 1) {printf("Connect:%d\n",nconn);}
	 	
		if (verbose == 1) {printf("Assigning Structure\n\tAll Atom\n");}
	
		// Array qui comprend tous les connects
	
		int **connect_bla=(int **)malloc(nconn*sizeof(int *)); 
	  for(i=0;i<nconn;i++) { connect_bla[i]=(int *)malloc(6*sizeof(int));}
	  
	  assign_connect(wholeinit,connect_bla);

		// Assign tous les atoms
	
		struct pdb_atom strc_init[all];
		build_all_strc(wholeinit,strc_init); // Retourne le nombre de Node
	
		if (verbose == 1) {printf("	Atom:%d\n",all);}
	//	check_lig(strc_init,connect_h,nconn,all);
		
		// Centre le BS
		print_vector(trans);
		translate_strc(strc_init,all,trans);
		write_strc("init_center.pdb",strc_init,all,1.00);
		
	}
	
	// Load the whole structure
	
	if (strcmp(wholetarg,"UNDEF") != 0) {
		printf("wholetarg:%s\n",wholetarg);	// Première strucutre
 		float tocenter[3];
 		build_all_strc(check_name,strc_all_t);
 		for(i=0;i<3;++i) {tocenter[i] = 0.0;}
 		for(i=0;i<all_t;++i) {
 			tocenter[0] -= strc_all_t[i].x_cord/all_t;
 			tocenter[1] -= strc_all_t[i].y_cord/all_t;
 			tocenter[2] -= strc_all_t[i].z_cord/all_t;
 		}
 		for(j=0;j<3;++j) {gsl_vector_set(trans,j,tocenter[j]);}
 		
	 	all = count_atom(wholetarg);
	 	nconn = count_connect(wholetarg);
	 	
	 	if (verbose == 1) {printf("Connect:%d\n",nconn);}
	 	
		if (verbose == 1) {printf("Assigning Structure\n\tAll Atom\n");}
	
		// Array qui comprend tous les connects
	
		int **connect_bla2=(int **)malloc(nconn*sizeof(int *)); 
	  for(i=0;i<nconn;i++) { connect_bla2[i]=(int *)malloc(6*sizeof(int));}
	  
	  assign_connect(wholetarg,connect_bla2);

		// Assign tous les atoms
	
		struct pdb_atom strc_init[all];
		build_all_strc(wholetarg,strc_init); // Retourne le nombre de Node
	
		if (verbose == 1) {printf("	Atom:%d\n",all);}
	//	check_lig(strc_init,connect_h,nconn,all);
		
		// Centre le BS
		print_vector(trans);
		translate_strc(strc_init,all,trans);
		
		rotate_matrix(rota,bgen[0],xaxis); // In radian !
		rotate_all(rota,strc_init,all);
		rotate_matrix(rota,bgen[1],yaxis); // In radian !
		rotate_all(rota,strc_init,all);
		rotate_matrix(rota,bgen[2],zaxis); // In radian !
		rotate_all(rota,strc_init,all);
		
		write_strc("target_center.pdb",strc_init,all,1.00);
		
	}
	
	
	
	
	
	
	
	
	gsl_matrix_free(rota);
}
