#include "STeM.h"
float write_strc_type(char filename[100], struct pdb_atom *newstrc,int nb_atom,float factor,int type);
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
	float temperature = 100;
	char liginit[500] = "UNDEF";
	char ligtarg[500] = "UNDEF";
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
 		if (strcmp("-temp",argv[i]) == 0) {float temp;sscanf(argv[i+1],"%f",&temp);temperature = temp;}
 		
 		if (strcmp("-li",argv[i]) == 0) {strcpy(liginit,argv[i+1]);}
		if (strcmp("-lt",argv[i]) == 0) {strcpy(ligtarg,argv[i+1]);}
 	
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
	
	// Load ligand Iinit
	int liginit_atom = count_atom(liginit);
	struct pdb_atom liginit_strc[liginit_atom];
	if (liginit_atom != 0) {
		build_all_strc(liginit,liginit_strc);
 	}
 	
 	// Load ligand Iinit
	int ligtarg_atom = count_atom(ligtarg);
	struct pdb_atom ligtarg_strc[ligtarg_atom];
	struct pdb_atom ligtarg_strct[ligtarg_atom];
	if (ligtarg_atom != 0) {
		build_all_strc(ligtarg,ligtarg_strc);
 	}
 	printf("Ligand atom:%d et %d\n",liginit_atom,ligtarg_atom);
 	assign_atom_type(strc_all,all);
 	assign_atom_type(strc_all_t,all_t);
 	
 	printf("I found %d Anisou with %d atom !\n",load_anisou(strc_all,file_name,all),all);
  printf("I found %d Anisou with %d atom !\n",load_anisou(strc_all_t,check_name,all_t),all_t);
	

	gsl_matrix *rota = gsl_matrix_alloc(3,3);
	gsl_vector *trans = gsl_vector_alloc(3);
	
	center_strc(all_t,strc_all_t,trans);
	if (ligtarg_atom != 0) {
		for(i = 0;i<ligtarg_atom;++i) {
			ligtarg_strc[i].x_cord -= gsl_vector_get(trans,0);
			ligtarg_strc[i].y_cord -= gsl_vector_get(trans,1);
			ligtarg_strc[i].z_cord -= gsl_vector_get(trans,2);
		}
	}
	center_strc(all,strc_all,trans);
	if (liginit_atom != 0) {
		for(i = 0;i<liginit_atom;++i) {
			liginit_strc[i].x_cord -= gsl_vector_get(trans,0);
			liginit_strc[i].y_cord -= gsl_vector_get(trans,1);
			liginit_strc[i].z_cord -= gsl_vector_get(trans,2);
		}
	}
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

	
	// Optimize the best super
	
	float range[6];
	for(i=0;i<3;++i) {range[i] = 3.1416;range[i+3] = 20;}
	
	
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
	
		// Copy strc
	
		copy_strc(tstrc,strc_all_t,all_t);
		if (ligtarg_atom != 0) {
			copy_strc(ligtarg_strct,ligtarg_strc,ligtarg_atom);
		}
		// Apply rotation
		int change = 0;
		if (i != 0) {
			while(change == 0) {
				for(j=0;j<6;++j) {
					gene[j] = qgen[j];
					if(ran2(&seed) > 0.90) {
						gene[j] += ran2(&seed)-0.5;
						
						++change;
					}
					if(ran2(&seed) > 0.95) {
						gene[j] = ran2(&seed)*range[j];
						if(j>2) {gene[j]-=10;}
						++change;
					}
				}
			}
		}
		rotate_matrix(rota,gene[0],xaxis); // In radian !
		rotate_all(rota,tstrc,all_t);
		if (ligtarg_atom != 0) {rotate_all(rota,ligtarg_strct,ligtarg_atom);}
		rotate_matrix(rota,gene[1],yaxis); // In radian !
		rotate_all(rota,tstrc,all_t);
		if (ligtarg_atom != 0) {rotate_all(rota,ligtarg_strct,ligtarg_atom);}
		rotate_matrix(rota,gene[2],zaxis); // In radian !
		rotate_all(rota,tstrc,all_t);
		if (ligtarg_atom != 0) {rotate_all(rota,ligtarg_strct,ligtarg_atom);}
		for(j=0;j<3;++j) {gsl_vector_set(trans,j,gene[j+3]);}
		translate_strc(tstrc,all_t,trans);
		if (ligtarg_atom != 0) {translate_strc(ligtarg_strct,ligtarg_atom,trans);}
		// Compare minimum distance
		
		// Comapre prob density
		float dens_sum = 0;
		float avg = 0;
		for(j=0;j<all_t;++j) {
			int type = tstrc[j].type-1;
			if (type < 0 || type > 7) {printf("I didn't have a type:\n");continue;}
			int k;
			float max = 20*20;
			float max_dens = 0;
			float dens_avg;
			for(k=0;k<twom[type];++k) {
				
				float dist = pow(tstrc[j].x_cord-strc_all[two[type][k]].x_cord,2)+pow(tstrc[j].y_cord-strc_all[two[type][k]].y_cord,2)+pow(tstrc[j].z_cord-strc_all[two[type][k]].z_cord,2);
			
				if (dist < max) {max = dist;}

				if (dist < 10*10) {
					
					if (conj_prob_init(&tstrc[j], &strc_all[two[type][k]], incov12,delr,&conj_dens12) != -1 && dens_flag == 1) {
						float dens = density_prob(incov12, delr, conj_dens12,pos);
					//	float prob =proxim_prob(incov12, delr, conj_dens12, 0.00, 4, 50);
						//if (i == 0) {printf("ID:%d %d %f %.12f \n",j,two[type][k],sqrt(dist),dens);}
						dens_sum += dens;
						if (dens > max_dens) {max_dens = dens;}
					}
				}
			
			}
			avg += (max);
			dens_avg += (max_dens);
		//	printf("Type:%d J:%d Dist:%f Dens:%f\n",type,j,sqrt(max),dens_sum);
			
		}
		float ligrmsd = 0;
		if (ligtarg_atom != 0) {
			int o;
			for (o=0;o<ligtarg_atom;++o) {
				ligrmsd += pow(ligtarg_strct[o].x_cord-liginit_strc[o].x_cord,2)+pow(ligtarg_strct[o].y_cord-liginit_strc[o].y_cord,2)+pow(ligtarg_strct[o].z_cord-liginit_strc[o].z_cord,2);
			}
			ligrmsd = sqrt(ligrmsd/ligtarg_atom);
		}

		
		/*if (i % 1000 == 0) {printf("IT:%d Best:%.12f Dens sum:%.12f Dist:%.6f Ligrmsd:%.6f Gene:",i,best,dens_sum,avg/all,ligrmsd);
		for(j=0;j<6;++j) {printf("%.6f ",gene[j]);}printf("\n");
		}*/
		//temperature = 1.0/(1.0+i);
		if (dens_flag == 1) {
			float prob = pow(2.71,-(old-dens_sum)/temperature);
			if (dens_sum > old || prob > ran2(&seed)) {
				printf("IT:%d Old:%g New:%g Prob:%.4f Ligrmsd:%.2f Gene:%.2f %.2f %.2f %.2f %.2f %.2f\n",i,old,dens_sum,prob,ligrmsd,gene[0],gene[1],gene[2],gene[3],gene[4],gene[5]);
				for(j=0;j<6;++j) {qgen[j] = gene[j];}
				old = dens_sum;		
			} 
			if (dens_sum > best) {
				for(j=0;j<6;++j) {bgen[j] = gene[j];}
				best = dens_sum;		
			}
		} else {
			float prob = pow(2.71,-(old-avg)/temperature);
			if (-avg > old) {
				printf("IT:%d Old:%g New:%g Prob:%.4f Ligrmsd:%.2f\n",i,old,dens_sum,prob,ligrmsd);
				for(j=0;j<6;++j) {qgen[j] = gene[j];}
				old = -avg;		
			}
			if (-avg > best) {
				for(j=0;j<6;++j) {bgen[j] = gene[j];}
				best = -avg;		
			}
		}
	}
	// MC finished
	
	
	// Apply rotation
	

		
	rotate_matrix(rota,bgen[0],xaxis); // In radian !
	rotate_all(rota,strc_all_t,all_t);
	rotate_matrix(rota,bgen[1],yaxis); // In radian !
	rotate_all(rota,strc_all_t,all_t);
	rotate_matrix(rota,bgen[2],zaxis); // In radian !
	rotate_all(rota,strc_all_t,all_t);
	
	for(j=0;j<3;++j) {gsl_vector_set(trans,j,bgen[j+3]);}
	translate_strc(strc_all_t,all_t,trans);
	


	

	printf("My best = %.12f\n",best);
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

float write_strc_type(char filename[100], struct pdb_atom *newstrc,int nb_atom,float factor,int type){
 	FILE *out_file;
 	int k;
 	out_file = fopen(filename,"w");
 	float avg = 0;
 	float tot = 1;
 	for (k = 0; k<nb_atom;k++) {
 		if(newstrc[k].type != type) {continue;}
 		if (newstrc[k].atom_type == 1) {fprintf(out_file,"ATOM  ");}
	 	if (newstrc[k].atom_type == 2) {fprintf(out_file,"HETATM");}
	 	if (newstrc[k].atom_type == 3) {fprintf(out_file,"HETATM");}
 		fprintf(out_file,"%5.d %s%s %s%4.d%12.3f%8.3f%8.3f  1.00  %2.2f\n",
 			newstrc[k].atom_number,
 			newstrc[k].atom_prot_type,
 			newstrc[k].res_type,
 			newstrc[k].chain,
 			newstrc[k].res_number,
 			newstrc[k].x_cord,
 			newstrc[k].y_cord,
 			newstrc[k].z_cord,
 			newstrc[k].b_factor
 			
 			);
	 		if (newstrc[k].main_vars[0] > 0) {
	 			// Rebuild matrix
	 			
	 			gsl_matrix *temp = gsl_matrix_alloc(3,3);
	 			gsl_matrix *temp2 = gsl_matrix_alloc(3,3);
	 			gsl_matrix *evec = gsl_matrix_alloc(3,3);
	 			gsl_matrix *tevec = gsl_matrix_alloc(3,3);
	 			gsl_matrix *eval = gsl_matrix_alloc(3,3);
	 			gsl_matrix_set_all(temp,0);
	 			gsl_matrix_set_all(temp2,0);
	 			gsl_matrix_set_all(eval,0);
	 			
	 			
	 			int i,j;

				for(i=0;i<3;++i) {
	 				for(j=0;j<3;++j) {
	 					gsl_matrix_set(evec,i,j,newstrc[k].global_evecs[i][j]);
	 					gsl_matrix_set(tevec,i,j,newstrc[k].global_evecs[j][i]);	
		 				}
		 				gsl_matrix_set(eval,i,i,newstrc[k].main_vars[i]);
		 		}
				multiplie_matrix(evec,3,3,eval,3,3,temp2);
				multiplie_matrix(temp2,3,3,tevec,3,3,temp);
				

	 			fprintf(out_file,"ANISOU");
	 			fprintf(out_file,"%5.d %s%s %s%4.d  %7d%7d%7d%7d%7d%7d\n",
		 			newstrc[k].atom_number,
		 			newstrc[k].atom_prot_type,
		 			newstrc[k].res_type,
		 			newstrc[k].chain,
		 			newstrc[k].res_number,
		 			(int) (gsl_matrix_get(temp,0,0)*factor),
		 			(int) (gsl_matrix_get(temp,1,1)*factor),
		 			(int) (gsl_matrix_get(temp,2,2)*factor),
		 			(int) (gsl_matrix_get(temp,0,1)*factor),
		 			(int) (gsl_matrix_get(temp,0,2)*factor),
		 			(int) (gsl_matrix_get(temp,1,2)*factor)
	 				
	 			);
	 			tot += 6.0;
	 			avg += gsl_matrix_get(temp,0,0)*factor;
	 			avg += gsl_matrix_get(temp,1,1)*factor;
	 			avg += gsl_matrix_get(temp,2,2)*factor;
	 			avg += gsl_matrix_get(temp,0,1)*factor;
	 			avg += gsl_matrix_get(temp,0,2)*factor;
	 			avg += gsl_matrix_get(temp,1,2)*factor;
	 			
	 			gsl_matrix_free(temp);
	 			gsl_matrix_free(temp2);
	 			gsl_matrix_free(eval);
	 			gsl_matrix_free(evec);
	 			gsl_matrix_free(tevec);
	 		}
 	}
 	 fclose(out_file);
 	 return(avg/tot);
 }
