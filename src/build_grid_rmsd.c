#include "STeM.h"
 
 // Construit une grille en se basant sur le RMSD comparativement à la forme initiale -> utiliser pour le docking
 
 
 void find_rmsd_step(struct pdb_atom *strc,gsl_matrix *m,int mode,int atom,float step,double *ampli);
int build_grid_rmsd(double **grid,double resolution,int nbr_mode,int atom,gsl_matrix *matrix,struct pdb_atom *strc,float limit,int mode, int now, double *actual);
int test_point_rmsd(struct pdb_atom *strc, double *actual, float limit, int atom,int mode, int nb_mode, gsl_matrix *matrix,double **grid,int total);
double find_zero(float a,float b,float c,int what);
int build_grid_math(double **grid,int nb_mode,float step,float limit,float resolution,struct pdb_atom *strc,int atom,gsl_matrix *matrix);
void mode_switching(gsl_vector *eval,gsl_matrix *evec,int atom1,int atom2,int *mode_list,int nm);
// Programme


 int main(int argc, char *argv[]) {

	int all; /*Nombre d'atomes dans pdb*/
	int atom; /*Nombre de carbone CA*/
 	int help_flag = 1;
 	char file_name[100];
 	char eigen_name[100] = "eigen.dat";
 	char grid_name[100] = "grid_motion.pdb";
 	char out_name[100] = "grid.dat";
 	char eigen_grid[100] = "eigen_grid.dat";
 	int verbose = 0;
	int nbr_mode = 2;
 	float resolution = 0.5;
	int mode = 6;
	int i;
	float limit = 0.01;
	int print = 0;
	int lig = 0;
	int nconn;
	int mode_list[20];
	int modef = 0;
 	for (i = 1;i < argc;i++) {
 		if (strcmp("-i",argv[i]) == 0) {strcpy(file_name,argv[i+1]);--help_flag;}
 		if (strcmp("-o",argv[i]) == 0) {strcpy(out_name,argv[i+1]);}
 		if (strcmp("-ieig",argv[i]) == 0) {strcpy(eigen_name,argv[i+1]);}
 		if (strcmp("-g",argv[i]) == 0) {strcpy(eigen_grid,argv[i+1]);}
 		if (strcmp("-h",argv[i]) == 0) {help_flag = 1;}
 		if (strcmp("-v",argv[i]) == 0) {verbose = 1;}
 		if (strcmp("-m",argv[i]) == 0) {int temp;sscanf(argv[i+1],"%d",&temp);mode = temp-1;} // Mode -1 pour éviter tout mélange
 		if (strcmp("-nm",argv[i]) == 0) {int temp;sscanf(argv[i+1],"%d",&temp);nbr_mode = temp;}
 		if (strcmp("-step",argv[i]) == 0) {float temp;sscanf(argv[i+1],"%f",&temp);resolution = temp;}
 		if (strcmp("-md",argv[i]) == 0) {float temp;sscanf(argv[i+1],"%f",&temp);limit = temp;}
 		if (strcmp("-lig",argv[i]) == 0) {lig = 1;}       
 		if (strcmp("-p",argv[i]) == 0) {print = 1;strcpy(grid_name,argv[i+1]);}
 		if (strcmp("-ml",argv[i]) == 0) {
 			printf("I detect a list of modes L:");
 			nbr_mode = 0;
 			++modef;
 			while(1) {
 				int temp;
 				if (i+nbr_mode+1 > argc-1) {break;}
 				int match = sscanf(argv[i+nbr_mode+1],"%d",&temp);
 				if (match == 0) {break;}
 				++nbr_mode;
 				printf("%d ",temp);
 				mode_list[nbr_mode-1] = temp-1;
 			}
 			printf("\nI have now %d modes (-nm)\n",nbr_mode);
 		}
 		
 		
 	}

 	// Printing Help
 	
 	if (help_flag == 0) { } else {
 		printf(
 			"****************************Help Section\n-i\tFile Input (PDB)\n-o\tOutput Name Grid\n-ieig\tFile Name Eigen\n-v\tVerbose\n-step\tMinimum RMSD between each points\n-m\tMode\n-md\tMax RMSD\n-nm\tNombre de mode\nlig\tTient compte du ligand\n-p\tprint la grid\n****************************\n");
 		return(0); 
 	} 
	
 	//****************************************************
 	//*																									 *
 	//*Build a structure contaning information on the pdb*
 	//*																									 *
 	//****************************************************
 	
 	all = count_atom(file_name);
 	nconn = count_connect(file_name);
 	
 	if (verbose == 1) {printf("Connect:%d\n",nconn);}
 	
	if (verbose == 1) {printf("Assigning Structure\n\tAll Atom\n");}
	
	// Array qui comprend tous les connects
	
	int **connect_h=(int **)malloc(nconn*sizeof(int *)); 
    for(i=0;i<nconn;i++) { connect_h[i]=(int *)malloc(7*sizeof(int));}
    if (verbose == 1) {printf("Assigning Connect\n");}
    
    assign_connect(file_name,connect_h);
	
	// Assign tous les atoms
	
	struct pdb_atom strc_all[all];
	atom = build_all_strc(file_name,strc_all); // Retourne le nombre de Node
	
	if (verbose == 1) {printf("	Node:%d\n	Atom:%d\n",atom,all);}

	check_lig(strc_all,connect_h,nconn,all);

	// Assign les Nodes
	
	if (verbose == 1) {printf("	CA Structure\n");}
	
	if (verbose == 1) {printf("	Node:%d\n",atom);}
	struct pdb_atom strc_node[atom];
	atom = build_cord_CA(strc_all, strc_node,all,lig,connect_h,nconn);
	if (verbose == 1) {printf("	Assign Node:%d\n",atom);}
	
	 //***************************************************
 	//*													*
 	//* Load eigenvector et eigenvalue					*
 	//*													*
 	//***************************************************
 	
 	if (verbose == 1) {printf("Loading Eigenvector\n");}
 	
 	gsl_vector *eval = gsl_vector_alloc(3*atom);
	gsl_matrix *evec= gsl_matrix_alloc (3*atom,3*atom);
	load_eigen(eval,evec,eigen_name,3*atom);
	
	//***************************************************
	//*													*
	//*Load Eigenvector and Eigenvalue and creater image
	//*													*
	//***************************************************
 	
 	double **grid=(double **)malloc(3000000*sizeof(double *));
 	for(i=0;i<(3000000);i++) { grid[i]=(double *)malloc(nbr_mode*sizeof(double));}
 	
 	// Building grid
 	
 	if (verbose == 1) {printf("\tCreating Grid\n");}
 	
 	

 	double step;
 	
 	

 		
 	find_rmsd_step(strc_node,evec,mode,atom,resolution,&step);
 	step = step*step;

	printf("Step:%f Eigen:%f\n",step,gsl_vector_get(eval,mode));
	
	// Generating grid
	// Ce qui va contenir la grid et tester les points
// 	double actual[nbr_mode];
//for(i=0;i<nbr_mode;++i) {actual[i] = 0;}
 	
//	int nb_grid = build_grid_rmsd(grid,step,nbr_mode,atom,evec,strc_node,limit,mode,0,actual);
 	
 	int nb_grid = build_grid_math(grid,nbr_mode,step,limit,resolution,strc_node,atom,evec) ;
 	printf("Points:%d\n",nb_grid);
 	// Assigning and writing
 	
 	gsl_matrix *grid_m = gsl_matrix_alloc (nb_grid,nbr_mode);
	assignArray(grid_m,grid,nb_grid,nbr_mode);
 	
 	if (verbose == 1) {printf("Writting Things\n");}	 
	if (verbose == 1) {printf("	Write Grid\n");}	   	
	write_grid_mat(out_name,grid_m,nb_grid,nbr_mode);
	if (verbose == 1) {printf("	Write Eigen\n");}
	if (modef != 0) {
		mode_switching(eval,evec,3*atom,3*atom,mode_list,nbr_mode);
		mode = 6;
	}  	
	write_eigen_mat(eigen_grid, evec,3*atom,3*atom,mode,nbr_mode);
	//write_matrix_pdb("grid_matrix.pdb", grid_m,nb_grid,nbr_mode);
	if (verbose == 1) {printf("	Write Motion\n");}	
	if (print == 1) {print_grid_motion(strc_all,evec,mode,all,grid_name,grid_m,nb_grid,nbr_mode,lig);}
 	
 	
	// Freeing
	
	if (verbose == 1) {printf("Freeing\n");}	

	gsl_vector_free(eval);
	gsl_matrix_free(evec);
	
		
	
	
	for(i=0;i<(200000);i++) { free(grid[i]);}
	free(grid);
}

void find_rmsd_step(struct pdb_atom *strc,gsl_matrix *m,int mode,int atom,float step,double *ampli) {
		
	// En temps normal, le RMSD varie de la fonction suivante: RMSD(x) = Ax^2... Il faut trouver le A
	int i;
	
	int temp_align[atom];
	
	for (i=0;i<atom;++i) { temp_align[i] = i;}
	
	
	// Crée une structure temporaire
	struct pdb_atom temp_strc[atom];
	copy_strc(temp_strc, strc, atom);
	
	// Applique l'amplitude
	apply_eigen(temp_strc,atom,m,7,1.0);
		
	
	
	
	
	
	
	*ampli = sqrt(rmsd_no(strc,temp_strc,atom, temp_align));

}

 int build_grid_rmsd(double **grid,double resolution,int nbr_mode,int atom,gsl_matrix *matrix,struct pdb_atom *strc,float limit,int mode, int now, double *actual) {
	
	int i,l=now;
	double passing[nbr_mode]; // Ce qui va etre teste
	for(i=0;i<nbr_mode;++i) {passing[i] = actual[i];}
		
	int ok = test_point_rmsd(strc, actual, limit, atom,mode, nbr_mode,matrix,grid,l);
	
	// Retourne dans le itératif	
	if (ok == 0) {return(l);}
	
	// Le point est accepté
	if (ok == 1) {
		for(i=0;i<nbr_mode;++i) {
			grid[l][i] = actual[i];
		}
		++l;
	}
	
	// Test dans une direction
	for (i=0;i<nbr_mode;++i) {
		// Faut changer le resolution
		
		double temp_reso = find_zero(1,2*passing[i],-1/resolution,1);
		passing[i] += temp_reso;
		l = build_grid_rmsd(grid,resolution,nbr_mode,atom,matrix,strc,limit,mode,l,passing);
	}
	
	// Reste les passing... parce qu'on avait mod avant
	for(i=0;i<nbr_mode;++i) {passing[i] = actual[i];}
	printf("Will be testing the other direction\n");
	for (i=0;i<nbr_mode;++i) {printf("Passing[%d] = %f\n",i,passing[i]);}
	// Test dans l'autre direction
	for (i=0;i<nbr_mode;++i) {
		
		double temp_reso = find_zero(1,2*passing[i],1/resolution,1);
		printf("Temp reso:%f\n",temp_reso);
		passing[i] += temp_reso;
		l = build_grid_rmsd(grid,resolution,nbr_mode,atom,matrix,strc,limit,mode,l,passing);
	}
	
	// Return... this is the end
	return(l);
}

int test_point_rmsd(struct pdb_atom *strc, double *actual, float limit, int atom,int mode, int nb_mode, gsl_matrix *m,double **grid,int total) {
 	
 	int i,j;
 	double epsilon = 0.001;
	// Regarde si le point a déjà été testé
	int temp = 0;
	
	// Pour chaque point dans la grid
	for (i=0;i<total;++i) {
		temp = 0; // Compte si le point[mode] est deja la
		for (j=0;j<nb_mode;++j) {
			if (grid[i][j] > actual[j] - epsilon && grid[i][j] < actual[j] + epsilon) {
				++temp;
			} else {
				// Match pas passe au prochain point
				break;
			}
		}
		
		if (temp == nb_mode) {return(0);} // Le points à déjà été testé
		
	}
	
 	// Crée une structure temporaire
	struct pdb_atom temp_strc[atom];
	copy_strc(temp_strc, strc, atom);
	
	int temp_align[atom];	
	
	for (i=0;i<atom;++i) { temp_align[i] = i;}
	// Applique l'amplitude
	for (i=0;i<nb_mode;++i) {
		printf("%.3f ",actual[i]);
		apply_eigen(temp_strc,atom,m,mode+i,actual[i]);
	}
 	printf("RMSD:%.5f\n",rmsd_no(strc,temp_strc,atom, temp_align));
 	if (rmsd_no(strc,temp_strc,atom, temp_align) > limit) {return(0);} else {return(1);}
 	
}

double find_zero(float a,float b,float c,int what) {
	if (what == 1) {
		return((-b + sqrt(b*b-4*a*c))/(2*a));
	} else {
		return((-b - sqrt(b*b-4*a*c))/(2*a));
	}
	return(0);
}

int build_grid_math(double **grid,int nb_mode,float step,float limit,float resolution,struct pdb_atom *strc,int atom,gsl_matrix *matrix) {
	
	// Construit une box et teste les points
	// Sa va être rapide parce que c est mathématiques !
	// y(x,y,z) = ax²+ay²+az²
	
	int i,j;
	
	float actual[nb_mode];
	int now[nb_mode];
	float grid_lim = sqrt(limit*limit/step);
	for (i=0;i<nb_mode;++i) {
		actual[i] = 0;
		now[i] = 0;
	}
	int exit = 1;
	int count = 0;
	int ok;
	float rmsd;
/*	int align[atom];
	for (i=0;i<atom;++i) {
		align[i] = i;
		
	}*/
	
	
	// Trouve toutes les bonnes amplitudes a appliquer avant de simuler...vraiment pas bon, mais bon
	int lenght = limit/resolution+10;
	float ampli[lenght];
	ampli[0] = 0;
	i = 0;
	while(exit) {
		rmsd = 0;
		printf("ampli[%d] = %.3f\n",i,ampli[i]);
		rmsd += step*ampli[i]*ampli[i];
		float temp_resol = (sqrt(rmsd)+resolution)*(sqrt(rmsd)+resolution)-rmsd;
		ampli[i+1] = ampli[i];
		ampli[i+1] += find_zero(1,2*ampli[i],-temp_resol/step,1);
		++i;
		if (rmsd > limit*limit) {exit = 0;}
	}
	exit = 1;
	while(exit) {
		
		
		
		// Calcul le RMSD de facon math
		rmsd = 0;
		
	//	struct pdb_atom temp_strc[atom];
	//	copy_strc(temp_strc, strc,  atom);
		
		for(i=0;i<nb_mode;++i) {
		//	printf("	Actual[%d] = %f\n",i,actual[i]);
		//	apply_eigen(temp_strc,atom,matrix,6+i,actual[i]);
			rmsd += step*actual[i]*actual[i];
		}
		
		
		
		
	//	for (i=0;i<nb_mode;++i) {printf("%.3f ",actual[i]);}printf("RMSD:%.3f\n",rmsd);
	
		// Test si dépasse la limite
		if (rmsd > limit*limit+0.1) {ok = 0;} else {ok = 1;}
		
		// Si c est vrai, copier dans la grid
		if (ok == 1) {
		//	printf("Predicted rms:%.4f Real rmsd:%.4f Real rms:%.4f\n",rmsd,sqrt(rmsd_no(temp_strc,strc,atom,align)),rmsd_no(temp_strc,strc,atom,align));
			for(i=0;i<nb_mode;++i) {
				grid[count][i] = actual[i];
			}
			++count;
		}
		
		// Modifie le points a tester
	//	float temp_resol = (sqrt(rmsd)+resolution)*(sqrt(rmsd)+resolution)-rmsd;
	/*	printf("I have a rms of %.4f and so my rmsd = %.4f\n",rmsd,sqrt(rmsd));
		printf("I want a RMSD of %.4f sor my rms need to be %.4f\n",sqrt(rmsd)+resolution,(sqrt(rmsd)+resolution)*(sqrt(rmsd)+resolution));
		printf("I Have a temp_step %.4f\n",temp_resol);
		printf("\n");*/
		++now[0];
		actual[0] = ampli[now[0]];
		
		for (i=0;i<nb_mode-1;++i) {
			if (actual[i] > grid_lim) {
				now[i] = 0;
				actual[i] = ampli[now[i]];
				++now[i+1];
				actual[i+1] = ampli[now[i+1]];
			}
		}
		if (actual[nb_mode-1] > grid_lim) {printf("actual[%d] = %f > %f\n",nb_mode-1,actual[nb_mode-1] ,grid_lim);exit=0;}
	//	if (count > 3) {exit = 0;}
		
	}

	
	// On a juste les points positifs de la grille, il faut mettre les points neg
	// ex (1,1) = (1,1),(-1,1),(1,-1),(-1,-1)

    printf("I make mirror image --> %d\n",count);

	int sign[nb_mode];
	int old_count = count;
	for (i=0;i<nb_mode;++i) {sign[i] = -1;}
	
	while(1) {
		
		int temp = 0;
		for (i=0;i<nb_mode;++i) {
			if (sign[i] == 1) {++temp;}
			
		}
		if (temp == nb_mode) {break;}
	//	for (i=0;i<nb_mode;++i) {printf("%2d ",sign[i]);} printf("\n");
		// Modifie
		
		for (j=1;j<old_count;++j) {
			for(i=0;i<nb_mode;++i) {
					grid[count][i] = grid[j][i]*sign[i];
			}
			++count;
			if (count % 10000 == 0) {printf("I have %d points and counting\n",count);}
		}
		
		
		sign[0] += 2;
		for (i=0;i<nb_mode-1;++i) {
			if (sign[i] > 1) {
				sign[i] = -1;
				sign[i+1] += 2;
			}
		}
		
		// Checks exit
		
		
		
		if (sign[nb_mode-1] > 1) {break;}
		
	}
	
	return(count);
	
	
}


void mode_switching(gsl_vector *eval,gsl_matrix *evec,int atom1,int atom2,int *mode_list,int nm) {
	// Function qui remplace les mode 7+ par la liste de mode

	gsl_vector *teval = gsl_vector_alloc(3*atom1);
	gsl_matrix *tevec= gsl_matrix_alloc (3*atom1,3*atom2);
	int i,j;
	
	// Set Eval in temp
	
	for (i = 0;i<atom1;++i) {
		gsl_vector_set(teval,i,gsl_vector_get(eval,i));
	}
	for (i = 0;i<atom1;++i) {
		for (j = 0;j<atom2;++j) {
			gsl_matrix_set(tevec,i,j,gsl_matrix_get(evec,i,j));
		}
	}
	
	// Replace in eval
	
	for (i=6;i<6+nm;++i) {
		gsl_vector_set(eval,i,gsl_vector_get(teval,mode_list[i-6]));
	}
	
	// Replace in evec
	
	for (i=6;i<6+nm;++i) {
		for (j = 0;j<atom2;++j) {
			gsl_matrix_set(evec,i,j,gsl_matrix_get(tevec,mode_list[i-6],j));
		}
	}
	


}
