#include "STeM.h"

int load_list(char filename[100],int **con);

int main(int argc, char *argv[]) {
	int all; /*Nombre d'atomes dans pdb*/
	int atom; /*Nombre de carbone CA*/
	int all_t; /*Nombre d'atomes dans pdb*/
	int atom_t; /*Nombre de carbone CA*/
 	int help_flag = 1;
 	char file_name[500];
 	char check_name[500];
 	int verbose = 0;
	int i,j;
	int lig = 0;
	int nconn;
	char eigen_name[500] = "eigen.dat";
	char eigen_name_two[500] = "eigen.dat";
	char init_list[500] = "init_list.dat";
	char targ_list[500] = "targ_list.dat";
	float scale = 1.0;
 	for (i = 1;i < argc;i++) {
 		if (strcmp("-i",argv[i]) == 0) {strcpy(file_name,argv[i+1]);--help_flag;}
 		if (strcmp("-h",argv[i]) == 0) {help_flag = 1;}
 		if (strcmp("-v",argv[i]) == 0) {verbose = 1;}
		if (strcmp("-ieig",argv[i]) == 0) {strcpy(eigen_name,argv[i+1]);}
		if (strcmp("-teig",argv[i]) == 0) {strcpy(eigen_name_two,argv[i+1]);}
		
		if (strcmp("-ilist",argv[i]) == 0) {strcpy(init_list,argv[i+1]);}
		if (strcmp("-tlist",argv[i]) == 0) {strcpy(targ_list,argv[i+1]);}
		
 		if (strcmp("-scale",argv[i]) == 0) {float temp;sscanf(argv[i+1],"%f",&temp);scale = temp;}
 		if (strcmp("-t",argv[i]) == 0) {strcpy(check_name,argv[i+1]);help_flag = 0;}
 		

 		
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
	atom = build_all_strc(file_name,strc_all); // Retourne le nombre de Node
	if (atom > 800) {printf("Too much node.... To fix, ask vincent.frappier@usherbrooke.ca\n");return(1);}
	if (verbose == 1) {printf("	Atom:%d\n",all);}
	check_lig(strc_all,connect_h,nconn,all);
	
	// Assign les Nodes
	
	if (verbose == 1) {printf("	CA Structure\n");}
	
	if (verbose == 1) {printf("	Node:%d\n",atom);}
	struct pdb_atom strc_node[atom];
	atom = build_cord_CA(strc_all, strc_node,all,lig,connect_h,nconn);
	
	if (verbose == 1) {printf("	Assign Node:%d\n",atom);}
	
	// Free Connect
		
	//for(i=0;i<nconn;i++) {printf("I:%d\n",i);free(connect_h[i]);}
	//free(connect_h);
	
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
	atom_t = build_all_strc(check_name,strc_all_t); // Retourne le nombre de Node
	if (atom_t > 800) {printf("Too much node.... To fix, ask vincent.frappier@usherbrooke.ca\n");return(1);}
	if (verbose == 1) {printf("	Atom:%d\n",all_t);}
	check_lig(strc_all_t,connect_t,nconn,all_t);
	
	// Assign les Nodes
	
	if (verbose == 1) {printf("	CA Structure\n");}
	
	if (verbose == 1) {printf("	Node:%d\n",atom_t);}
	struct pdb_atom strc_node_t[atom_t];

	atom_t = build_cord_CA(strc_all_t, strc_node_t,all_t,lig,connect_t,nconn);
	if (verbose == 1) {printf("	Assign Node:%d\n",atom_t);}
 
 	int align[atom];
 	int score = node_align(strc_node,atom,strc_node_t,atom_t,align);
 	printf("RMSD:%8.5f Score: %d/%d\n",sqrt(rmsd_no(strc_node,strc_node_t,atom, align)),score,atom);
 	if ((float)score/(float)atom < 0.8 && (float)score/(float)atom_t < 0.8 ) {
 		printf("Low Score... Will try an homemade alignement !!!\n");
 		score = node_align_low(strc_node,atom,strc_node_t,atom_t,align);
 		printf("RMSD:%8.5f Score: %d/%d\n",sqrt(rmsd_no(strc_node,strc_node_t,atom, align)),score,atom);
 	}
 	
 	
 	// **********************************************
 	// Load Eigen de la strc init
 	// **********************************************
 	//***************************************************
	//*                                                                                                     *
	//* Load eigenvector et eigenvalue									*
	//*                                                                                                     *
	//***************************************************
	
	if (verbose == 1) {printf("Loading Eigenvector One\n");}
	
	gsl_vector *eval = gsl_vector_alloc(3*atom);
	
	gsl_matrix *evec= gsl_matrix_alloc (3*atom,3*atom);
	
	load_eigen(eval,evec,eigen_name,3*atom);
	
	if (verbose == 1) {printf("Loading Eigenvector Two\n");}
	
	gsl_vector *eval_two = gsl_vector_alloc(3*atom_t);
	
	gsl_matrix *evec_two= gsl_matrix_alloc (3*atom_t,3*atom_t);
	
	load_eigen(eval_two,evec_two,eigen_name_two,3*atom_t);
	
	// Scaling eval
	
	gsl_vector_scale(eval,scale);
	gsl_vector_scale(eval_two,scale);
	
	if (verbose == 1) {printf("Building Cross-Correlation Two\n");}
	
	gsl_matrix *k_totinv_two = gsl_matrix_alloc(atom_t*3, atom_t*3); /* Déclare et crée une matrice qui va être le pseudo inverse */
	
	gsl_matrix_set_all(k_totinv_two, 0);
	
	k_cov_inv_matrix_stem(k_totinv_two,atom_t,eval_two,evec_two,6,atom_t*3-6); /* Génère une matrice contenant les superéléments diagonaux de la pseudo-inverse. */
	
	// Il va falloir crer un objet qui va etre constitué de toutes les align à comparer
	
	// Load list 1
	if (verbose == 1) {printf("Loading list init:%s\n",init_list);}
	int **pair_one=(int **)malloc(1000*sizeof(int *)); 
  for(i=0;i<1000;i++) { pair_one[i]=(int *)malloc(10*sizeof(int));}
  
  int num_list_init = load_list(init_list,pair_one);
  
  if (verbose == 1) {printf("Loading list targ:%s\n",targ_list);}
  int **pair_two=(int **)malloc(1000*sizeof(int *)); 
  for(i=0;i<1000;i++) { pair_two[i]=(int *)malloc(10*sizeof(int));}
	
	int num_list_targ =load_list(targ_list,pair_two);
	
	
	int paira;
	for(paira = 0;paira<num_list_init;++paira) {
		int pairb;
		for(pairb = 0;pairb<num_list_targ;++pairb) {
		
			
			// Test manuel, va changer l'align
			int k;
			for(k=0;k<atom;++k) {align[k] = -1;} // Reset align
			// Assign pair A
			for(k=0;k<6;++k) {
				if (pair_one[paira][k] != -1) {
					align[pair_one[paira][k]] = pair_two[pairb][k];
					
				}
			}

			score = 0;
			for(k=0;k<atom;++k) {if(align[k] != -1) {++score;}}
		 	if (verbose == 1) {printf("\tRotating Eigenvector One\n");}
			
			// Rotationne les eigenvector
	
			float myRmsd = rmsd_yes_eigen(strc_node,strc_node_t,atom, align,strc_all,all,evec);
			//printf("\tRMSD = %f\n",myRmsd);
			if (verbose == 1) {printf("\tBuilding Cross-Correlation One\n");}

			// Transforme en matrice de correlation
	
			gsl_matrix *k_totinv = gsl_matrix_alloc(atom*3, atom*3); /* Déclare et crée une matrice qui va être le pseudo inverse */
	
			gsl_matrix_set_all(k_totinv, 0);
	
			k_cov_inv_matrix_stem(k_totinv,atom,eval,evec,6,atom*3-6); /* Génère une matrice contenant les superéléments diagonaux de la pseudo-inverse. */

	
	
	
	
	
			// Il faut convolutionner les deux matrices de covariance (pour les atomes correspondants !)

			if (verbose == 1) {printf("\tAdding Two Cross-Correlation Matrix\n");}	
			// Variable importante pour mon inversion
			gsl_matrix *incov12 = gsl_matrix_alloc(score*3,score*3);
			double conj_dens12 = 1;
			gsl_matrix *cov12 = gsl_matrix_alloc(score*3,score*3);
			gsl_matrix_set_all(cov12,0);
			// Build mix des deux cov, il faut les additioner
	
			int kept_one[score];
			int kept_two[score];
			int count = 0;
	
			gsl_vector *delr = gsl_vector_alloc(score*3); // Differene entre deux positions
			gsl_vector *pos = gsl_vector_alloc(score*3); // Le vecteur de différence qu'on va vouloir evaluer
	
			gsl_vector_set_all(pos,0);
	
			for(i=0;i<atom;++i) {
				if (align[i] != -1) {
					kept_one[count] = i;
					kept_two[count] = align[i];
					gsl_vector_set(delr,count*3+0,strc_node[i].x_cord-strc_node_t[align[i]].x_cord);
					gsl_vector_set(delr,count*3+1,strc_node[i].y_cord-strc_node_t[align[i]].y_cord);
					gsl_vector_set(delr,count*3+2,strc_node[i].z_cord-strc_node_t[align[i]].z_cord);
					++count;
				}
			}
	
			/*print_matrix(k_totinv);printf("\n\n");
	
			print_matrix(k_totinv_two);printf("\n\n");*/
	
			//gsl_vector_scale(delr,-1);
			for(i = 0;i<score;++i) {
				for(j = 0;j<score;++j) {
					int k,l;
					for (l = 0;l<3;++l) {for (k = 0;k<3;++k) {
					gsl_matrix_set(cov12,i*3+k,j*3+l,
						gsl_matrix_get(k_totinv,kept_one[i]*3+k,kept_one[j]*3+l)+gsl_matrix_get(k_totinv_two,kept_two[i]*3+k,kept_two[j]*3+l)
					);
					}}
	
				}
	
			}

	
			// Fonction qui build ma nouvelle matrice de covariance pour ma conjugé
	
		/*	print_matrix(cov12);
			printf("\n\n");*/
	
			// On va essayer de la diagonalyser
	
			gsl_vector *eval_cov = gsl_vector_alloc(3*score); 
			gsl_vector_set_all(eval_cov,0);
			gsl_matrix *evec_cov = gsl_matrix_alloc (3*score,3*score); 
			gsl_matrix_set_all(evec_cov,0);
			diagonalyse_matrix(cov12,3*score,eval_cov,evec_cov);
			/*print_vector(eval_cov);
	
			print_matrix(evec_cov);*/

		const double PI = 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421;
			k_cov_inv_matrix_stem(incov12,score,eval_cov,evec_cov,0,3*score);
	
	
	
			conj_dens12= -0.5*score*3*log(2*PI);
			for(i=0;i<score*3;++i) {
				 //if (i < 30) {printf("I:%d -> %g Log -> %g\n",i, gsl_vector_get (eval_cov, i),log(gsl_vector_get (eval_cov, i)));}
				 if  (gsl_vector_get (eval_cov, i) < 0.000001) {continue;}
		
				 conj_dens12 -= 0.5*log(gsl_vector_get (eval_cov, i));
			}
			// Re-verse

		/*	conj_prob_init_n(incov12 ,&conj_dens12,cov12,score*3);
	
			print_matrix(incov12);
			printf("\n\n");*/
	
			//printf("\tConj_dens = %g\n",conj_dens12);
			float dens = density_prob_n(incov12, delr, conj_dens12,pos,score*3);
			//	printf("Pair:%d et %d\n",paira,pairb);
			printf("%d %d %f %g\n",paira,pairb,myRmsd,dens);
		
			// Free some stuff
			gsl_matrix_free(k_totinv);
			gsl_matrix_free(evec_cov);		
			gsl_matrix_free(cov12);
			gsl_matrix_free(incov12);
		
			gsl_vector_free(eval_cov);
			gsl_vector_free(pos);
			gsl_vector_free(delr);
		}
	}
}


int load_list(char filename[100],int **con) {
	FILE *file;
 	file = fopen(filename,"r");
 	if (file == NULL) {return(-1);}

 	char line[100];
 	int temp1,temp2,temp3,temp4,temp5,temp6;
 	int j=-1;
 	while(fgets(line,100,file)) {
 		
 			++j;
 			temp1 =-1;
 			temp2 =-1;
 			temp3 =-1;
 			temp4 =-1;
 			temp5 =-1;
 			temp6 =-1;
 			//printf("1:%d\t2:%d\t3:%d\t4:%d\n",temp1,temp2,temp3,temp4);
 			sscanf(line, "%d %d %d %d %d %d",&temp1,&temp2,&temp3,&temp4,&temp5,&temp6);
 			con[j][0] =temp1;
 			con[j][1] =temp2;
 			con[j][2] =temp3;
 			con[j][3] =temp4;
 			con[j][4] =temp5;
 			con[j][5] =temp6;
 			//printf("J:%d Value:%d %d -- %d %d %d %d %d %d\n",j,temp6,con[j][5],temp1,temp2,temp3,temp4,temp5,temp6);
 			//printf("%s",line);
 		
	}
	
	//printf("FCLOSE:%d\n",fclose(file));
	fclose(file);
	return(j+1);



}
