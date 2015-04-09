#include "STeM.h"

void fit_svd_fold(struct pdb_atom *init,struct pdb_atom *targ,int atom,int all,int atom_t,int all_t,struct pdb_atom *all_init,struct pdb_atom *all_targ,gsl_matrix *evec, int *align, int nb_mode, int mode,gsl_vector *eval,gsl_vector *answer) {
	int i,j,l,k,m;
 	float newstrc,old;
	double amp[nb_mode];
	float quad[3];
	float b,a;
	int alit = 0;// Value dans le vecteur diff
	int ngila[atom_t]; for (i=0;i<atom_t;++i) {ngila[i] = -1;  }
	for (i=0;i<atom;++i) {
		if (align[i] !=-1 ){ ngila[align[i]] = i; }
	}

	center_yes(init,targ,atom,atom_t, align); // Targ rotate, donc pas d'impact sur vecteurs
	rmsd_yes(targ,init,atom_t, ngila,all_targ,all_t);

	
 	gsl_vector *diff     = gsl_vector_alloc(3*atom);
 	
 	// Minimization
 	old = rmsd_no(init,targ,atom,align);
	// Construit vector des différences

 	for(i=0;i<atom;++i) {
 		j = align[i];
 		if (j == -1) {continue;}
 		//printf("Res: %s :: %s Num: %d :: %d\n",init[i].res_type,targ[j].res_type,init[i].res_number,targ[j].res_number);
 		gsl_vector_set(diff,alit*3+0,init[i].x_cord - targ[j].x_cord);
		gsl_vector_set(diff,alit*3+1,init[i].y_cord - targ[j].y_cord);
		gsl_vector_set(diff,alit*3+2,init[i].z_cord - targ[j].z_cord);
	//	printf("I:%d J:%d Alit:%d Init:(%f,%f,%f) Targ:(%f,%f,%f) Diff:(%f,%f,%f)\n",i,j,alit,init[i].x_cord,init[i].y_cord,init[i].z_cord,targ[j].x_cord,targ[j].y_cord,targ[j].z_cord,gsl_vector_get(diff,alit*3+0),gsl_vector_get(diff,alit*3+1),gsl_vector_get(diff,alit*3+2));
		++alit;
		
 	}
 	//printf("Alit:%d\n",alit);
 	gsl_vector *B = gsl_vector_alloc((alit)*3);
 	for(i=0;i<alit;++i) {
 		//printf("I:%4d Vector:(%8.5f,%8.5f,%8.5f)\n",i,gsl_vector_get(diff,i*3+0),gsl_vector_get(diff,i*3+1),gsl_vector_get(diff,i*3+2));
 		gsl_vector_set(B,i*3+0,gsl_vector_get(diff,i*3+0));
 		gsl_vector_set(B,i*3+1,gsl_vector_get(diff,i*3+1));
 		gsl_vector_set(B,i*3+2,gsl_vector_get(diff,i*3+2));
 	}

	
	for(j=0;j<nb_mode;++j) {
		amp[j] = 0;
	}
	
	// On veut ||Ax-b||² soit minimiser
	// http://www.ces.clemson.edu/~petersj/Agents/MatLabNA/MatLabNA004.html
	// A: Matrice des eigenvectors
	// B: Le vecteur de différence
	// x: Vecteur des amplitudes a appliquer sur chacun des Eigenvectors.... CE QU ON VEUT TROUVER
	// On décompose A = USV
	// Nouveau vecteur C = U^T*b
	// Y = S(diag)/C... Y = S(i,i)/C(i)
	// x = Y*V
	// e = Ax-b.... RMSD = ||e||
	
	// Go !!!
	// Construit matrix A
	gsl_matrix *A = gsl_matrix_alloc(3*alit,nb_mode);
	
	for (i=0;i<nb_mode;++i) {
		l = -1;
		for (j=0;j<atom;++j) {		
			if (align[j] == -1) {continue;}
			++l;
			//printf("Vector pour %d node: (%f,%f,%f) et eigen (%f,%f,%f)\n",j,gsl_vector_get(B,l*3),gsl_vector_get(B,l*3+1),gsl_vector_get(B,l*3+2),gsl_matrix_get(evec,j*3+0,i+mode-1),gsl_matrix_get(evec,j*3+1,i+mode-1),gsl_matrix_get(evec,j*3+2,i+mode-1));
			//printf("I:%d L:%d J:%d :: %d Max:%d\n",i,l,j,align[j],l*3+2);
			gsl_matrix_set(A,l*3+0,i,gsl_matrix_get(evec,j*3+0,i+mode-1));
			gsl_matrix_set(A,l*3+1,i,gsl_matrix_get(evec,j*3+1,i+mode-1));
			gsl_matrix_set(A,l*3+2,i,gsl_matrix_get(evec,j*3+2,i+mode-1));
			
		/*	gsl_matrix_set(A,l*3+0,i,gsl_matrix_get(evec,i+mode-1,j*3+0));
			gsl_matrix_set(A,l*3+1,i,gsl_matrix_get(evec,i+mode-1,j*3+1));
			gsl_matrix_set(A,l*3+2,i,gsl_matrix_get(evec,i+mode-1,j*3+2));*/
			
		}
	
	}
	
	// Décompose A
	
	
	gsl_matrix *mat_U = gsl_matrix_alloc(nb_mode,nb_mode);
	gsl_vector *vec_S = gsl_vector_alloc(nb_mode);
	gsl_vector *vec_w = gsl_vector_alloc(nb_mode);
	gsl_vector *result = gsl_vector_alloc(nb_mode);
	//printf("Decomposition\n");
	gsl_linalg_SV_decomp(A, mat_U,vec_S,vec_w);
	//printf("SVD\n");
	gsl_linalg_SV_solve (A, mat_U, vec_S,B, result);
	
	for(i=0;i<nb_mode;++i) {
		
		amp[i] = -gsl_vector_get(result,i);
		//if (overlap(atom,mode+i-1,evec,diff,align) < 0.01) {amp[i] = 0.0;}
		gsl_vector_set(answer,i,amp[i]);
	}
		float energy = 0;
	/*for(j=0;j<nb_mode;++j) {
 		old = overlap(atom,mode+j-1,evec,diff,align);	
 		apply_eigen(init,atom,evec,mode+j-1,amp[j]);
		apply_eigen(all_init,all,evec,mode+j-1,amp[j]);
		newstrc = rmsd_no(init,targ,atom,align);
		energy += 0.5*pow(gsl_vector_get(eval,j+mode-1),2)*pow(amp[j],2);
		printf("J:%4d %10.6f %10.6f %10.6f %12.6f\n",j+mode,amp[j],old,sqrt(newstrc),energy);
	}*/
	
	
 	gsl_vector_free(diff);
}

int main(int argc, char *argv[]) {

	int all; /*Nombre d'atomes dans pdb*/
	int atom; /*Nombre de carbone CA*/
	int all_t; /*Nombre d'atomes dans pdb*/
	int atom_t; /*Nombre de carbone CA*/
 	int help_flag = 1;
 	int i,j;
 	char file_name[500];
 	int lig = 0;
 	char out_name[500] = "templaate.dat";
 	char matrix_name[500] = "interaction_m.dat";
 	int verbose = 0;
	int super_node = 1;
	int print_flag = 0;
	int nconn;
	float vinit = 0.0001; // Valeur de base 
	float epsilon = 0.36;					// Epsilon			
	float bond_factor = 100*epsilon;		// Facteur pour poid des bond strechcing
	float angle_factor = 20*epsilon;		// Facteur pour poid des angles
	double K_phi1 = epsilon;				// Facteurs pour angles dièdres
	double K_phi3 = 0.5*epsilon;
	long seed; seed = time_seed();			// Random
	float kp_factor = 1;					// acteur pour poid des angles dièdres
	int it = 0;
	float init_templaate = 1;
	char check_name[500];
	int encom_flag = 0;
	float cutoff = 18.00;
	int number_mode = 2;
 	for (i = 1;i < argc;i++) {
 		if (strcmp("-i",argv[i]) == 0) {strcpy(file_name,argv[i+1]);--help_flag;}
 		if (strcmp("-o",argv[i]) == 0) {strcpy(out_name,argv[i+1]);}
 		if (strcmp("-h",argv[i]) == 0) {++help_flag;}
 		if (strcmp("-v",argv[i]) == 0) {verbose = 1;}
 		if (strcmp("-p",argv[i]) == 0) {print_flag = 1;}
 		if (strcmp("-sp",argv[i]) == 0) {super_node = 3;}
 		if (strcmp("-spc",argv[i]) == 0) {super_node = 4;}
 		if (strcmp("-lig",argv[i]) == 0) {lig= 1;}
 		if (strcmp("-m",argv[i]) == 0) {strcpy(matrix_name,argv[i+1]);}
 		if (strcmp("-init",argv[i]) == 0) {float temp;sscanf(argv[i+1],"%f",&temp);vinit = temp;}
 		if (strcmp("-t",argv[i]) == 0) {float temp;sscanf(argv[i+1],"%f",&temp);init_templaate = temp;}
 		if (strcmp("-kr",argv[i]) == 0) {float temp;sscanf(argv[i+1],"%f",&temp);bond_factor = temp * bond_factor;}
 		if (strcmp("-kt",argv[i]) == 0) {float temp;sscanf(argv[i+1],"%f",&temp);angle_factor = temp * angle_factor;}
 		if (strcmp("-kpf",argv[i]) == 0) {float temp;sscanf(argv[i+1],"%f",&temp); kp_factor = temp;}
 		if (strcmp("-nm",argv[i]) == 0) {int temp;sscanf(argv[i+1],"%d",&temp); number_mode = temp;}
 		if (strcmp("-target",argv[i]) == 0) {strcpy(check_name,argv[i+1]);help_flag = 0;}    
 	}
 	
 	//***************************************************
 	//*													*
 	//*Build a structure contaning information on the pdb
 	//*													*
 	//***************************************************
 	
 	all = count_atom(file_name);
 	nconn = count_connect(file_name);
 	
 	if (verbose == 1) {printf("Connect:%d\n",nconn);}
 	
	if (verbose == 1) {printf("Assigning Structure\n\tAll Atom\n");}
	
	// Array qui comprend tous les connects
	
	int **connect_h=(int **)malloc(nconn*sizeof(int *)); 
    for(i=0;i<nconn;i++) { connect_h[i]=(int *)malloc(7*sizeof(int));}
    
    assign_connect(file_name,connect_h);
	
	// Assign tous les atoms
	
	struct pdb_atom strc_all[all*2];
	atom = build_all_strc(file_name,strc_all); // Retourne le nombre de Node
	
	if (verbose == 1) {printf("	Node:%d\n	Atom:%d\n",atom,all);}

	check_lig(strc_all,connect_h,nconn,all);

	// Assign les Nodes
	
	if (verbose == 1) {printf("	CA Structure\n");}
	
	//atom = count_atom_CA_n(strc_all,all,super_node,lig);
	if (verbose == 1) {printf("	Node:%d\n",atom);}
	struct pdb_atom strc_node[atom*2];
	atom = build_cord_CA(strc_all, strc_node,all,lig,connect_h,nconn);
	
	if (verbose == 1) {printf("	Assign Node:%d\n",atom);}
	if (atom > 800) {printf("To much node !!! I QUIT BECAUSE VINCE TOLD ME\n");return(0);}
	//write_strc("structure.pdb",strc_all,all);
 	
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
 	if (score/atom < 0.8) {
 		printf("Low Score... Will try an homemade alignement !!!\n");
 		score = node_align_low(strc_node,atom,strc_node_t,atom_t,align);
 		printf("RMSD:%8.5f Score: %d/%d\n",sqrt(rmsd_no(strc_node,strc_node_t,atom, align)),score,atom);
 	}
 	
 	
 	gsl_vector *eval = gsl_vector_alloc(3*atom); /*Déclare un vector qui va contenir les eigenvalue */	
	gsl_matrix *evec= gsl_matrix_alloc (3*atom,3*atom); /*Déclare une matrice qui va contenir les eigenvector correspondant à eval*/
 	
 	 double **hessian=(double **)malloc(3*atom*sizeof(double *)); // Matrix of the Hessian 1 2 3 (bond, angle, dihedral)
  	 for(i=0;i<3*atom;i++) { hessian[i]=(double *)malloc(3*atom*sizeof(double));}
  	 for(i=0;i<3*atom;i++)for(j=0;j<(3*atom);j++){hessian[i][j]=0;}
  	      
	 gsl_matrix *h_matrix = gsl_matrix_alloc(3*atom, 3*atom);	/*Déclare une matrice hessian matrix de 3n par 3n*/
	 gsl_matrix *templaate = gsl_matrix_alloc(3*atom, 3*atom);	/*Déclare une matrice hessian matrix de 3n par 3n*/
	 if (number_mode == -1) {number_mode = 3*atom -10;}
 	gsl_vector *answer = gsl_vector_alloc(number_mode);
 	FILE *file;
 	file = fopen("movie.pdb","w");
 	for (it = 0; it < 20;++it) {
 		//***************************************************
 		//*													*
 		//*Build Hessian									*
 		//*													*
 		//***************************************************

 		for(i=0;i<3*atom;i++) {for(j=0;j<3*atom;j++) {hessian[i][j] = 0;}}
   	 	
		
		gsl_matrix_set_all(templaate,1);
		//printf("Build Matrix\n");
		build_enm(strc_node,hessian,atom,epsilon,templaate,cutoff,0);
		assignArray(h_matrix,hessian,3*atom,3*atom);
		//printf("Diago\n");
		gsl_vector *eval = gsl_vector_alloc(3*atom); /*Déclare un vector qui va contenir les eigenvalue */	
		gsl_matrix *evec= gsl_matrix_alloc (3*atom,3*atom); /*Déclare une matrice qui va contenir les eigenvector correspondant à eval*/
		diagonalyse_matrix(h_matrix,(3)*atom,eval,evec);
		//write_matrix("hessian.dat",h_matrix,3*atom,3*atom);
		write_eigen("eigen.dat",evec,eval,3*atom);
 		
 		fit_svd_fold(strc_node,strc_node_t,atom,all, atom_t, all_t,strc_all,strc_all_t, evec, align, number_mode, 7,eval,answer);
 		float energ = 0;
 		for (j=0;j<number_mode;++j) {
 			energ += gsl_vector_get(answer,j)*gsl_vector_get(answer,j)*gsl_vector_get(eval,j+7-1)*gsl_vector_get(eval,j+7-1);
 			//printf("ENERG:%f += %f * %f\n",energ,gsl_vector_get(answer,j),gsl_vector_get(eval,j+7-1));
 		}
 		
 		float factor = sqrt(energ/100);
 		printf("FACTOR:%f\n",factor);
 		if (factor < 1) {factor = 1;}
 		energ = 0;
 		for (j=0;j<number_mode;++j) {
 			energ += gsl_vector_get(answer,j)*gsl_vector_get(answer,j)*gsl_vector_get(eval,j+7-1)*gsl_vector_get(eval,j+7-1)/factor/factor;
 		//	printf("ENERG:%f += %f * %f / %f\n",energ,gsl_vector_get(answer,j),gsl_vector_get(eval,j+7-1));
 		}
 		
 		for (i=0;i<20;++i) {
 			
 			
 			for (j=0;j<number_mode;++j) {
 				if (gsl_vector_get(answer,j) < 0.001) {continue;}
				apply_eigen(strc_node,atom,evec,7+j-1,gsl_vector_get(answer,j)/20/factor);
				apply_eigen(strc_all,all,evec,7+j-1,gsl_vector_get(answer,j)/20/factor);
				//write_movie(file,strc_all ,all,(it)*20+i*10+j+1);
			}
			write_movie(file,strc_all ,all,(it)*20+i+1);
		}
		
 		printf("IT: %4d RMSD:%8.5f Score: %d/%d Energ:%2f\n",it,sqrt(rmsd_no(strc_node,strc_node_t,atom, align)),score,atom,energ);
 		gsl_vector_free(eval);
 		gsl_matrix_free(evec);
 	
 	}
 	
 	
 	
 	
 	
 	
 	
 	
 	
 	
 	
 	
 	return(1);
}
