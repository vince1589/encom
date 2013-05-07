#include "STeM.h"

int main(int argc, char *argv[]) {

	int all; 								// Nombre d'atoms dans le PDB
	int atom; 								// Nombre de nodes
	float epsilon = 0.36;					// Epsilon			
	float bond_factor = 100*epsilon;		// Facteur pour poid des bond strechcing
	float angle_factor = 20*epsilon;		// Facteur pour poid des angles
	double K_phi1 = epsilon;				// Facteurs pour angles dièdres
	double K_phi3 = 0.5*epsilon;
	long seed; seed = time_seed();			// Random
 	int help_flag = 1;						// Help Flag
	int i,j;								// Dummy Counter
 	char file_name[100];					// Nom du PDB
 	char templaate_name[100];				// Nom du fichier du templaate
 	char eigen_name[100] = "eigen.dat";		// Nom pour output des eigen	
 	float init_templaate = 1;
 	int verbose = 0;						// Loud
 	int templaate_flag = 0;					// Flag pour savoir si tempalte ou non
	int b_factor_flag = 0;					// Output b_factor flag
	int lig = 0;							// Flag pour savoir si tient compte ligand
	int nconn;								// Nombre de connect
	float kp_factor = 1;					// acteur pour poid des angles dièdres
	int enm = 0;							//flag pour enm instead of STeM
	float tot_factor = 1;
	int energy = 0;							// Flag pour printer l'energy
	float cutoff = 18.0;
	int weight_factor = 0;
 	for (i = 1;i < argc;i++) {
 		if (strcmp("-i",argv[i]) == 0) {strcpy(file_name,argv[i+1]);--help_flag;}
 		if (strcmp("-o",argv[i]) == 0) {strcpy(eigen_name,argv[i+1]);}
 		if (strcmp("-h",argv[i]) == 0) {help_flag = 1;}
 		if (strcmp("-w",argv[i]) == 0) {weight_factor = 1;}
 		if (strcmp("-v",argv[i]) == 0) {verbose = 1;}
 		if (strcmp("-t",argv[i]) == 0) {float temp;sscanf(argv[i+1],"%f",&temp);init_templaate = temp;}  
 		if (strcmp("-lt",argv[i]) == 0) {strcpy(templaate_name,argv[i+1]);templaate_flag = -1;}
 		if (strcmp("-kr",argv[i]) == 0) {float temp;sscanf(argv[i+1],"%f",&temp);bond_factor = temp * bond_factor;}
 		if (strcmp("-kt",argv[i]) == 0) {float temp;sscanf(argv[i+1],"%f",&temp);angle_factor = temp * angle_factor;}
 		if (strcmp("-kpf",argv[i]) == 0) {float temp;sscanf(argv[i+1],"%f",&temp); kp_factor = temp;}  
 		if (strcmp("-b",argv[i]) == 0) {b_factor_flag = 1;}
 		if (strcmp("-lig",argv[i]) == 0) {lig= 1;} 
 		if (strcmp("-enm",argv[i]) == 0) {enm= 1;}  
 		if (strcmp("-ee",argv[i]) == 0) {energy= 1;}  
 		if (strcmp("-tot",argv[i]) == 0) {float temp;sscanf(argv[i+1],"%f",&temp); tot_factor = temp;}
 		if (strcmp("-cut",argv[i]) == 0) {float temp;sscanf(argv[i+1],"%f",&temp); cutoff = temp;}  
 		
 	}

 	printf("\n********************\nFile:%s\nEpsilon:%f\nKr:\t%f\nKt:\t%f\nK phi1:\t%f\nK phi3:\t%f\n********************\n",file_name,epsilon,bond_factor,angle_factor,K_phi1,K_phi3);

 	if (help_flag == 0) { } else {
 		printf("****************************\nHelp Section\n-i\tFile Input (PDB)\n-oh\tOutput Name Hessian\n-oeig\tFile Name Eigen\n-v\tVerbose\n-w\tWeight Vector\n-t\tInitial Value of templaate (negative value for random)\n\tIf Load templaate, multiply the templaate\n-lt\tLoad input templaate\n-sp\tSuper Node Mode (CA, C et N)\n-kt\tPoid de l'angle entre les nodes (1)\n-kr\tPoid de la distance entre les nodes (1)\n-lig\tTient compte duligand (sauf HOH)\n-kp\tPoid de l'angles dièdres\n-kpf\tFacteur entre \n****************************\n");
 		return(0); 
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
	
	struct pdb_atom strc_all[all];
	atom = build_all_strc(file_name,strc_all); // Retourne le nombre de Node
	
	if (verbose == 1) {printf("	Node:%d\n	Atom:%d\n",atom,all);}

	check_lig(strc_all,connect_h,nconn,all);

	// Assign les Nodes
	
	if (verbose == 1) {printf("	CA Structure\n");}
	
	//atom = count_atom_CA_n(strc_all,all,super_node,lig);
	if (verbose == 1) {printf("	Node:%d\n",atom);}
	struct pdb_atom strc_node[atom];
	atom = build_cord_CA(strc_all, strc_node,all,lig,connect_h,nconn);
	if (verbose == 1) {printf("	Assign Node:%d\n",atom);}
	
	//write_strc("structure.pdb",strc_all,all);
	
	// Pour les besoins... on limite à 800 atomes
	
	if (atom > 800) {
		printf("Too much atom, if you want to remove the limite.... vincent.frappier@usherbrooke.ca\n");
		return(1);
	}
	
	//***************************************************
 	//*													*
 	//*Build Hessian									*
 	//*													*
 	//***************************************************


    double **hessian=(double **)malloc(3*atom*sizeof(double *)); // Matrix of the Hessian 1 2 3 (bond, angle, dihedral)
    for(i=0;i<3*atom;i++) { hessian[i]=(double *)malloc(3*atom*sizeof(double));}
    for(i=0;i<3*atom;i++)for(j=0;j<(3*atom);j++){hessian[i][j]=0;}
        
	gsl_matrix *h_matrix = gsl_matrix_alloc(3*atom, 3*atom);	/*Déclare une matrice hessian matrix de 3n par 3n*/
	gsl_matrix *templaate = gsl_matrix_alloc(3*atom, 3*atom);	/*Déclare une matrice hessian matrix de 3n par 3n*/
	
	//Setting templaate
	
	if (templaate_flag == -1) {
		printf("Loading templaate... Help flag:%d\n",help_flag);
		load_matrix(templaate,templaate_name);
		//for (i=0;i<atom;++i) {printf("%10.5f ",gsl_matrix_get(templaate,0,i));}printf("\n");
		gsl_matrix_scale (templaate, init_templaate);
	} else {
		if (init_templaate >= 0) {
			gsl_matrix_set_all (templaate, init_templaate/tot_factor);
		} else {
			for (i=0;i<atom;++i) {
				for (j=0;j<atom;++j) {
					float random =ran2(&seed);
					if (abs(i-j) < 4) {random = 0;}
					gsl_matrix_set(templaate,i,j,random);
					gsl_matrix_set(templaate,j,i,random);
				}
			} // Set random number if negatif number
		}
	}
	//write_matrix("out_templaate", templaate,atom);
	//Building Matrix
		
	if (verbose == 1) {printf("Building Hessian\n");}
	
	if (verbose == 1) {printf("\tCovalent Bond Potential\n");}		
	build_1st_matrix(strc_node,hessian,atom,bond_factor/tot_factor);
	
	if (verbose == 1) {printf("\tAngle Potential\n");}	
	build_2_matrix(strc_node,hessian,atom,angle_factor/tot_factor);
	
	if (verbose == 1) {printf("\tDihedral Potential\n");}	
	build_3_matrix(strc_node, hessian,atom,K_phi1/2+K_phi3*9/2,kp_factor/tot_factor);
	
	if (verbose == 1) {printf("\tNon Local Interaction Potential\n");}	
	build_4h_matrix(strc_node,hessian,atom,epsilon,templaate);
	
	if (enm == 1) {
		 for(i=0;i<3*atom;i++)for(j=0;j<(3*atom);j++){hessian[i][j]=0;}
		 build_enm(strc_node,hessian,atom,epsilon,templaate,cutoff);
	}
	
	if (verbose == 1) {printf("\tAssigning Array\n");}	
	assignArray(h_matrix,hessian,3*atom,3*atom);
	if (weight_factor == 1) {mass_weight_hessian(h_matrix,atom,strc_node);}
//	write_matrix("hessian.dat",h_matrix,3*atom,3*atom);
	//***************************************************
	//*													*
	//*Diagonalyse the matrix							*
	//*													*
	//***************************************************
	

	if (verbose == 1) {printf("Diagonalizing Hessian\n");}
	gsl_vector *eval = gsl_vector_alloc(3*atom); /*Déclare un vector qui va contenir les eigenvalue */	
	gsl_matrix *evec= gsl_matrix_alloc (3*atom,3*atom); /*Déclare une matrice qui va contenir les eigenvector correspondant à eval*/
	diagonalyse_matrix(h_matrix,3*atom,eval,evec); /*Diagonalyse la matrice*/
	if (weight_factor == 1) {adjust_weight_evec(evec,atom,strc_node);}
	write_eigen(eigen_name,evec,eval,3*atom);
	if (energy == 1) {printf("Energy:%10.5f\n",calc_energy(atom,eval));}
	
	//***************************************************
	//*													*
	//*Inverse matrix									*
	//*													*
	//***************************************************
	 if (b_factor_flag == 1) {
		 if (verbose == 1) {printf("Inversing Matrix\n");}
		gsl_matrix *k_inverse = gsl_matrix_alloc(atom, atom); /*Déclare et crée une matrice qui va être le pseudo inverse*/
		k_inverse_matrix_stem(k_inverse,atom,eval,evec);
	
		//write_matrix("b_factor.dat", k_inverse,atom,atom);

		printf("Correlation:%f\n",correlate(k_inverse,strc_node, atom));
		/*for (i=0;i<all;++i) {
			strc_all[i].b_factor = gsl_matrix_get(k_inverse,strc_all[i].node,strc_all[i].node);
		}*/
		//write_strc("b_factor.pdb",strc_all,all);
		//if (super_node != 3) {write_strc_b("b_factor.pdb",strc_all,all,k_inverse,super_node);}
		gsl_matrix_free(k_inverse);
		
	}
	
	// Freeing
	
	gsl_matrix_free(h_matrix);
	gsl_matrix_free(templaate);
	gsl_vector_free(eval);
	gsl_matrix_free(evec);
	
	for(i=0;i<3*atom;i++){free(hessian[i]);}
	free(hessian);
	for(i=0;i<nconn;i++) {free(connect_h[i]);}
	free(connect_h);
}


