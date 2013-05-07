#include "STeM.h"

int main(int argc, char *argv[]) {
	int all; /*Nombre d'atomes dans pdb*/
	int atom; /*Nombre de carbone CA*/
	int all_t; /*Nombre d'atomes dans pdb*/
	int atom_t; /*Nombre de carbone CA*/
	int i,j;
	int help_flag = 1;
 	char file_name[500];
 	char check_name[500];
 	int lig = 0;
 	char matrix_name[500] = "interaction_m.dat";
 	char param_name[500] = "param_m.dat";
 	int verbose = 0;
	int nconn;
	float vinit = 0.0001; // Valeur de base
	float epsilon = 0.36;					// Epsilon			
	float bond_factor = 100*epsilon;		// Facteur pour poid des bond strechcing
	float angle_factor = 20*epsilon;		// Facteur pour poid des angles
	double K_phi1 = epsilon;				// Facteurs pour angles dièdres
	double K_phi3 = 0.5*epsilon;
	float init_templaate = 1;
	float kp_factor = 1;					// acteur pour poid des angles dièdres
	int weight_factor = 0;
	int nbr_mode = 2;
	int mode = 7;
	float cutoff = 18.0;
 	for (i = 1;i < argc;i++) {
 		if (strcmp("-i",argv[i]) == 0) {strcpy(file_name,argv[i+1]);--help_flag;}
 		if (strcmp("-h",argv[i]) == 0) {++help_flag;}
 		if (strcmp("-v",argv[i]) == 0) {verbose = 1;}
 		if (strcmp("-lig",argv[i]) == 0) {lig= 1;}
 		if (strcmp("-targ",argv[i]) == 0) {strcpy(check_name,argv[i+1]);help_flag = 0;}
 		if (strcmp("-mode",argv[i]) == 0) {int temp;sscanf(argv[i+1],"%d",&temp);mode = temp;}
 		if (strcmp("-nm",argv[i]) == 0) {int temp;sscanf(argv[i+1],"%d",&temp);nbr_mode = temp;}
 		
 		
 		
 		
 	}
 	
 	//***************************************************
 	//*													*
 	//*Build a structure contaning information on the pdb
 	//*													*
 	//***************************************************
        if (verbose == 1) {printf("First file:%s\n",file_name);} 	
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
	
	if (verbose == 1) {printf("	Node:%d\n",atom);}
	struct pdb_atom strc_node[atom];
	atom = build_cord_CA(strc_all, strc_node,all,lig,connect_h,nconn);
	if (verbose == 1) {printf("	Assign Node:%d\n",atom);}
		
	// Pour les besoins... on limite à 800 atomes
	
	if (atom > 2000) {
		printf("Too much atom, if you want to remove the limite.... vincent.frappier@usherbrooke.ca\n");
		return(1);
	}
	if (verbose == 1) {printf("Build Second Structure\n");}
	//Construit la structure a comparer
 	nconn = 0;
 	if (verbose == 1) {printf("Sec file:%s\n",check_name);}
 	all_t = count_atom(check_name);
 	nconn = count_connect(check_name);
 	
 	if (verbose == 1) {printf("Connect:%d\n",nconn);}
 	
	if (verbose == 1) {printf("Assigning Structure\n\tAll Atom\n");}
	
	// Array qui comprend tous les connects
	
	int **connect_t=(int **)malloc(nconn*sizeof(int *)); 
    for(i=0;i<nconn;i++) { connect_t[i]=(int *)malloc(7*sizeof(int));}
    
    assign_connect(check_name,connect_t);
	
	// Assign tous les atoms
	
	struct pdb_atom strc_all_t[all_t];
	atom_t = build_all_strc(check_name,strc_all_t); // Retourne le nombre de Node
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
 	if (verbose == 1) {printf("Score: %d/%d\n",score,atom);}
	//if (atom_t != atom) {printf("Not the same number of Node... Terminating\n");return(0);}
	
	if ( (float) score/ (float) atom < 0.8) {
 		printf("Low Score... Will try an homemade alignement !!!\n");
 		score = node_align_low(strc_node,atom,strc_node_t,atom_t,align);
 		
 	}
	printf("RMSD:%8.5f Score: %d/%d\n",sqrt(rmsd_no(strc_node,strc_node_t,atom, align)),score,atom);
	
	//***************************************************
 	//*													*
 	//*Build Hessian									*
 	//*													*
 	//***************************************************


    double **hessian=(double **)malloc(3*atom*sizeof(double *)); // Matrix of the Hessian 1 2 3 (bond, angle, dihedral)
    for(i=0;i<3*atom;i++) { hessian[i]=(double *)malloc(3*atom*sizeof(double));}
    for(i=0;i<3*atom;i++)for(j=0;j<(3*atom);j++){hessian[i][j]=0;}
        
	gsl_matrix *h_matrix = gsl_matrix_alloc(3*atom, 3*atom);	/*Déclare une matrice hessian matrix de 3n par 3n*/
 	
 	//***************************************************
 	//*													*
 	//*Build template									*
 	//*													*
 	//***************************************************
	
	if (verbose == 1) {printf("Building Hessian\n");}
	
	if (verbose == 1) {printf("\tCovalent Bond Potential\n");}		
	
	if (verbose == 1) {printf("\tAngle Potential\n");}	
	
	if (verbose == 1) {printf("\tDihedral Potential\n");}	
	
	if (verbose == 1) {printf("\tNon Local Interaction Potential\n");}	
gsl_matrix *templaate = gsl_matrix_alloc(3*atom, 3*atom);
gsl_matrix_set_all(templaate,1);
	build_enm(strc_node,hessian,atom,epsilon,templaate,cutoff);	
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
	if (verbose == 1) {
		printf("First eigenvalues\n");
		for (i=0;i<10;++i) {printf("I:%d %.10f\n",i,gsl_vector_get(eval,i));}
	}
 	
 	
 	fit_svd(strc_node,strc_node_t,atom,all,atom_t,all_t,strc_all,strc_all_t,evec,align,nbr_mode,mode,eval);
 	
 	gsl_matrix_free(templaate);
 	
 	gsl_matrix_free(h_matrix);
	gsl_vector_free(eval);
	gsl_matrix_free(evec);
 	
 	free(connect_h);
	free(connect_t);
 	
 	
 	
 	
 	return(1);
 	
}
