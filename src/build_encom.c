#include "STeM.h"

// Fonciton qui construit les eigen

int main(int argc, char *argv[]) {
	int all; /*Nombre d'atomes dans pdb*/
	int atom; /*Nombre de carbone CA*/

	int i,j;
	int help_flag = 1;
 	char file_name[500];

 	int lig = 0;
 	char matrix_name[500] = "interaction_m.dat";
 	char param_name[500] = "param_m.dat";
 	char inputname[500] ="none";
 	int verbose = 0;
	int nconn;
	float vinit = 4.85489969815978; // Valeur de base
	float epsilon = 0.36;					// Epsilon			
	float bond_factor = 100*epsilon*459.32638757416;		// Facteur pour poid des bond strechcing
	float angle_factor = 20*epsilon*89.0845206090432;		// Facteur pour poid des angles
	double K_phi1 = epsilon;				// Facteurs pour angles dièdres
	double K_phi3 = 0.5*epsilon;
	float init_templaate = 11.1469110430119;
	float kp_factor = 41644.7785340819;					// Facteur pour poid des angles dièdres
	int weight_factor = 0;
	char eigen_name[500] = "eigen.dat";	
	float temperature = 310;
 	for (i = 1;i < argc;i++) {
 		if (strcmp("-i",argv[i]) == 0) {strcpy(file_name,argv[i+1]);--help_flag;}
 		if (strcmp("-inp",argv[i]) == 0) {strcpy(inputname,argv[i+1]);--help_flag;}
 		if (strcmp("-h",argv[i]) == 0) {++help_flag;}
 		if (strcmp("-v",argv[i]) == 0) {verbose = 1;}
 		if (strcmp("-lig",argv[i]) == 0) {lig= 1;}
 		if (strcmp("-m",argv[i]) == 0) {strcpy(matrix_name,argv[i+1]);}
 		if (strcmp("-param",argv[i]) == 0) {strcpy(param_name,argv[i+1]);}
 		if (strcmp("-init",argv[i]) == 0) {float temp;sscanf(argv[i+1],"%f",&temp);vinit = temp;}
 		if (strcmp("-t",argv[i]) == 0) {float temp;sscanf(argv[i+1],"%f",&temp);init_templaate = temp;} 
 		if (strcmp("-kr",argv[i]) == 0) {float temp;sscanf(argv[i+1],"%f",&temp);bond_factor = 100*epsilon * bond_factor;}
 		if (strcmp("-kt",argv[i]) == 0) {float temp;sscanf(argv[i+1],"%f",&temp);angle_factor = 20*epsilon * angle_factor;}
 		if (strcmp("-kpf",argv[i]) == 0) {float temp;sscanf(argv[i+1],"%f",&temp); kp_factor = temp;}
 		if (strcmp("-temp",argv[i]) == 0) {float temp;sscanf(argv[i+1],"%f",&temp); temperature = temp;}
 		if (strcmp("-w",argv[i]) == 0) {weight_factor = 1;}
 		if (strcmp("-o",argv[i]) == 0) {strcpy(eigen_name,argv[i+1]);}
 	}
 	
 	printf("\n********************\nFile:%s\nEpsilon:%f\nKr:\t%f\nKt:\t%f\nK phi1:\t%f\nK phi3:\t%f\n********************\n",file_name,epsilon,bond_factor,angle_factor,K_phi1,K_phi3);
 	
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
	
	if (verbose == 1) {printf("	Node:%d\n",atom);}
	struct pdb_atom strc_node[atom];
	atom = build_cord_CA(strc_all, strc_node,all,lig,connect_h,nconn);
	if (verbose == 1) {printf("	Assign Node:%d\n",atom);}
		
	// Pour les besoins... on limite à 800 atomes
	
	if (atom > 2000) {
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
 	
 	//***************************************************
 	//*													*
 	//*Build template									*
 	//*													*
 	//***************************************************
	
	assign_atom_type(strc_all, all);
	if (strcmp(inputname,"none") == 0) {} else {assign_lig_type(strc_all, all, inputname);}
	gsl_matrix *vcon = gsl_matrix_alloc(all,all);
	gsl_matrix *inter_m = gsl_matrix_alloc(8,8);
	gsl_matrix *templaate = gsl_matrix_alloc(atom*3, atom*3);
	gsl_matrix_set_all(templaate,vinit);
	gsl_matrix_set_all(vcon,0);
	
	if (verbose == 1) {printf("Do Vcon !!!\n");}
	
	vcon_file_dom(strc_all,vcon,all);
	
	
	
	if (verbose == 1) {printf("Reading Interaction Matrix %s\n",matrix_name);}
	load_matrix(inter_m,matrix_name);
	//write_matrix("vcon_vince.dat", vcon,all,all);
	if (verbose == 1) {printf("Building templaate\n");}
	all_interaction(strc_all,all, atom, templaate,lig,vcon,inter_m,strc_node);
	gsl_matrix_scale (templaate, init_templaate);
	
	//write_matrix("template_diff.dat", templaate,atom,atom);
	
	if (verbose == 1) {printf("Building Hessian\n");}
	
	if (verbose == 1) {printf("\tCovalent Bond Potential\n");}		
	build_1st_matrix(strc_node,hessian,atom,bond_factor);
	
	if (verbose == 1) {printf("\tAngle Potential\n");}	
	build_2_matrix(strc_node,hessian,atom,angle_factor);
	
	if (verbose == 1) {printf("\tDihedral Potential\n");}	
	build_3_matrix(strc_node, hessian,atom,K_phi1/2+K_phi3*9/2,kp_factor);
	
	if (verbose == 1) {printf("\tNon Local Interaction Potential\n");}	
	build_4h_matrix(strc_node,hessian,atom,epsilon,templaate);
	
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
 	write_eigen(eigen_name,evec,eval,3*atom);


		 	
 	float ener = calc_energy(atom,eval,temperature);
 	printf("Energy:%10.8f\n",ener);
 	printf("Energy/node:%10.8f\n",ener/(float)(atom*3));

 	gsl_matrix_free(templaate);
	gsl_matrix_free(inter_m);
	gsl_matrix_free(vcon);
 	
 	gsl_matrix_free(h_matrix);
	gsl_vector_free(eval);
	gsl_matrix_free(evec);
 	
 	free(connect_h);

 	
 	
 	
 	
 	return(1);
 	
}

