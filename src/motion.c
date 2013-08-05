#include "STeM.h"





int main(int argc, char *argv[]) {

	int all; /*Nombre d'atomes dans pdb*/
	int atom; /*Nombre de carbone CA*/
 	int help_flag = 1;
 	char file_name[100];
 	char eigen_name[100] = "eigen.dat";
 	char out_name[100] = "motion.pdb";
 	int verbose = 0;
	int mode = 7;
	int i;
	float max_dist = 1.05;
	float max_val = 1,min_val = -1;
	int nconn;
	int lig = 0;
	int torsion = 0;
 	for (i = 1;i < argc;i++) {
 		if (strcmp("-i",argv[i]) == 0) {strcpy(file_name,argv[i+1]);--help_flag;}
 		if (strcmp("-o",argv[i]) == 0) {strcpy(out_name,argv[i+1]);}
 		if (strcmp("-ieig",argv[i]) == 0) {strcpy(eigen_name,argv[i+1]);}
 		if (strcmp("-h",argv[i]) == 0) {help_flag = 1;}
 		if (strcmp("-v",argv[i]) == 0) {verbose = 1;}
 		if (strcmp("-md",argv[i]) == 0) {float temp;sscanf(argv[i+1],"%f",&temp);max_dist = temp;}
 		if (strcmp("-m",argv[i]) == 0) {int temp;sscanf(argv[i+1],"%d",&temp);mode = temp;}

 		if (strcmp("-lig",argv[i]) == 0) {lig= 1;}
 		if (strcmp("-angle",argv[i]) == 0) {torsion = 1;}     
 	}
 	
 	if (help_flag == 0) { } else {
 		printf("****************************\nHelp Section\n-i\tFile Input (PDB)\n-o\tOutput Name Motion\n-ieig\tFile Name Eigen\n-v\tVerbose\n-sp\tSuper Node Mode (CA, N, C)\n-m\tMode\n-nm\tNombre de mode\n-lig\tTient compte duligand (sauf HOH)\n****************************\n");
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
	
	
	
	//***************************************************
 	//*													*
 	//* Load eigenvector et eigenvalue					*
 	//*													*
 	//***************************************************
 	
 	if (verbose == 1) {printf("Loading Eigenvector\n");}
 	
 	gsl_vector *eval = gsl_vector_alloc((3-torsion)*atom+10);
	gsl_matrix *evec= gsl_matrix_alloc ((3-torsion)*atom+10,(3-torsion)*atom+10);
	load_eigen_motion(eval,evec,eigen_name,(3-torsion)*atom,mode);
	
	//***************************************************
 	//*													*
 	//* Print Motion									*
 	//*													*
 	//***************************************************
 	if (verbose == 1) {printf("Creating Motion\n");}
 	if (verbose == 1) {printf("\tFinding Max and Min movement\n");}

	if (torsion == 0) {
	 	find_max_ampli(strc_node,evec,mode,atom,max_dist,&max_val,&min_val);
	 	
	 	max_val = max_dist;
	 	min_val = -max_dist;
	 	printf("Max:%f\tMin:%f\n",max_val,min_val);
	 	if (verbose == 1) {printf("Printing Motion\n");}
	 	print_image_motion(strc_all,evec,mode,max_val,min_val,all,out_name,lig);
	 } else {
	 	print_image_torsion(strc_all,evec,mode,max_val,min_val,all,out_name,atom,strc_node);
	 
	 }
 	

	gsl_vector_free(eval);
	gsl_matrix_free(evec);
	
}
