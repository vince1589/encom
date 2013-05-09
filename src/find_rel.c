#include "STeM.h"

int main(int argc, char *argv[]) {

	int all; /*Nombre d'atomes dans pdb*/
	int atom; /*Nombre de carbone CA*/

 	int help_flag = 1;
 	char file_name[500];
 	char eigen_name[100] = "eigen.dat";
 	int verbose = 0;
	int i,j,k;
	int nb_mode = 2;
	int mode = 6;
	int lig = 1;
	int nconn;

	float ligalign = 5; // Flag/valeur pour aligner seulement les résidus dans un cutoff du ligand, 0, one le fait pas... > 0... le cutoff
 	for (i = 1;i < argc;i++) {
 		if (strcmp("-i",argv[i]) == 0) {strcpy(file_name,argv[i+1]);--help_flag;}
 		if (strcmp("-h",argv[i]) == 0) {help_flag = 1;}
 		if (strcmp("-v",argv[i]) == 0) {verbose = 1;}
 		if (strcmp("-lig",argv[i]) == 0) {lig= 1;}  		
 		if (strcmp("-m",argv[i]) == 0) {int temp;sscanf(argv[i+1],"%d",&temp);mode = temp;}
 		if (strcmp("-nm",argv[i]) == 0) {int temp;sscanf(argv[i+1],"%d",&temp);nb_mode = temp;}
 		if (strcmp("-ieig",argv[i]) == 0) {strcpy(eigen_name,argv[i+1]);}
 		if (strcmp("-ligc",argv[i]) == 0) {float temp;sscanf(argv[i+1],"%f",&temp);ligalign = temp;}
 		
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
 	if (verbose == 1) {printf("First file:%s\n",file_name);}
 	if (verbose == 1) {printf("Connect:%d\n",nconn);}
 	
	if (verbose == 1) {printf("Assigning Structure\n\tAll Atom\n");}
	
	// Array qui comprend tous les connects
	
	int **connect_h=(int **)malloc(nconn*sizeof(int *)); 
    for(i=0;i<nconn;i++) { connect_h[i]=(int *)malloc(7*sizeof(int));}
   
    assign_connect(file_name,connect_h);
	 printf("HERE\n");
	// Assign tous les atoms
	
	struct pdb_atom strc_all[all];
	atom = build_all_strc(file_name,strc_all); // Retourne le nombre de Node
	
	if (verbose == 1) {printf("	Atom:%d\n",all);}
	check_lig(strc_all,connect_h,nconn,all);
	
	// Assign les Nodes
	
	if (verbose == 1) {printf("	CA Structure\n");}
	
	if (verbose == 1) {printf("	Node:%d\n",atom);}
	struct pdb_atom strc_node[atom];
	atom = build_cord_CA(strc_all, strc_node,all,lig,connect_h,nconn);
	
	if (verbose == 1) {printf("	Assign Node:%d\n",atom);}

	

	struct pdb_atom strc_node_t[atom];
	copy_strc(strc_node_t, strc_node, atom);
	
 
 	int align[atom];
 	int score = node_align(strc_node,atom,strc_node_t,atom,align);
 	
 	if (verbose == 1) {printf("Score: %d/%d\n",score,atom);}
		
	if (score/atom < 0.8) {
 		printf("Low Score... Will try an homemade alignement !!!\n");
 		score = node_align_low(strc_node,atom,strc_node_t,atom,align);
 		
 	}
	printf("RMSD:%8.5f Score: %d/%d\n",sqrt(rmsd_no(strc_node,strc_node_t,atom, align)),score,atom);	
	if (ligalign > 0) {

		score = node_align_lig(strc_node,atom,strc_node_t,atom,align,strc_all,all,strc_all,all,ligalign);
		
		printf("RMSD:%8.5f Score: %d/%d\n",sqrt(rmsd_no(strc_node,strc_node_t,atom, align)),score,atom);
	}
	//***************************************************
 	//*													*
 	//* Load eigenvector et eigenvalue					*
 	//*													*
 	//***************************************************
 	
 	if (verbose == 1) {printf("Loading Eigenvector\n");}
 	
 	gsl_vector *eval = gsl_vector_alloc(3*atom+10);
	gsl_matrix *evec = gsl_matrix_alloc (3*atom+10,3*atom+10);
 	
	
	load_eigen(eval,evec,eigen_name,3*atom);
	
	
	
	
	for (k=0;k<nb_mode;++k) {
		copy_strc(strc_node_t, strc_node, atom);
		apply_eigen(strc_node_t,atom,evec,mode+k-1,1.0);
		printf("Mode:%d RMSD:%.4f\n",mode+k,sqrt(rmsd_no(strc_node,strc_node_t,atom, align)));
		
	
	}
	
    
   free(connect_h);

	return(0);
 }
 
 

