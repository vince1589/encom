#include "STeM.h"

float simple_vector_rmsd(float a,float b,float c,float x,float y,float z);
float simple_vector_lenght(float a,float b,float c);
int main(int argc, char *argv[]) {

	int all; /*Nombre d'atomes dans pdb*/
	int atom; /*Nombre de carbone CA*/
	int all_t; /*Nombre d'atomes dans pdb*/
	int atom_t; /*Nombre de carbone CA*/
 	int help_flag = 1;
 	char file_name[500];
 	char check_name[500];
 	char eigen_name[100] = "eigen.dat";
 	int verbose = 0;
	int i,j,k;
	int nb_mode = 2;
	int mode = 6;
	int lig = 1;
	int nconn;
	int torsion = 0;
	float ligalign = 5; // Flag/valeur pour aligner seulement les résidus dans un cutoff du ligand, 0, one le fait pas... > 0... le cutoff
 	for (i = 1;i < argc;i++) {
 		if (strcmp("-i",argv[i]) == 0) {strcpy(file_name,argv[i+1]);--help_flag;}
 		if (strcmp("-h",argv[i]) == 0) {help_flag = 1;}
 		if (strcmp("-v",argv[i]) == 0) {verbose = 1;}
 		if (strcmp("-lig",argv[i]) == 0) {lig= 1;} 
 		if (strcmp("-t",argv[i]) == 0) {strcpy(check_name,argv[i+1]);help_flag = 0;}
 		if (strcmp("-m",argv[i]) == 0) {int temp;sscanf(argv[i+1],"%d",&temp);mode = temp;}
 		if (strcmp("-nm",argv[i]) == 0) {int temp;sscanf(argv[i+1],"%d",&temp);nb_mode = temp;}
 		if (strcmp("-ieig",argv[i]) == 0) {strcpy(eigen_name,argv[i+1]);}
 		if (strcmp("-ligc",argv[i]) == 0) {float temp;sscanf(argv[i+1],"%f",&temp);ligalign = temp;}
 		if (strcmp("-angle",argv[i]) == 0) {torsion = 1;}
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

	// Free Connect
		
	//for(i=0;i<nconn;i++) {printf("I:%d\n",i);free(connect_h[i]);}
	//free(connect_h);
	
	//Construit la structure a comparer
 	nconn = 0;
 	all_t = count_atom(check_name);
 	nconn = count_connect(check_name);
 	if (verbose == 1) {printf("Sec file:%s\n",check_name);}
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
	
	if (score/atom < 0.8) {
 		printf("Low Score... Will try an homemade alignement !!!\n");
 		score = node_align_low(strc_node,atom,strc_node_t,atom_t,align);
 		
 	}
	printf("RMSD:%8.5f Score: %d/%d\n",sqrt(rmsd_no(strc_node,strc_node_t,atom, align)),score,atom);	
	if (ligalign > 0) {

		score = node_align_lig(strc_node,atom,strc_node_t,atom_t,align,strc_all,all,strc_all_t,all_t,ligalign);
		
		printf("RMSD:%8.5f Score: %d/%d\n",sqrt(rmsd_no(strc_node,strc_node_t,atom, align)),score,atom);
	}
	//***************************************************
 	//*													*
 	//* Load eigenvector et eigenvalue					*
 	//*													*
 	//***************************************************
 	
 	if (verbose == 1) {printf("Loading Eigenvector\n");}
 	
 	gsl_vector *eval = gsl_vector_alloc((3-torsion)*atom+10);
	gsl_matrix *evec = gsl_matrix_alloc ((3-torsion)*atom+10,(3-torsion)*atom+10);
 	
	
	load_eigen(eval,evec,eigen_name,3*atom);
	
	// Mode 6
	for (k=0;k<nb_mode;++k) {
	float sum = 0;
	float longueur = 0;
	for(i=0;i<atom;++i) {
		if (align[i] == -1) {continue;}
		longueur += simple_vector_lenght(gsl_matrix_get(evec,3*i,mode+k),gsl_matrix_get(evec,3*i+1,mode+k),gsl_matrix_get(evec,3*i+2,mode+k));
		for(j=i+1;j<atom;++j) {
			if (align[j] == -1) {continue;}
			//printf("I:%d J:%d\n",i,j);
			//printf("I:%d %.3f %.3f %.3f\n",i,gsl_matrix_get(evec,3*i,mode),gsl_matrix_get(evec,3*i+1,mode),gsl_matrix_get(evec,3*i+2,mode));
			float rmsd = simple_vector_rmsd(gsl_matrix_get(evec,3*i,mode+k),gsl_matrix_get(evec,3*i+1,mode+k),gsl_matrix_get(evec,3*i+2,mode+k),gsl_matrix_get(evec,3*j,mode+k),gsl_matrix_get(evec,3*j+1,mode+k),gsl_matrix_get(evec,3*j+2,mode+k));
			//printf("rmsd:%f\n",rmsd);
			sum += rmsd;
	
		}	
	}
	printf("Mode:%d RMSD:%.3f Longueur:%.3f\n",mode+k+1,sum,longueur);
	}
	
    
    free(connect_h);
	free(connect_t);
	return(0);
 }
 
 
 float simple_vector_rmsd(float a,float b,float c,float x,float y,float z) {
 	
 	
 	float temp = pow(a-x,2)+pow(b-y,2)+pow(c-z,2);
 	
 	
 
 	return(temp);
 }
 
 float simple_vector_lenght(float a,float b,float c) {
 
 	float temp = pow(a,2)+pow(b,2)+pow(c,2);
 	
 	
 
 	return(temp);
 
 }
