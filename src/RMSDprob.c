#include "STeM.h"

// Programme la probabilité qu'une proteine est la même conformation (pairwise sequence alignement)

int main(int argc, char *argv[]) {
	int all; /*Nombre d'atomes dans pdb*/

	int all_t; /*Nombre d'atomes dans pdb*/

 	int help_flag = 1;
 	char file_name[500];
 	char check_name[500];
 	char out_name[500];
 	int verbose = 0;
	int i;
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
 	
 	if (verbose == 1) {printf("File I:%s\n",file_name);}
 	
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
 	
 	if (verbose == 1) {printf("File T:%s\n",check_name);}
 	
 	if (verbose == 1) {printf("Connect:%d\n",nconn);}
 	
	if (verbose == 1) {printf("Assigning Structure\n\tAll Atom\n");}
	
	// Array qui comprend tous les connects
	
	int **connect_t=(int **)malloc(nconn*sizeof(int *)); 
    for(i=0;i<nconn;i++) { connect_t[i]=(int *)malloc(6*sizeof(int));}

  assign_connect(check_name,connect_t);
	
	// Assign tous les atoms
	
	struct pdb_atom strc_all_t[all_t];
	build_all_strc(check_name,strc_all_t); // Retourne le nombre de Node
	
	
	// Va faire l'alignement de strc
	
	assign_atom_type(strc_all,all);
 	assign_atom_type(strc_all_t,all_t);

  printf("I found %d Anisou with %d atom !\n",load_anisou(strc_all_t,check_name,all_t),all_t);
	
	
}
