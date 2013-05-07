#include "STeM.h"

int main(int argc, char *argv[]) {

	int all; /*Nombre d'atomes dans pdb*/
	int atom; /*Nombre de carbone CA*/
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
 		
 	}
 	 	 	
 	if (help_flag == 0) { } else {
 		printf("****************************\nHelp Section\n-i\tFile Input (PDB)\n-init\tValeur de base de la matrice\n-m\tMatrix Interaction File (8X8)\n-o\tOutput templaate\n-v\tVerbose\n****************************\n");
 		return(0); 
 	} 

 	//***************************************************
 	//*													*
 	//*Build a structure contaning information on the pdb
 	//*													*
 	//***************************************************
 	
 	all = count_atom(file_name);
 	nconn = count_connect(file_name);
	if (verbose == 1) {printf("Filename:%s\n",file_name);} 	
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
	
	//***************************************************
 	//*													*
 	//*Build template									*
 	//*													*
 	//***************************************************
	
	assign_atom_type(strc_all, all);
	gsl_matrix *vcon = gsl_matrix_alloc(all,all);
	gsl_matrix *inter_m = gsl_matrix_alloc(8,8);
	gsl_matrix *templaate = gsl_matrix_alloc(atom, atom);
	gsl_matrix_set_all(templaate,vinit);
	gsl_matrix_set_all(vcon,0);
	
	if (verbose == 1) {printf("Do Vcon !!!\n");}
	
	vcon_file_dom(strc_all,vcon,all);
	
	
	
	if (verbose == 1) {printf("Reading Interaction Matrix %s\n",matrix_name);}
	load_matrix(inter_m,matrix_name);
	//write_matrix("vcon_vince.dat", vcon,all,all);
	if (verbose == 1) {printf("Building templaate\n");}
	all_interaction(strc_all,all, atom, templaate,lig,vcon,inter_m,strc_node);
	if (verbose == 1) {printf("Writing Things\n");}
//	write_matrix(out_name, templaate,atom,atom);
	if (print_flag == 1) {print_templaate(strc_node,atom,templaate,"templaate.pdb",1,9999);}
	//print_templaate(pdb_CA_atom,atom,templaate,"templaate_neg.pdb",-999,0);
	printf("Success\n");
	gsl_matrix_free(templaate);
	gsl_matrix_free(inter_m);
	gsl_matrix_free(vcon);
	
}
