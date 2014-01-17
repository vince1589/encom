#include "STeM.h"
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>



int main(int argc, char *argv[])
{
	int all; /*Nombre d'atomes dans pdb*/
	int atom; /*Nombre de carbone CA*/
	int help_flag = 0;
	
	
	char file_name[500];
	char eigen_name[500] = "eigen.dat";
	char out_name[500] = "out.pdb";
	int nm = -1;
	int verbose = 0;
	
	
	int i;
	
	int nconn;
	int lig = 0;
	double beta = 0.0000001;
	float factor = 1.0;
	for (i = 1;i < argc;i++)
	{
		if (strcmp("-i",argv[i]) == 0) {strcpy(file_name,argv[i+1]);--help_flag;}
		if (strcmp("-o",argv[i]) == 0) {strcpy(out_name,argv[i+1]);}
		if (strcmp("-ieig",argv[i]) == 0) {strcpy(eigen_name,argv[i+1]);}
		if (strcmp("-h",argv[i]) == 0) {help_flag = 1;}
		if (strcmp("-v",argv[i]) == 0) {verbose = 1;}
		if (strcmp("-b",argv[i]) == 0) {float temp;sscanf(argv[i+1],"%f",&temp);beta = temp;}
		if (strcmp("-nm",argv[i]) == 0) {int temp;sscanf(argv[i+1],"%d",&temp);nm = temp;}
		
		
	}
	
	if (help_flag == 1)
	{
		printf("****************************\nHelp Section\n-i\tFile Input (PDB)\n-o\tOutput Name Motion\n-ieig\tFile Name Eigen\n-v\tVerbose\n-sp\tSuper Node Mode (CA, N, C)\n-m\tMode\n-nm\tNombre de mode\n-lig\tTient compte duligand (sauf HOH)\n-prox\tAffiche la probabilite que deux CA puissent interagir\n-prev\tAffiche les directions principales du mouvement de chaque CA ponderees par leur ecart-type\n****************************\n");
		
		return(0);
	}
	
	//***************************************************
	//*                                                                                                     *
	//*Build a structure contaning information on the pdb
	//*                                                                                                     *
	//***************************************************
	
	all = count_atom(file_name);
	
	nconn = count_connect(file_name);
	
	if (verbose == 1) {printf("Connect:%d\n",nconn);}
	
	if (verbose == 1) {printf("Assigning Structure\n\tAll Atom\n");}
	
	// Array qui comprend tous les connects
	
	int **connect_h=(int **)malloc(nconn*sizeof(int *));
	
	for(i=0;i<nconn;i++) { connect_h[i] = (int *)malloc(7*sizeof(int)); }
	
	assign_connect(file_name,connect_h);
	
	// Assign tous les atoms
	
	struct pdb_atom strc_all[all];
	
	atom = build_all_strc(file_name,strc_all); // Retourne le nombre de Node
	
	if (verbose == 1) { printf("    Node:%d\n       Atom:%d\n",atom,all); }
	
	check_lig(strc_all,connect_h,nconn,all);
	
	// Assign les Nodes
	
	if (verbose == 1) { printf("    CA Structure\n"); }
	
	//atom = count_atom_CA_n(strc_all,all,super_node,lig);
	if (verbose == 1) { printf("    Node:%d\n",atom); }
	
	struct pdb_atom strc_node[atom];
	
	atom = build_cord_CA(strc_all, strc_node,all,lig,connect_h,nconn);
	
	if (verbose == 1) { printf("    Assign Node:%d\n",atom); }
	
	
	//***************************************************
	//*                                                                                                     *
	//* Load eigenvector et eigenvalue									*
	//*                                                                                                     *
	//***************************************************
	
	if (verbose == 1) {printf("Loading Eigenvector\n");}
	
	gsl_vector *eval = gsl_vector_alloc(3*atom);
	
	gsl_matrix *evec= gsl_matrix_alloc (3*atom,3*atom);
	
	load_eigen(eval,evec,eigen_name,3*atom);
	
	gsl_matrix *k_totinv = gsl_matrix_alloc(atom*3, atom*3); /* Déclare et crée une matrice qui va être le pseudo inverse */
	
	k_cov_inv_matrix_stem(k_totinv,atom,eval,evec,6,atom*3-6); /* Génère une matrice contenant les superéléments diagonaux de la pseudo-inverse. */
	
	printf("Generating gaussians\n");
	
	double_gauss(strc_node, k_totinv, atom, beta);
	
	gsl_vector_free(eval);
	gsl_matrix_free(evec);
	
	assign_anisou_all(strc_node,atom, strc_all, all, 1);
	
	write_anisou_file(out_name, strc_all, all, 1);
	
	free(connect_h);
	
	return(1);
	
}