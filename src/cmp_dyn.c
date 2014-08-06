#include "STeM.h"

int main(int argc, char *argv[])
{
	int all; /*Nombre d'atomes dans pdb*/
	int atom; /*Nombre de carbone CA*/
	int help_flag = 0;
	
	char file_name[500];
	char check_name[500];
	int verbose = 0;
	int mode = 6;
	
	float vinit = 1; // Valeur de base
	float bond_factor = 1;		// Facteur pour poid des bond strechcing
	float angle_factor = 1;		// Facteur pour poid des angles
	double K_phi1 = 1;				// Facteurs pour angles dièdres
	double K_phi3 = 0.5;
	float init_templaate = 1;
	float kp_factor = 1;					// Facteur pour poid des angles dièdres
	char inputname[500] ="none";
	char matrix_name[500];
	
	int i, j, k, l;
	
	int nconn;
	int lig = 0;
	int lig_t = 0;
	float factor = 1.0;
	for (i = 1;i < argc;i++)
	{
		if (strcmp("-i",argv[i]) == 0) {strcpy(file_name,argv[i+1]);--help_flag;}
		if (strcmp("-t",argv[i]) == 0) {strcpy(check_name,argv[i+1]);--help_flag;}
		
		if (strcmp("-init",argv[i]) == 0) {float temp;sscanf(argv[i+1],"%f",&temp);vinit = temp;}
 		if (strcmp("-kr",argv[i]) == 0) {float temp;sscanf(argv[i+1],"%f",&temp);bond_factor = temp;}
 		if (strcmp("-kt",argv[i]) == 0) {float temp;sscanf(argv[i+1],"%f",&temp);angle_factor = temp;}
 		if (strcmp("-kpf",argv[i]) == 0) {float temp;sscanf(argv[i+1],"%f",&temp); kp_factor = temp;}
 		
 		if (strcmp("-m",argv[i]) == 0) {strcpy(matrix_name,argv[i+1]);help_flag = 0;}
		
		if (strcmp("-lig",argv[i]) == 0) {lig = 1;}
		if (strcmp("-ligt",argv[i]) == 0) {lig_t = 1;}
		
		if (strcmp("-h",argv[i]) == 0) {help_flag = 1;}
		if (strcmp("-v",argv[i]) == 0) {verbose = 1;}
	}
	
	if (help_flag == 1)
	{
		printf("****************************\nHelp Section\n-i\tFile Input (PDB)\n-o\tOutput Name Motion\n-ieig\tFile Name Eigen\n-v\tVerbose\n-sp\tSuper Node Mode (CA, N, C)\n-m\tMode\n-nm\tNombre de mode\n-lig\tTient compte duligand (sauf HOH)\n-prox\tAffiche la probabilite que deux CA puissent interagir\n-prev\tAffiche les directions principales du mouvement de chaque CA ponderees par leur ecart-type\n****************************\n");
		
		return(0);
	}
	
	//***************************************************
 	//*													*
 	//*Builds a structure contaning information on the initial pdb structure
 	//*													*
 	//***************************************************
 	
 	all = count_atom(file_name);
	
 	nconn = count_connect(file_name);
 	
 	if (verbose == 1) {printf("Connect:%d\n",nconn);}
 	
	if (verbose == 1) {printf("Assigning Structure\n\tAll Atom\n");}
	
	// Array with all connects
	
	int **connect_h=(int **)malloc(nconn*sizeof(int *)); 
	
	for(i=0;i<nconn;i++) { connect_h[i]=(int *)malloc(6*sizeof(int));}
	
	assign_connect(file_name,connect_h);
	
	// Assigns all the atoms
	
	struct pdb_atom strc_all[all];
	
	atom = build_all_strc(file_name,strc_all); // Retourne le nombre de Node
	
	if (atom > 800) {printf("Too much nodes .... To fix, ask vincent.frappier@usherbrooke.ca\n");return(1);}
	
	if (verbose == 1) {printf("	Atom:%d\n",all);}
	
	check_lig(strc_all,connect_h,nconn,all);
	
	// Assigns all Nodes
	
	if (verbose == 1) {printf("	CA Structure\n");}
	
	if (verbose == 1) {printf("	Node:%d\n",atom);}
	
	struct pdb_atom strc_node[atom];
	
	atom = build_cord_CA(strc_all, strc_node,all,lig,connect_h,nconn);
	
	if (verbose == 1) {printf("	Assign Node:%d\n",atom);}
	
	// Free Connect
		
	//for(i=0;i<nconn;i++) {printf("I:%d\n",i);free(connect_h[i]);}
	//free(connect_h);
	
	printf("Check 1\n");
	
	//***************************************************
	//*													*
	//*Builds a structure contaning information on the target pdb structure
	//*													*
	//***************************************************
	
 	nconn = 0;
	
 	int all_t = count_atom(check_name);
	
 	nconn = count_connect(check_name);
 	
 	if (verbose == 1) {printf("Connect:%d\n",nconn);}
 	
	if (verbose == 1) {printf("Assigning Structure\n\tAll Atom\n");}
	
	// Array with all connects
	
	int **connect_t=(int **)malloc(nconn*sizeof(int *));
	
	for(i=0;i<nconn;i++) { connect_t[i]=(int *)malloc(6*sizeof(int));}
	
	assign_connect(check_name,connect_t);
	
	// Assigns all the atoms
	
	struct pdb_atom strc_all_t[all_t];
	
	int atom_t = build_all_strc(check_name,strc_all_t); // Retourne le nombre de Node
	
	if (atom_t > 800) {printf("Too much node.... To fix, ask vincent.frappier@usherbrooke.ca\n");return(1);}
	
	if (verbose == 1) {printf("	Atom:%d\n",all_t);}
	
	check_lig(strc_all_t,connect_t,nconn,all_t);
	
	// Assigns all Nodes
	
	if (verbose == 1) {printf("	CA Structure\n");}
	
	if (verbose == 1) {printf("	Node:%d\n",atom_t);}
	
	struct pdb_atom strc_node_t[atom_t];

	atom_t = build_cord_CA(strc_all_t, strc_node_t,all_t,lig_t,connect_t,nconn);
	
	if (verbose == 1) {printf("	Assign Node:%d\n",atom_t);}
	
	printf("Check 2\n");
	
	//***************************************************
	//*													*
	//*Aligns both structures										*
	//*													*
	//***************************************************
	
 	int align[atom];
	
 	int score = node_align(strc_node,atom,strc_node_t,atom_t,align);
	
 	printf("RMSD:%8.5f Score: %d/%d\n",sqrt(rmsd_no(strc_node,strc_node_t,atom, align)),score,atom);
	
	if ((float)score/(float)atom < 0.8)
	{
		printf("Low Score... Will try an homemade alignement !!!\n");
		
		score = node_align_onechain(strc_node,atom,strc_node_t,atom_t,align);
		
		printf("RMSD:%8.5f Score: %d/%d\n",sqrt(rmsd_no(strc_node,strc_node_t,atom, align)),score,atom);
	}
	
 	if ((float)score/(float)atom < 0.8)
	{
 		printf("Low Score... Will try an homemade alignement !!!\n");
		
 		score = node_align_low(strc_node,atom,strc_node_t,atom_t,align);
		
 		printf("RMSD:%8.5f Score: %d/%d\n",sqrt(rmsd_no(strc_node,strc_node_t,atom, align)),score,atom);
 	}
	
	printf("RMSD:%8.5f Score: %d/%d\n",sqrt(rmsd_yes(strc_node,strc_node_t,atom, align,strc_all,all)),score,atom);
	
	printf("Check 4\n");
	
	//***************************************************
	//*													*
	//*Build hessian matrices										*
	//*													*
	//***************************************************
	
	double **hessian=(double **)malloc(3*atom*sizeof(double *)); // Matrix of the Hessian 1 2 3 (bond, angle, dihedral)
	for(i=0;i<3*atom;i++) { hessian[i]=(double *)malloc(3*atom*sizeof(double));}
	for(i=0;i<3*atom;i++)for(j=0;j<(3*atom);j++){hessian[i][j]=0;}
	
	gsl_matrix *hess = gsl_matrix_alloc(3*atom,3*atom);
	gsl_matrix_set_all(hess, 0);
	gsl_matrix *hess_t = gsl_matrix_alloc(3*atom_t,3*atom_t);
	gsl_matrix_set_all(hess_t, 0);
	
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
	
	if (verbose == 1) {printf("Building Hessian\n");}
	
	if (verbose == 1) {printf("\tCovalent Bond Potential\n");}		
	build_1st_matrix(strc_node,hessian,atom,bond_factor);
	
	if (verbose == 1) {printf("\tAngle Potential\n");}	
	build_2_matrix(strc_node,hessian,atom,angle_factor);
	
	if (verbose == 1) {printf("\tDihedral Potential\n");}	
	build_3_matrix(strc_node, hessian,atom,K_phi1/2+K_phi3*9/2,kp_factor);
	
	if (verbose == 1) {printf("\tNon Local Interaction Potential\n");}	
	build_4h_matrix(strc_node,hessian,atom,1.0,templaate);
	
	if (verbose == 1) {printf("\tAssigning Array\n");}	
	assignArray(hess,hessian,3*atom,3*atom);
	
	gsl_matrix_free(vcon);
	gsl_matrix_free(templaate);
	
	double **hessian_t=(double **)malloc(3*atom_t*sizeof(double *)); // Matrix of the Hessian 1 2 3 (bond, angle, dihedral)
	for(i=0;i<3*atom_t;i++) { hessian_t[i]=(double *)malloc(3*atom_t*sizeof(double));}
	for(i=0;i<3*atom_t;i++)for(j=0;j<(3*atom_t);j++){hessian_t[i][j]=0;}
	
	assign_atom_type(strc_all_t, all_t);
	
	if (strcmp(inputname,"none") == 0) {} else {assign_lig_type(strc_all_t, all_t, inputname);}
	
	gsl_matrix *vcon_t = gsl_matrix_alloc(all_t,all_t);
	gsl_matrix *templaate_t = gsl_matrix_alloc(atom_t*3, atom_t*3);
	gsl_matrix_set_all(templaate_t,vinit);
	gsl_matrix_set_all(vcon_t,0);
	
	if (verbose == 1) {printf("Do Vcon !!!\n");}
	
	vcon_file_dom(strc_all_t,vcon_t,all_t);
	
	//write_matrix("vcon_vince.dat", vcon,all,all);
	if (verbose == 1) {printf("Building templaate\n");}
	all_interaction(strc_all_t,all_t, atom_t, templaate_t,lig,vcon_t,inter_m,strc_node_t);
	
	gsl_matrix_scale (templaate_t, init_templaate);
	
	if (verbose == 1) {printf("Building Hessian\n");}
	
	if (verbose == 1) {printf("\tCovalent Bond Potential\n");}
	build_1st_matrix(strc_node_t,hessian_t,atom_t,bond_factor);
	
	if (verbose == 1) {printf("\tAngle Potential\n");}
	build_2_matrix(strc_node_t,hessian_t,atom_t,angle_factor);
	
	if (verbose == 1) {printf("\tDihedral Potential\n");}	
	build_3_matrix(strc_node_t, hessian_t,atom_t,K_phi1/2+K_phi3*9/2,kp_factor);
	
	if (verbose == 1) {printf("\tNon Local Interaction Potential\n");}	
	build_4h_matrix(strc_node_t,hessian_t,atom_t,1.0,templaate_t);
	
	if (verbose == 1) {printf("\tAssigning Array\n");}
	assignArray(hess_t,hessian_t,3*atom_t,3*atom_t);
	
	gsl_matrix_free(vcon_t);
	gsl_matrix_free(templaate_t);
	
	printf("Check 5\n");
	
	//***************************************************
	//*													*
	//*Build mini evec matrices										*
	//*													*
	//***************************************************
	
	gsl_vector *eval = gsl_vector_alloc(3*atom);
	gsl_matrix *evec = gsl_matrix_alloc (3*atom,3*atom);
	
	gsl_vector *eval_t = gsl_vector_alloc(3*atom_t);
	gsl_matrix *evec_t = gsl_matrix_alloc (3*atom_t,3*atom_t);
	
	diagonalyse_matrix(hess,3*atom,eval,evec);
	diagonalyse_matrix(hess_t,3*atom_t,eval_t,evec_t);
	
	printf("Check 5.1\n");
	
	gsl_matrix_free(hess);
	gsl_matrix_free(hess_t);
	
	gsl_vector *mini_eval = gsl_vector_alloc(3*score);
	gsl_matrix *mini_evec = gsl_matrix_alloc (3*score,3*score);
	
	gsl_vector *mini_eval_t = gsl_vector_alloc(3*score);
	gsl_matrix *mini_evec_t = gsl_matrix_alloc (3*score,3*score);
	
	int sup_line = 0;
	
	int sup_to_node[score];
	
	printf("%1i\n", score);
	
	for(j = 0; j < score; j++)
	{
		//printf("OK\n");
		
		for(k = 0; k < 3; k++)
		{
			//printf("OK\n");
			
			gsl_vector_set(mini_eval, 3*j + k, gsl_vector_get(eval, 3*j + k));
			gsl_vector_set(mini_eval_t, 3*j + k, gsl_vector_get(eval_t, 3*j + k));
			
			sup_line = 0;
			
			for(i = 0; i < atom; i++)
			{
				//printf("OK\n");
				
				if(align[i] == -1) {continue;}
				
				//printf("OK\n");
				
				sup_to_node[sup_line] = i;
				
				//printf("OK\n");
				
				for(l = 0; l < 3; l++)
				{
					//printf("OK\n");
					//printf("%1i\t%1i\n", 3*sup_line + k, 3*j + l);
					
					gsl_matrix_set(mini_evec, 3*sup_line + k, 3*j + l, gsl_matrix_get(evec, 3*i + k, 3*j + l));
					gsl_matrix_set(mini_evec_t, 3*sup_line + k, 3*j + l, gsl_matrix_get(evec_t, 3*align[i] + k, 3*j + l));
				}
				
				sup_line++;
			}
		}
	}
	
	gsl_matrix_free(evec);
	gsl_matrix_free(evec_t);
	
	gsl_vector_free(eval);
	gsl_vector_free(eval_t);
	
	printf("Check 6\n");
	
	//***************************************************
	//*													*
	//*Build mini covariance matrices									*
	//*													*
	//***************************************************
	
	gsl_matrix *mini_cov = gsl_matrix_alloc(3*score, 3*score);
	gsl_matrix_set_all(mini_cov, 0);
	gsl_matrix *mini_cov_t = gsl_matrix_alloc(3*score, 3*score);
	gsl_matrix_set_all(mini_cov_t, 0);
	
	for(i = 0; i < 3*score; i++)
	{
		if(gsl_vector_get(mini_eval, i) > 0.0000001)
		{
			gsl_vector_set(mini_eval, i, 1.0 / gsl_vector_get(mini_eval, i));
			
			for(j = 0; j < 3*score; j++)
			{
				for(k = 0; k < 3*score; k++)
				{
					gsl_matrix_set(mini_cov, j, k, gsl_matrix_get(mini_cov, j, k) + gsl_matrix_get(mini_evec, j, i)*gsl_matrix_get(mini_evec, k, i)*gsl_vector_get(mini_eval, i));
				}
			}
		}
	}
	
	for(i = 0; i < 3*score; i++)
	{
		if(gsl_vector_get(mini_eval_t, i) > 0.0000001)
		{
			gsl_vector_set(mini_eval_t, i, 1.0 / gsl_vector_get(mini_eval_t, i));
			
			for(j = 0; j < 3*score; j++)
			{
				for(k = 0; k < 3*score; k++)
				{
					gsl_matrix_set(mini_cov_t, j, k, gsl_matrix_get(mini_cov_t, j, k) + gsl_matrix_get(mini_evec_t, j, i)*gsl_matrix_get(mini_evec_t, k, i)*gsl_vector_get(mini_eval_t, i));
				}
			}
		}
	}
	
	gsl_matrix *mini_cov_tmp = gsl_matrix_alloc(3*score, 3*score);
	gsl_matrix *mini_cov_tmp_t = gsl_matrix_alloc(3*score, 3*score);
	
	for(i = 0; i < 3*score; i++)
	{
		for(j = 0; j < 3*score; j++)
		{
			gsl_matrix_set(mini_cov_tmp, i, j, gsl_matrix_get(mini_cov, i, j));
			gsl_matrix_set(mini_cov_tmp_t, i, j, gsl_matrix_get(mini_cov_t, i, j));
		}
	}
	
	diagonalyse_matrix(mini_cov_tmp,3*score,mini_eval,mini_evec);
	diagonalyse_matrix(mini_cov_tmp_t,3*score,mini_eval_t,mini_evec_t);
	
	gsl_matrix_free(mini_cov_tmp);
	gsl_matrix_free(mini_cov_tmp_t);
	
	printf("Check 7\n");
	
	printf("Dynamic distance : %1.10f\n", cmp_gauss(mini_cov, mini_eval, mini_cov_t, mini_eval_t, 3*score));
	
	gsl_matrix_free(mini_cov);
	gsl_matrix_free(mini_evec);
	gsl_vector_free(mini_eval);
	gsl_matrix_free(mini_cov_t);
	gsl_matrix_free(mini_evec_t);
	gsl_vector_free(mini_eval_t);
	
	
	free(connect_h);
	
	return(1);
}
