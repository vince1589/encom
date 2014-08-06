#include "STeM.h"

int main(int argc, char *argv[])
{
	int all; /*Number of atoms in the initial PDB*/
	int atom; /*Number of initial CAs*/
	int all_t; /*Number of atoms in the target PDB*/
	int atom_t; /*Number of target CAs*/
	
 	int help_flag = 1;
 	char file_name[500];
 	char check_name[500];
 	char matrix_name[500];
	int print_inter = 0;
	int print_inter_t = 0;
 	char out_name[500];
	char out_name_t[500];
	char out_movie[500];
	int print_movie = 0;
 	int verbose = 0;
	float vinit = 1; // Valeur de base
	float bond_factor = 1;		// Facteur pour poid des bond strechcing
	float angle_factor = 1;		// Facteur pour poid des angles
	double K_phi1 = 1;				// Facteurs pour angles dièdres
	double K_phi3 = 0.5;
	float init_templaate = 1;
	float kp_factor = 1;					// Facteur pour poid des angles dièdres
	char inputname[500] ="none";
	double beta = 0.000005;
	int morph = 0;
	char morph_name[500];
	
	int change_density = 0;
	
	int i,j,k,l;
	
	int lig = 0;
	int ligt = 0;
	int nconn;
	int print_flag = 0;
	float ligalign = 0; // Flag/valeur pour aligner seulement les résidus dans un cutoff du ligand, 0, one le fait pas... > 0... le cutoff
 	for (i = 1;i < argc;i++) {
 		if (strcmp("-i",argv[i]) == 0) {strcpy(file_name,argv[i+1]);--help_flag;}
 		
 		if (strcmp("-h",argv[i]) == 0) {help_flag = 1;}
 		if (strcmp("-v",argv[i]) == 0) {verbose = 1;}
 		if (strcmp("-lig",argv[i]) == 0) {lig= 1;}
 		if (strcmp("-ligt",argv[i]) == 0) {ligt= 1;}
 		
 		if (strcmp("-init",argv[i]) == 0) {float temp;sscanf(argv[i+1],"%f",&temp);vinit = temp;}
 		if (strcmp("-kr",argv[i]) == 0) {float temp;sscanf(argv[i+1],"%f",&temp);bond_factor = temp;}
 		if (strcmp("-kt",argv[i]) == 0) {float temp;sscanf(argv[i+1],"%f",&temp);angle_factor = temp;}
 		if (strcmp("-kpf",argv[i]) == 0) {float temp;sscanf(argv[i+1],"%f",&temp); kp_factor = temp;}
 		
 		if (strcmp("-b",argv[i]) == 0) {float temp;sscanf(argv[i+1],"%f",&temp);beta = temp;}
 		
 		if (strcmp("-t",argv[i]) == 0) {strcpy(check_name,argv[i+1]);help_flag = 0;}
 		
 		if (strcmp("-m",argv[i]) == 0) {strcpy(matrix_name,argv[i+1]);help_flag = 0;}
 		
 		if (strcmp("-o",argv[i]) == 0) {strcpy(out_name,argv[i+1]); print_inter = 1;}
 		if (strcmp("-ot",argv[i]) == 0) {strcpy(out_name_t,argv[i+1]); print_inter_t = 1;}
 		if (strcmp("-om",argv[i]) == 0) {strcpy(out_movie,argv[i+1]); print_movie = 1;}
 		
 		if (strcmp("-ligc",argv[i]) == 0) {float temp;sscanf(argv[i+1],"%f",&temp);ligalign = temp;}
 		
 		if (strcmp("-conf",argv[i]) == 0) {change_density = 1;}
 		
 		if (strcmp("-morph",argv[i]) == 0) {strcpy(morph_name,argv[i+1]); morph = 1;}
 	}
	
 	if (help_flag == 1)
	{
 		printf("****************************\nHelp Section\n-i\tFile Input (PDB)\n-v\tVerbose\n-w\tWeight Vector\n-t\tInitial Value of template (negative value for random)\n\tIf Load Template, multiply the template\n-lt\tLoad input template\n-sp\tSuper Node Mode (CA, N, C)\n-kt\tPoid de l'angle entre les nodes (1)\n-kr\tPoid de la distance entre les nodes (1)\n-f\tFile to fit\n****************************\n");
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
	
	//***************************************************
	//*													*
	//*Builds a structure contaning information on the target pdb structure
	//*													*
	//***************************************************
	
 	nconn = 0;
	
 	all_t = count_atom(check_name);
	
 	nconn = count_connect(check_name);
 	
 	if (verbose == 1) {printf("Connect:%d\n",nconn);}
 	
	if (verbose == 1) {printf("Assigning Structure\n\tAll Atom\n");}
	
	// Array with all connects
	
	int **connect_t=(int **)malloc(nconn*sizeof(int *));
	
	for(i=0;i<nconn;i++) { connect_t[i]=(int *)malloc(6*sizeof(int));}
	
	assign_connect(check_name,connect_t);
	
	// Assigns all the atoms
	
	struct pdb_atom strc_all_t[all_t];
	
	atom_t = build_all_strc(check_name,strc_all_t); // Retourne le nombre de Node
	
	if (atom_t > 800) {printf("Too much node.... To fix, ask vincent.frappier@usherbrooke.ca\n");return(1);}
	
	if (verbose == 1) {printf("	Atom:%d\n",all_t);}
	
	check_lig(strc_all_t,connect_t,nconn,all_t);
	
	// Assigns all Nodes
	
	if (verbose == 1) {printf("	CA Structure\n");}
	
	if (verbose == 1) {printf("	Node:%d\n",atom_t);}
	
	struct pdb_atom strc_node_t[atom_t];

	atom_t = build_cord_CA(strc_all_t, strc_node_t,all_t,ligt,connect_t,nconn);
	
	if (verbose == 1) {printf("	Assign Node:%d\n",atom_t);}
	
	//***************************************************
	//*													*
	//*Aligns both structures
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
 	
 	if (ligalign > 0)
	{
		score = node_align_lig(strc_node,atom,strc_node_t,atom_t,align,strc_all,all,strc_all_t,all_t,ligalign);
		
		printf("RMSD:%8.5f Score: %d/%d\n",sqrt(rmsd_no(strc_node,strc_node_t,atom, align)),score,atom);
	}
	
	printf("RMSD:%8.5f Score: %d/%d\n",sqrt(rmsd_yes(strc_node,strc_node_t,atom, align,strc_all,all)),score,atom);
	
	// Build hessians
	
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
	
	printf("Check 1\n");
	
	// Build mini-hessians (calculated mini_hess object is equal to mini_hess + mini_hess_t
	
	gsl_matrix *mini_hess = gsl_matrix_alloc(3*score, 3*score);
	gsl_matrix_set_all(mini_hess, 0);
	gsl_matrix *mini_hess_i = gsl_matrix_alloc(3*score, 3*score);
	gsl_matrix_set_all(mini_hess_i, 0);
	gsl_matrix *mini_hess_t = gsl_matrix_alloc(3*score, 3*score);
	gsl_matrix_set_all(mini_hess_t, 0);
	
	int sup_line = 0;
	
	int sup_to_node[score];
	
	for(i = 0; i < atom; i++)
	{
		if(align[i] == -1) {continue;}
		
		sup_to_node[sup_line] = i;
		
		int sup_col = 0;
		
		for(j = 0; j < atom; j++)
		{
			if(align[j] == -1) {continue;}
			
			for(k = 0; k < 3; k++)
			{
				for(l = 0; l < 3; l++)
				{
					gsl_matrix_set(mini_hess_i, 3*sup_line + k, 3*sup_col + l, gsl_matrix_get(hess, 3*i + k, 3*j + l));
					gsl_matrix_set(mini_hess, 3*sup_line + k, 3*sup_col + l, gsl_matrix_get(hess, 3*i + k, 3*j + l) + gsl_matrix_get(hess_t, 3*align[i] + k, 3*align[j] + l));
					gsl_matrix_set(mini_hess_t, 3*sup_line + k, 3*sup_col + l, gsl_matrix_get(hess_t, 3*align[i] + k, 3*align[j] + l));
				}
			}
			
			sup_col ++;
		}
		
		sup_line++;
	}
	
	gsl_matrix_free(hess);
	gsl_matrix_free(hess_t);
	
	printf("Check 2\n");
	
	// Invert mini_hess
	
	gsl_vector *eval2 = gsl_vector_alloc(3*score);
	
	gsl_matrix *evec2 = gsl_matrix_alloc (3*score,3*score);
	
	diagonalyse_matrix(mini_hess, 3*score, eval2, evec2);
	
	gsl_matrix_set_all(mini_hess, 0);
	
	for(i = 0; i < 3*score; i++)
	{
		if(gsl_vector_get(eval2, i) > 0.00001)
		{
			for(j = 0; j < 3*score; j++)
			{
				for(k = 0; k < 3*score; k++)
				{
					gsl_matrix_set(mini_hess, j, k, gsl_matrix_get(mini_hess, j, k) + gsl_matrix_get(evec2, j, i)*gsl_matrix_get(evec2, k, i)/gsl_vector_get(eval2, i));
				}
			}
		}
	}
	
	gsl_matrix_free(evec2);
	gsl_vector_free(eval2);
	
	printf("Check 3\n");
	
	// Evaluate delta-conf of init to most stable transitionnal conformer and store delta-conf of init to target.
	
	gsl_vector *del_conf = gsl_vector_alloc(3*score);
	gsl_vector_set_all(del_conf, 0);
	gsl_vector *copy_conf = gsl_vector_alloc(3*score);
	
	for(i = 0; i < score; i++)
	{
		gsl_vector_set(copy_conf, 3*i, strc_node_t[align[i]].x_cord - strc_node[i].x_cord);
		gsl_vector_set(copy_conf, 3*i + 1, strc_node_t[align[i]].y_cord - strc_node[i].y_cord);
		gsl_vector_set(copy_conf, 3*i + 2, strc_node_t[align[i]].z_cord - strc_node[i].z_cord);
	}
	
	for(i = 0; i < 3*score; i++)
	{
		for(j = 0; j < 3*score; j++)
		{
			gsl_vector_set(del_conf, i, gsl_vector_get(del_conf, i) + gsl_matrix_get(mini_hess_t, i, j)*gsl_vector_get(copy_conf, j));
		}
	}
	
	for(i = 0; i < 3*score; i++)
	{
		gsl_vector_set(copy_conf, i, gsl_vector_get(del_conf, i));
	}
	
	gsl_vector_set_all(del_conf, 0);
	
	for(i = 0; i < 3*score; i++)
	{
		for(j = 0; j < 3*score; j++)
		{
			gsl_vector_set(del_conf, i, gsl_vector_get(del_conf, i) + gsl_matrix_get(mini_hess, i, j)*gsl_vector_get(copy_conf, j));
		}
	}
	
	for(i = 0; i < score; i++)
	{
		gsl_vector_set(copy_conf, 3*i, strc_node_t[align[i]].x_cord - strc_node[i].x_cord);
		gsl_vector_set(copy_conf, 3*i + 1, strc_node_t[align[i]].y_cord - strc_node[i].y_cord);
		gsl_vector_set(copy_conf, 3*i + 2, strc_node_t[align[i]].z_cord - strc_node[i].z_cord);
	}
	
	printf("Check 4\n");
	
	// Translate all the nodes (and all its atoms) of both structures to delta-conf of most stable conformer.
	
	for(i = 0; i < score; i++)
	{
		for(j = 0; j < all; j++)
		{
			if(strc_all[j].node == sup_to_node[i])
			{
				strc_all[j].x_cord += gsl_vector_get(del_conf, 3*i);
				strc_all[j].y_cord += gsl_vector_get(del_conf, 3*i + 1);
				strc_all[j].z_cord += gsl_vector_get(del_conf, 3*i + 2);
			}
		}
		
		for(j = 0; j < all_t; j++)
		{
			if(strc_all_t[j].node == align[sup_to_node[i]])
			{
				strc_all_t[j].x_cord += gsl_vector_get(del_conf, 3*i) - gsl_vector_get(copy_conf, 3*i);
				strc_all_t[j].y_cord += gsl_vector_get(del_conf, 3*i + 1) - gsl_vector_get(copy_conf, 3*i + 1);
				strc_all_t[j].z_cord += gsl_vector_get(del_conf, 3*i + 2) - gsl_vector_get(copy_conf, 3*i + 2);
			}
		}
	}
	
	gsl_matrix_free(mini_hess);
	
	printf("Check 5\n");
	
	int aligned[atom_t];
	
	for(i = 0; i < atom_t; i++)
	{
		aligned[i] = -1;
	}
	
	for(i = 0; i < atom; i++)
	{
		if(align[i] != -1)
		{
			aligned[align[i]] = 1;
		}
	}
	
	// Print transition from init
	
	if(print_inter == 1)
	{
		FILE *out_file;
		
		out_file = fopen(out_name, "w");
		
		for (i = 0; i < all; i++)
		{
			if(align[strc_all[i].node] != -1)
			{
				if (strc_all[i].atom_type == 1) {fprintf(out_file,"ATOM  ");}
				if (strc_all[i].atom_type == 2) {fprintf(out_file,"HETATM");}
				if (strc_all[i].atom_type == 3) {fprintf(out_file,"HETATM");}
				fprintf(out_file,"%5.d %s%s %s%4.d%12.3f%8.3f%8.3f  1.00  %2.2f\n",
					strc_all[i].atom_number,
					strc_all[i].atom_prot_type,
					strc_all[i].res_type,
					strc_all[i].chain,
					strc_all[i].res_number,
					strc_all[i].x_cord,
					strc_all[i].y_cord,
					strc_all[i].z_cord,
					strc_all[i].b_factor
					);
			}
		}
		
		fclose(out_file);
	}
	
	printf("Check 6\n");
	
	// Print transition from target
	
	if(print_inter_t == 1)
	{
		FILE *out_file_t;
		
		out_file_t = fopen(out_name_t, "w");
		
		for (i = 0; i < all_t; i++)
		{
			if(aligned[strc_all_t[i].node] != -1)
			{
				if (strc_all_t[i].atom_type == 1) {fprintf(out_file_t,"ATOM  ");}
				if (strc_all_t[i].atom_type == 2) {fprintf(out_file_t,"HETATM");}
				if (strc_all_t[i].atom_type == 3) {fprintf(out_file_t,"HETATM");}
				fprintf(out_file_t,"%5.d %s%s %s%4.d%12.3f%8.3f%8.3f  1.00  %2.2f\n",
					strc_all_t[i].atom_number,
					strc_all_t[i].atom_prot_type,
					strc_all_t[i].res_type,
					strc_all_t[i].chain,
					strc_all_t[i].res_number,
					strc_all_t[i].x_cord,
					strc_all_t[i].y_cord,
					strc_all_t[i].z_cord,
					strc_all_t[i].b_factor
					);
			}
		}
		
		fclose(out_file_t);
	}
	
	printf("Check 7\n");
	
	// Translate all the nodes (and all its atoms) of init structure back to its original conformation.
	
	for(i = 0; i < score; i++)
	{
		for(j = 0; j < all; j++)
		{
			if(strc_all[j].node == sup_to_node[i])
			{
				strc_all[j].x_cord -= gsl_vector_get(del_conf, 3*i);
				strc_all[j].y_cord -= gsl_vector_get(del_conf, 3*i + 1);
				strc_all[j].z_cord -= gsl_vector_get(del_conf, 3*i + 2);
			}
		}
	}
	
	// Print transition from init to target passing by inter
	
	if(print_movie == 1)
	{
		FILE *out_file_m;
		
		out_file_m = fopen(out_movie, "w");
		
		for(j = 0; j < 30; j++)
		{
			fprintf(out_file_m, "Model %1i\n", j + 1);
			
			for (i = 0; i < all; i++)
			{
				if(align[strc_all[i].node] != -1)
				{
					if (strc_all[i].atom_type == 1) {fprintf(out_file_m,"ATOM  ");}
					if (strc_all[i].atom_type == 2) {fprintf(out_file_m,"HETATM");}
					if (strc_all[i].atom_type == 3) {fprintf(out_file_m,"HETATM");}
					fprintf(out_file_m,"%5.d %s%s %s%4.d%12.3f%8.3f%8.3f  1.00  %2.2f\n",
						strc_all[i].atom_number,
						strc_all[i].atom_prot_type,
						strc_all[i].res_type,
						strc_all[i].chain,
						strc_all[i].res_number,
						strc_all[i].x_cord,
						strc_all[i].y_cord,
						strc_all[i].z_cord,
						strc_all[i].b_factor
						);
				}
			}
			
			fprintf(out_file_m, "TER\nENDMDL\n\n");
			
			for(i = 0; i < score; i++)
			{
				for(k = 0; k < all; k++)
				{
					if(strc_all[k].node == sup_to_node[i])
					{
						strc_all[k].x_cord += gsl_vector_get(del_conf, 3*i)/30.0;
						strc_all[k].y_cord += gsl_vector_get(del_conf, 3*i + 1)/30.0;
						strc_all[k].z_cord += gsl_vector_get(del_conf, 3*i + 2)/30.0;
					}
				}
			}
		}
		
		for(j = 0; j < 31; j++)
		{
			fprintf(out_file_m, "Model %1i\n", j + 31);
			
			for (i = 0; i < all_t; i++)
			{
				if(aligned[strc_all_t[i].node] != -1)
				{
					if (strc_all_t[i].atom_type == 1) {fprintf(out_file_m,"ATOM  ");}
					if (strc_all_t[i].atom_type == 2) {fprintf(out_file_m,"HETATM");}
					if (strc_all_t[i].atom_type == 3) {fprintf(out_file_m,"HETATM");}
					fprintf(out_file_m,"%5.d %s%s %s%4.d%12.3f%8.3f%8.3f  1.00  %2.2f\n",
						strc_all_t[i].atom_number,
						strc_all_t[i].atom_prot_type,
						strc_all_t[i].res_type,
						strc_all_t[i].chain,
						strc_all_t[i].res_number,
						strc_all_t[i].x_cord,
						strc_all_t[i].y_cord,
						strc_all_t[i].z_cord,
						strc_all_t[i].b_factor
						);
				}
			}
			
			fprintf(out_file_m, "TER\nENDMDL\n\n");
			
			for(i = 0; i < score; i++)
			{
				for(k = 0; k < all_t; k++)
				{
					if(strc_all_t[k].node == align[sup_to_node[i]])
					{
						strc_all_t[k].x_cord -= (gsl_vector_get(del_conf, 3*i) - gsl_vector_get(copy_conf, 3*i))/30.0;
						strc_all_t[k].y_cord -= (gsl_vector_get(del_conf, 3*i + 1) - gsl_vector_get(copy_conf, 3*i + 1))/30.0;
						strc_all_t[k].z_cord -= (gsl_vector_get(del_conf, 3*i + 2) - gsl_vector_get(copy_conf, 3*i + 2))/30.0;
					}
				}
			}
		}
		
		fclose(out_file_m);
	}
	
	printf("Check 8\n");
	
	// Translate all the nodes (and all its atoms) of init structure back to its original conformation if print_movie == 1
	
	if(print_movie == 1)
	{
		for(i = 0; i < score; i++)
		{
			for(j = 0; j < all; j++)
			{
				if(strc_all[j].node == sup_to_node[i])
				{
					strc_all[j].x_cord -= gsl_vector_get(del_conf, 3*i);
					strc_all[j].y_cord -= gsl_vector_get(del_conf, 3*i + 1);
					strc_all[j].z_cord -= gsl_vector_get(del_conf, 3*i + 2);
				}
			}
		}
	}
	
	printf("Check 8.1\n");
	
	// Print transition from init to target by morphing
	
	if(morph == 1)
	{
		FILE *out_file_morph;
		
		out_file_morph = fopen(morph_name, "w");
		
		printf("Check 8.2\n");
		
		for(j = 0; j < 60; j++)
		{
			fprintf(out_file_morph, "Model %1i\n", j + 1);
			
			for (i = 0; i < all; i++)
			{
				if(align[strc_all[i].node] != -1)
				{
					if (strc_all[i].atom_type == 1) {fprintf(out_file_morph,"ATOM  ");}
					if (strc_all[i].atom_type == 2) {fprintf(out_file_morph,"HETATM");}
					if (strc_all[i].atom_type == 3) {fprintf(out_file_morph,"HETATM");}
					fprintf(out_file_morph,"%5.d %s%s %s%4.d%12.3f%8.3f%8.3f  1.00  %2.2f\n",
						strc_all[i].atom_number,
						strc_all[i].atom_prot_type,
						strc_all[i].res_type,
						strc_all[i].chain,
						strc_all[i].res_number,
						strc_all[i].x_cord,
						strc_all[i].y_cord,
						strc_all[i].z_cord,
						strc_all[i].b_factor
						);
				}
			}
			
			fprintf(out_file_morph, "TER\nENDMDL\n\n");
			
			printf("Check 8.3\n");
			
			for(i = 0; i < score; i++)
			{
				for(k = 0; k < all; k++)
				{
					if(strc_all[k].node == sup_to_node[i])
					{
						strc_all[k].x_cord += gsl_vector_get(copy_conf, 3*i)/60.0;
						strc_all[k].y_cord += gsl_vector_get(copy_conf, 3*i + 1)/60.0;
						strc_all[k].z_cord += gsl_vector_get(copy_conf, 3*i + 2)/60.0;
					}
				}
			}
			
			printf("Check 8.4\n");
		}
		
		fclose(out_file_morph);
	}
	
	printf("Check 9\n");
	
	// Evaluate theoretical delta-G
	
	gsl_vector *dummy_conf = gsl_vector_alloc(3*score);
	gsl_vector_set_all(dummy_conf, 0);
	gsl_vector *dummy_conf_t = gsl_vector_alloc(3*score);
	gsl_vector_set_all(dummy_conf_t, 0);
	
	for(i = 0; i < 3*score; i++)
	{
		for(j = 0; j < 3*score; j++)
		{
			gsl_vector_set(dummy_conf, i, gsl_vector_get(dummy_conf, i) + gsl_vector_get(del_conf, j)*gsl_matrix_get(mini_hess_i, j, i));
			
			gsl_vector_set(dummy_conf_t, i, gsl_vector_get(dummy_conf_t, i) + (gsl_vector_get(del_conf, j) - gsl_vector_get(copy_conf, j))*gsl_matrix_get(mini_hess_t, j, i));
		}
	}
	
	double energy_ic = 0.0;
	double energy_tc = 0.0;
	
	for(i = 0; i < 3*score; i++)
	{
		energy_ic += gsl_vector_get(dummy_conf, i)*gsl_vector_get(del_conf, i);
		
		energy_tc += gsl_vector_get(dummy_conf_t, i)*(gsl_vector_get(del_conf, i) - gsl_vector_get(copy_conf, i));
	}
	
	gsl_vector_free(dummy_conf_t);
	gsl_vector_free(del_conf);
	
	printf("Check 10\n");
	
	if(change_density == 1)
	{
		gsl_vector_set_all(dummy_conf, 0);
		
		gsl_vector *eval = gsl_vector_alloc(3*score);
		gsl_vector_set_all(eval, 0);
		gsl_vector *eval_t = gsl_vector_alloc(3*score);
		gsl_vector_set_all(eval_t, 0);
		
		gsl_matrix *evec = gsl_matrix_alloc (3*score,3*score);
		gsl_matrix_set_all(evec, 0);
		gsl_matrix *evec_t = gsl_matrix_alloc (3*score,3*score);
		gsl_matrix_set_all(evec_t, 0);
		
		diagonalyse_matrix(mini_hess_i, 3*score, eval, evec);
		diagonalyse_matrix(mini_hess_t, 3*score, eval_t, evec_t);
		
		gsl_matrix_free(mini_hess_i);
		gsl_matrix_free(mini_hess_t);
		
		gsl_matrix *comb_hess = gsl_matrix_alloc(3*score, 3*score);
		gsl_matrix_set_all(comb_hess, 0);
		
		double entro = 0.0;
		
		double entro_t = 0.0;
		
		for(i = 0; i < 3*score; i++)
		{
			if(gsl_vector_get(eval, i) > 0.00001)
			{
				for(j = 0; j < 3*score; j++)
				{
					for(k = 0; k < 3*score; k++)
					{
						gsl_matrix_set(comb_hess, j, k, gsl_matrix_get(comb_hess, j, k) + gsl_matrix_get(evec, j, i)*gsl_matrix_get(evec, k, i)/gsl_vector_get(eval, i));
					}
				}
				
				entro += log(3.141592653589793238462643383279) - log(beta * gsl_vector_get(eval, i));
			}
			
			if(gsl_vector_get(eval_t, i) > 0.00001)
			{
				for(j = 0; j < 3*score; j++)
				{
					for(k = 0; k < 3*score; k++)
					{
						gsl_matrix_set(comb_hess, j, k, gsl_matrix_get(comb_hess, j, k) + gsl_matrix_get(evec_t, j, i)*gsl_matrix_get(evec_t, k, i)/gsl_vector_get(eval_t, i));
					}
				}
				
				entro_t += log(3.141592653589793238462643383279) - log(beta * gsl_vector_get(eval_t, i));
			}
		}
		
		gsl_vector_free(eval_t);
		gsl_matrix_free(evec_t);
		
		gsl_vector_set_all(eval, 0);
		gsl_matrix_set_all(evec, 0);
		
		diagonalyse_matrix(comb_hess, 3*score, eval, evec);
		gsl_matrix_set_all(comb_hess, 0);
		
// 		double comb_det = 1.0;
		
		for(i = 0; i < 3*score; i++)
		{
			if(gsl_vector_get(eval, i) > 0.00001)
			{
				for(j = 0; j < 3*score; j++)
				{
					for(k = 0; k < 3*score; k++)
					{
						gsl_matrix_set(comb_hess, j, k, gsl_matrix_get(comb_hess, j, k) + gsl_matrix_get(evec, j, i)*gsl_matrix_get(evec, k, i)/gsl_vector_get(eval, i));
					}
				}
				
// 				comb_det *= beta / (3.141592653589793238462643383279 * gsl_vector_get(eval, i));
			}
		}
		
		for(i = 0; i < 3*score; i++)
		{
			for(j = 0; j < 3*score; j++)
			{
				gsl_vector_set(dummy_conf, i, gsl_vector_get(dummy_conf, i) + gsl_vector_get(copy_conf, j)*gsl_matrix_get(comb_hess, j, i));
			}
		}
		
		double energy_conf = 0.0;
		
		for(i = 0; i < 3*score; i++)
		{
			energy_conf += gsl_vector_get(dummy_conf, i)*gsl_vector_get(copy_conf, i);
		}
		
		double dummy_energy = 0.0;
		
		double prob_conf = 0.0;
		
		dummy_energy = -1.0 * beta * energy_conf;
		
		prob_conf = exp(dummy_energy);
		
		printf("Weighted probability density of exact conformational change at beta = %1.10f : %1.100f\n", beta, prob_conf);
		
		printf("Delta-S (target - init) : %1.10f\n", entro_t - entro);
		
		printf("Delta-H (target - init) : %1.10f\n", energy_ic - energy_tc);
		
		dummy_energy = -310.25 * (entro_t - entro) + energy_ic - energy_tc;
		
		printf("Delta-G (target - init) : %1.10f\n", dummy_energy);
		
		double K_eq = exp(-1.0 * dummy_energy / (310.25 * 8.3145));
		
		printf("Equilibrium constant (target / init) : %1.10f\n", K_eq);
	}
	else
	{
		printf("Energy from init to inter : %1.10f\nEnergy from target to inter : %1.10f\nDelta-H (target - init) : %1.10f\n", energy_ic, energy_tc, energy_ic - energy_tc);
	}
	
	gsl_vector_free(copy_conf);
	gsl_vector_free(dummy_conf);
}
