#include "STeM.h"




int main(int argc, char *argv[])
{
	int all; /*Number of atoms in the initial PDB*/
	int atom; /*Number of initial CAs*/
	int all_t; /*Number of atoms in the target PDB*/
	int atom_t; /*Number of target CAs*/
	
 	int help_flag = 1;
 	char file_name[500];
	char file_eigen[500];
 	char check_name[500];
	char check_eigen[500];
	int print_inter = 0;
	int print_inter_t = 0;
 	char out_name[500];
	char out_name_t[500];
	char out_movie[500];
	int print_movie = 0;
 	int verbose = 0;
	
	int i,j,k,l;
	
	int lig = 0;
	int nconn;
	int print_flag = 0;
	float ligalign = 0; // Flag/valeur pour aligner seulement les résidus dans un cutoff du ligand, 0, one le fait pas... > 0... le cutoff
 	for (i = 1;i < argc;i++) {
 		if (strcmp("-i",argv[i]) == 0) {strcpy(file_name,argv[i+1]);--help_flag;}
 		if (strcmp("-ieig",argv[i]) == 0) {strcpy(file_eigen,argv[i+1]);}
 		
 		if (strcmp("-h",argv[i]) == 0) {help_flag = 1;}
 		if (strcmp("-v",argv[i]) == 0) {verbose = 1;}
 		if (strcmp("-lig",argv[i]) == 0) {lig= 1;} 
 		
 		if (strcmp("-t",argv[i]) == 0) {strcpy(check_name,argv[i+1]);help_flag = 0;}
 		if (strcmp("-teig",argv[i]) == 0) {strcpy(check_eigen,argv[i+1]);}
 		
 		if (strcmp("-o",argv[i]) == 0) {strcpy(out_name,argv[i+1]); print_inter = 1;}
 		if (strcmp("-ot",argv[i]) == 0) {strcpy(out_name_t,argv[i+1]); print_inter_t = 1;}
 		if (strcmp("-om",argv[i]) == 0) {strcpy(out_movie,argv[i+1]); print_movie = 1;}
 		
 		if (strcmp("-ligc",argv[i]) == 0) {float temp;sscanf(argv[i+1],"%f",&temp);ligalign = temp;}
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

	atom_t = build_cord_CA(strc_all_t, strc_node_t,all_t,lig,connect_t,nconn);
	
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
		
 		score = node_align_low(strc_node,atom,strc_node_t,atom_t,align);
		
 		printf("RMSD:%8.5f Score: %d/%d\n",sqrt(rmsd_no(strc_node,strc_node_t,atom, align)),score,atom);
 	}
 	
 	if (ligalign > 0)
	{
		score = node_align_lig(strc_node,atom,strc_node_t,atom_t,align,strc_all,all,strc_all_t,all_t,ligalign);
		
		printf("RMSD:%8.5f Score: %d/%d\n",sqrt(rmsd_no(strc_node,strc_node_t,atom, align)),score,atom);
	}
	
	printf("RMSD:%8.5f Score: %d/%d\n",sqrt(rmsd_yes(strc_node,strc_node_t,atom, align,strc_all,all)),score,atom);
	
	// Load eigenvalues and eigenvectors of each structure
	
	gsl_vector *eval = gsl_vector_alloc(3*atom);
	gsl_vector *eval_t = gsl_vector_alloc(3*atom_t);
	
	gsl_matrix *evec = gsl_matrix_alloc (3*atom,3*atom);
	gsl_matrix *evec_t = gsl_matrix_alloc (3*atom_t,3*atom_t);
	
	load_eigen(eval,evec,file_eigen,3*atom);
	load_eigen(eval_t,evec_t,check_eigen,3*atom_t);
	
	printf("Check 0\n");
	
	// Rebuild hessians
	
	gsl_matrix *hess = gsl_matrix_alloc(3*atom,3*atom);
	gsl_matrix_set_all(hess, 0);
	gsl_matrix *hess_t = gsl_matrix_alloc(3*atom_t,3*atom_t);
	gsl_matrix_set_all(hess_t, 0);
	
	for(i = 0; i < 3*atom; i++)
	{
		if(gsl_vector_get(eval, i) < 0.00001) {printf("Skip node %1i of init : %1.5f\n", i, gsl_vector_get(eval, i)); continue;}
		
		for(j = 0; j < 3*atom; j++)
		{
			for(k = 0; k < 3*atom; k++)
			{
				gsl_matrix_set(hess, j, k, gsl_matrix_get(hess, j, k) + gsl_matrix_get(evec, j, i)*gsl_matrix_get(evec, k, i)*gsl_vector_get(eval, i));
			}
		}
	}
	
	for(i = 0; i < 3*atom_t; i++)
	{
		if(gsl_vector_get(eval_t, i) < 0.00001) {printf("Skip node %1i of target : %1.5f\n", i, gsl_vector_get(eval_t, i)); continue;}
		
		for(j = 0; j < 3*atom_t; j++)
		{
			for(k = 0; k < 3*atom_t; k++)
			{
				gsl_matrix_set(hess_t, j, k, gsl_matrix_get(hess_t, j, k) + gsl_matrix_get(evec_t, j, i)*gsl_matrix_get(evec_t, k, i)*gsl_vector_get(eval_t, i));
			}
		}
	}
	
	gsl_vector_free(eval);
	gsl_matrix_free(evec);
	gsl_vector_free(eval_t);
	gsl_matrix_free(evec_t);
	
	printf("Check 1\n");
	
	// Build mini-hessians (calculated mini_hess object is equal to mini_hess + mini_hess_t
	
	gsl_matrix *mini_hess = gsl_matrix_alloc(3*score, 3*score);
	gsl_matrix_set_all(mini_hess, 10);
	gsl_matrix *mini_hess_i = gsl_matrix_alloc(3*score, 3*score);
	gsl_matrix_set_all(mini_hess_i, 10);
	gsl_matrix *mini_hess_t = gsl_matrix_alloc(3*score, 3*score);
	gsl_matrix_set_all(mini_hess_t, 10);
	
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
	
	gsl_vector_free(del_conf);
	gsl_vector_free(copy_conf);
	gsl_vector_free(dummy_conf);
	gsl_vector_free(dummy_conf_t);
	gsl_matrix_free(mini_hess_i);
	gsl_matrix_free(mini_hess_t);
	
	printf("Energy from init to inter : %1.10f\nEnergy from target to inter : %1.10f\nDelta-E (init - target) : %1.10f\n", energy_ic, energy_tc, energy_ic - energy_tc);
}