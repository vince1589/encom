#include "STeM.h"




int main(int argc, char *argv[])
{
	int all; /*Number of atoms in the initial PDB*/
	int atom; /*Number of initial CAs*/
	int all_t; /*Number of atoms in the initial PDB*/
	int atom_t; /*Number of target CAs*/
	
 	int help_flag = 1;
 	char file_name[500];
 	char check_name[500];
 	char out_name[500];
 	int verbose = 0;
	
	int i,j,k;
	
	int lig = 0;
	int nconn;
	int print_flag = 0;
	float ligalign = 0; // Flag/valeur pour aligner seulement les résidus dans un cutoff du ligand, 0, one le fait pas... > 0... le cutoff
 	for (i = 1;i < argc;i++) {
 		if (strcmp("-i",argv[i]) == 0) {strcpy(file_name,argv[i+1]);--help_flag;}
 		if (strcmp("-h",argv[i]) == 0) {help_flag = 1;}
 		if (strcmp("-v",argv[i]) == 0) {verbose = 1;}

 		if (strcmp("-lig",argv[i]) == 0) {lig= 1;} 
 		if (strcmp("-t",argv[i]) == 0) {strcpy(check_name,argv[i+1]);help_flag = 0;}
 		if (strcmp("-o",argv[i]) == 0) {strcpy(out_name,argv[i+1]);++print_flag;}
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
	
	// Loads the ANISOUs in the initial node structure;
	
	printf("I found %d ANISOU\n",load_anisou(strc_node,file_name,atom));
	
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
	
	// Loads the ANISOUs in the target node structure
	
	printf("I found %d ANISOU\n",load_anisou(strc_node_t,check_name,atom_t));
	
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
	
// 	gsl_vector *iso_distances = gsl_vector_alloc(atom);
	
	gsl_vector *distances = gsl_vector_alloc(atom);
	
	gsl_matrix *aniso = gsl_matrix_alloc(atom, 4); // First anisotropy of initial, First anisotropy of target, Second anisotropy of initial, Second anisotropy of target
	
// 	gsl_matrix *angles = gsl_matrix_alloc(atom, 3); // cos(theta), cos(phi), cos(psi)
	
	double start = -4.0;
	
	double increment = -0.1;
	
	double scaler;
	
	double min_avg90 = 10000.0;
	
	double min_sd90;
	
	double min_med90 = 10000.0;
	
	double min_avgall = 10000.0;
	
	double min_sdall;
	
	double min_medall = 10000.0;
	
	double min_scaler;
	
	double min_scaler_med;
	
	double min_scaler90;
	
	double min_scaler90_med;
	
	int rep;
	
	gsl_vector *min_dist = gsl_vector_alloc(atom);
	
	gsl_matrix *min_aniso = gsl_matrix_alloc(atom, 4);
	
	for(rep = 0; rep < 81; rep++)
	{
		increment += 0.1;
		
		scaler = pow(10.0, start + increment);
		
// 		printf("%1.10f\n", scaler);
		
		for(i = 0; i < atom; ++i)
		{
			int mat = align[i];
			if (mat == -1) {continue;}
			
			for(j = 0; j < 3; j++)
			{
				if (strc_node[i].main_vars[0] <= 0) {continue;}
				if (strc_node_t[mat].main_vars[0] <= 0) {continue;}
			}
			
			gsl_matrix *anisou1 = gsl_matrix_alloc(3,3);
			gsl_vector *vars1 = gsl_vector_alloc(3);
			
			gsl_matrix *anisou2 = gsl_matrix_alloc(3,3);
			gsl_vector *vars2 = gsl_vector_alloc(3);
			
	// 		printf("Cov1 :\n");
			
			for(j = 0; j < 3; j++)
			{
				for(k = 0; k < 3; k++)
				{
	// 				printf("%6.10f\t", strc_node[i].anisou[j][k]);
					gsl_matrix_set(anisou1, j, k, strc_node[i].covar[j][k]);
					gsl_matrix_set(anisou2, j, k, strc_node_t[mat].covar[j][k] * scaler);
				}
				
	// 			printf("%6.10f\n", strc_node[i].anisou_vars[j]);
				gsl_vector_set(vars1, j, strc_node[i].main_vars[j]);
				gsl_vector_set(vars2, j, strc_node_t[mat].main_vars[j] * scaler);
			}
			
			gsl_vector_set(distances, i, cmp_gauss(anisou1, vars1, anisou2, vars2));
			
			gsl_matrix_set(aniso, i, 0, strc_node[i].main_vars[2] / strc_node[i].main_vars[0]);
			gsl_matrix_set(aniso, i, 2, strc_node[i].main_vars[2] / strc_node[i].main_vars[1]);
			gsl_matrix_set(aniso, i, 1, strc_node_t[mat].main_vars[2] / strc_node_t[mat].main_vars[0]);
			gsl_matrix_set(aniso, i, 3, strc_node_t[mat].main_vars[2] / strc_node_t[mat].main_vars[1]);
			
	// 		printf("%1i\t%6.10f\n", i, gsl_vector_get(distances, i));
		}
		
		gsl_permutation *perms = gsl_permutation_alloc(atom);
		
		gsl_sort_vector_index(perms, distances);
		
		double mean_all;
		
		double sd_all;
		
		double sd90;
		
		double avg90;
		
		double med_all;
		
		double med_90;
		
		int dum_atom = atom;
		
		for(i = 0; i < 2; i++)
		{
			if(i == 1)
			{
				dum_atom = 9 * dum_atom / 10;
			}
			
			double mean = 0.0;
			
			for(j = 0; j < dum_atom; j++)
			{
				mean += gsl_vector_get(distances,j);
			}
			
			mean /= dum_atom;
			
			double sd = 0.0;
			
			for(j = 0; j < dum_atom; j++)
			{
				sd += (gsl_vector_get(distances, j) - mean) * (gsl_vector_get(distances, j) - mean);
			}
			
			sd /= dum_atom;
			
			if(i == 0)
			{
				mean_all = mean; sd_all = sd;
				
				med_all = gsl_vector_get(distances, gsl_permutation_get(perms, dum_atom / 2));
			}
			else
			{
				avg90 = mean; sd90 = sd;
				
				med_90 = gsl_vector_get(distances, gsl_permutation_get(perms, dum_atom / 2));
			}
			
			/*
			if(i == 0) {printf("Mean value, total : %1.10f\nS.D. total : %1.10f\n", mean, sd);}
			else {printf("Mean value, 90pc : %1.10f\nS.D. 90pc : %1.10f\n", mean, sd);}
			*/
		}
		
		if(avg90 < min_avg90)
		{
// 			printf("%1.10f\t%1.10f\n", avg90, min_avg90);
			
			for(i = 0; i < atom; i++)
			{
				gsl_vector_set(min_dist, i, gsl_vector_get(distances, gsl_permutation_get(perms, i)));
				
				gsl_matrix_set(min_aniso, i, 0, gsl_matrix_get(aniso, gsl_permutation_get(perms, i), 0));
				gsl_matrix_set(min_aniso, i, 1, gsl_matrix_get(aniso, gsl_permutation_get(perms, i), 1));
				gsl_matrix_set(min_aniso, i, 2, gsl_matrix_get(aniso, gsl_permutation_get(perms, i), 2));
				gsl_matrix_set(min_aniso, i, 3, gsl_matrix_get(aniso, gsl_permutation_get(perms, i), 3));
			}
			
			min_avg90 = avg90;
			
			min_sd90 = sd90;
			
			min_scaler90 = scaler;
		}
		
		if(mean_all < min_avgall)
		{
			min_avgall = mean_all;
			
			min_sdall = sd_all;
			
			min_scaler = scaler;
		}
		
		if(med_all < min_medall)
		{
			min_medall = med_all;
			
			min_scaler_med = scaler;
		}
		
		if(med_90 < min_med90)
		{
			min_med90 = med_90;
			
			min_scaler90_med = scaler;
		}
		
		gsl_permutation_free(perms);
	}
	
	if(min_avg90 == 10000.0)
	{
		min_avg90 = 0.0;
		min_sd90 = 0.0;
		min_scaler90 = 0.0;
	}
	
	if(min_med90 == 10000.0)
	{
		min_med90 = 0.0;
		min_scaler90_med = 0.0;
	}
	
	if(min_avgall == 10000.0)
	{
		min_avgall = 0.0;
		min_sdall = 0.0;
		min_scaler = 0.0;
	}
	
	if(min_medall == 10000.0)
	{
		min_medall = 0.0;
		min_scaler_med = 0.0;
	}
	
	printf("Min avg value, total : %1.10f\nS.D. total : %1.10f\nMedian, total : %1.10f\n", min_avgall, min_sdall, min_medall);
	printf("Min scaler, all : %1.10f\n", min_scaler);
	printf("Min scaler, median, all : %1.10f\n", min_scaler_med);
	printf("Min avg value, 90pc : %1.10f\nS.D. 90pc : %1.10f\nMedian, 90pc : %1.10f\n", min_avg90, min_sd90, min_med90);
	printf("Min scaler, 90pc : %1.10f\n", min_scaler90);
	printf("Min scaler, median, 90pc : %1.10f\n", min_scaler90_med);
	
		
	printf("Data :\nNode\tDist\tAniso_1I\tAniso_1T\tAniso_2I\tAniso_2T\n");
	
	for(i = 0; i < atom; i++)
	{
		printf("%1i\t%1.5f\t%3.4f\t%3.4f\t%3.4f\t%3.4f\n", i,
			gsl_vector_get(min_dist, i),
			gsl_matrix_get(min_aniso, i, 0),
			gsl_matrix_get(min_aniso, i, 1),
			gsl_matrix_get(min_aniso, i, 2),
			gsl_matrix_get(min_aniso, i, 3)
			);
	}
}


