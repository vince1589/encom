#include "STeM.h"

void avg_dist(struct pdb_atom *init,int node,struct pdb_atom *targ,int *align,float scale,float *ret,int print_flag);

void assign_bfactor(struct pdb_atom *strc, int atom,struct pdb_atom *strc_all,int all) {
	int i;
	float max = 0;
	float min = 100000;
	float sum = 0;
	float count = 0;
	for(i=0;i< atom ;++i) {
		double val = strc[i].b_factor;
		
		sum += val;
		count += 1.0;
		if (val < min) {min = val;}
		if (val > max) {max =val;}
	}
	sum /= count;
	float ratio = (max-min)/90.0;
	
	printf("Average:%.10f Min:%.4f Max:%.4f\n",sum,min,max);
	for(i=0;i<all;++i) {
		
		strc_all[i].b_factor = (strc[strc_all[i].node].b_factor-min)/ratio+5.0;
	
	}

}

int main(int argc, char *argv[])
{
	int all; /*Number of atoms in the initial PDB*/
	int atom; /*Number of initial CAs*/
	int all_t; /*Number of atoms in the initial PDB*/
	int atom_t; /*Number of target CAs*/
	
 	int help_flag = 1;
 	char file_name[500];
 	char check_name[500];
 	char out_name[500] = "UNDEF";
 	int verbose = 0;
	
	int i;
	
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
	
 	if ((float)score/(float)atom > 0.0)
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
	float retu[4];
	
	avg_dist(strc_node,atom,strc_node_t,align,1.00,retu,1);
	
	if (strcmp(out_name,"UNDEF") != 0) {
		avg_dist(strc_node,atom,strc_node_t,align,1.00,retu,0);
		assign_bfactor(strc_node,atom,strc_all,all);
		write_strc(out_name,strc_all,all,1.0);
	}

}

void avg_dist(struct pdb_atom *init,int node,struct pdb_atom *targ,int *align,float scaler,float *ret,int print_flag) {
	
	ret[0] = 0;
	ret[1] = -1;
	ret[2] = 0;
	ret[3] = -1;
	
	gsl_vector *distances = gsl_vector_alloc(node);
	gsl_vector_set_all(distances, -1);
	
	gsl_matrix *aniso = gsl_matrix_alloc(node, 4);
	
	int i,j,k;	
			
	for(i = 0; i < node; ++i)	{
		int mat = align[i];
		init[i].b_factor = 0;
		for(j = 0; j < 3; j++) {
			if (init[i].main_vars[0] <= 0) {mat = -1; continue;}
			if (targ[mat].main_vars[0] <= 0) {mat = -1; continue;}
		}

		if (mat == -1) {continue;}

		gsl_matrix *anisou1 = gsl_matrix_alloc(3,3);
		gsl_vector *vars1 = gsl_vector_alloc(3);

		gsl_matrix *anisou2 = gsl_matrix_alloc(3,3);
		gsl_vector *vars2 = gsl_vector_alloc(3);

		// 		printf("Cov1 :\n");

		for(j = 0; j < 3; j++)	{
			for(k = 0; k < 3; k++) {

				gsl_matrix_set(anisou1, j, k, init[i].covar[j][k]);
				gsl_matrix_set(anisou2, j, k, targ[mat].covar[j][k] * scaler);
			}
	

			gsl_vector_set(vars1, j, init[i].main_vars[j]);
			gsl_vector_set(vars2, j, targ[mat].main_vars[j] * scaler);
		}



		gsl_vector_set(distances, i, cmp_gauss(anisou1, vars1, anisou2, vars2));



		gsl_matrix_set(aniso, i, 0, init[i].main_vars[2] / init[i].main_vars[0]);
		gsl_matrix_set(aniso, i, 2, init[i].main_vars[2] / init[i].main_vars[1]);
		gsl_matrix_set(aniso, i, 1, targ[mat].main_vars[2] / targ[mat].main_vars[0]);
		gsl_matrix_set(aniso, i, 3, targ[mat].main_vars[2] / targ[mat].main_vars[1]);
		init[i].b_factor = gsl_vector_get(distances,i);
		if (print_flag != 0) {
			printf("I:%4d %4.10f %s %d %s %s %d %s\n",i,gsl_vector_get(distances,i),init[i].res_type,init[i].res_number,init[i].chain,targ[mat].res_type,targ[mat].res_number,targ[mat].chain);
		}
	}
	gsl_sort_vector(distances);

	int count[2];
	count[0] = 0;
	count[1] = 0;
	for(i = 0; i < node; ++i)	{
		double val = gsl_vector_get(distances,i);
		if ( val < -0.000001) {continue;}
		count[0] += 1;
		if ((ret[1] < -0.00001) && ((float) i / (float) node > 0.5)) {
			ret[1] = val;
		}
		ret[0] += val;
		if ((float) i / (float) node < 0.9) {
			count[1] += 1;
			ret[2] += val;
			if ((ret[3] < -0.00001) && ((float) i / (float) node > 0.45)) {
				ret[3] = val;
			}
		}
	}

	ret[0] /= (float) count[0];
	ret[2] /= (float) count[1];
	printf("AVG_all:%.10f Median_all:%.10f AVG_90:%.10f Median_90:%.10f\n",ret[0],ret[1],ret[2],ret[3]);

}

