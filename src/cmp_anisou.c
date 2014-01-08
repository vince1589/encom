#include "STeM.h"




int main(int argc, char *argv[]) {
	int all; /*Nombre d'atomes dans pdb*/
	int atom; /*Nombre de carbone CA*/
	int all_t; /*Nombre d'atomes dans pdb*/
	int atom_t; /*Nombre de carbone CA*/
 	int help_flag = 1;
 	char file_name[500];
 	char check_name[500];
 	char out_name[500];
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
 	
 	if (verbose == 1) {printf("Connect:%d\n",nconn);}
 	
	if (verbose == 1) {printf("Assigning Structure\n\tAll Atom\n");}
	
	// Array qui comprend tous les connects
	
	int **connect_h=(int **)malloc(nconn*sizeof(int *)); 
    for(i=0;i<nconn;i++) { connect_h[i]=(int *)malloc(6*sizeof(int));}
    
    assign_connect(file_name,connect_h);

	// Assign tous les atoms
	
	struct pdb_atom strc_all[all];
	atom = build_all_strc(file_name,strc_all); // Retourne le nombre de Node
	if (atom > 800) {printf("Too much node.... To fix, ask vincent.frappier@usherbrooke.ca\n");return(1);}
	if (verbose == 1) {printf("	Atom:%d\n",all);}
	check_lig(strc_all,connect_h,nconn,all);
	
	// Assign les Nodes
	
	if (verbose == 1) {printf("	CA Structure\n");}
	
	if (verbose == 1) {printf("	Node:%d\n",atom);}
	struct pdb_atom strc_node[atom];
	atom = build_cord_CA(strc_all, strc_node,all,lig,connect_h,nconn);
	
	if (verbose == 1) {printf("	Assign Node:%d\n",atom);}
	printf("I found %d ANISOU\n",load_anisou(strc_node,file_name,atom));
	// Free Connect
		
	//for(i=0;i<nconn;i++) {printf("I:%d\n",i);free(connect_h[i]);}
	//free(connect_h);
	
	//Construit la structure a comparer
 	nconn = 0;
 	all_t = count_atom(check_name);
 	nconn = count_connect(check_name);
 	
 	if (verbose == 1) {printf("Connect:%d\n",nconn);}
 	
	if (verbose == 1) {printf("Assigning Structure\n\tAll Atom\n");}
	
	// Array qui comprend tous les connects
	
	int **connect_t=(int **)malloc(nconn*sizeof(int *)); 
    for(i=0;i<nconn;i++) { connect_t[i]=(int *)malloc(6*sizeof(int));}

    assign_connect(check_name,connect_t);
	
	// Assign tous les atoms
	
	struct pdb_atom strc_all_t[all_t];
	atom_t = build_all_strc(check_name,strc_all_t); // Retourne le nombre de Node
	if (atom_t > 800) {printf("Too much node.... To fix, ask vincent.frappier@usherbrooke.ca\n");return(1);}
	if (verbose == 1) {printf("	Atom:%d\n",all_t);}
	check_lig(strc_all_t,connect_t,nconn,all_t);
	
	// Assign les Nodes
	
	if (verbose == 1) {printf("	CA Structure\n");}
	
	if (verbose == 1) {printf("	Node:%d\n",atom_t);}
	struct pdb_atom strc_node_t[atom_t];

	atom_t = build_cord_CA(strc_all_t, strc_node_t,all_t,lig,connect_t,nconn);
	if (verbose == 1) {printf("	Assign Node:%d\n",atom_t);}
 printf("I found %d ANISOU\n",load_anisou(strc_node_t,check_name,atom_t));
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
	outlier_bfact(strc_node,atom,strc_node_t,atom_t,align,5);
	
	int it = 10;
	int k;
	for(k=0;k<it;++k) {
	
		double conj_dens12;
		gsl_vector *delr = gsl_vector_alloc(3);
		gsl_matrix *incov12 = gsl_matrix_alloc(3,3);
		gsl_vector *pos = gsl_vector_alloc(3);
		gsl_vector_set_all(pos,0);
		float prob_sum = 0;
		for(i=0;i<atom;++i) {
			int mat = align[i];
			if (mat == -1) {continue;}
			if (strc_node[i].main_vars[0] < 0) {continue;}
			if (strc_node_t[mat].main_vars[0] < 0) {continue;}
			if (conj_prob_init(&strc_node[i], &strc_node_t[mat], incov12,delr,&conj_dens12) == -1) {continue;}
			float dens = density_prob(incov12, delr, conj_dens12,pos);

		
			float prob = proxim_prob(incov12, delr, conj_dens12, 0.0, 4.0, 30);
			prob_sum += prob;
			/*printf("%s %d %s %s %d %s %f %f %f %f\n",strc_node[i].res_type,	strc_node[i].res_number,strc_node[i].chain,strc_node_t[mat].res_type,	strc_node_t[mat].res_number,strc_node_t[mat].chain,prob*1000000,dens*1000000,
				sqrt(pow(strc_node[i].x_cord-strc_node_t[mat].x_cord,2)+pow(strc_node[i].y_cord-strc_node_t[mat].y_cord,2)+pow(strc_node[i].z_cord-strc_node_t[mat].z_cord,2)),
				sqrt((pow(strc_node[i].main_vars[0]-strc_node_t[mat].main_vars[0],2)+pow(strc_node[i].main_vars[1]-strc_node_t[mat].main_vars[1],2)+pow(strc_node[i].main_vars[2]-strc_node_t[mat].main_vars[2],2))/3)
		
			);*/
			
		
		
		
		}
		printf("Sprob:%f\n",prob_sum);
	}
	
}


