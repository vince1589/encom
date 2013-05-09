#include "STeM.h"
void overlap_prob(struct pdb_atom *init,struct pdb_atom *targ,int atom,int atom_t,gsl_matrix *evec, int *align,gsl_vector *eval,struct pdb_atom *all_targ,int all_t);


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
	int i;
	int nbr_mode = 2;
	int mode = 7;
	int lig = 0;
	int nconn;

	float ligalign = 0; // Flag/valeur pour aligner seulement les résidus dans un cutoff du ligand, 0, one le fait pas... > 0... le cutoff
 	for (i = 1;i < argc;i++) {
 		if (strcmp("-i",argv[i]) == 0) {strcpy(file_name,argv[i+1]);--help_flag;}
 		if (strcmp("-h",argv[i]) == 0) {help_flag = 1;}
 		if (strcmp("-v",argv[i]) == 0) {verbose = 1;}
 		if (strcmp("-lig",argv[i]) == 0) {lig= 1;} 
 		if (strcmp("-t",argv[i]) == 0) {strcpy(check_name,argv[i+1]);help_flag = 0;}
 		if (strcmp("-m",argv[i]) == 0) {int temp;sscanf(argv[i+1],"%d",&temp);mode = temp;}
 		if (strcmp("-nm",argv[i]) == 0) {int temp;sscanf(argv[i+1],"%d",&temp);nbr_mode = temp;}
 		if (strcmp("-ieig",argv[i]) == 0) {strcpy(eigen_name,argv[i+1]);}
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
 	
 	gsl_vector *eval = gsl_vector_alloc(3*atom+10);
	gsl_matrix *evec = gsl_matrix_alloc (3*atom+10,3*atom+10);
 	
	
	load_eigen(eval,evec,eigen_name,3*atom);
	
	
	// Fit Using Eigenvalue
	if (nbr_mode == -1) {nbr_mode = 3*atom-7;}
	
	// Build vector of difference

		
	overlap_prob(strc_node,strc_node_t,atom,atom_t,evec,align,eval,strc_all_t,all_t);
	

   
  free(connect_h);
	free(connect_t);
	return(0);
 }
 
 void overlap_prob(struct pdb_atom *init,struct pdb_atom *targ,int atom,int atom_t,gsl_matrix *evec, int *align,gsl_vector *eval,struct pdb_atom *all_targ,int all_t) {
	int i,j;
 	int nb_mode = atom*3-6;
 	//nb_mode = 10;
 	int mode = 7;
		
	int alit = 0;// Value dans le vecteur diff
	
	int ngila[atom_t]; for (i=0;i<atom_t;++i) {ngila[i] = -1;  }
	for (i=0;i<atom;++i) {
		if (align[i] !=-1 ){ ngila[align[i]] = i; }
	}

	center_yes(init,targ,atom,atom_t, align); // Targ rotate, donc pas d'impact sur vecteurs
	rmsd_yes(targ,init,atom_t, ngila,all_targ,all_t);

	
 	gsl_vector *diff     = gsl_vector_alloc(3*atom);
 	
 	// Minimization
 
	// Construit vector des différences

 	for(i=0;i<atom;++i) {
 		j = align[i];
 		if (j == -1) {continue;}
 		//printf("Res: %s :: %s Num: %d :: %d\n",init[i].res_type,targ[j].res_type,init[i].res_number,targ[j].res_number);
 		gsl_vector_set(diff,alit*3+0,init[i].x_cord - targ[j].x_cord);
		gsl_vector_set(diff,alit*3+1,init[i].y_cord - targ[j].y_cord);
		gsl_vector_set(diff,alit*3+2,init[i].z_cord - targ[j].z_cord);
	//	printf("I:%d J:%d Alit:%d Init:(%f,%f,%f) Targ:(%f,%f,%f) Diff:(%f,%f,%f)\n",i,j,alit,init[i].x_cord,init[i].y_cord,init[i].z_cord,targ[j].x_cord,targ[j].y_cord,targ[j].z_cord,gsl_vector_get(diff,alit*3+0),gsl_vector_get(diff,alit*3+1),gsl_vector_get(diff,alit*3+2));
		++alit;
		
 	}
 	printf("Alit:%d\n",alit);
 	gsl_vector *B = gsl_vector_alloc((alit)*3);
 	for(i=0;i<alit;++i) {
 		//printf("I:%4d Vector:(%8.5f,%8.5f,%8.5f)\n",i,gsl_vector_get(diff,i*3+0),gsl_vector_get(diff,i*3+1),gsl_vector_get(diff,i*3+2));
 		gsl_vector_set(B,i*3+0,gsl_vector_get(diff,i*3+0));
 		gsl_vector_set(B,i*3+1,gsl_vector_get(diff,i*3+1));
 		gsl_vector_set(B,i*3+2,gsl_vector_get(diff,i*3+2));
 	}


	
	printf("Calculate probability of mode\n");
	double repar_func = 0;
	double beta = 10.0;
	for(j=0;j<nb_mode;++j) {
		//printf("%.4f / %.4f\n",-gsl_vector_get(eval,mode+j-1),beta);
 		repar_func += pow(2.71828,-gsl_vector_get(eval,mode+j-1)/beta);
 	
 	}
 	printf("Reparation function:%g\n",repar_func);
	
 	printf("ASSIGNING MODE\n");
	double overlap_proability = 0.0000000000;
	double total_prob = 0.00000000;
 	for(j=0;j<nb_mode;++j) {
 		float over = overlap(atom,mode+j-1,evec,diff,align);		
		//printf("J:%4d %10.6f %10.6f %10.6f\n",j+mode,gsl_vector_get(eval,mode+j-1),over,pow(2.71828,-gsl_vector_get(eval,mode+j-1)/beta)/repar_func);
		overlap_proability += over*pow(2.71828,-gsl_vector_get(eval,mode+j-1)/beta)/repar_func;
		total_prob += pow(2.71828,-gsl_vector_get(eval,mode+j-1)/beta)/repar_func;
	}
 	printf("Probability:%.10g\n",overlap_proability);
// 	printf("Tot:%.10g\n",total_prob);
 	gsl_vector_free(diff);

 }
