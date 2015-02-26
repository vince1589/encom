#include "STeM.h"

void inv_mat_align(gsl_matrix *m,int nb_atom, gsl_vector *evl,gsl_matrix *evc,int mode,int nm,int *align,gsl_matrix *done);
int bb_sc(int node,struct pdb_atom *strc,int all);

int main(int argc, char *argv[]) {
	int all; /*Nombre d'atomes dans pdb*/
	int atom; /*Nombre de carbone CA*/
	int all_t; /*Nombre d'atomes dans pdb*/
	int atom_t; /*Nombre de carbone CA*/
 	int help_flag = 1;
 	char file_name[500];
 	char check_name[500];
 	int verbose = 0;
	int i,j;
	int lig = 0;
	int nconn;
	char eigen_name[500] = "eigen.dat";
	char eigen_name_two[500] = "eigen.dat";
	char init_list[500] = "init_list.dat";
	char targ_list[500] = "targ_list.dat";
	float rmsd_cutoff = 4.0;
	int constraint_flag = 1;
	int  ipos = -1;
	char posname[100] = "UNK";
	int posnum = -1;
	int npair = 2;
	int long_range = 0;
	int design = 0;
	int no_bb = 0;
	int non_local = 0;
	int max_pair = -1;
 	for (i = 1;i < argc;i++) {
 		if (strcmp("-i",argv[i]) == 0) {strcpy(file_name,argv[i+1]);--help_flag;}
 		if (strcmp("-h",argv[i]) == 0) {help_flag = 1;}
 		if (strcmp("-lr",argv[i]) == 0) {long_range = 1;}
 		if (strcmp("-nl",argv[i]) == 0) {non_local = 1;}
 		if (strcmp("-design",argv[i]) == 0) {design = 1;}
 		if (strcmp("-no_bb",argv[i]) == 0) {no_bb = 1;}
 		if (strcmp("-no_const",argv[i]) == 0) {constraint_flag = 0;}
 		if (strcmp("-v",argv[i]) == 0) {verbose = 1;}
		if (strcmp("-ieig",argv[i]) == 0) {strcpy(eigen_name,argv[i+1]);}
		if (strcmp("-teig",argv[i]) == 0) {strcpy(eigen_name_two,argv[i+1]);}
		if (strcmp("-pos",argv[i]) == 0) {strcpy(posname,argv[i+1]);int temp;sscanf(argv[i+2],"%d",&temp); posnum = temp;}
 		if (strcmp("-max_rmsd",argv[i]) == 0) {float temp;sscanf(argv[i+1],"%f",&temp);rmsd_cutoff = temp;}
 		if (strcmp("-size",argv[i]) == 0) {int temp;sscanf(argv[i+1],"%d",&temp); npair = temp-1;}
 		if (strcmp("-max_pair",argv[i]) == 0) {int temp;sscanf(argv[i+1],"%d",&temp); max_pair = temp;}
		if (strcmp("-t",argv[i]) == 0) {strcpy(check_name,argv[i+1]);help_flag = 0;}
 		
 	}
	 	
 	if (help_flag == 1) {
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
	if (atom > 5000) {printf("Too much node.... To fix, ask vincent.frappier@usherbrooke.ca\n");return(1);}
	if (verbose == 1) {printf("	Atom:%d\n",all);}
	check_lig(strc_all,connect_h,nconn,all);
	
	// Assign les Nodes
	
	if (verbose == 1) {printf("	CA Structure\n");}
	
	if (verbose == 1) {printf("	Node:%d\n",atom);}
	struct pdb_atom strc_node[atom];
	atom = build_cord_CA(strc_all, strc_node,all,lig,connect_h,nconn);
	
	if (verbose == 1) {printf("	Assign Node:%d\n",atom);}
	
	
	// Calcul surface en contact
	gsl_matrix *vcon = gsl_matrix_alloc(all,all);
	gsl_matrix_set_all(vcon,0);
	gsl_matrix *contact = gsl_matrix_alloc(atom, atom);
	gsl_matrix_set_all(contact,0);

	vcon_file_dom(strc_all,vcon,all);
	all_interaction_leStatium(strc_all,all, atom, contact,lig,vcon,strc_node);
	gsl_matrix_free(vcon);
	
	
		// Maintenant loader strc numero deux !!!!
	
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
	if (atom_t > 5000) {printf("Too much node.... To fix, ask vincent.frappier@usherbrooke.ca\n");return(1);}
	if (verbose == 1) {printf("	Atom:%d\n",all_t);}
	check_lig(strc_all_t,connect_t,nconn,all_t);
	
	// Assign les Nodes
	
	if (verbose == 1) {printf("	CA Structure\n");}
	
	if (verbose == 1) {printf("	Node:%d\n",atom_t);}
	struct pdb_atom strc_node_t[atom_t];

	atom_t = build_cord_CA(strc_all_t, strc_node_t,all_t,lig,connect_t,nconn);
	if (verbose == 1) {printf("	Assign Node:%d\n",atom_t);}
	
	// On va chercher les ALL dans la structure initial
	
	for(i=0;i<atom;++i) {
		if (strncmp(strc_node[i].res_type,"ALL",3)) {
			printf("I:%d\n",i);
		}
	
	}
	
	
	
	
	
	// All Scale
	char allAA[21][4];
	
	strcpy(allAA[0],"ALA");
	strcpy(allAA[1],"LEU");
	strcpy(allAA[2],"ILE");
	strcpy(allAA[3],"VAL");
	strcpy(allAA[4],"GLY");
	strcpy(allAA[5],"TRP");
	strcpy(allAA[6],"ARG");
	strcpy(allAA[7],"GLN");
	strcpy(allAA[8],"MET");
	strcpy(allAA[9],"CYS");
	strcpy(allAA[10],"SER");
	strcpy(allAA[11],"PHE");
	strcpy(allAA[12],"GLU");
	strcpy(allAA[13],"TYR");
	strcpy(allAA[14],"ASP");
	strcpy(allAA[15],"THR");
	strcpy(allAA[16],"LYS");
	strcpy(allAA[17],"ASN");
	strcpy(allAA[18],"HIS");
	strcpy(allAA[19],"PRO");
	

}

int bb_sc(int node,struct pdb_atom *strc,int all) {
	int max = all;
	int low = 0;
	int k;
	int i;
	// On va chercher inteligment la node
	for (k = 0 ;k<all;++k) {
		 i = low+(max-low)/2;
		//printf("I:%d Node:%d\n",i,strc[i].node);
		if (strc[i].node < node) {low = i;continue;}
		if (strc[i].node > node) {max = i;continue;}
		if (strc[i].node == node) {break;}
	}
	int bb = 0;
	int to = 0;
	for (k = i-5 ;k<i+5;++k) {
		if (k < 0) {continue;}
		if (k > all-1) {continue;}
		if (strc[k].node != node) {continue;}
		++to;
		if (strncmp(strc[k].atom_prot_type," N ",3) == 0) {++bb;}
		if (strncmp(strc[k].atom_prot_type," CA",3) == 0) {++bb;}
		if (strncmp(strc[k].atom_prot_type," C ",3) == 0) {++bb;}
		if (strncmp(strc[k].atom_prot_type," O ",3) == 0) {++bb;}
		//printf("I:%d ->%s<- %d %d\n",k,strc[k].atom_prot_type,to,bb);
	}
	if (bb == 4 && to == 4) {return(1);} else {return(0);}

}

void inv_mat_align(gsl_matrix *m,int nb_atom, gsl_vector *evl,gsl_matrix *evc,int mode,int nm,int *align,gsl_matrix *done) {
	//gsl_matrix *buffer = gsl_matrix_alloc(nb_atom, nb_atom); /*Matrix buffer a additionner*/
	//gsl_matrix_set_all (m, 0);
	int i,j,k;
	 
	for (i=0;i<nb_atom*3;++i)	{
		if (align[i/3] == -1) {continue;}
		//printf("I:%d / 3 = %d -> align -> %d\n",i, i/3,align[i/3]);
		for (j=0;j<nb_atom*3;++j) {
			if (align[j/3] == -1) {continue;}
			if (gsl_matrix_get(done, i,j) > 0) {continue;}
			//printf("I:%d J:%d\n",i,j);
			for (k=mode;k<mode+nm;++k) {
				if (k > int (evl->size-1)) {break;}
				if  (gsl_vector_get (evl, k) < 0.000001) {
					continue;
				}
				gsl_matrix_set(m, i,j,
					gsl_matrix_get(evc,i,k)*gsl_matrix_get(evc, j, k)/gsl_vector_get(evl, k) 
					+ gsl_matrix_get(m, i,j)
				);
				

			}
			gsl_matrix_set(done, i,j,1);
		}
	}
}

