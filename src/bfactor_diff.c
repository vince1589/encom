#include "STeM.h"
void assign_bfactor(struct pdb_atom *strc, gsl_matrix *m,int atom,int lig) {
	int i;
	float max = 0;
	float min = 100000;
	float sum = 0;
	float count = 0;
	for(i=0;i< (int) m->size1 ;++i) {
		sum += gsl_matrix_get(m,i,i);
		count += 1.0;
		if (gsl_matrix_get(m,i,i) < min) {min = gsl_matrix_get(m,i,i);}
		if (gsl_matrix_get(m,i,i) > max) {max = gsl_matrix_get(m,i,i);}
	}
	sum /= count;
	float ratio = (max-min)/90.0;
	
	printf("Average:%.10f Min:%.4f Max:%.4f\n",sum,min,max);
	for(i=0;i<atom;++i) {
		int node = strc[i].node;
		strc[i].b_factor = (gsl_matrix_get(m,node,node)-min)/ratio+5.0;
	}

}

int count_eigen(char filename[100]) {
	int i = 0;
	int j,max=0;
	FILE *file;
 	file = fopen(filename,"r");
 	if (file == NULL) {return(1);}
 	char line[100];
	
 	while(fgets(line,100,file)) {
 	
 		
 		if (strncmp("2th",line,3) == 0) {break;}
 		sscanf(line,"%d\t%d",&i,&j);
 		if (i>max) {max=i;}
	}
	
	fclose(file);
	
	
	
	return(max);
}


int main(int argc, char *argv[]) {
	int i;
	char eigen_name[500] = "eigen.dat";
	char pca_name[500] = "pca_eigen.dat";
	char out_name[500] = "b_factor.pdb";
	int atom[2];
	char file_name[500];
	int lig = 0;
	int mode = 7;
	int nm = 10;
	for (i = 1;i < argc;i++) {
 		if (strcmp("-ieig1",argv[i]) == 0) {strcpy(eigen_name,argv[i+1]);}
		if (strcmp("-ieig2",argv[i]) == 0)  {strcpy(pca_name,argv[i+1]);}
		if (strcmp("-o",argv[i]) == 0)  {strcpy(out_name,argv[i+1]);}
		if (strcmp("-i",argv[i]) == 0) {strcpy(file_name,argv[i+1]);}
		if (strcmp("-lig",argv[i]) == 0) {lig= 1;}
		if (strcmp("-m",argv[i]) == 0)  {sscanf(argv[i+1],"%d",&mode);}
		if (strcmp("-nm",argv[i]) == 0)  {sscanf(argv[i+1],"%d",&nm);}
		
	}
	
	
	//***************************************************
 	//*													*
 	//*Build a structure contaning information on the pdb
 	//*													*
 	//***************************************************
 	
 	int all;
 	int nconn;
 	all = count_atom(file_name);
 	nconn = count_connect(file_name);
 	int verbose = 1;
 	if (verbose == 1) {printf("Connect:%d\n",nconn);}
 	
	if (verbose == 1) {printf("Assigning Structure\n\tAll Atom\n");}
	
	// Array qui comprend tous les connects
	
	int **connect_h=(int **)malloc(nconn*sizeof(int *)); 
    for(i=0;i<nconn;i++) { connect_h[i]=(int *)malloc(7*sizeof(int));}
    
    assign_connect(file_name,connect_h);
	
	// Assign tous les atoms
	int atom_strc;
	struct pdb_atom strc_all[all];
	atom_strc = build_all_strc(file_name,strc_all); // Retourne le nombre de Node
	
	if (verbose == 1) {printf("	Node:%d\n	Atom:%d\n",atom_strc,all);}

	check_lig(strc_all,connect_h,nconn,all);

	// Assign les Nodes
	
	if (verbose == 1) {printf("	CA Structure\n");}
	
	//atom = count_atom_CA_n(strc_all,all,super_node,lig);
	if (verbose == 1) {printf("	Node:%d\n",atom_strc);}
	struct pdb_atom strc_node[atom_strc];
	atom_strc = build_cord_CA(strc_all, strc_node,all,lig,connect_h,nconn);
	if (verbose == 1) {printf("	Assign Node:%d\n",atom_strc);}
	
	
	// Load eigen
	
	printf("My Eigen1:%s    My Eigen2:%s\n",eigen_name,pca_name);
	
	atom[0] = count_eigen(eigen_name);
	
	gsl_vector *eval = gsl_vector_alloc(atom[0]);
	gsl_matrix *evec = gsl_matrix_alloc (atom[0],atom[0]);
	
	load_eigen(eval,evec,eigen_name,atom[0]);
	
	atom[1] = count_eigen(pca_name);
	
	gsl_vector *evalpca = gsl_vector_alloc(atom[1]);
	gsl_matrix *evecpca = gsl_matrix_alloc (atom[1],atom[1]);
	
	load_eigen(evalpca,evecpca,pca_name,atom[1]);
	
	printf("Atom:   %d::%d\n",atom[0],atom[1]);
	if (atom[0] != atom[1]) {printf("I exited, car pas même nombre de Ca entre les eigenvecoctors.... la première structure de build_pca, est la strc de référence pour comparer\n");return(1);}
	
	// Inverse eigen
	
	printf("Inversing Matrix 1\n");
	gsl_matrix *k_inverse = gsl_matrix_alloc(atom[0]/3, atom[0]/3); /*Déclare et crée une matrice qui va être le pseudo inverse*/
	k_inverse_matrix_stem(k_inverse,atom[0]/3,eval,evec,mode,nm);
	
	printf("Inversing Matrix 2\n");
	gsl_matrix *k_inverse_pca = gsl_matrix_alloc(atom[1]/3, atom[1]/3); /*Déclare et crée une matrice qui va être le pseudo inverse*/
	k_inverse_matrix_stem(k_inverse_pca,atom[1]/3,evalpca,evecpca,mode,nm);
	
	for(i=0;i<atom[0]/3;++i) {
		float diff = gsl_matrix_get(k_inverse,i,i)-	gsl_matrix_get(k_inverse_pca,i,i);
		
		//if (sqrt(diff*diff) > 0.0001) {
		printf("I: %d Diff: %f -> %s %d %s\n",i,diff*1000000000.0,strc_node[i].res_type,strc_node[i].res_number,strc_node[i].chain);
		//}
		//gsl_matrix_set(k_inverse,i,i,diff*10000);
	}
	int j;
	char cc_name[50] = "cross_correlation_diff.dat";
	FILE *file;	
 	file = fopen(cc_name,"w");
	for(i=0;i<atom[0]/3;++i) {
		for(j=0;j<atom[0]/3;++j) {
			float diff = gsl_matrix_get(k_inverse,i,j)-	gsl_matrix_get(k_inverse_pca,i,j);
		
			//if (sqrt(diff*diff) > 0.0001) {
			fprintf(file,"%d %d %f\n",i,j,diff*1000000000.0);
			//}

		}
	}
	fclose(file);
	assign_bfactor(strc_all,k_inverse,all,lig);
	printf("I write strc\n");
	write_strc(out_name, strc_all,all,1.0);
	
	
	return(1);
}
