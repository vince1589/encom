#include "STeM.h"

void calculate_entropie(gsl_vector *entro,gsl_vector *evalt,gsl_matrix *evec,struct pdb_atom *strc,int atom) {
	const double E = 2.718281828;
	const double PI = 3.14159265359;
	const double beta = 0.00000000001;
	int i,j;
	int nm = atom*3-6;
	int node = 20;
	//nm = 4;
	// Find Q
	gsl_vector *eval = gsl_vector_alloc(nm);
	double q = 1.0000;
	double e_prod = 1.000000;
	for (i=6;i<nm+6;++i) {
		gsl_vector_set(eval,i-6,gsl_vector_get(evalt,i));
		q *= sqrt(PI/(beta*pow(gsl_vector_get(eval,i-6),2)));
		
		e_prod *= pow(gsl_vector_get(eval,i-6),2);
	}
	printf("Q:%f\n",q);
	
	
	for (node = 0;node<atom;++node) {
		printf("Node:%d\n",node);
		gsl_matrix *GJ = gsl_matrix_alloc(3,nm+3);
		gsl_matrix_set_all(GJ,0);
		
		
		// Set Rx,Ry,Rz
		for(i=0;i<3;++i)
		{
			gsl_matrix_set(GJ,i,nm+i,1);
		}
	
		// Set evec
	
		for (j = 0;j<nm;++j)
		{
			for(i=0;i<3;++i)
			{
				gsl_matrix_set(GJ,i,j,gsl_matrix_get(evec,i+node*3,j));
			}
		}

		double div1 = gsl_matrix_get(GJ,0,0);
		double div2 = gsl_matrix_get(GJ,1,0);
		double div3 = gsl_matrix_get(GJ,2,0);
	
		for (i = 0;i<nm+3;++i)
		{
			gsl_matrix_set(GJ,1,i,gsl_matrix_get(GJ,1,i)-gsl_matrix_get(GJ,0,i)/div1*div2);		
			gsl_matrix_set(GJ,2,i,gsl_matrix_get(GJ,2,i)-gsl_matrix_get(GJ,0,i)/div1*div3);
		}
	
		// Print matrix
		// printf("\nElim_1\n");
	
		// for (j = 0 ; j<3;++j)
		// {
			// for (i = 0;i<nm+3;++i)
			// {
			// 	printf("%6.3f ",gsl_matrix_get(GJ,j,i));
			// }
		
			// printf("\n");
		// }
	
	
		// Elimine A dans ligne 1 et 3
	
		div1 = gsl_matrix_get(GJ,0,1);
		div2 = gsl_matrix_get(GJ,1,1);
		div3 = gsl_matrix_get(GJ,2,1);
	
		for (i = 0;i<nm+3;++i)
		{
			gsl_matrix_set(GJ,0,i,gsl_matrix_get(GJ,0,i)-gsl_matrix_get(GJ,1,i)/div2*div1);
		
			gsl_matrix_set(GJ,2,i,gsl_matrix_get(GJ,2,i)-gsl_matrix_get(GJ,1,i)/div2*div3);
		}
	
		// Print matrix
		// printf("\nElim_2\n");
		// for (j = 0 ; j<3;++j)
		// {
			// for (i = 0;i<nm+3;++i)
			// {
			// 	printf("%6.3f ",gsl_matrix_get(GJ,j,i));
			// }
		
			// printf("\n");
		// }
	
			// Elimine A dans ligne 1 et 2
	
		div1 = gsl_matrix_get(GJ,0,2);
		div2 = gsl_matrix_get(GJ,1,2);
		div3 = gsl_matrix_get(GJ,2,2);
	
		for (i = 0;i<nm+3;++i)
		{
			gsl_matrix_set(GJ,0,i,gsl_matrix_get(GJ,0,i)-gsl_matrix_get(GJ,2,i)/div3*div1);
		
			gsl_matrix_set(GJ,1,i,gsl_matrix_get(GJ,1,i)-gsl_matrix_get(GJ,2,i)/div3*div2);
		}
	
		// Print matrix
		/* printf("\nElim_3\n");
		 for (j = 0 ; j<3;++j)
		 {
			 for (i = 0;i<nm+3;++i)
			 {
			 	printf("%6.3f ",gsl_matrix_get(GJ,j,i));
			 }
		
			 printf("\n");
		 }*/
		
		double det123 = gsl_matrix_get(GJ,0,0)*gsl_matrix_get(GJ,1,1)*gsl_matrix_get(GJ,2,2);
	
		// Set 012 diagonal to 1
	
		div1 = gsl_matrix_get(GJ,0,0);
		div2 = gsl_matrix_get(GJ,1,1);
		div3 = gsl_matrix_get(GJ,2,2);
	
		for (i = 0;i<nm+3;++i)
		{
			gsl_matrix_set(GJ,0,i,gsl_matrix_get(GJ,0,i)/div1);
			gsl_matrix_set(GJ,1,i,gsl_matrix_get(GJ,1,i)/div2);
			gsl_matrix_set(GJ,2,i,gsl_matrix_get(GJ,2,i)/div3);	
		}
	
		// Print matrix
		/* printf("\nNormal\n");
		 for (j = 0 ; j<3;++j)
		 {
			 for (i = 0;i<nm+3;++i)
			 {
				 printf("%6.3f ",gsl_matrix_get(GJ,j,i));
			 }
		
			 printf("\n");
		 }
	
		 printf("\n");*/
		// Visualisation de l'exposant

		// Matrice A de covariance de A_3,A_4,..A_n,Rx,Ry,Rz)
	
		// Il faut avoir une forme -1/2 E Aij xi xj, donc il faut multipler par 2 à chaque fois (parce que j ai deja -beta)
	
		gsl_matrix *K = gsl_matrix_alloc(nm,nm);
		gsl_matrix_set_all(K,0);
	
		for (i=0;i<nm;++i)
		{
			for (j=0;j<nm;++j)
			{
				gsl_matrix_set(K, i, j, beta * ( pow(gsl_vector_get(eval,0),2) * gsl_matrix_get(GJ,0,i+3) * gsl_matrix_get(GJ,0,j+3) + pow(gsl_vector_get(eval,1),2) * gsl_matrix_get(GJ,1,i+3) * gsl_matrix_get(GJ,1,j+3) + pow(gsl_vector_get(eval,2),2) * gsl_matrix_get(GJ,2,i+3) * gsl_matrix_get(GJ,2,j+3) ));
			}
		}
	
		gsl_matrix *A = gsl_matrix_alloc(nm-3,nm-3);
		gsl_matrix_set_all(A,0);
	
		for (i=0;i<nm-3;++i)
		{
			for (j=0;j<nm-3;++j)
			{
				gsl_matrix_set(A, i, j, 2 * gsl_matrix_get(K,i,j));
			}
		}
	
		for (i=0;i<nm-3;++i)
		{
			gsl_matrix_set(A, i, i, gsl_matrix_get(A,i,i) + 2 * beta * pow(gsl_vector_get(eval,i+3),2));
		}
	
		/* printf("\nMatrice A\n");
		 for (i=0;i<nm-3;++i)
		 {
			 for (j=0;j<nm-3;++j)
			 {
				 printf("%6.10f ",gsl_matrix_get(A,i,j));
			 }
		
			 printf("\n");
		 }*/
	
		// Matrix inversion (IA)
		int s;
		gsl_matrix *IA = gsl_matrix_alloc(nm-3,nm-3);
		gsl_matrix_set_all(IA,0);
		gsl_permutation *perm = gsl_permutation_alloc (nm-3);
	
		// Make LU decomposition of matrix m
		gsl_linalg_LU_decomp (A, perm, &s);

		// Invert the matrix m
		gsl_linalg_LU_invert (A, perm, IA);
		
		gsl_matrix_scale(A,10.0);
		int stemp;
		gsl_permutation *permtemp = gsl_permutation_alloc (nm-3);
		
		gsl_linalg_LU_decomp (A, permtemp, &stemp);
		
		// Get derterminant
		double detA = gsl_linalg_LU_det(A, stemp);
	
/*		 printf("\nMatrice I A\n");
		 for (i=0;i<nm-3;++i)
		 {
			 for(j=0;j<nm-3;++j)
			 {
			 	printf("%6.10f ",gsl_matrix_get(IA,i,j));
			 }
		
			 printf("\n");
		 }*/
	
		double KXX = 0;
		double KYY = 0;
		double KZZ = 0;
		double KXY = 0;
		double KXZ = 0;
		double KYZ = 0;
	
		for (i=0;i<nm-3;++i)
		{
			for (j=0;j<nm-3;++j)
			{
				KXX += 2 * gsl_matrix_get(K,i,nm-3) * gsl_matrix_get(K,j,nm-3)* gsl_matrix_get(IA,i,j);
				KYY += 2 * gsl_matrix_get(K,i,nm-2) * gsl_matrix_get(K,j,nm-2)* gsl_matrix_get(IA,i,j);
				KZZ += 2 * gsl_matrix_get(K,i,nm-1) * gsl_matrix_get(K,j,nm-1) * gsl_matrix_get(IA,i,j);
				KXY += 2 * (gsl_matrix_get(K,i,nm-3) * gsl_matrix_get(K,j,nm-2) + gsl_matrix_get(K,i,nm-2) * gsl_matrix_get(K,j,nm-3)) * gsl_matrix_get(IA,i,j);
				KXZ += 2 * (gsl_matrix_get(K,i,nm-3) * gsl_matrix_get(K,j,nm-1) + gsl_matrix_get(K,i,nm-1) * gsl_matrix_get(K,j,nm-3)) * gsl_matrix_get(IA,i,j);
				KYZ += 2 * (gsl_matrix_get(K,i,nm-2) * gsl_matrix_get(K,j,nm-1) + gsl_matrix_get(K,i,nm-1) * gsl_matrix_get(K,j,nm-2)) * gsl_matrix_get(IA,i,j);
			}
		}
	
		KXX -= gsl_matrix_get(K,nm-3,nm-3);
		KYY -= gsl_matrix_get(K,nm-2,nm-2);
		KZZ -= gsl_matrix_get(K,nm-1,nm-1);
		KXY -= gsl_matrix_get(K,nm-3,nm-2) + gsl_matrix_get(K,nm-2,nm-3);
		KXZ -= gsl_matrix_get(K,nm-3,nm-1) + gsl_matrix_get(K,nm-1,nm-3);
		KYZ -= gsl_matrix_get(K,nm-2,nm-1) + gsl_matrix_get(K,nm-1,nm-2);
	
		KXX *= -1;
		KYY *= -1;
		KZZ *= -1;
		KXY *= -1;
		KXZ *= -1;
		KYZ *= -1;
	
		if (det123 < 0) { det123 *= -1; }
	
		double Damping_Factor = sqrt( pow(2,nm-3) * pow(beta,nm) / (pow(PI,3) * detA*pow(1.0/10.0,nm-3)) ) / det123;
	
		double ConfEnt = 0.5 *(3*log(2*PI) + 3 - log((pow(Damping_Factor,2)*pow(2*PI,3))));
	
		for(i = 0; i < nm; ++i)
		{
			Damping_Factor *= gsl_vector_get(eval,i);
		}
		gsl_vector_set(entro,node,ConfEnt);
		printf("\nConstants : \nq = %6.10f\ndet123 = %6.10f\ndetA = %6.10f\nDamping factor = %6.10f\nKXX = %6.10f\nKYY = %6.10f\nKZZ = %6.10f\nKXY = %6.10f\nKXZ = %6.10f\nKYZ = %6.10f\n\nDifferential entropy = %6.10f\n\n", q, det123, detA, Damping_Factor, KXX, KYY, KZZ, KXY, KXZ, KYZ, ConfEnt);
		break;
	}
}


int main(int argc, char *argv[]) {

	int all; /*Nombre d'atomes dans pdb*/
	int atom; /*Nombre de carbone CA*/
 	int help_flag = 1;
 	char file_name[500];
 	char eigen_name[500] = "eigen.dat";
 	char out_name[500] = "b_factor.pdb";
 	int verbose = 0;

	int i;

	int nconn;
	int lig = 0;

 	for (i = 1;i < argc;i++) {
 		if (strcmp("-i",argv[i]) == 0) {strcpy(file_name,argv[i+1]);--help_flag;}
 		if (strcmp("-o",argv[i]) == 0) {strcpy(out_name,argv[i+1]);}
 		if (strcmp("-ieig",argv[i]) == 0) {strcpy(eigen_name,argv[i+1]);}
 		if (strcmp("-h",argv[i]) == 0) {help_flag = 1;}
 		if (strcmp("-v",argv[i]) == 0) {verbose = 1;}
 		if (strcmp("-lig",argv[i]) == 0) {lig= 1;}    
 	}
 	
 	if (help_flag == 0) { } else {
 		printf("****************************\nHelp Section\n-i\tFile Input (PDB)\n-o\tOutput Name Motion\n-ieig\tFile Name Eigen\n-v\tVerbose\n-sp\tSuper Node Mode (CA, N, C)\n-m\tMode\n-nm\tNombre de mode\n-lig\tTient compte duligand (sauf HOH)\n****************************\n");
 		return(0); 
 	} 

 	
 	//***************************************************
 	//*													*
 	//*Build a structure contaning information on the pdb
 	//*													*
 	//***************************************************
 	
 	all = count_atom(file_name);
 	nconn = count_connect(file_name);
 	
 	if (verbose == 1) {printf("Connect:%d\n",nconn);}
 	
	if (verbose == 1) {printf("Assigning Structure\n\tAll Atom\n");}
	
	// Array qui comprend tous les connects
	
	int **connect_h=(int **)malloc(nconn*sizeof(int *)); 
    for(i=0;i<nconn;i++) { connect_h[i]=(int *)malloc(7*sizeof(int));}
    
    assign_connect(file_name,connect_h);
	
	// Assign tous les atoms
	
	struct pdb_atom strc_all[all];
	atom = build_all_strc(file_name,strc_all); // Retourne le nombre de Node
	
	if (verbose == 1) {printf("	Node:%d\n	Atom:%d\n",atom,all);}

	check_lig(strc_all,connect_h,nconn,all);

	// Assign les Nodes
	
	if (verbose == 1) {printf("	CA Structure\n");}
	
	//atom = count_atom_CA_n(strc_all,all,super_node,lig);
	if (verbose == 1) {printf("	Node:%d\n",atom);}
	struct pdb_atom strc_node[atom];
	atom = build_cord_CA(strc_all, strc_node,all,lig,connect_h,nconn);
	if (verbose == 1) {printf("	Assign Node:%d\n",atom);}
	
	
	
	//***************************************************
 	//*													*
 	//* Load eigenvector et eigenvalue					*
 	//*													*
 	//***************************************************
 	
 	if (verbose == 1) {printf("Loading Eigenvector\n");}
 	
 	gsl_vector *eval = gsl_vector_alloc(3*atom);
	gsl_matrix *evec= gsl_matrix_alloc (3*atom,3*atom);
	load_eigen(eval,evec,eigen_name,3*atom);
	
	gsl_vector *entro = gsl_vector_alloc(atom);
	
	calculate_entropie(entro,eval,evec,strc_node,atom);


	
}
