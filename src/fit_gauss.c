#include "STeM.h"
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>

/*
void calculate_entropie(gsl_vector *entro,gsl_vector *evalt,gsl_matrix *evec,struct pdb_atom *strc,int atom, int zmodes)
{
	const double PI = 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421;
	const double beta = 0.0000001;
	int i,j;
	
	int nm = atom*3-zmodes;
	
	//int nm = 50;
	
	int node = 20;
	//nm = 4;
	// Find Q
	gsl_vector *eval = gsl_vector_alloc(nm);
	
	for (i=zmodes;i<nm+zmodes;++i)
	{
		gsl_vector_set(eval,i-6,gsl_vector_get(evalt,i));
	}
	
	(/*)
	double evaldiv = 0;
	
	// double evaldiv = gsl_vector_get(eval, (int) (nm-1)*(41+nm/170)/100);
	
	
	for(i = 0; i < nm; ++i)
	{
		evaldiv += gsl_vector_get(eval,i);
	}
	
	evaldiv /= 10*nm;
	
	
	for(i = 0; i < nm; ++i)
	{
		gsl_vector_set(eval,i,gsl_vector_get(eval,i)/evaldiv);
	}
	(end/*)
	
	for (node = 0;node<atom;++node)
	{
		gsl_matrix *GJ = gsl_matrix_alloc(3,nm+3);
		gsl_matrix_set_all(GJ,0);
		
		
		// Set Rx,Ry,Rz
		for(i=0;i<3;++i)
		{
			gsl_matrix_set(GJ,i,nm+i,1);
		}
	
		// Set evec
	
		for (j = 6;j<nm+6;++j)
		{
			for(i=0;i<3;++i)
			{
				gsl_matrix_set(GJ,i,j - 6,gsl_matrix_get(evec,i+node*3,j));
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
		 }(end/*)
		
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
	
		 printf("\n");(end/*)
		// Visualisation de l'exposant

		// Matrice K de transformation energetique
	
		gsl_matrix *K = gsl_matrix_alloc(nm,nm);
		gsl_matrix_set_all(K,0);
	
		for (i=0;i<nm;++i)
		{
			for (j=0;j<nm;++j)
			{
				gsl_matrix_set(K, i, j, beta * ( pow(gsl_vector_get(eval,0),2) * gsl_matrix_get(GJ,0,i+3) * gsl_matrix_get(GJ,0,j+3) + pow(gsl_vector_get(eval,1),2) * gsl_matrix_get(GJ,1,i+3) * gsl_matrix_get(GJ,1,j+3) + pow(gsl_vector_get(eval,2),2) * gsl_matrix_get(GJ,2,i+3) * gsl_matrix_get(GJ,2,j+3) ));
			}
		}
		
		gsl_matrix_free(GJ);
		
		// Matrice A d'energie-correlation (des covariances)
		
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
	
		/*
		printf("\nMatrice A\n");
		for (i=0;i<nm-3;++i)
		{
			for (j=0;j<nm-3;++j)
			{
				 printf("%6.10f ",gsl_matrix_get(A,i,j));
			}
			
			printf("\n");
		 }
		 (end/*)
		
		// make Cholesky decomposition of matrix A
		
		gsl_linalg_cholesky_decomp(A);
		
		// get cholesky matrix diagonal (product of diagonal elements will give sqrt(detA))
		
		gsl_vector *predet = gsl_vector_alloc(nm-3);
		
		for(i = 0; i < nm-3; ++i)
		{
			gsl_vector_set(predet, i, gsl_matrix_get(A,i,i));
		}
		
		gsl_sort_vector(predet);
		
		// invert matrix A
		
		gsl_linalg_cholesky_invert(A);
		
		/*
		 printf("\nMatrice I A\n");
		 for (i=0;i<nm-3;++i)
		 {
			 for(j=0;j<nm-3;++j)
			 {
			 	printf("%6.10f ",gsl_matrix_get(IA,i,j));
			 }
		
			 printf("\n");
		 }
		 (end/*)
		
		// get XYZ covariances
		
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
				KXX += 2 * gsl_matrix_get(K,i,nm-3) * gsl_matrix_get(K,j,nm-3)* gsl_matrix_get(A,i,j);
				KYY += 2 * gsl_matrix_get(K,i,nm-2) * gsl_matrix_get(K,j,nm-2)* gsl_matrix_get(A,i,j);
				KZZ += 2 * gsl_matrix_get(K,i,nm-1) * gsl_matrix_get(K,j,nm-1) * gsl_matrix_get(A,i,j);
				KXY += 2 * (gsl_matrix_get(K,i,nm-3) * gsl_matrix_get(K,j,nm-2) + gsl_matrix_get(K,i,nm-2) * gsl_matrix_get(K,j,nm-3)) * gsl_matrix_get(A,i,j);
				KXZ += 2 * (gsl_matrix_get(K,i,nm-3) * gsl_matrix_get(K,j,nm-1) + gsl_matrix_get(K,i,nm-1) * gsl_matrix_get(K,j,nm-3)) * gsl_matrix_get(A,i,j);
				KYZ += 2 * (gsl_matrix_get(K,i,nm-2) * gsl_matrix_get(K,j,nm-1) + gsl_matrix_get(K,i,nm-1) * gsl_matrix_get(K,j,nm-2)) * gsl_matrix_get(A,i,j);
			}
		}
		
		gsl_matrix_free(A);
		
		KXX -= gsl_matrix_get(K,nm-3,nm-3);
		KYY -= gsl_matrix_get(K,nm-2,nm-2);
		KZZ -= gsl_matrix_get(K,nm-1,nm-1);
		KXY -= gsl_matrix_get(K,nm-3,nm-2) + gsl_matrix_get(K,nm-2,nm-3);
		KXZ -= gsl_matrix_get(K,nm-3,nm-1) + gsl_matrix_get(K,nm-1,nm-3);
		KYZ -= gsl_matrix_get(K,nm-2,nm-1) + gsl_matrix_get(K,nm-1,nm-2);
		
		gsl_matrix_free(K);
		
		KXX *= -1;
		KYY *= -1;
		KZZ *= -1;
		KXY *= -1;
		KXZ *= -1;
		KYZ *= -1;
	
		if (det123 < 0) { det123 *= -1; }
	
		long double Dens_Fac = sqrt( pow(beta/PI,3) ) / det123;
		
		// multiply Dens_Fac with (product of eigenvalues)/sqrt(detA)
		
		for(i = 3; i < nm; ++i)
		{
			Dens_Fac *= sqrt(2*beta)*gsl_vector_get(eval,i)/gsl_vector_get(predet,i-3);
		}
		
		Dens_Fac *= gsl_vector_get(eval,0)*gsl_vector_get(eval,1)*gsl_vector_get(eval,2);
		
		// Calculate differential entropy
		
		double ConfEnt = 1.5 - log(Dens_Fac);
		
		gsl_vector_set(entro,node,ConfEnt);
		
		printf("Node:%d\n",node);
		
		//printf("\nConstants : \nq = %6.10f\ndet123 = %6.10f\ndetA = %6.10f\nDamping factor = %6.50Lf\nKXX = %6.10f\nKYY = %6.10f\nKZZ = %6.10f\nKXY = %6.10f\nKXZ = %6.10f\nKYZ = %6.10f\n\nDifferential entropy = %20.20Lf\n\n", q, det123, detA, Damping_Factor, KXX, KYY, KZZ, KXY, KXZ, KYZ, ConfEnt);
		printf("\nDifferential entropy = %6.10f\t\tDensity factor = %6.20Lf\n\nKXX = %6.10f\t\tKYY = %6.10f\t\tKZZ = %6.10f\t\tKXY = %6.10f\t\tKXZ = %6.10f\t\tKYZ = %6.10f\n\n", ConfEnt, Dens_Fac, KXX, KYY, KZZ, KXY, KXZ, KYZ);
	}
}
*/

int main(int argc, char *argv[])
{
	int all; /*Nombre d'atomes dans pdb*/
	int atom; /*Nombre de carbone CA*/
	int help_flag = 1;
	int all_t; /*Nombre d'atomes dans pdb*/
	int atom_t; /*Nombre de carbone CA*/
	char file_name[500];
	char eigen_name[500] = "eigen.dat";
	char out_name[500] = "b_factor.pdb";
	char check_name[500];
	int verbose = 0;
	float ligalign = 0;
	
	int i,j,k;
	
	int nconn;
	int lig = 0;
	
	for (i = 1;i < argc;i++) {
		if (strcmp("-i",argv[i]) == 0) {strcpy(file_name,argv[i+1]);--help_flag;}
		if (strcmp("-o",argv[i]) == 0) {strcpy(out_name,argv[i+1]);}
		if (strcmp("-ieig",argv[i]) == 0) {strcpy(eigen_name,argv[i+1]);}
		if (strcmp("-h",argv[i]) == 0) {help_flag = 1;}
		if (strcmp("-v",argv[i]) == 0) {verbose = 1;}
		if (strcmp("-lig",argv[i]) == 0) {lig= 1;}
		if (strcmp("-t",argv[i]) == 0) {strcpy(check_name,argv[i+1]);help_flag = 0;}
		if (strcmp("-ligc",argv[i]) == 0) {float temp;sscanf(argv[i+1],"%f",&temp);ligalign = temp;}
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
 	
 	gsl_vector *eval = gsl_vector_alloc(3*atom);
	gsl_matrix *evec= gsl_matrix_alloc (3*atom,3*atom);
	load_eigen(eval,evec,eigen_name,3*atom);
	
	gen_gauss(strc_node,evec,eval,atom,0.0000001);

	// Correlation
	
	double lbfactmes = 0;
	double bfactmes = 0;
	double diffentmes = 0;
	
	gsl_vector *bfacs = gsl_vector_alloc(atom);
	gsl_vector *lbfacs = gsl_vector_alloc(atom);
	gsl_vector *entro = gsl_vector_alloc(atom);
	
	for(k=0;k<atom;++k)
	{
		lbfactmes += log(strc_node[k].b_factor);
		gsl_vector_set(lbfacs, k, log(strc_node[k].b_factor));
		bfactmes += strc_node[k].b_factor;
		gsl_vector_set(bfacs, k, strc_node[k].b_factor);
		diffentmes += strc_node[k].entro;
		gsl_vector_set(entro, k, strc_node[k].entro);
	}
	
	bfactmes /= atom;
	lbfactmes /= atom;
	diffentmes /= atom;
	
	for(k=0;k<atom;++k)
	{
		gsl_vector_set(bfacs, k, gsl_vector_get(bfacs,k) - bfactmes);
		gsl_vector_set(lbfacs, k, gsl_vector_get(lbfacs,k) - lbfactmes);
		gsl_vector_set(entro, k, gsl_vector_get(entro,k) - diffentmes);
	}
	
	double correl = 0;
	double lcorrel = 0;
	bfactmes = 0;
	lbfactmes = 0;
	diffentmes = 0;
	
	// printf("Differential entropy vs b-factors (b,S) :\n");
	
	for(k=0;k<atom;++k)
	{
		correl += gsl_vector_get(entro,k)*gsl_vector_get(bfacs,k);
		lcorrel += gsl_vector_get(entro,k)*gsl_vector_get(lbfacs,k);
		bfactmes += gsl_vector_get(bfacs,k)*gsl_vector_get(bfacs,k);
		lbfactmes += gsl_vector_get(lbfacs,k)*gsl_vector_get(lbfacs,k);
		diffentmes += gsl_vector_get(entro,k)*gsl_vector_get(entro,k);
		
		// printf("{%6.10f,%6.10f},",gsl_vector_get(bfacs,k),gsl_vector_get(entro2,k));
	}
	
	// printf("\nDifferential entropy vs ln(b-factor)s (ln(b),S) :\n");
	/*
	for(k=0;k<atom;++k)
	{
		printf("{%6.10f,%6.10f},",gsl_vector_get(lbfacs,k),gsl_vector_get(entro2,k));
	}
	*/
	
	correl *= 1/(sqrt(diffentmes*bfactmes));
	
	lcorrel *= 1/(sqrt(diffentmes*lbfactmes));
	
	printf("Differential entropy to b-factor correlation : %6.10f\nDifferential entropy to ln(b-factor) correlation : %6.10f\n", correl,lcorrel);
	
	/*
	for(i = 1; i < atom; i++)
	{
		for(j = 0; j < i; j++)
		{
			printf("\nAtoms %3i and %3i :\n", i, j);
			
			over_prob(&strc_node[i], &strc_node[j]);
		}
	}
	*/
	
	for(i = 1; i < atom; i++)
	{
		printf("\nProbability of overlap between atoms %3i and %3i : %6.10f\n", i, i-1, over_prob(&strc_node[i], &strc_node[i-1]));
	}
	
	/*
	printf("Test :\n\n[");
	
	for(i = 0; i < atom; i++)
	{
		/*
		printf("Weighted eigenvector matrix for node %4i :\n", i);
		
		for(j = 0; j < 3; j++)
		{
			for(k = 0; k < 3; k++)
			{
				printf("%6.10f\t", strc_node[i].global_evecs[j][k]);
			}
			
			printf("\n");
		}
		
		printf("\n");
		
		
		printf("[[%6.10f,%6.10f,%6.10f],[%6.10f,%6.10f,%6.10f]],", strc_node[i].x_cord, strc_node[i].y_cord, strc_node[i].z_cord, strc_node[i].global_evecs[0][0], strc_node[i].global_evecs[1][0], strc_node[i].global_evecs[2][0]);
		printf("[[%6.10f,%6.10f,%6.10f],[%6.10f,%6.10f,%6.10f]],", strc_node[i].x_cord, strc_node[i].y_cord, strc_node[i].z_cord, strc_node[i].global_evecs[0][1], strc_node[i].global_evecs[1][1], strc_node[i].global_evecs[2][1]);
		printf("[[%6.10f,%6.10f,%6.10f],[%6.10f,%6.10f,%6.10f]],", strc_node[i].x_cord, strc_node[i].y_cord, strc_node[i].z_cord, strc_node[i].global_evecs[0][2], strc_node[i].global_evecs[1][2], strc_node[i].global_evecs[2][2]);
	}
	
	printf("]\n\n");
	*/
	int alit = 0;// Value dans le vecteur diff
	int ngila[atom_t]; for (i=0;i<atom_t;++i) {ngila[i] = -1;  }
	for (i=0;i<atom;++i) {
		if (align[i] !=-1 ){ ngila[align[i]] = i; }
	}

	center_yes(strc_node,strc_node_t,atom,atom_t, align); // Targ rotate, donc pas d'impact sur vecteurs
	rmsd_yes(strc_node_t,strc_node,atom_t, ngila,strc_all_t,all_t);
	
 	gsl_vector *diff     = gsl_vector_alloc(3*atom);
 	
 	// Minimization
	// Construit vector des différences

 	for(i=0;i<atom;++i) {
 		j = align[i];
 		if (j == -1) {continue;}
 		//printf("Res: %s :: %s Num: %d :: %d\n",strc_node[i].res_type,strc_node_t[j].res_type,strc_node[i].res_number,strc_node_t[j].res_number);
 		gsl_vector_set(diff,alit*3+0,strc_node[i].x_cord - strc_node_t[j].x_cord);
		gsl_vector_set(diff,alit*3+1,strc_node[i].y_cord - strc_node_t[j].y_cord);
		gsl_vector_set(diff,alit*3+2,strc_node[i].z_cord - strc_node_t[j].z_cord);
		//printf("I:%d J:%d Alit:%d Init:(%f,%f,%f) Targ:(%f,%f,%f) Diff:(%f,%f,%f)\n",i,j,alit,strc_node[i].x_cord,strc_node[i].y_cord,strc_node[i].z_cord,strc_node_t[j].x_cord,strc_node_t[j].y_cord,strc_node_t[j].z_cord,gsl_vector_get(diff,alit*3+0),gsl_vector_get(diff,alit*3+1),gsl_vector_get(diff,alit*3+2));
		++alit;
		
 	}
 	printf("Alit:%d\n",alit);
	
	write_strc("init.pdb",strc_node,atom);
	write_strc("target.pdb",strc_node_t,atom_t);
}


