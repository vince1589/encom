#include "STeM.h"
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>

int main(int argc, char *argv[])
{
	int all; /*Nombre d'atomes dans pdb*/
	int atom; /*Nombre de carbone CA*/
	int help_flag = 0;
	int printev = 0;
	int prox = 0;
	int all_t; /*Nombre d'atomes dans pdb*/
	int atom_t; /*Nombre de carbone CA*/
	char file_name[500];
	char eigen_name[500] = "eigen.dat";
	char out_name[500] = "b_factor.pdb";
	char check_name[500];
	int verbose = 0;
	float ligalign = 0;
	
	int i,j,k,l;
	
	int nconn;
	int lig = 0;
	
	for (i = 1;i < argc;i++)
	{
		if (strcmp("-i",argv[i]) == 0) {strcpy(file_name,argv[i+1]);--help_flag;}
		if (strcmp("-o",argv[i]) == 0) {strcpy(out_name,argv[i+1]);}
		if (strcmp("-ieig",argv[i]) == 0) {strcpy(eigen_name,argv[i+1]);}
		if (strcmp("-h",argv[i]) == 0) {help_flag = 1;}
		if (strcmp("-v",argv[i]) == 0) {verbose = 1;}
		if (strcmp("-lig",argv[i]) == 0) {lig= 1;}
		if (strcmp("-t",argv[i]) == 0) {strcpy(check_name,argv[i+1]);help_flag = 0;}
		if (strcmp("-ligc",argv[i]) == 0) {float temp;sscanf(argv[i+1],"%f",&temp);ligalign = temp;}
		if (strcmp("-prox", argv[i]) == 0) { prox = 1; }
		if (strcmp("-prev", argv[i]) == 0) { printev = 1; }
	}
	
	if (help_flag == 1)
	{
		printf("****************************\nHelp Section\n-i\tFile Input (PDB)\n-o\tOutput Name Motion\n-ieig\tFile Name Eigen\n-v\tVerbose\n-sp\tSuper Node Mode (CA, N, C)\n-m\tMode\n-nm\tNombre de mode\n-lig\tTient compte duligand (sauf HOH)\n-prox\tAffiche la probabilite que deux CA puissent interagir\n-prev\tAffiche les directions principales du mouvement de chaque CA ponderees par leur ecart-type\n****************************\n");
		
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
	
	for(i=0;i<nconn;i++) { connect_h[i] = (int *)malloc(7*sizeof(int)); }
	
	assign_connect(file_name,connect_h);
	
	// Assign tous les atoms
	
	struct pdb_atom strc_all[all];
	
	atom = build_all_strc(file_name,strc_all); // Retourne le nombre de Node
	
	if (verbose == 1) { printf("	Node:%d\n	Atom:%d\n",atom,all); }
	
	check_lig(strc_all,connect_h,nconn,all);
	
	// Assign les Nodes
	
	if (verbose == 1) { printf("	CA Structure\n"); }
	
	//atom = count_atom_CA_n(strc_all,all,super_node,lig);
	if (verbose == 1) { printf("	Node:%d\n",atom); }
	
	struct pdb_atom strc_node[atom];
	
	atom = build_cord_CA(strc_all, strc_node,all,lig,connect_h,nconn);
	
	if (verbose == 1) { printf("	Assign Node:%d\n",atom); }
	
	/*
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
	
	if (verbose == 1) { printf("Score: %d/%d\n",score,atom); }
	
	//if (atom_t != atom) { printf("Not the same number of Node... Terminating\n");return(0); }
	
	if (score/atom < 0.8)
	{
		printf("Low Score... Will try an homemade alignement !!!\n");
		
		score = node_align_low(strc_node,atom,strc_node_t,atom_t,align);
	}
	
	printf("RMSD:%8.5f Score: %d/%d\n",sqrt(rmsd_no(strc_node,strc_node_t,atom, align)),score,atom);
	
	if (ligalign > 0)
	{
		score = node_align_lig(strc_node,atom,strc_node_t,atom_t,align,strc_all,all,strc_all_t,all_t,ligalign);
		
		printf("RMSD:%8.5f Score: %d/%d\n",sqrt(rmsd_no(strc_node,strc_node_t,atom, align)),score,atom);
	}
	*/
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
	
	for(k = 0; k < atom; k++)
	{
		bfactmes += strc_node[k].b_factor;
		lbfactmes += log(strc_node[k].b_factor);
		diffentmes += strc_node[k].entro;
		
		gsl_vector_set(lbfacs, k, log(strc_node[k].b_factor));
		gsl_vector_set(bfacs, k, strc_node[k].b_factor);
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
	
	printf("Differential entropy to b-factor correlation : %6.10f\nDifferential entropy to ln(b-factor) correlation : %6.10f\n", correl, lcorrel);
	
	if(prox == 1)
	{
		// Display the probability of all nodes to be from minrad to maxrad angstroms apart.
		
		double minrad = 3.2;
		
		double maxrad = 5.0; // Must be larger than minrad
		
		for(i = 1; i < atom; i++)
		{
			for(j = 0; j < i; j++)
			{
				// Initialize the objects needed for joint probability
				
				gsl_matrix *conj_var_inv = gsl_matrix_alloc(3,3);
				gsl_matrix_set_all(conj_var_inv, 0);
				
				gsl_vector *delta_pos = gsl_vector_alloc(3);
				gsl_vector_set_all (delta_pos, 0);
				
				double conj_dens = 0;
				
				conj_prob_init(&strc_node[i], &strc_node[j], conj_var_inv, delta_pos, &conj_dens);
				
				/*
				printf("Inverse conjugate covariance matrix :\n", i);
			
				for(k = 0; k < 3; k++)
				{
					for(int l = 0; l < 3; l++)
					{
						printf("%6.10f\t", gsl_matrix_get(conj_var_inv, k, l));
					}
					
					printf("\n");
				}
				
				printf("\n");
				*/
				
				// printf("\nConjugate density : %6.10f\n", conj_dens);
				
				// Calculates the maximum probability density of the evaluation range
				
				gsl_vector *evpos = gsl_vector_alloc(3);
				gsl_vector_set_all(evpos, 0);
				
				gsl_vector_set(evpos, 0, gsl_vector_get(delta_pos, 0));
				gsl_vector_set(evpos, 1, gsl_vector_get(delta_pos, 1));
				gsl_vector_set(evpos, 2, gsl_vector_get(delta_pos, 2));
				
				double delta_rad = sqrt(pow(gsl_vector_get(evpos, 0), 2) + pow(gsl_vector_get(evpos, 1), 2) + pow(gsl_vector_get(evpos, 2), 2));
				
				if(minrad < delta_rad)
				{
					if(maxrad < delta_rad)
					{
						gsl_vector_scale(evpos, 1 - maxrad / delta_rad);
					}
					else
					{
						gsl_vector_set_all(evpos, 0);
					}
				}
				else
				{
					gsl_vector_scale(evpos, 1 - minrad / delta_rad);
				}
				
				double argexp = 0;
				
				for(k = 0; k < 3; k++)
				{
					for(l = 0; l < 3; l++)
					{
						argexp += gsl_vector_get(evpos, k)*gsl_vector_get(evpos, l)*gsl_matrix_get(conj_var_inv, k, l);
					}
				}
				
				gsl_vector_free(evpos);
				
				argexp *= -0.5;
				
				double max_est = conj_dens*exp(argexp);
				
				// Evaluates the probability of nodes i and j to be from minrad to maxrad angstroms apart if max_est >= 0.0000000001
				
				double proxprob = 0;
				
				if(max_est >= 0.0000000001)
				{
					proxprob = proxim_prob(conj_var_inv, delta_pos, conj_dens, minrad, maxrad, 70);
				}
				
				gsl_vector_free(delta_pos);
				gsl_matrix_free(conj_var_inv);
				
				// Displays the probability density
				
				if(proxprob > 0)
				{
					printf("\nAtoms %2i and %2i :\t", i, j);
					
					printf("Distance : %6.11f\t\t", sqrt(pow(strc_node[j].x_cord - strc_node[i].x_cord, 2) + pow(strc_node[j].y_cord - strc_node[i].y_cord, 2) + pow(strc_node[j].z_cord - strc_node[i].z_cord, 2)));
					
					printf("P(%1.2f < delta-r < %1.2f) = %6.10f", minrad, maxrad, proxprob);
				}
			}
		}
		
		printf("\n");
	}
	
	/*
	for(i = 1; i < atom; i++)
	{
		printf("\nDistance between atoms %3i and %3i : %6.10f\t\t", i, i-1, sqrt(pow(strc_node[i].x_cord - strc_node[i-1].x_cord, 2) + pow(strc_node[i].y_cord - strc_node[i-1].y_cord, 2) + pow(strc_node[i].z_cord - strc_node[i-1].z_cord, 2)));
		printf("Probability of overlap : %6.10f\n", over_prob(&strc_node[i], &strc_node[i-1]));
	}
	*/
	
	// Rotate the protein
	
	/*
	gsl_matrix *rotest = gsl_matrix_alloc(3,3);
	gsl_matrix_set_all(rotest, 0);
	
	gsl_vector *rotaxis = gsl_vector_alloc(3);
	
	gsl_vector_set(rotaxis, 0, 1);
	gsl_vector_set(rotaxis, 1, 1);
	gsl_vector_set(rotaxis, 2, 1);
	
	rotate_matrix(rotest, 2.5, rotaxis);
	
	rotate_all(rotest, strc_node, atom);
	
	for(i = 1; i < atom; i++)
	{
		printf("\nDistance between atoms %3i and %3i : %6.10f\t\t", i, i-1, sqrt(pow(strc_node[i].x_cord - strc_node[i-1].x_cord, 2) + pow(strc_node[i].y_cord - strc_node[i-1].y_cord, 2) + pow(strc_node[i].z_cord - strc_node[i-1].z_cord, 2)));
		printf("Probability of overlap : %6.10f\n", over_prob(&strc_node[i], &strc_node[i-1]));
	}
	*/
	
	if(printev == 1)
	{
		printf("Eigenvectors with position :\n\n[");
		
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
			*/
			
			
			
			printf("[[%6.10f,%6.10f,%6.10f],[%6.10f,%6.10f,%6.10f]],", strc_node[i].x_cord, strc_node[i].y_cord, strc_node[i].z_cord, strc_node[i].global_evecs[0][0] * sqrt(strc_node[i].main_vars[0]), strc_node[i].global_evecs[1][0] * sqrt(strc_node[i].main_vars[0]), strc_node[i].global_evecs[2][0] * sqrt(strc_node[i].main_vars[0]));
			printf("[[%6.10f,%6.10f,%6.10f],[%6.10f,%6.10f,%6.10f]],", strc_node[i].x_cord, strc_node[i].y_cord, strc_node[i].z_cord, strc_node[i].global_evecs[0][1] * sqrt(strc_node[i].main_vars[1]), strc_node[i].global_evecs[1][1] * sqrt(strc_node[i].main_vars[1]), strc_node[i].global_evecs[2][1] * sqrt(strc_node[i].main_vars[1]));
			printf("[[%6.10f,%6.10f,%6.10f],[%6.10f,%6.10f,%6.10f]],", strc_node[i].x_cord, strc_node[i].y_cord, strc_node[i].z_cord, strc_node[i].global_evecs[0][2] * sqrt(strc_node[i].main_vars[2]), strc_node[i].global_evecs[1][2] * sqrt(strc_node[i].main_vars[2]), strc_node[i].global_evecs[2][2] * sqrt(strc_node[i].main_vars[2]));
		}
		
		printf("\b]\n\n");
	}
	
	/*
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
	*/
}


