#include "STeM.h"
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>

void outlier_bfact(struct pdb_atom *strc, int atom,int next);

int main(int argc, char *argv[])
{
	int all; /*Nombre d'atomes dans pdb*/
	int atom; /*Nombre de carbone CA*/
	int help_flag = 0;


	char file_name[500];
	char eigen_name[500] = "eigen.dat";


	int verbose = 0;

	
	int i;
	
	int nconn;
	int lig = 0;
	double beta = 0.000001;
	for (i = 1;i < argc;i++)
	{
		if (strcmp("-i",argv[i]) == 0) {strcpy(file_name,argv[i+1]);--help_flag;}
		if (strcmp("-ieig",argv[i]) == 0) {strcpy(eigen_name,argv[i+1]);}
		if (strcmp("-h",argv[i]) == 0) {help_flag = 1;}
		if (strcmp("-v",argv[i]) == 0) {verbose = 1;}
		if (strcmp("-b",argv[i]) == 0) {float temp;sscanf(argv[i+1],"%f",&temp);beta = temp;}


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
	
	if (verbose == 1) { printf("	Loading Anisou\n"); }
	printf("I found %d ANISOU\n",load_anisou(strc_node,file_name,atom));
	//***************************************************
	//*													*
	//* Load eigenvector et eigenvalue					*
	//*													*
	//***************************************************
	
	if (verbose == 1) {printf("Loading Eigenvector\n");}
	
	gsl_vector *eval = gsl_vector_alloc(3*atom);
	
	gsl_matrix *evec= gsl_matrix_alloc (3*atom,3*atom);
	
	load_eigen(eval,evec,eigen_name,3*atom);
	if (verbose == 1) {printf("Gen Gauss with beta = %f\n",beta);}
	gen_gauss(strc_node,evec,eval,atom,beta);

	printf("Eval correlation\n");
	outlier_bfact(strc_node,atom,atom/10);
	

	// General direction vector comparison
	
	float esens[3];
	float psens[3];
	
	
	
	gsl_vector * gsens_evec = gsl_vector_alloc(atom*3);
	gsl_vector * gsens_exp = gsl_vector_alloc(atom*3);
	gsl_vector_set_all(gsens_evec,0);
	for(i=0;i<atom;++i) {
		int l;
		float max = 0;
		for (l=0;l<3;++l) {
			esens[l] = 1.0;
			psens[l] = -1.0;
		}
		// For each atom we want to try all direction of 
	
		
		for (l=0;l<9;++l) {
			gsl_vector *esum = gsl_vector_alloc(3);
			gsl_vector *psum = gsl_vector_alloc(3);
			gsl_vector_set_all(esum,0);
			gsl_vector_set_all(psum,0);
			int k;
			for (k=0;k<3;++k) {
				gsl_vector_set(psum,0,gsl_vector_get(psum,0)+strc_node[i].main_vars[k]*strc_node[i].global_evecs[k][0]*psens[k]);
				gsl_vector_set(psum,1,gsl_vector_get(psum,1)+strc_node[i].main_vars[k]*strc_node[i].global_evecs[k][1]*psens[k]);
				gsl_vector_set(psum,2,gsl_vector_get(psum,2)+strc_node[i].main_vars[k]*strc_node[i].global_evecs[k][2]*psens[k]);
			
				gsl_vector_set(esum,0,gsl_vector_get(esum,0)+strc_node[i].main_vars[k]*strc_node[i].avec[k][0]*esens[k]);
				gsl_vector_set(esum,1,gsl_vector_get(esum,1)+strc_node[i].main_vars[k]*strc_node[i].avec[k][1]*esens[k]);
				gsl_vector_set(esum,2,gsl_vector_get(esum,2)+strc_node[i].main_vars[k]*strc_node[i].avec[k][2]*esens[k]);
			
			}
				float elen = vector_lenght(esum,3);
				float plen = vector_lenght(psum,3);

				float a = 0.0;
				float b = 0.0;
				float c = 0.0; 		

			 	for (k=0;k<3;++k) {
			 		float x = gsl_vector_get(esum,k)/elen;
			 		float z = gsl_vector_get(psum,k)/plen;
			 		a += x*z;
			 		b += x*x;
			 		c += z*z;
				}

			 	a = sqrt(a*a);
			 
			 	//printf("I:%d %f %f %f Over:%f\n",i,psens[0],psens[1],psens[2],(a/sqrt(b*c)));
			 	if (a/sqrt(b*c) > max) {
			 		for (k = 0;k<3;++k) {
			 			gsl_vector_set(gsens_evec,i*3+k,gsl_vector_get(psum,k));
			 			gsl_vector_set(gsens_exp,i*3+k,gsl_vector_get(esum,k));
			 		}
			 		max = a/sqrt(b*c);
			 	}
			 	psens[0] += 2.0;
			 	if (psens[0] > 1.0) {
			 		psens[0] = -1.0;
			 		psens[1] += 2.0;
			 	}
			 	if (psens[1] > 1.0) {
			 		psens[1] = -1.0;
			 		psens[2] += 2.0;
			 	}
		 	}
		 	//printf("I:%d Max:%f\n",i,max);
	}
	// Petit hack pour utiliser le bootstrap, on va mettre les vecteur dans les evals
	for(i=0;i<atom;++i) {
		int k;
		for(k=0;k<3;++k) {
			strc_node[i].aval[k] = gsl_vector_get(gsens_exp,i*3+k);
			strc_node[i].main_vars[k] = gsl_vector_get(gsens_evec,i*3+k);
		}
	}
	printf("Evec Corr\n");
	outlier_bfact(strc_node,atom,atom/10);	
}

void outlier_bfact(struct pdb_atom *strc, int atom,int next) {
	int it;
	int it_max = atom;
	int k;
	int i;
	int remove[next];
	for (i = 0;i<next;++i) {
		remove[i] = -1;
	}
	int l = 0;
	for (l=0;l<next;++l) {
		int mind = -1;
		float max = -2.0;
		for (it=0;it<it_max;++it) {
			
			
			float avg[2];
	
			avg[0] = 0.0;
			avg[1] = 0.0;
	
			float besl[4];
			besl[0] = 0.0;
			besl[1] = 0.0;

			for(i = 0;i<atom;++i) {
				if (i == it) {continue;}
				int nf = 0;
				for (k = 0;k<next;++k) {
					if (i == remove[k]) {++nf;}
				}
				if (nf != 0) {continue;}
				for (k = 0 ;k<3;++k) {

					avg[0] += strc[i].main_vars[k]/(atom*3);
					avg[1] += strc[i].aval[k]/(atom*3);
					besl[0] += strc[i].aval[k]*strc[i].aval[k];
					besl[1] += strc[i].aval[k]*strc[i].main_vars[k];
				}
			}
			float cor[3];
			cor[0] = 0;
			cor[1] = 0;
			cor[2] = 0;
		
			for(i = 0;i<atom;++i) {
				int k;
				if (i == it) {continue;}
				int nf = 0;
				for (k = 0;k<next;++k) {
					if (i == remove[k]) {++nf;}
				}
				if (nf != 0) {continue;}
				for (k = 0 ;k<3;++k) {
					cor[0] += (strc[i].main_vars[k]-avg[0])*(strc[i].aval[k]-avg[1]);
					cor[1] += (strc[i].main_vars[k]-avg[0])*(strc[i].main_vars[k]-avg[0]);
					cor[2] += (strc[i].aval[k]-avg[1])*(strc[i].aval[k]-avg[1]);
				}
		
			}
	
			// Find y = ax + b
	
	
			if (max < cor[0]/sqrt(cor[1])/sqrt(cor[2])) {
				max = cor[0]/sqrt(cor[1])/sqrt(cor[2]);
				mind = it;
			//	printf("IT:%d Eval Cor:%f  ",it,cor[0]/sqrt(cor[1])/sqrt(cor[2]));
		
			//	float slope = (besl[1] - avg[1]*atom*3 * avg[0]) / (besl[0] - avg[1]*atom*3 * avg[1]);
			//	float inter = avg[0] - slope * avg[1];
			//	printf("A = %f B = %f Xm = %f Ym = %f\n",slope,inter	,avg[1],avg[0]);
			}
		}
		printf("Max:%f Mind:%d\n",max,mind);
		remove[l] = mind;
	}
	

}
