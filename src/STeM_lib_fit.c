#include "STeM.h"
#include <gsl/gsl_matrix.h>
#include <math.h>
#include <gsl/gsl_eigen.h>
#include <time.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>
 
 
int node_align(struct pdb_atom *strc,int atom,struct pdb_atom *strc_t,int atom_t, int *align) {
	int i,j,l;
	int score;


	for(i=0;i<atom;++i) {
		align[i] = -1;
	}
		

		
	//Retourne un score
	score = 0;
	

	for(i=0;i<atom;++i) {
		for(j=0;j<atom_t;++j) {
			
			if (
				(strcmp(strc[i].res_type,strc_t[j].res_type) == 0) && 
				(strc[i].res_number == strc_t[j].res_number) && 
				(strcmp(strc[i].chain,strc_t[j].chain) == 0)) {
				
				

				if (align[i] == -1) {++score;align[i] = j;}
				
				
			}
		}
	}
	
	// Look if node match two times
	score = 0;
	for(i=0;i<atom;++i) {
		for(j=0;j<atom;++j) {
			if (i == j) {continue;}
			if (align[i] == align[j]) {align[j] = -1;}
		}
		if (align[i] != -1) {++score;}
	
	}
	for(i=0;i<atom;++i) {
			if (align[i]==-1) {
				//printf("%3d :: %3d   %s::%s  %3d::%3d\n",i,align[i],strc[i].res_type,"xxx",strc[i].res_number,0);
				continue;
			}
			//printf("%3d :: %3d   %s::%s  %3d::%3d\n",i,align[i],strc[i].res_type,strc_t[align[i]].res_type,strc[i].res_number,strc_t[align[i]].res_number);
			j = align[i];
			//printf("-%s- -%d- %s :: -%s- -%d- %s\n",strc[i].res_type,strc[i].res_number,strc[i].chain,strc_t[j].res_type,strc_t[j].res_number,strc_t[j].chain);
		}
	
	return(score);
	

}

int node_align_lig(struct pdb_atom *strc,int atom,struct pdb_atom *strc_t,int atom_t, int *align,struct pdb_atom *strc_all,int all,struct pdb_atom *strc_all_t,int all_t,float cutoff) {
	
	int i,j,k = 0;
	int score = 0;
	// Fonction qui ne garde que les aa aligne autour du ligand de la TARG
	
	float coord[200][3]; // Array qui va comprendre les coordone du ligand
	
	// On trouve les ligands et coordone storer dans coord
	
	for (i=0;i<atom_t;++i) {
			if (strc_t[i].atom_type == 3) {
				printf("Lig:%d %s %d %s\n",i,strc_t[i].res_type,strc_t[i].res_number,strc_t[i].chain);
				for (j=0;j<all_t;++j) {
					// Parce que je suis pas certains que les type des strc_all sont bien defini
					// On cherche les nodes correspondante pour trouver coord des lig, pas meilleur pratique
					
					//printf("%d == ? %d\n",strc_all_t[j].node , strc_t[i].node);
					
					if (strc_all_t[j].node == strc_t[i].node) {
					//	printf("K:%d  (%.3f,%.3f,%.3f)\n",k, strc_all_t[j].x_cord, strc_all_t[j].y_cord, strc_all_t[j].z_cord);
						coord[k][0] = strc_all_t[j].x_cord;
						coord[k][1] = strc_all_t[j].y_cord;
						coord[k][2] = strc_all_t[j].z_cord;
						++k;
					}
				}
			}
	} 
	
	// On trouve les node proche du ligand dans TARG
	int kept_node[atom_t]; // Node qui match
	int l = 0;
	for (i=0;i<all_t;++i) {
		//printf("I:%d K:%d\n",i,k);
		
		// Regarde si on a pas deja ajoute node
		int flag = 0;
		for (j=0;j<l;++j) {
			//printf("J:%d I:%d %d ==? %d\n",j,i,kept_node[j] , strc_all_t[i].node);
			if (kept_node[j] == strc_all_t[i].node) {++flag;break;}
		}
		if (flag != 0) {continue;}
		// Ajoute node de targe si dans cutoff du ligand
		for (j=0;j<k;++j) {
			float dist = (coord[j][0] - strc_all_t[i].x_cord)*(coord[j][0] - strc_all_t[i].x_cord);
			dist += (coord[j][1] - strc_all_t[i].y_cord)*(coord[j][1] - strc_all_t[i].y_cord);
			dist += (coord[j][2] - strc_all_t[i].z_cord)*(coord[j][2] - strc_all_t[i].z_cord);
			/*printf("K:%d Dist = %.3f (%.3f,%.3f,%.3f) (%.3f,%.3f,%.3f)\n",k,dist,coord[j][0],coord[j][1],coord[j][2],strc_all_t[i].x_cord,strc_all_t[i].y_cord,strc_all_t[i].z_cord);
			printf("%f < %f\n",dist,cutoff*cutoff);*/
			if (dist < cutoff*cutoff) {
					//printf("I add I:%d L:%d Node:%d %s %d %s\n",i,l,strc_all_t[i].node,strc_all_t[i].res_type,strc_all_t[i].res_number,strc_all_t[i].chain);
					kept_node[l] = strc_all_t[i].node;
					++l;
					break;
			}
		}
	
	}
	
	
	// On regarde pour dans align pour garder seulement ceux qui ont matche
	
	for (i=0;i<atom;++i) {
		int flag = 0;
		if (align[i] == -1) {continue;}
		for (j=0;j<l;++j) {
			if (strc_t[align[i]].node == kept_node[j]) {
				//printf("%d :: %d Init: %s%d%s :: %s%d%s :Targ\n",i,align[i],strc[i].res_type,strc[i].res_number,strc[i].chain,strc_t[align[i]].res_type,strc_t[align[i]].res_number,strc_t[align[i]].chain);
				++flag;
			}
		}
		if (flag == 0) {
			align[i] = -1;
		} else {
			++score;
		}
	}
	
	return(score);

}

int next_cons(struct pdb_atom *strc,int atom,struct pdb_atom *strc_t,int atom_t,int i,int j) {
	int k = 1;
	int cons = 0;
	while(1) {
		
		if (strcmp(strc[i+k].res_type,strc_t[i+j+k].res_type) == 0) {++cons;++k;} else {break;}
	
	}
	//if (cons > 2) {printf("I:%d %d Cons:%d\n",i,j,cons);}
	return(cons);
}

int node_align_low(struct pdb_atom *strc,int atom,struct pdb_atom *strc_t,int atom_t, int *align) {
	int i,j,l;
	int score;
	int num;
	int cons;
	//printf("Atom:%d Atom_t:%d\n",atom,atom_t);
	for(i=0;i<atom;++i) {
		align[i] = -1;
	}
		

		
		//Retourne un score
		score = 0;
		for(j=-atom_t;j<atom_t;++j) {
			cons = 0;
			
			for(i=0;i<atom;++i) {
				if (i+j < 0) {continue;}
				if (i+j > atom_t-1) {continue;}
				
				// Élimine les Ligands de l'alignement
				
				if (strc[i].atom_type == 3) {continue;}
				if (strc_t[i+j].atom_type == 3) {continue;}
				
				// On veut des clusters de séquences identiques	
				
			//	if (strcmp(strc[i].res_type,strc_t[i+j].res_type) == 0 || ( 3 < next_cons(strc,atom,strc_t,atom_t,i,j))) {
				if (strcmp(strc[i].res_type,strc_t[i+j].res_type) == 0) {
					++cons;
					
					if (cons > 8) {

						align[i] = i+j;
						if (cons == 9) {
							for (l=0;l<9;++l) {
								align[i-l] = i+j-l;
							} 
						}
					}
				} else {cons = 0;}
				
			}

		}
		
		// Aligne si single point mutation ou si deux mutations de suites
		
		for(i=0;i<atom;++i) {
			if (align[i] == -1) {
			
				// Cas ou une mutation autour d'une sequence qui fit (ABCDEF est égale a ABGDEF)
				if (i > 0 && i<atom-1 && align[i-1] != -1 && align[i+1] != -1 && align[i+1] - align[i-1] == 2) {
					align[i] = align[i+1] -1;
				}
				// Cas ou deux mutations consécutives
				if (i > 0 && i<atom-2 && align[i-1] != -1 && align[i+2] != -1 && align[i+2] - align[i-1] == 3 && align[i+1] == -1) {
					align[i] = align[i+2] -2;
				}
				
							
			}
		}
		
		
		// Look if node match two times
		score = 0;
		for(i=0;i<atom;++i) {
			for(j=0;j<atom;++j) {
				if (i == j) {continue;}
				if (align[i] == align[j]) {align[j] = -1;}
			}
			if (align[i] != -1) {++score;}
	
		}
		
		
	for(i=0;i<atom;++i) {
			if (align[i]==-1) {
				//printf("%3d :: %3d   %s::%s  %3d::%3d\n",i,align[i],strc[i].res_type,"xxx",strc[i].res_number,0);
				continue;
			}
			//printf("%3d :: %3d   %s::%s  %3d::%3d\n",i,align[i],strc[i].res_type,strc_t[align[i]].res_type,strc[i].res_number,strc_t[align[i]].res_number);
			j = align[i];
			//printf("-%s- -%d- %s :: -%s- -%d- %s\n",strc[i].res_type,strc[i].res_number,strc[i].chain,strc_t[j].res_type,strc_t[j].res_number,strc_t[j].chain);
		}
		
		
		
		
		return(score);

}
 
void rotate_all(gsl_matrix *rota,struct pdb_atom *all_init,int all)
{
	int i,j,k,l;
	
	double temp_x, temp_y, temp_z;
	
	for(i=0;i<all;++i)
	{
		temp_x = all_init[i].x_cord;
		temp_y = all_init[i].y_cord;
		temp_z = all_init[i].z_cord;
		
		all_init[i].x_cord = gsl_matrix_get(rota,0,0) * temp_x + gsl_matrix_get(rota,1,0) * temp_y + gsl_matrix_get(rota,2,0) * temp_z;
		all_init[i].y_cord = gsl_matrix_get(rota,0,1) * temp_x + gsl_matrix_get(rota,1,1) * temp_y + gsl_matrix_get(rota,2,1) * temp_z;
		all_init[i].z_cord = gsl_matrix_get(rota,0,2) * temp_x + gsl_matrix_get(rota,1,2) * temp_y + gsl_matrix_get(rota,2,2) * temp_z;
		
		gsl_matrix *neweigen = gsl_matrix_alloc(3,3);
		
		for(j = 0; j < 3; j++)
		{
			for(k = 0; k < 3; k++)
			{
				double eigenjk = 0;
				
				for(l = 0; l < 3; l++)
				{
					eigenjk += gsl_matrix_get(rota,l,j)*all_init[i].global_evecs[l][k];
				}
				
				gsl_matrix_set(neweigen, j, k, eigenjk);
			}
		}
		
		for(j = 0; j < 3; j++)
		{
			for(k = 0; k < 3; k++)
			{
				all_init[i].global_evecs[j][k] = gsl_matrix_get(neweigen, j, k);
			}
		}
		
		gsl_matrix_free(neweigen);
	}

}

 
 void write_movie(FILE *out_file, struct pdb_atom *newstrc,int nb_atom,int model){
 	int k;
 	fprintf(out_file,"MODEL	%d\n",model);
 	
 	for (k = 0; k<nb_atom;k++) {
 		if (newstrc[k].atom_type == 1) {fprintf(out_file,"ATOM  ");}
	 	if (newstrc[k].atom_type == 2) {fprintf(out_file,"HETATM");}
	 	if (newstrc[k].atom_type == 3) {fprintf(out_file,"HETATM");}

	 		fprintf(out_file,"%5.d %s%s %s%4.d%12.3f%8.3f%8.3f  1.00  %2.2f\n",
	 			newstrc[k].atom_number,
	 			newstrc[k].atom_prot_type,
	 			newstrc[k].res_type,
	 			newstrc[k].chain,
	 			newstrc[k].res_number,
	 			newstrc[k].x_cord,
	 			newstrc[k].y_cord,
	 			newstrc[k].z_cord,
	 			newstrc[k].b_factor
	 			
	 			);
 	}
 	fprintf(out_file,"ENDMDL\n");
 }

double gsl_matrix_Det3D(gsl_matrix *M){
  int i,j;
  double det;

  //  guide to indexes: 0=x, 1=y, 2=z
  //  00 01 02
  //  10 11 12
  //  20 21 22

  det = gsl_matrix_get(M,0,0)*gsl_matrix_get(M,1,1)*gsl_matrix_get(M,2,2)
    + gsl_matrix_get(M,0,1)*gsl_matrix_get(M,1,2)*gsl_matrix_get(M,2,0)
    + gsl_matrix_get(M,0,2)*gsl_matrix_get(M,2,1)*gsl_matrix_get(M,1,0)
    - gsl_matrix_get(M,0,2)*gsl_matrix_get(M,1,1)*gsl_matrix_get(M,2,0)
    - gsl_matrix_get(M,1,2)*gsl_matrix_get(M,2,1)*gsl_matrix_get(M,0,0)
    - gsl_matrix_get(M,2,2)*gsl_matrix_get(M,1,0)*gsl_matrix_get(M,0,1);

  return det;
}

 void multiplie_matrix(gsl_matrix *a,int a_one,int a_two,gsl_matrix *b,int b_one,int b_two,gsl_matrix *c) {
 	int i,j,k;
 	double val;
 	// Ligne fois collone
 	if (a_two != b_one) {printf("CANNOT DO\n");}
 	for(i=0;i<a_one;++i) {
 		for(j=0;j<b_two;++j) {
 		val = 0.0;
 		for(k=0;k<a_two;++k) {
 			//printf("K:%d\n",k);
 			val += gsl_matrix_get(a,i,k)*gsl_matrix_get(b,k,j);
 		}
 		//printf("I:%d	J:%d	Val:%f\n",i,j,val);
 		gsl_matrix_set(c,i,j,val);
 		}
 	}
 	
 }

 float rmsd_no(struct pdb_atom *init,struct pdb_atom *targ,int atom, int *align) {
 	// Fonction qui suimperpose deux structures en utilisant single value decomposition
 	
 	int i,j,k; // Dummy
 	int t_atom;	
 	double cen_init[3],cen_targ[3];
 	double rmsd;
	int d;
	int temp;
	//printf("IN FUNCTION\n");
 	gsl_matrix *init_v = gsl_matrix_alloc(atom,3);
 	gsl_matrix *targ_v = gsl_matrix_alloc(3,atom);
 	gsl_matrix *corr   = gsl_matrix_alloc(3,3);
 	gsl_vector *vec_s  = gsl_vector_alloc(3);
 	gsl_vector *vec_w  = gsl_vector_alloc(3);
 	gsl_matrix *mat_v  = gsl_matrix_alloc(3,3);
 	gsl_matrix *rota   = gsl_matrix_alloc(3,3);
	
	// Construit structure qui comprend seulement les atomes correpondants
	
	struct pdb_atom t_init[atom];
	struct pdb_atom t_targ[atom];
	k=-1;
	for(i=0;i<atom;++i) {
 		d = align[i];
 		//printf("D:%d	%d\n",d,align[i]);
 		if(d == -1) {continue;}
 		++k;
 		t_init[k].x_cord = init[i].x_cord;
 		t_init[k].y_cord = init[i].y_cord;
 		t_init[k].z_cord = init[i].z_cord;
 		
 		t_targ[k].x_cord = targ[d].x_cord;
 		t_targ[k].y_cord = targ[d].y_cord;
 		t_targ[k].z_cord = targ[d].z_cord;
 	}
 	t_atom =k+1;
 	// Centre les deux structures
 	if (k == -1) {return(-1);}
 	for(i=0;i<3;++i) {
 		cen_init[i] = 0;
 		cen_targ[i] = 0;
 	}
 	rmsd = 0;
 	for(i=0;i<t_atom;++i) {
 		
 		cen_init[0] += t_init[i].x_cord;
 		cen_init[1] += t_init[i].y_cord;
 		cen_init[2] += t_init[i].z_cord;
 		
 		cen_targ[0] += t_targ[i].x_cord;
 		cen_targ[1] += t_targ[i].y_cord;
 		cen_targ[2] += t_targ[i].z_cord;
 		rmsd += ((t_init[i].x_cord-t_targ[i].x_cord)*(t_init[i].x_cord-t_targ[i].x_cord)+
 				 (t_init[i].y_cord-t_targ[i].y_cord)*(t_init[i].y_cord-t_targ[i].y_cord)+
 				 (t_init[i].z_cord-t_targ[i].z_cord)*(t_init[i].z_cord-t_targ[i].z_cord));
 	}
 	
 	//printf("	Init RMSD:%f	",rmsd/t_atom);
 	for(i=0;i<3;++i) {
 		cen_init[i] /= t_atom;
 		cen_targ[i] /= t_atom;
 	}

	// Translationne et crée des matrice à multiplier

 	rmsd = 0;
 	for(i=0;i<t_atom;++i) {
 	
 		t_init[i].x_cord -= cen_init[0];
 		t_init[i].y_cord -= cen_init[1];
 		t_init[i].z_cord -= cen_init[2];
 		
 		gsl_matrix_set(init_v,i,0,t_init[i].x_cord);
 		gsl_matrix_set(init_v,i,1,t_init[i].y_cord);
 		gsl_matrix_set(init_v,i,2,t_init[i].z_cord);
 		
 		t_targ[i].x_cord -= cen_targ[0];
 		t_targ[i].y_cord -= cen_targ[1];
 		t_targ[i].z_cord -= cen_targ[2];
 		

 		gsl_matrix_set(targ_v,0,i,t_targ[i].x_cord);
 		gsl_matrix_set(targ_v,1,i,t_targ[i].y_cord);
 		gsl_matrix_set(targ_v,2,i,t_targ[i].z_cord);
 
 		rmsd += ((t_init[i].x_cord-t_targ[i].x_cord)*(t_init[i].x_cord-t_targ[i].x_cord)+
 				 (t_init[i].y_cord-t_targ[i].y_cord)*(t_init[i].y_cord-t_targ[i].y_cord)+
 				 (t_init[i].z_cord-t_targ[i].z_cord)*(t_init[i].z_cord-t_targ[i].z_cord));
 	}
 	//printf("	Trans RMSD:%f\n",rmsd/t_atom);
 	
 	// Multiply les deux vecteurs qui comprennent les valeurs
 	 	
 	multiplie_matrix(targ_v,3,t_atom,init_v,t_atom,3,corr);

 	gsl_linalg_SV_decomp (corr, mat_v,vec_s,vec_w);
	//gsl_linalg_SV_decomp_jacobi (corr, mat_v, vec_s);
 	
 	// Multipli les matrix
 	
	gsl_matrix_transpose (corr);
	
	multiplie_matrix(mat_v,3,3,corr,3,3,rota);

 	rotate_all(rota,t_init,t_atom);
 	
 	rmsd = 0;
 	float max_displacement = 0;
 	for(i=0;i<t_atom;++i) {
 		/*printf("I:%d %8.4f\n",i,((t_init[i].x_cord-t_targ[i].x_cord)*(t_init[i].x_cord-t_targ[i].x_cord)+
 				 (t_init[i].y_cord-t_targ[i].y_cord)*(t_init[i].y_cord-t_targ[i].y_cord)+
 				 (t_init[i].z_cord-t_targ[i].z_cord)*(t_init[i].z_cord-t_targ[i].z_cord)));*/
 		float dist = ((t_init[i].x_cord-t_targ[i].x_cord)*(t_init[i].x_cord-t_targ[i].x_cord)+
 				 (t_init[i].y_cord-t_targ[i].y_cord)*(t_init[i].y_cord-t_targ[i].y_cord)+
 				 (t_init[i].z_cord-t_targ[i].z_cord)*(t_init[i].z_cord-t_targ[i].z_cord));
 		rmsd += dist;
 				 if (dist > max_displacement) {max_displacement = dist;}
 	}
 	//printf("Max_displacement:%.4f\n",sqrt(max_displacement));
 	//printf("	Rota RMSD:%f = %f/%d\n",rmsd/t_atom,rmsd,t_atom);
 	
 	gsl_matrix_free(corr);
 	gsl_matrix_free(init_v);
 	gsl_matrix_free(targ_v);
 	gsl_matrix_free(mat_v);
    gsl_vector_free(vec_s);
  	gsl_vector_free(vec_w);
 	gsl_matrix_free(rota);

 	
 	return(rmsd/t_atom);
 }
 
 float rmsd_yes(struct pdb_atom *init,struct pdb_atom *targ,int atom, int *align,struct pdb_atom *all_init,int all) {
 	// Fonction qui suimperpose deux structures en utilisant single value decomposition
 	// Centre les deux structures autour des atomes dans align et rotate init et init all pour fitter !
 	
 	int i,j,k; // Dummy
 	int t_atom;	
 	double cen_init[3],cen_targ[3];
 	double rmsd;
	int d;
	int temp;
	
 	gsl_matrix *init_v = gsl_matrix_alloc(atom,3);
 	gsl_matrix *targ_v = gsl_matrix_alloc(3,atom);
 	gsl_matrix *corr   = gsl_matrix_alloc(3,3);
 	gsl_vector *vec_s  = gsl_vector_alloc(3);
 	gsl_vector *vec_w  = gsl_vector_alloc(3);
 	gsl_matrix *mat_v  = gsl_matrix_alloc(3,3);
 	gsl_matrix *rota   = gsl_matrix_alloc(3,3);
	
	// Construit structure qui comprend seulement les atomes aligné
	
	struct pdb_atom t_init[atom];
	struct pdb_atom t_targ[atom];
	k=-1;
	for(i=0;i<atom;++i) {
 		d = align[i];
 		//printf("D:%d	%d\n",d,align[i]);
 		if(d == -1) {continue;}

 		//printf("D:%d	%d\n",i,align[i]);
 		//printf("(%8.5f,%8.5f,%8.5f) (%8.5f,%8.5f,%8.5f)\n",init[i].x_cord,init[i].y_cord,init[i].z_cord,targ[d].x_cord,targ[d].y_cord,targ[d].z_cord);
 		++k;
 		t_init[k].x_cord = init[i].x_cord;
 		t_init[k].y_cord = init[i].y_cord;
 		t_init[k].z_cord = init[i].z_cord;
 		
 		t_targ[k].x_cord = targ[d].x_cord;
 		t_targ[k].y_cord = targ[d].y_cord;
 		t_targ[k].z_cord = targ[d].z_cord;
 	}
 	t_atom =k+1;
 	// Centre les deux structures
 	if (k == -1) {return(-1);}
 	for(i=0;i<3;++i) {
 		cen_init[i] = 0;
 		cen_targ[i] = 0;
 	}
 	rmsd = 0;
 	// Trouve le centre de masse pour chaque structure
 	for(i=0;i<t_atom;++i) {
 		
 		cen_init[0] += t_init[i].x_cord;
 		cen_init[1] += t_init[i].y_cord;
 		cen_init[2] += t_init[i].z_cord;
 		
 		cen_targ[0] += t_targ[i].x_cord;
 		cen_targ[1] += t_targ[i].y_cord;
 		cen_targ[2] += t_targ[i].z_cord;
 		rmsd += ((t_init[i].x_cord-t_targ[i].x_cord)*(t_init[i].x_cord-t_targ[i].x_cord)+
 				 (t_init[i].y_cord-t_targ[i].y_cord)*(t_init[i].y_cord-t_targ[i].y_cord)+
 				 (t_init[i].z_cord-t_targ[i].z_cord)*(t_init[i].z_cord-t_targ[i].z_cord));
 	}
 	
 	
 	for(i=0;i<3;++i) {
 		cen_init[i] /= t_atom;
 		cen_targ[i] /= t_atom;
 	}

	// Translationne et crée des matrice à multiplier (on retournera au centre tantot)
		
 	rmsd = 0;
 	for(i=0;i<t_atom;++i) {
 	 		

 		
 		t_init[i].x_cord -= cen_init[0];
 		t_init[i].y_cord -= cen_init[1];
 		t_init[i].z_cord -= cen_init[2];
 		
 		
 		t_targ[i].x_cord -= cen_targ[0];
 		t_targ[i].y_cord -= cen_targ[1];
 		t_targ[i].z_cord -= cen_targ[2];

 		gsl_matrix_set(init_v,i,0,t_init[i].x_cord);
 		gsl_matrix_set(init_v,i,1,t_init[i].y_cord);
 		gsl_matrix_set(init_v,i,2,t_init[i].z_cord);

 		gsl_matrix_set(targ_v,0,i,t_targ[i].x_cord);
 		gsl_matrix_set(targ_v,1,i,t_targ[i].y_cord);
 		gsl_matrix_set(targ_v,2,i,t_targ[i].z_cord);
 
 		rmsd += ((t_init[i].x_cord-t_targ[i].x_cord)*(t_init[i].x_cord-t_targ[i].x_cord)+
 				 (t_init[i].y_cord-t_targ[i].y_cord)*(t_init[i].y_cord-t_targ[i].y_cord)+
 				 (t_init[i].z_cord-t_targ[i].z_cord)*(t_init[i].z_cord-t_targ[i].z_cord));
 	}
 	//printf("	Trans RMSD:%f\n",rmsd/t_atom);
 	
 	// Multiply les deux vecteurs qui comprennent les valeurs
 	 	
 	multiplie_matrix(targ_v,3,t_atom,init_v,t_atom,3,corr);


 	gsl_linalg_SV_decomp (corr, mat_v,vec_s,vec_w);
	//gsl_linalg_SV_decomp_jacobi (corr, mat_v, vec_s);
 	
 	// Multiplie les matrix
 	
	gsl_matrix_transpose (corr);
	
	multiplie_matrix(mat_v,3,3,corr,3,3,rota);
	//print_matrix(rota);
	
	// Centre avant de rotationer, les pdb qui n'ont pas été copier dans t_init
	
	for(i=0;i<atom;++i) {
		init[i].x_cord -= cen_init[0];
 		init[i].y_cord -= cen_init[1];
 		init[i].z_cord -= cen_init[2];
	}
	
	for(i=0;i<all;++i) {
		all_init[i].x_cord -= cen_init[0];
 		all_init[i].y_cord -= cen_init[1];
 		all_init[i].z_cord -= cen_init[2];
	}
	
 	rotate_all(rota,init,atom);
 	rotate_all(rota,t_init,t_atom);
 	rotate_all(rota,all_init,all);
 	
 	
 	rmsd = 0;
 	float max_displacement = 0;
 	for(i=0;i<t_atom;++i) {
 		/*printf("I:%d %8.4f\n",i,((t_init[i].x_cord-t_targ[i].x_cord)*(t_init[i].x_cord-t_targ[i].x_cord)+
 				 (t_init[i].y_cord-t_targ[i].y_cord)*(t_init[i].y_cord-t_targ[i].y_cord)+
 				 (t_init[i].z_cord-t_targ[i].z_cord)*(t_init[i].z_cord-t_targ[i].z_cord)));*/
 				float	dist = ((t_init[i].x_cord-t_targ[i].x_cord)*(t_init[i].x_cord-t_targ[i].x_cord)+
 				 (t_init[i].y_cord-t_targ[i].y_cord)*(t_init[i].y_cord-t_targ[i].y_cord)+
 				 (t_init[i].z_cord-t_targ[i].z_cord)*(t_init[i].z_cord-t_targ[i].z_cord));
 					rmsd += dist;
 				 if (dist > max_displacement) {max_displacement = dist;}
 	}
// 	printf("Max_displacement:%.4f\n",sqrt(max_displacement));
 	//printf("	Rota RMSD:%f\n",rmsd/t_atom);
 	
 	
 	for(i=0;i<all;++i) {
 		/*all_init[i].x_cord -= cen_init[0];
 		all_init[i].y_cord -= cen_init[1];
 		all_init[i].z_cord -= cen_init[2];*/
 		
 		all_init[i].x_cord += cen_targ[0];
 		all_init[i].y_cord += cen_targ[1];
 		all_init[i].z_cord += cen_targ[2];
 	}
 	
	for(i=0;i<atom;++i) {
		/*init[i].x_cord -= cen_init[0];
 		init[i].y_cord -= cen_init[1];
 		init[i].z_cord -= cen_init[2];*/
	
 		init[i].x_cord += cen_targ[0];
 		init[i].y_cord += cen_targ[1];
 		init[i].z_cord += cen_targ[2];
 	}
 	
 	
 	gsl_matrix_free(corr);
 	gsl_matrix_free(init_v);
 	gsl_matrix_free(targ_v);
 	gsl_matrix_free(mat_v);
    gsl_vector_free(vec_s);
  	gsl_vector_free(vec_w);
 	gsl_matrix_free(rota);

 	
 	return(rmsd/t_atom);
 }
 
 void apply_eigen(struct pdb_atom *strc,int atom,gsl_matrix *m,int mode,float amplitude) {
 	int k;
 //	amplitude = 0;
 	//printf("In Function\t");
 	//printf("Ampli Eigen:%8.5f\n",amplitude);
 	for(k=0;k<atom;++k) {
 		//printf("K:%d %f %f %f\n",k,strc[k].x_cord,strc[k].y_cord,strc[k].z_cord);
 		//printf("K:%d Mode:%d	Node:%d	Place:%d	Amplit:%f Atom type:%s Res:%s\n",k,mode,strc[k].node,strc[k].node*3,amplitude,strc[k].atom_prot_type,strc[k].res_type);
 		if (strc[k].node < 0 || strc[k].node*3+2 > m->size1) {
 			//printf("K:%d Mode:%d	Node:%d	Place:%d	Amplit:%f Atom type:%s Res:%s\n",k,mode,strc[k].node,strc[k].node*3,amplitude,strc[k].atom_prot_type,strc[k].res_type);
 			continue;
 		}
 		if (strc[k].atom_type == 3) {continue;}
	 	strc[k].x_cord += (amplitude*gsl_matrix_get(m,strc[k].node*3,mode));
		strc[k].y_cord += (amplitude*gsl_matrix_get(m,strc[k].node*3+1,mode));
		strc[k].z_cord += (amplitude*gsl_matrix_get(m,strc[k].node*3+2,mode));
		//printf("K:%d %f %f %f\n",k,strc[k].x_cord,strc[k].y_cord,strc[k].z_cord);
	}
 }
 
 float vector_lenght(gsl_vector *d,int atom) {
 	int i;
 	float a = 0;
 	for(i=0;i<atom;++i) {
 		a += gsl_vector_get(d,i)*gsl_vector_get(d,i);
 	}
 	return(sqrt(a));
 }
 
 float overlap(int atom, int mode,gsl_matrix *m,gsl_vector *d, int *align) {
 	int i;
 	float a = 0,b = 0,c = 0;
 	float x,z;
 	float lenght = vector_lenght(d,atom);
 	int j = 0;
 	int k;
 	for(i=0;i<atom;++i) {
 		if (align[i] == -1) {continue;}
 	//	printf("I:%d (%8.5f,%8.5f,%8.5f) (%8.5f,%8.5f,%8.5f)\n",i,gsl_vector_get(d,j*3+0)/lenght,gsl_vector_get(d,j*3+1)/lenght,gsl_vector_get(d,j*3+2)/lenght,gsl_matrix_get(m,i*3+0,mode),gsl_matrix_get(m,i*3+1,mode),gsl_matrix_get(m,i*3+2,mode));
 		for (k=0;k<3;++k) {
	 		x = gsl_vector_get(d,j*3+k)/lenght;
	 		z = gsl_matrix_get(m,i*3+k,mode);
	 		a += x*z;
	 		b += x*x;
	 		c += z*z;
	 	}
 		++j;
 	}
 	a = sqrt(a*a);
 	//printf("A:%f	B:%f	C:%f\n",a,b,c);
 	return(a/sqrt(b*c));
 }

 void fit(struct pdb_atom *init,struct pdb_atom *targ,int atom,int all,struct pdb_atom *all_init,struct pdb_atom *all_targ,gsl_matrix *eval, int *align) {
 	


	int i,j,l,k;
 	long seed; 
 	int mc; // Nombre de Step de Monte Carlo
 	float newstrc,old = 9999999;
 	int mode;
 	double random;
 	float step;
 	//FILE *file;
 	//file = fopen("motion_align.pdb","w");
 	struct pdb_atom store_init[atom];
 	struct pdb_atom work_init[atom];
	//seed = time_seed();
 	gsl_matrix *rota     = gsl_matrix_alloc(3,3);
 	gsl_vector *diff     = gsl_vector_alloc(3*atom);
 	
 	// Store Init dans store init
 	
 	for(i=0;i<atom;++i) {
 		store_init[i].x_cord = init[i].x_cord;
		store_init[i].y_cord = init[i].y_cord;
		store_init[i].z_cord = init[i].z_cord;
		store_init[i].res_number = init[i].res_number;
		store_init[i].atom_number = init[i].atom_number;
		store_init[i].node = init[i].node;
		store_init[i].atom_type = init[i].atom_type;
		strcpy(store_init[i].atom_prot_type,init[i].atom_prot_type);
		strcpy(store_init[i].res_type,init[i].res_type);
		strcpy(store_init[i].chain,init[i].chain);
		
 	}
 	
 	// Minimization
 	old = rmsd_no(store_init,targ,atom,align);

 	// Construit vector des différences
 	
 	for(i=0;i<atom;++i) {
 		gsl_vector_set(diff,i*3+0,store_init[i].x_cord - targ[i].x_cord);
		gsl_vector_set(diff,i*3+1,store_init[i].y_cord - targ[i].y_cord);
		gsl_vector_set(diff,i*3+2,store_init[i].z_cord - targ[i].z_cord); 	
 	}
 	printf("INIT:%f\n",old);
	mode = 6;
	int nb_mode = atom*3-6;
	double amp[nb_mode];
	double namp[nb_mode];
	int rand_node;
	
	
	
	for(i=0;i<nb_mode;++i) {
		amp[i]  = 0;
		namp[i] = 0;
	}
	

	for(j = 0;j<nb_mode;++j) {
		step = 1;
	 	for (i = 0;i < 200;++i) {
	 		for(l=0;l<atom;++l) {
	 			store_init[l].x_cord = init[l].x_cord;
				store_init[l].y_cord = init[l].y_cord;
				store_init[l].z_cord = init[l].z_cord;
	 		}
	 		
	 		for(l=0;l<nb_mode;++l) {
	 			apply_eigen(store_init,atom,eval,mode+l,namp[l]);
	 		}
		 	
	 		newstrc = rmsd_no(store_init,targ,atom,align);
	 		if (i == 0) {old = newstrc+1;}
	 		if (step < 0.0001 && step > -0.0001) {amp[j] = namp[j];printf("%3d Mode:%4d Amp:%7.4f	Step:%7.4f	RMS:%10.9f\n",i,mode+j+1,namp[j],step,newstrc);break;}
	 		//
	 		if (i > nb_mode * atom *100) {break;}
			if (old < newstrc) {step = -step/10;}
			//if (old == newstrc && i > 4) {amp[j] = namp[j];break;}
	 		old = newstrc;
	 		namp[j] += step;
		}
	}
	// Monte Carlo Approch
	for(i=0;i<nb_mode*1000;++i) {
		amp[i] = namp[i];
	}
	printf("RMS before MC:%f\n",old);
 	for(i=0;i<50*1;++i) {
 		
 		for(l=0;l<atom;++l) {
 			store_init[l].x_cord = init[l].x_cord;
			store_init[l].y_cord = init[l].y_cord;
			store_init[l].z_cord = init[l].z_cord;
 		}
 		
 		for(j=0;j<nb_mode;++j) {
 			if (namp[j] == 0) {continue;}
	 		apply_eigen(store_init,atom,eval,mode+j,namp[j]);
	 	}
	 	
 		newstrc = rmsd_no(store_init,targ,atom,align);
 		
 		if (newstrc < old) {
 			//printf("I:%4d Old:%10.9f	newstrc:%10.9f\n",i,old,newstrc);
 			
 			for(j=0;j<nb_mode;++j) {
 				//printf("%5.3f	",amp[j]);
 				amp[j] = namp[j];
 			}
 			
 			old = newstrc;
 		} else {
 			for(j=0;j<nb_mode;++j) {
 				namp[j] = amp[j];
 			}
 		}
 		
 		// Modify namp
 		
 		rand_node = ran2(&seed)*nb_mode;
 		
 		namp[rand_node] += ran2(&seed)*0.5 - 0.25;
		 		
 	}
 	printf("RMS after  MC:%f\n",newstrc);
 	newstrc = 0;
 	for(j=0;j<nb_mode;++j) {
 		old = overlap(atom,mode+j,eval,diff,align);	
 		apply_eigen(init,atom,eval,mode+j,amp[j]);
		apply_eigen(all_init,all,eval,mode+j,amp[j]);
		newstrc = rmsd_no(init,targ,atom,align);
		printf("J:%4d %10.6f	%10.6f	%10.6f\n",j+mode+1,amp[j],old,newstrc);
	}
 	printf("END:%f\n",rmsd_no(init,targ,atom,align));
 	//write_strc("target.pdb",all_targ,all);
 	//write_strc("final.pdb",all_init,all);
 	//fclose(file);
 	gsl_vector_free(diff);
 }
 
  void fit_math(struct pdb_atom *init,struct pdb_atom *targ,int atom,int all,struct pdb_atom *all_init,struct pdb_atom *all_targ,gsl_matrix *eval, int *align) {
 	


	int i,j,l,k;
 	long seed; //seed = time_seed();
 	int mc; // Nombre de Step de Monte Carlo
 	float newstrc,old = 9999999;
 	int mode;
 	double random;
 	float step;
 		mode = 6;
	int nb_mode = 2;
 	//FILE *file;
 	//file = fopen("motion_align.pdb","w");
 	struct pdb_atom store_init[atom];
 	struct pdb_atom work_init[atom];

 	gsl_vector *diff     = gsl_vector_alloc(3*atom);
 	gsl_matrix *T   = gsl_matrix_alloc(3*atom,1);
 	gsl_matrix *M   = gsl_matrix_alloc(3*atom,nb_mode);
 	gsl_matrix *Mt  = gsl_matrix_alloc(nb_mode,3*atom);
 	gsl_matrix *LU  = gsl_matrix_alloc(nb_mode,nb_mode);
 	gsl_matrix *MtT = gsl_matrix_alloc(nb_mode,1);
 	gsl_matrix *inv = gsl_matrix_alloc(nb_mode,nb_mode);
 	gsl_matrix *B   = gsl_matrix_alloc(nb_mode,1);
 	gsl_permutation *p = gsl_permutation_alloc (nb_mode);
 	int s;
 	
 	// Store Init dans store init
 	
 	for(i=0;i<atom;++i) {
 		store_init[i].x_cord = init[i].x_cord;
		store_init[i].y_cord = init[i].y_cord;
		store_init[i].z_cord = init[i].z_cord;
		store_init[i].res_number = init[i].res_number;
		store_init[i].atom_number = init[i].atom_number;
		store_init[i].node = init[i].node;
		store_init[i].atom_type = init[i].atom_type;
		strcpy(store_init[i].atom_prot_type,init[i].atom_prot_type);
		strcpy(store_init[i].res_type,init[i].res_type);
		strcpy(store_init[i].chain,init[i].chain);
		
 	}
 	
 	// Minimization
 	old = rmsd_no(store_init,targ,atom,align);

 	// Construit vector des différences
 	
 	for(i=0;i<atom;++i) {
 		gsl_vector_set(diff,i*3+0,store_init[i].x_cord - targ[i].x_cord);
		gsl_vector_set(diff,i*3+1,store_init[i].y_cord - targ[i].y_cord);
		gsl_vector_set(diff,i*3+2,store_init[i].z_cord - targ[i].z_cord); 	
 	}
 	printf("INIT:%f\n",old);

	double amp[nb_mode];
	
	// Aproche Mathématique
	
	// Construct a 3N vector thaht show difference between structure

	
	for(i=0;i<atom;++i) {
		gsl_matrix_set(T,i*3+0,0,init[i].x_cord-targ[i].x_cord);
		gsl_matrix_set(T,i*3+1,0,init[i].y_cord-targ[i].y_cord);
		gsl_matrix_set(T,i*3+2,0,init[i].z_cord-targ[i].z_cord);
	}
	
	// Construct a 3N x Mode matrix of the eigenvector

	
	for(i=0;i<atom;++i) {
		for(l=0;l<nb_mode;++l) {
			gsl_matrix_set(M,i,l,gsl_matrix_get(eval,i,mode+l));
			gsl_matrix_set(Mt,l,i,gsl_matrix_get(eval,i,mode+l));
		}
	}
	
	// Solve B = (Mt*M)^(-1)*Mt*T
	
	multiplie_matrix(Mt,nb_mode,3*atom,M,3*atom,nb_mode,LU);
	multiplie_matrix(Mt,nb_mode,3*atom,T,3*atom,1,MtT);
	
	// Inverse Mt*M = LU
	
	gsl_linalg_LU_decomp (LU, p,&s);
	gsl_linalg_LU_invert (LU, p, inv);
	
	// Mutlupli Inv*MtT
	
	multiplie_matrix(inv,nb_mode,nb_mode,MtT,nb_mode,1,B);
	
	for(i=0;i<nb_mode;++i) {
		amp[j] = gsl_matrix_get(B,i,0);
	}
 	
 	newstrc = 0;
 	for(j=0;j<nb_mode;++j) {
 		old = overlap(atom,mode+j,eval,diff,align);	
 		apply_eigen(init,atom,eval,mode+j,amp[j]);
		apply_eigen(all_init,all,eval,mode+j,amp[j]);
		newstrc = rmsd_no(init,targ,atom,align);
		printf("J:%4d %10.6f	%10.6f	%10.6f\n",j+mode+1,amp[j],old,newstrc);
	}
 	printf("END:%f\n",rmsd_no(init,targ,atom,align));
 	//write_strc("target.pdb",all_targ,all);
 	//write_strc("final.pdb",all_init,all);
 	//fclose(file);
 	gsl_vector_free(diff);
 }
 
  void fit_vince(struct pdb_atom *init,struct pdb_atom *targ,int atom,int all,struct pdb_atom *all_init,struct pdb_atom *all_targ,gsl_matrix *evec, int *align, int nb_mode, int mode,gsl_vector *eval) {
	int i,j,l,k,m;
 	float newstrc,old;
	double amp[nb_mode];
	float quad[3];
	float b,a;
	
 	struct pdb_atom store_init[atom];

 	gsl_vector *diff     = gsl_vector_alloc(3*atom);
 	
 	// Store Init dans store init
 	
 	for(i=0;i<atom;++i) {
 		store_init[i].x_cord = init[i].x_cord;
		store_init[i].y_cord = init[i].y_cord;
		store_init[i].z_cord = init[i].z_cord;
		store_init[i].node = init[i].node;
		store_init[i].atom_type = init[i].atom_type;
		
		store_init[i].res_number = init[i].res_number;
		store_init[i].atom_number = init[i].atom_number;
		store_init[i].node = init[i].node;
		store_init[i].atom_type = init[i].atom_type;
		strcpy(store_init[i].atom_prot_type,init[i].atom_prot_type);
		strcpy(store_init[i].res_type,init[i].res_type);
		strcpy(store_init[i].chain,init[i].chain);
 	}
 	
 	// Minimization
 	old = rmsd_no(store_init,targ,atom,align);
	// Construit vector des différences et on va storer dans eigenvector 6
 	
 	for(i=0;i<atom;++i) {
 		gsl_vector_set(diff,i*3+0,store_init[i].x_cord - targ[i].x_cord);
		gsl_vector_set(diff,i*3+1,store_init[i].y_cord - targ[i].y_cord);
		gsl_vector_set(diff,i*3+2,store_init[i].z_cord - targ[i].z_cord);
		//gsl_matrix_set(evec,	
 	}
 	printf("INIT:%f\n",sqrt(old));
	
	for(j=0;j<nb_mode;++j) {
		amp[j] = 0;
	}
	//if(nb_mode > 3*atom) {nb_mode = 3*atom;}

		for(j = 0;j<nb_mode;++j) {
		 	for (i = -1;i < 2;++i) {
		 		for(l=0;l<atom;++l) {
		 			store_init[l].x_cord = init[l].x_cord;
					store_init[l].y_cord = init[l].y_cord;
					store_init[l].z_cord = init[l].z_cord;
		 		}
			
				for(k=0;k<nb_mode;++k) {
					if (k == j) {continue;}
					apply_eigen(store_init,atom,evec,mode+k-1,amp[k]);
				}
			
		 		apply_eigen(store_init,atom,evec,mode+j-1,i);
		 		//write_strc("no_strc.pdb",store_init,atom);	
		 		newstrc = rmsd_no(store_init,targ,atom,align);
		 		//printf("Mode:%4d I:%2d RMS:%10.6f\n",mode+j,i,newstrc);
		 		quad[i+1] = newstrc;

			}
			b = (quad[0] - quad[2])/2;		
			a = quad[0]/2 + quad[2]/2 - quad[1];

			amp[j] = b/(2*a);
			if (amp[j] > 50) {amp[j] = 0;}
			if (amp[j] < -50) {amp[j] = 0;}
			if ((2*a) == 0) {amp[j] = 0;}
			//printf("Amp:%f\t",amp[j]);
			for(l=0;l<atom;++l) {
		 		store_init[l].x_cord = init[l].x_cord;
				store_init[l].y_cord = init[l].y_cord;
				store_init[l].z_cord = init[l].z_cord;
		 	}
		 	for(k=0;k<nb_mode;++k) {
			//	if (k == j) {continue;}
				apply_eigen(store_init,atom,evec,mode+k-1,amp[k]);
			}
			newstrc = rmsd_no(store_init,targ,atom,align);
		 	//printf("Mode:%4d RMS:%10.6f\n",j+mode,newstrc);
			
		}
	
	
 	newstrc = 0;
 	float impro = sqrt(old);
 	old = 0;
 	printf("ASSIGNING MODE\n");
 	float energy = 0;
 	
 	for(j=0;j<nb_mode;++j) {
 		old = overlap(atom,mode+j-1,evec,diff,align);	
 		apply_eigen(init,atom,evec,mode+j-1,amp[j]);
		apply_eigen(all_init,all,evec,mode+j-1,amp[j]);
		newstrc = rmsd_no(init,targ,atom,align);
		energy += 0.5*pow(gsl_vector_get(eval,j+mode-1),2)*pow(amp[j],2);
		printf("J:%4d %10.6f %10.6f %10.6f %12.6f %10.6f %10.6f\n",j+mode,amp[j],old,sqrt(newstrc),energy,gsl_vector_get(eval,j+mode-1),impro-(newstrc));
		impro = sqrt(newstrc);
	}
 	printf("END:%f\n",sqrt(rmsd_no(init,targ,atom,align)));

 	gsl_vector_free(diff);
 }
 
 void nrg_rmsd(struct pdb_atom *init,int atom,gsl_matrix *evec, int *align, int nb_mode, int mode,gsl_vector *eval) {
	
	int i,j,l,k;
	float rmsd;
	struct pdb_atom store_init[atom];
	
	// Store Init dans store init
 	
 	for(i=0;i<atom;++i) {
 		store_init[i].x_cord = init[i].x_cord;
		store_init[i].y_cord = init[i].y_cord;
		store_init[i].z_cord = init[i].z_cord;
		store_init[i].node = init[i].node;
		store_init[i].atom_type = init[i].atom_type;
		
		store_init[i].res_number = init[i].res_number;
		store_init[i].atom_number = init[i].atom_number;
		store_init[i].node = init[i].node;
		store_init[i].atom_type = init[i].atom_type;
		strcpy(store_init[i].atom_prot_type,init[i].atom_prot_type);
		strcpy(store_init[i].res_type,init[i].res_type);
		strcpy(store_init[i].chain,init[i].chain);
 	}
	
	for(j = 0;j<nb_mode;++j) {
		for(l=0;l<atom;++l) {
			store_init[l].x_cord = init[l].x_cord;
			store_init[l].y_cord = init[l].y_cord;
			store_init[l].z_cord = init[l].z_cord;
		}
		apply_eigen(store_init,atom,evec,mode+j,1);
		rmsd = rmsd_no(store_init,init,atom,align);
		printf("%4d %10.5f\n",j+mode,gsl_vector_get(eval,j+mode)/rmsd);
	}
}
 
 void center_yes(struct pdb_atom *init,struct pdb_atom *targ,int atom,int atom_t, int *align) {
 	
 	int i,j,k; // Dummy
 	int t_atom;	
 	double cen_init[3],cen_targ[3];
	int d;
	
	// Construit structure qui comprend seulement les atomes correpondants
	
	struct pdb_atom t_init[atom];
	struct pdb_atom t_targ[atom];
	k=-1;
	for(i=0;i<atom;++i) {
 		d = align[i];
 		//printf("D:%d	%d\n",i,align[i]);
 		if(d == -1) {continue;}
 		//printf("D:%d	%d\n",i,align[i]);
 		//printf("(%8.5f,%8.5f,%8.5f) (%8.5f,%8.5f,%8.5f)\n",init[i].x_cord,init[i].y_cord,init[i].z_cord,targ[d].x_cord,targ[d].y_cord,targ[d].z_cord);
 		++k;
 		t_init[k].x_cord = init[i].x_cord;
 		t_init[k].y_cord = init[i].y_cord;
 		t_init[k].z_cord = init[i].z_cord;
 		
 		t_targ[k].x_cord = targ[d].x_cord;
 		t_targ[k].y_cord = targ[d].y_cord;
 		t_targ[k].z_cord = targ[d].z_cord;
 	}
 	t_atom =k+1;
 	// Centre les deux structures
 	for(i=0;i<3;++i) {
 		cen_init[i] = 0;
 		cen_targ[i] = 0;
 	}
 	
 	for(i=0;i<t_atom;++i) {
 		
 		cen_init[0] += t_init[i].x_cord;
 		cen_init[1] += t_init[i].y_cord;
 		cen_init[2] += t_init[i].z_cord;
 		
 		cen_targ[0] += t_targ[i].x_cord;
 		cen_targ[1] += t_targ[i].y_cord;
 		cen_targ[2] += t_targ[i].z_cord;
 	}	
 	
 	//printf("	Init RMSD:%f	",rmsd/t_atom);
 	for(i=0;i<3;++i) {
 		cen_init[i] /= t_atom;
 		cen_targ[i] /= t_atom;
 	}

	// Translationne et crée des matrice à multiplier
	
	
	for(i=0;i<atom;++i) {
 		init[i].x_cord -= cen_init[0];
 		init[i].y_cord -= cen_init[1];
 		init[i].z_cord -= cen_init[2];
 	}
 	for(i=0;i<atom_t;++i) {
 		targ[i].x_cord -= cen_targ[0];
 		targ[i].y_cord -= cen_targ[1];
 		targ[i].z_cord -= cen_targ[2];
 	}
 
 
 }

double get_diedre(int node,struct pdb_atom *init,int atom,int ty) {
	// Function qui retourne l'angle diedre d'une node (-1) pour avant et 1 pour l'autre
	if (ty != 1 && ty != -1) {return(0);}
	if (node -1 < 0) {return(0);}

	int i;
	gsl_vector *a = gsl_vector_alloc(3);
	gsl_vector *b = gsl_vector_alloc(3);
	gsl_vector *c = gsl_vector_alloc(3);
	
	if (ty == 1) {
		assign_vector(init,node-1," C   ",atom,init,node," N   ",atom,a);
		assign_vector(init,node," N   ",atom,init,node," CA  ",atom,b);
		assign_vector(init,node," CA  ",atom,init,node," C   ",atom,c);
	}
	if (ty == -1) {
		
		assign_vector(init,node," N   ",atom,init,node," CA  ",atom,a);
		assign_vector(init,node," CA  ",atom,init,node," C   ",atom,b);
		assign_vector(init,node," C   ",atom,init,node+1," N   ",atom,c);
	}
	
	gsl_vector *ab = gsl_vector_alloc(3);
	gsl_vector *bc = gsl_vector_alloc(3);
	
	dot_product_v(a,b,ab);
	dot_product_v(b,c,bc);
	
	return(vector_angle_v(ab,bc));
	
	
	
	}
	
void apply_eigen_rot(struct pdb_atom *old,gsl_matrix *evec,int mode,int atom,int nbnode,struct pdb_atom *strc_node,float amplitude) {
	
	
	int k;
	gsl_vector *vect = gsl_vector_alloc(3);
	gsl_matrix *rota = gsl_matrix_alloc(3,3);
	
	 for (k=0;k<nbnode;++k) {	
			
			center_node(atom,old,k);
			center_node(nbnode,strc_node,k);
	 		
	 		assign_vector(old,k," CA  ",atom,old,k," C   ",atom,vect);
	 		rotate_matrix(rota,gsl_matrix_get(evec,k*2+0,mode-1)*amplitude,vect);
		 	rotate_phi(rota,old,atom,k);
		 	rotate_phi(rota,strc_node,nbnode,k);
		 	
		 	assign_vector(old,k," CA  ",atom,old,k," N   ",atom,vect);
			rotate_matrix(rota,gsl_matrix_get(evec,k*2+1,mode-1)*amplitude,vect);
		 	rotate_psy(rota,old,atom,k);
		 	rotate_psy(rota,strc_node,nbnode,k);
		 	
		 	
		 //	rmsd_yes(newnode,strc_node,nbnode, align,newstrc,atom);
		 	
		// 	write_strc("nod_after_rmsd.pdb",newnode,nbnode);
		 	//break;
	 	}











}

  void fit_mc_torsion(struct pdb_atom *init,struct pdb_atom *targ,int atom,int all,int atom_t,int all_t,struct pdb_atom *all_init,struct pdb_atom *all_targ,gsl_matrix *evec, int *align, int nb_mode, int mode,gsl_vector *eval) {
	int j,i;
	float amp[nb_mode];
	float amp_bck[nb_mode];
	float old = 9999.999;
	float now = 9999.999;
	int rand_node = 0;
	long seed; 
	seed = time_seed();
	struct pdb_atom temp_all[all];
	struct pdb_atom temp_node[atom];
	for (j =0 ;j<nb_mode;++j) {
		amp[j] = 0;
		amp_bck[j] = amp[j];
	}
	
	for (i=0;i<10000;++i) {
		copy_strc(temp_all, all_init, all);
		copy_strc(temp_node, init, atom);
		
		// Apply the amplitude on temp_strc
	
	 	for(j=0;j<nb_mode;++j) {
			
	 		apply_eigen_rot(temp_all,evec,mode+j,all,atom,temp_node,amp[j]);
		
		
		}
		now = rmsd_no(temp_node,targ,atom,align);
		
		
		if (now < old) {
			// Accept
			old = now;
		
			printf("IT:%5d ",i);
			for (j=0;j<nb_mode;++j) {
				printf("%5.2f ",amp[j]);
				amp_bck[j] = amp[j];
			}
			printf("RMSD:%8.5f\n",sqrt(old));
			
		} else {
			// Refuse
			for (j=0;j<nb_mode;++j) {
				amp[j] = amp_bck[j];
			}
			
		}
		
		// Modify ampli
		
		rand_node = int(ran2(&seed) * nb_mode);
		//printf("rand_node:%d\n",rand_node);
		amp[rand_node] += ran2(&seed)-0.5 ;
		
	}


 }





float fit_svd(struct pdb_atom *init,struct pdb_atom *targ,int atom,int all,int atom_t,int all_t,struct pdb_atom *all_init,struct pdb_atom *all_targ,gsl_matrix *evec, int *align, int nb_mode, int mode,gsl_vector *eval) {
	int i,j,l,k,m;
 	float newstrc,old;
	double amp[nb_mode];
	float quad[3];
	float b,a;
	int alit = 0;// Value dans le vecteur diff
	int ngila[atom_t]; for (i=0;i<atom_t;++i) {ngila[i] = -1;  }
	for (i=0;i<atom;++i) {
		if (align[i] !=-1 ){ ngila[align[i]] = i; }
	}

	center_yes(init,targ,atom,atom_t, align); // Targ rotate, donc pas d'impact sur vecteurs
	rmsd_yes(targ,init,atom_t, ngila,all_targ,all_t);

	
 	gsl_vector *diff     = gsl_vector_alloc(3*atom);
 	
 	// Minimization
 	old = rmsd_no(init,targ,atom,align);
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


	for(j=0;j<nb_mode;++j) {
		amp[j] = 0;
	}
	
	// On veut ||Ax-b||² soit minimiser
	// http://www.ces.clemson.edu/~petersj/Agents/MatLabNA/MatLabNA004.html
	// A: Matrice des eigenvectors
	// B: Le vecteur de différence
	// x: Vecteur des amplitudes a appliquer sur chacun des Eigenvectors.... CE QU ON VEUT TROUVER
	// On décompose A = USV
	// Nouveau vecteur C = U^T*b
	// Y = S(diag)/C... Y = S(i,i)/C(i)
	// x = Y*V
	// e = Ax-b.... RMSD = ||e||
	
	// Go !!!
	// Construit matrix A
	gsl_matrix *A = gsl_matrix_alloc(3*alit,nb_mode);
	//printf("gsl_matrix_alloc(%d,%d)\n",3*alit,nb_mode);

	for (i=0;i<nb_mode;++i) {
		l = -1;
		for (j=0;j<atom;++j) {		
			if (align[j] == -1) {continue;}
			++l;
			//printf("Vector pour %d node: (%f,%f,%f) et eigen (%f,%f,%f)\n",j,gsl_vector_get(B,l*3),gsl_vector_get(B,l*3+1),gsl_vector_get(B,l*3+2),gsl_matrix_get(evec,j*3+0,i+mode-1),gsl_matrix_get(evec,j*3+1,i+mode-1),gsl_matrix_get(evec,j*3+2,i+mode-1));
			//printf("I:%d L:%d J:%d :: %d Max:%d\n",i,l,j,align[j],l*3+2);
			//printf("Mode:%d blabla:%d %d\n",mode,j*3+2,i+mode-1);
			gsl_matrix_set(A,l*3+0,i,gsl_matrix_get(evec,j*3+0,i+mode-1));
			gsl_matrix_set(A,l*3+1,i,gsl_matrix_get(evec,j*3+1,i+mode-1));
			gsl_matrix_set(A,l*3+2,i,gsl_matrix_get(evec,j*3+2,i+mode-1));
			
		/*	gsl_matrix_set(A,l*3+0,i,gsl_matrix_get(evec,i+mode-1,j*3+0));
			gsl_matrix_set(A,l*3+1,i,gsl_matrix_get(evec,i+mode-1,j*3+1));
			gsl_matrix_set(A,l*3+2,i,gsl_matrix_get(evec,i+mode-1,j*3+2));*/
			
		}
	
	}
	
	// Décompose A
	
	
	gsl_matrix *mat_U = gsl_matrix_alloc(nb_mode,nb_mode);
	gsl_vector *vec_S = gsl_vector_alloc(nb_mode);
	gsl_vector *vec_w = gsl_vector_alloc(nb_mode);
	gsl_vector *result = gsl_vector_alloc(nb_mode);
	printf("Decomposition\n");
	gsl_linalg_SV_decomp(A, mat_U,vec_S,vec_w);
	printf("SVD\n");
	gsl_linalg_SV_solve (A, mat_U, vec_S,B, result);
	
	for(i=0;i<nb_mode;++i) {
		amp[i] = -gsl_vector_get(result,i);
	}
	
 	newstrc = 0;
 	old = 0;
 	printf("ASSIGNING MODE\n");
 	float energy = 0;
 	float best_overlap = 0;
 	int best_mode = 0;
 	for(j=0;j<nb_mode;++j) {
 		old = overlap(atom,mode+j-1,evec,diff,align);
 		if (old > best_overlap) {best_overlap = old;best_mode = j + mode;}
 		apply_eigen(init,atom,evec,mode+j-1,amp[j]);
		apply_eigen(all_init,all,evec,mode+j-1,amp[j]);
		newstrc = rmsd_no(init,targ,atom,align);
		energy += 0.5*pow(gsl_vector_get(eval,j+mode-1),2)*pow(amp[j],2);
		printf("J:%4d %10.6f %10.6f %10.6f %12.6f\n",j+mode,amp[j],old,sqrt(newstrc),energy);
	}
 	printf("END:%f\n",sqrt(rmsd_no(init,targ,atom,align)));
	printf("Best Overlap:%f Mode:%d\n",best_overlap,best_mode);
 	gsl_vector_free(diff);
 	return(energy);
 }

void gen_gauss(struct pdb_atom *init, gsl_matrix *evec, gsl_vector *eval, int atom, double beta)
{
	const double PI = 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421;
	
	int i, j, node;
	
	int nm = atom*3 - 6;
	
	for (node = 0;node<atom;++node)
	{
		// printf("Node %3i\n", node);
		gsl_matrix *GJ = gsl_matrix_alloc(3,nm+3);
		gsl_matrix_set_all(GJ,0);
		
		// Set Rx,Ry,Rz
		for(i=0;i<3;++i)
		{
			gsl_matrix_set(GJ,i,nm+i,1);
		}
		
		// Set evec
		
		// printf("Set GJ\n");
		
		for (j = 0;j<nm;++j)
		{
			for(i=0;i<3;++i)
			{
				gsl_matrix_set(GJ,i,j,gsl_matrix_get(evec,i+node*3,j+6));
			}
		}
		
		// printf("GJ is set\n");
		
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
		
		// printf("Calc det123\n");
		
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
		
		// printf("Set K\n");
		
		// Initialize kappa matrix (equal to kappa, except for the first n-3 term of the diagonal)
		
		gsl_matrix *K = gsl_matrix_alloc(nm,nm);
		gsl_matrix_set_all(K,0);
		
		for (i=0;i<nm;++i)
		{
			for (j=0;j<nm;++j)
			{
				gsl_matrix_set(K, i, j, beta * ( pow(gsl_vector_get(eval,6),2) * gsl_matrix_get(GJ,0,i+3) * gsl_matrix_get(GJ,0,j+3) + pow(gsl_vector_get(eval,7),2) * gsl_matrix_get(GJ,1,i+3) * gsl_matrix_get(GJ,1,j+3) + pow(gsl_vector_get(eval,8),2) * gsl_matrix_get(GJ,2,i+3) * gsl_matrix_get(GJ,2,j+3) ));
			}
		}
		
		gsl_matrix_free(GJ);
		
		// Initialize A covariance matrix (amplitude covariances, equal to kappa_ij with correct diagonal)
		
		// printf("Set A\n");
		
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
			gsl_matrix_set(A, i, i, gsl_matrix_get(A,i,i) + 2 * beta * pow(gsl_vector_get(eval,i+9),2));
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
		 */
		
		// Make Cholesky decomposition of matrix A
		
		gsl_linalg_cholesky_decomp(A);
		
		// Get cholesky matrix diagonal (product of diagonal elements will give sqrt(detA))
		
		// printf("Set predet\n");
		
		gsl_vector *predet = gsl_vector_alloc(nm-3);
		
		for(i = 0; i < nm-3; ++i)
		{
			gsl_vector_set(predet, i, gsl_matrix_get(A,i,i));
		}
		
		gsl_sort_vector(predet);
		
		// Invert matrix A
		
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
		 */
		
		// Get XYZ covariances
		
		double KXX = 0;
		double KYY = 0;
		double KZZ = 0;
		double KXY = 0;
		double KXZ = 0;
		double KYZ = 0;
		
		// This part corresponds to sum(sum(k_ai * k_bj * sig_ij, j = 1...nm-3), i = 1...nm-3), or xi_ab (and 2 * xi_ab for XY, XZ and YZ)
		
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
		
		// This part substracts kappa_ab to xi_ab (2* for XY, XZ and YZ)
		
		KXX -= gsl_matrix_get(K,nm-3,nm-3);
		KYY -= gsl_matrix_get(K,nm-2,nm-2);
		KZZ -= gsl_matrix_get(K,nm-1,nm-1);
		KXY -= gsl_matrix_get(K,nm-3,nm-2) + gsl_matrix_get(K,nm-2,nm-3);
		KXZ -= gsl_matrix_get(K,nm-3,nm-1) + gsl_matrix_get(K,nm-1,nm-3);
		KYZ -= gsl_matrix_get(K,nm-2,nm-1) + gsl_matrix_get(K,nm-1,nm-2);
		
		gsl_matrix_free(K);
		
		// Inverts the sign of xi_ab - kappa_ab (in the equation, the covariances are kappa_ab - xi_ab)
		
		KXX *= -1;
		KYY *= -1;
		KZZ *= -1;
		KXY *= -1;
		KXZ *= -1;
		KYZ *= -1;
		
		// Gets the absolute value of det123, if needed
		
		if (det123 < 0) { det123 *= -1; }
		
		// Initialize Dens_Fac (d in the equation)
		
		long double Dens_Fac = sqrt( pow(beta/PI,3) ) / det123;
		
		// multiply Dens_Fac with (product of eigenvalues)/sqrt(detA)
		
		for(i = 0; i < nm-3; ++i)
		{
			Dens_Fac *= sqrt(2*beta)*gsl_vector_get(eval,i+9)/gsl_vector_get(predet,i);
		}
		
		Dens_Fac *= gsl_vector_get(eval,6)*gsl_vector_get(eval,7)*gsl_vector_get(eval,8);
		
		gsl_vector_free(predet);
		
		/*
		gsl_vector_set(entro,node,ConfEnt);
		
		
		
		//printf("\nConstants : \nq = %6.10f\ndet123 = %6.10f\ndetA = %6.10f\nDamping factor = %6.50Lf\nKXX = %6.10f\nKYY = %6.10f\nKZZ = %6.10f\nKXY = %6.10f\nKXZ = %6.10f\nKYZ = %6.10f\n\nDifferential entropy = %20.20Lf\n\n", q, det123, detA, Damping_Factor, KXX, KYY, KZZ, KXY, KXZ, KYZ, ConfEnt);
		*/
		
		printf("Node:%d\n",node);
		
		printf("\nDifferential entropy = %6.10f\n\n", 1.5 - log(Dens_Fac));
		
		// Assign density factor and entropy
		
		init[node].entro = 1.5 - log(Dens_Fac);
		
		init[node].dens = Dens_Fac;
		
		// Initialize inverse covariance matrix
		
		gsl_matrix *incov1 = gsl_matrix_alloc(3,3);
		gsl_matrix_set_all(incov1,0);
		
		gsl_matrix_set(incov1, 0, 0, 2*KXX);
		gsl_matrix_set(incov1, 1, 1, 2*KYY);
		gsl_matrix_set(incov1, 2, 2, 2*KZZ);
		gsl_matrix_set(incov1, 0, 1, KXY); gsl_matrix_set(incov1, 1, 0, KXY);
		gsl_matrix_set(incov1, 0, 2, KXZ); gsl_matrix_set(incov1, 2, 0, KXZ);
		gsl_matrix_set(incov1, 1, 2, KYZ); gsl_matrix_set(incov1, 2, 1, KYZ);
		
		/*
		printf("Inverse covariance matrix :\n");
		
		for(i = 0; i < 3; i++)
		{
			for(j = 0; j < 3; j++)
			{
				printf("%6.10f\t", gsl_matrix_get(incov1, i, j));
			}
			
			printf("\n");
		}
		
		printf("\n");
		*/
		
		// Invert inverse covariance matrix
		
		gsl_linalg_cholesky_decomp(incov1);
		
		gsl_linalg_cholesky_invert(incov1);
		
		/*
		printf("Covariance matrix :\n");
		
		for(i = 0; i < 3; i++)
		{
			for(j = 0; j < 3; j++)
			{
				printf("%6.10f\t", gsl_matrix_get(incov1, i, j));
			}
			
			printf("\n");
		}
		
		printf("\n");
		*/
		
		// Set workspace to get eigenvectors and eigenvalues
		
		gsl_matrix *glevecs = gsl_matrix_alloc(3,3);
		gsl_matrix_set_all(glevecs,0);
		
		gsl_vector *mainvars = gsl_vector_alloc(3);
		gsl_vector_set_all(mainvars,0);
		
		gsl_eigen_symmv_workspace *solver = gsl_eigen_symmv_alloc(3);
		
		gsl_eigen_symmv(incov1, mainvars, glevecs, solver);
		
		gsl_matrix_free(incov1);
		gsl_eigen_symmv_free(solver);
		
		/*
		printf("Eigenvector matrix :\n");
		
		for(i = 0; i < 3; i++)
		{
			for(j = 0; j < 3; j++)
			{
				printf("%6.10f\t", gsl_matrix_get(glevecs, i, j));
			}
			
			printf("\n");
		}
		
		printf("\n");
		*/
		
		/*
		printf("Eigenvalues :\n");
		
		for(i = 0; i < 3; i++)
		{
			printf("%6.10f\t", gsl_vector_get(mainvars, i));
		}
		
		printf("\n\n");
		*/
		
		// Assign global eigenvectors (main axes) and global variances to node
		
		for(j = 0; j < 3; j++)
		{
			for(i = 0; i < 3; i++)
			{
				init[node].global_evecs[i][j] = gsl_matrix_get(glevecs, i, j);
			}
			
			init[node].main_vars[j] = gsl_vector_get(mainvars, j);
		}
		
		/*
		printf("Weighted eigenvector matrix :\n");
		
		for(i = 0; i < 3; i++)
		{
			for(j = 0; j < 3; j++)
			{
				printf("%6.10f\t", init[node].global_evecs[i][j]);
			}
			
			printf("\n");
		}
		
		printf("\n");
		*/
		
		// printf("Done\n");
	}
	// printf("EXIT\n");
}

void conj_prob_init(struct pdb_atom *atm1, struct pdb_atom *atm2, gsl_matrix *incov12, gsl_vector *delr, double *conj_dens12) // Builds inverse conjugate covariance matrix, conjugate density factor and delta-r for use with over_prob and proxim_prob; 
{
	const double PI = 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421;
	
	int i,j,k;
	
	// Build the covariance matrix
	
	gsl_matrix *cov12 = gsl_matrix_alloc(3,3);
	
	for(i = 0; i < 3; i++)
	{
		for(j = 0; j < 3; j++)
		{
			double covar12 = 0;
			
			for(k = 0; k < 3; k++)
			{
				covar12 += atm1->global_evecs[i][k]*atm1->global_evecs[j][k]*atm1->main_vars[k] + atm2->global_evecs[i][k]*atm2->global_evecs[j][k]*atm2->main_vars[k];
			}
			
			gsl_matrix_set(cov12, i, j, covar12);
		}
	}
	
	// Make a Cholesky decomposition of the matrix and uses the decomposition to divide conj_dens12 by sqrt(det(cov12))
	
	gsl_linalg_cholesky_decomp(cov12);
	
	*conj_dens12 = 1/sqrt(pow(2*PI,3));
	
	for(i = 0; i < 3; i++)
	{
		*conj_dens12 /= gsl_matrix_get(cov12, i, i);
	}
	
	// Invert the matrix and stores it in incov12
	
	gsl_linalg_cholesky_invert(cov12);
	
	for(i = 0; i < 3; i++)
	{
		for(j = 0; j < 3; j++)
		{
			gsl_matrix_set(incov12, i, j, gsl_matrix_get(cov12, i, j));
		}
	}
	
	// Store delta_x, delta_y and delta_z in delr
	
	gsl_vector_set(delr, 0, atm2->x_cord - atm1->x_cord);
	gsl_vector_set(delr, 1, atm2->y_cord - atm1->y_cord);
	gsl_vector_set(delr, 2, atm2->z_cord - atm1->z_cord);
}

double proxim_prob(gsl_matrix *incov12, gsl_vector *delr, double conj_dens12, double minrad, double maxrad, int nsteps) // -----> Must have run conj_prob_init beforehand. <-------
{
	int i,j,k,r,s;
	
	double overprob = 0;
	
	// Ajouter un if(minrad == 0)
	/*
	double evalrange = 2*maxrad;
	
	for(i = 0; i < nsteps; i++)
	{
		for(j = 0; j < nsteps; j++)
		{
			for(k = 0; k < nsteps; k++)
			{
				if(sqrt(pow((2*i + 1)*evalrange/(2 * nsteps) - evalrange/2, 2) + pow((2*j + 1)*evalrange/(2 * nsteps) - evalrange/2, 2) + pow((2*k + 1)*evalrange/(2 * nsteps) - evalrange/2, 2)) <= radius)
				{
					gsl_vector *evpos = gsl_vector_alloc(3);
					gsl_vector_set_all(evpos, 0);
					
					gsl_vector_set(evpos, 0, gsl_vector_get(delr, 0) - evalrange/2 + (2*i + 1)*evalrange/(2 * nsteps));
					gsl_vector_set(evpos, 1, gsl_vector_get(delr, 1) - evalrange/2 + (2*j + 1)*evalrange/(2 * nsteps));
					gsl_vector_set(evpos, 2, gsl_vector_get(delr, 2) - evalrange/2 + (2*k + 1)*evalrange/(2 * nsteps));
					
					double argexp = 0;
					
					for(r = 0; r < 3; r++)
					{
						for(s = 0; s < 3; s++)
						{
							argexp += gsl_vector_get(evpos, r)*gsl_vector_get(evpos, s)*gsl_matrix_get(incov12, r, s);
						}
					}
					
					gsl_vector_free(evpos);
					
					argexp *= -0.5;
					
					overprob += conj_dens12*exp(argexp)*pow(evalrange/nsteps, 3);
				}
			}
		}
	}
	
	return overprob;
	*/
	
	const double reduct = 0.7071067811865475244008443621048490392848359376884740365883398689953662392310535194251937671638207864; // sqrt(PI)/2
	
	
	//Exceptions
	if(minrad == maxrad)
	{
		return 0;
	}
	else
	{
		if(minrad > maxrad)
		{
			return 0;
		}
	}
	
	if(minrad == 0)
	{
		double evalrange = 2*maxrad;
		
		for(i = 0; i < nsteps; i++)
		{
			for(j = 0; j < nsteps; j++)
			{
				for(k = 0; k < nsteps; k++)
				{
					if(sqrt(pow((2*i + 1)*evalrange/(2 * nsteps) - evalrange/2, 2) + pow((2*j + 1)*evalrange/(2 * nsteps) - evalrange/2, 2) + pow((2*k + 1)*evalrange/(2 * nsteps) - evalrange/2, 2)) <= maxrad)
					{
						gsl_vector *evpos = gsl_vector_alloc(3);
						gsl_vector_set_all(evpos, 0);
						
						gsl_vector_set(evpos, 0, gsl_vector_get(delr, 0) - evalrange/2 + (2*i + 1)*evalrange/(2 * nsteps));
						gsl_vector_set(evpos, 1, gsl_vector_get(delr, 1) - evalrange/2 + (2*j + 1)*evalrange/(2 * nsteps));
						gsl_vector_set(evpos, 2, gsl_vector_get(delr, 2) - evalrange/2 + (2*k + 1)*evalrange/(2 * nsteps));
						
						double argexp = 0;
						
						for(r = 0; r < 3; r++)
						{
							for(s = 0; s < 3; s++)
							{
								argexp += gsl_vector_get(evpos, r)*gsl_vector_get(evpos, s)*gsl_matrix_get(incov12, r, s);
							}
						}
						
						gsl_vector_free(evpos);
						
						argexp *= -0.5;
						
						overprob += conj_dens12*exp(argexp)*pow(evalrange/nsteps, 3);
					}
				}
			}
		}
		
		return overprob;
	}
	else
	{
		// Generate lower bound of evaluation
		double minev = minrad * reduct;
		
		
		// Generate all volume ratios of all three zones of evaluation where evaluation points will be scattered
		double rat_1 = (maxrad - minev) * pow(maxrad, 2) / (pow(maxrad, 3) - pow(minev, 3));
		double rat_2 = minev * (maxrad - minev) * maxrad / (pow(maxrad, 3) - pow(minev, 3));
		double rat_3 = pow(minev, 2) * (maxrad - minev) / (pow(maxrad, 3) - pow(minev, 3));
		
		//Calculate probability in the first zone (z planks, with X element of [delta_x - maxrad, delta_x + maxrad], y element of [delta_y - maxrad, delta_y + maxrad] and z element of [delta_z - maxrad, delta_z - minev] U [delta_z + minev, delta_z + maxrad])
		
		// printf("Part 1 :\n");
		
		int nz = round(pow(rat_1 * pow(nsteps, 3) * pow(maxrad - minev, 2) / pow(maxrad, 2), 1/3.0));
		
		if(nz < 1)
		{
			nz = 1;
		}
		
		// printf("nz : %1i\t\t", nz);
		
		int ny = round(sqrt(rat_1 * pow(nsteps, 3) / nz));
		
		if(ny < 1)
		{
			ny = 1;
		}
		
		int nx = ny;
		
		// printf("nx and ny : %1i\n", nx);
		
		for(i = 0; i < nx; i++)
		{
			for(j = 0; j < ny; j++)
			{
				for(k = 0; k < nz; k++)
				{
					if(sqrt(pow((2 * i + 1) * maxrad / nx - maxrad, 2) + pow((2 * j + 1) * maxrad / ny - maxrad, 2) + pow(minev + (2 * k + 1) * (maxrad - minev) / (2 * nz), 2)) > minrad && sqrt(pow((2 * i + 1) * maxrad / nsteps - maxrad, 2) + pow((2 * i + 1) * maxrad / nsteps - maxrad, 2) + pow(minev + (2 * k + 1) * (maxrad - minev) / (2 * nz), 2)) < maxrad)
					{
						gsl_vector *evpos = gsl_vector_alloc(3);
						gsl_vector_set_all(evpos, 0);
						
						gsl_vector_set(evpos, 0, gsl_vector_get(delr, 0) - maxrad + (2 * i + 1) * maxrad / nx);
						gsl_vector_set(evpos, 1, gsl_vector_get(delr, 1) - maxrad + (2 * j + 1) * maxrad / ny);
						gsl_vector_set(evpos, 2, gsl_vector_get(delr, 2) + minev + (2 * k + 1) * (maxrad - minev) / (2 * nz));
						
						double argexp = 0;
						
						for(r = 0; r < 3; r++)
						{
							for(s = 0; s < 3; s++)
							{
								argexp += gsl_vector_get(evpos, r)*gsl_vector_get(evpos, s)*gsl_matrix_get(incov12, r, s);
							}
						}
						
						argexp *= -0.5;
						
						overprob += conj_dens12*exp(argexp) * pow(2 * maxrad / nx, 2) * (maxrad - minev) / nz;
						
						gsl_vector_set(evpos, 0, gsl_vector_get(delr, 0) - maxrad + (2 * i + 1) * maxrad / nx);
						gsl_vector_set(evpos, 1, gsl_vector_get(delr, 1) - maxrad + (2 * j + 1) * maxrad / ny);
						gsl_vector_set(evpos, 2, gsl_vector_get(delr, 2) - minev - (2 * k + 1) * (maxrad - minev) / (2 * nz));
						
						argexp = 0;
						
						for(r = 0; r < 3; r++)
						{
							for(s = 0; s < 3; s++)
							{
								argexp += gsl_vector_get(evpos, r)*gsl_vector_get(evpos, s)*gsl_matrix_get(incov12, r, s);
							}
						}
						
						gsl_vector_free(evpos);
						
						argexp *= -0.5;
						
						overprob += conj_dens12*exp(argexp) * pow(2 * maxrad / nx, 2) * (maxrad - minev) / nz;
					}
				}
			}
		}
		
		//Calculate probability in the second zone (y sides, with X element of [delta_x - maxrad, delta_x + maxrad], y element of [delta_y - maxrad, delta_y - minev] U [delta_y + minev, delta_y + maxrad] and z element of [delta_z - minev, delta_z + minev])
		
		// printf("Part 2 :\n");
		
		nz = round(pow(rat_2 * pow(nsteps, 3) * pow(minev, 2) / (maxrad * (maxrad - minev)), 1/3.0));
		
		if(nz < 1)
		{
			nz = 1;
		}
		
		// printf("nz : %1i\t\t", nz);
		
		ny = round(nz * (maxrad - minev) / minev);
		
		if(ny < 1)
		{
			ny = 1;
		}
		
		// printf("ny : %1i\t\t", ny);
		
		nx = round(nz * maxrad / minev);
		
		if(nx < 1)
		{
			nx = 1;
		}
		
		// printf("nx : %1i\n", nx);
		
		for(i = 0; i < nx; i++)
		{
			for(j = 0; j < ny; j++)
			{
				for(k = 0; k < nz; k++)
				{
					if(sqrt(pow((2 * i + 1) * maxrad / nx - maxrad, 2) + pow(minev + (2 * j + 1) * (maxrad - minev) / (2 * ny), 2) + pow((2 * k + 1) * minev / nz - minev, 2)) > minrad && sqrt(pow((2 * i + 1) * maxrad / nx - maxrad, 2) + pow(minev + (2 * j + 1) * (maxrad - minev) / (2 * ny), 2) + pow((2 * k + 1) * minev / nz - minev, 2)) < maxrad)
					{
						gsl_vector *evpos = gsl_vector_alloc(3);
						gsl_vector_set_all(evpos, 0);
						
						gsl_vector_set(evpos, 0, gsl_vector_get(delr, 0) - maxrad + (2 * i + 1) * maxrad / nx);
						gsl_vector_set(evpos, 1, gsl_vector_get(delr, 1) + minev + (2 * j + 1) * (maxrad - minev) / (2 * ny));
						gsl_vector_set(evpos, 2, gsl_vector_get(delr, 2) - minev + (2 * k + 1) * minev / nz);
						
						double argexp = 0;
						
						for(r = 0; r < 3; r++)
						{
							for(s = 0; s < 3; s++)
							{
								argexp += gsl_vector_get(evpos, r)*gsl_vector_get(evpos, s)*gsl_matrix_get(incov12, r, s);
							}
						}
						
						argexp *= -0.5;
						
						overprob += conj_dens12*exp(argexp) * 2 * minev / nz * (maxrad - minev) / ny * 2 * maxrad / nx;
						
						gsl_vector_set(evpos, 0, gsl_vector_get(delr, 0) - maxrad + (2 * i + 1) * maxrad / nx);
						gsl_vector_set(evpos, 1, gsl_vector_get(delr, 1) - minev - (2 * j + 1) * (maxrad - minev) / (2 * ny));
						gsl_vector_set(evpos, 2, gsl_vector_get(delr, 2) - minev + (2 * k + 1) * minev / nz);
						
						argexp = 0;
						
						for(r = 0; r < 3; r++)
						{
							for(s = 0; s < 3; s++)
							{
								argexp += gsl_vector_get(evpos, r)*gsl_vector_get(evpos, s)*gsl_matrix_get(incov12, r, s);
							}
						}
						
						gsl_vector_free(evpos);
						
						argexp *= -0.5;
						
						overprob += conj_dens12*exp(argexp) * 2 * minev / nz * (maxrad - minev) / ny * 2 * maxrad / nx;
					}
				}
			}
		}
		
		//Calculate probability in the second zone (x plugs, with X element of [delta_x - maxrad, delta_x - minev] U [delta_x + minev, delta_x + maxrad], y element of [delta_y - minev, delta_y + minev] and z element of [delta_z - minev, delta_z + minev])
		
		// printf("Part 3 :\n");
		
		nz = round(pow(rat_3 * pow(nsteps, 3) * minev / (maxrad - minev), 1/3.0));
		
		if(nz < 1)
		{
			nz = 1;
		}
		
		// printf("nz and ny : %1i\t\t", nz);
		
		ny = nz;
		
		nx = round(nz * maxrad / minev);
		
		if(nx < 1)
		{
			nx = 1;
		}
		
		// printf("nx : %1i\n", nx);
		
		for(i = 0; i < nx; i++)
		{
			for(j = 0; j < ny; j++)
			{
				for(k = 0; k < nz; k++)
				{
					if(sqrt(pow((2 * i + 1) * maxrad / nx - maxrad, 2) + pow(minev + (2 * j + 1) * (maxrad - minev) / (2 * ny), 2) + pow((2 * k + 1) * minev / nz - minev, 2)) > minrad && sqrt(pow((2 * i + 1) * maxrad / nx - maxrad, 2) + pow(minev + (2 * j + 1) * (maxrad - minev) / (2 * ny), 2) + pow((2 * k + 1) * minev / nz - minev, 2)) < maxrad)
					{
						gsl_vector *evpos = gsl_vector_alloc(3);
						gsl_vector_set_all(evpos, 0);
						
						gsl_vector_set(evpos, 0, gsl_vector_get(delr, 0) + minev + (2 * i + 1) * (maxrad - minev) / (2 * nx));
						gsl_vector_set(evpos, 1, gsl_vector_get(delr, 1) - minev + (2 * j + 1) * minev / ny);
						gsl_vector_set(evpos, 2, gsl_vector_get(delr, 2) - minev + (2 * k + 1) * minev / nz);
						
						double argexp = 0;
						
						for(r = 0; r < 3; r++)
						{
							for(s = 0; s < 3; s++)
							{
								argexp += gsl_vector_get(evpos, r)*gsl_vector_get(evpos, s)*gsl_matrix_get(incov12, r, s);
							}
						}
						
						argexp *= -0.5;
						
						overprob += conj_dens12 * exp(argexp) * pow(2 * minev / nz, 2) * (maxrad - minev) / nx;
						
						gsl_vector_set(evpos, 0, gsl_vector_get(delr, 0) - minev - (2 * i + 1) * (maxrad - minev) / (2 * nx));
						gsl_vector_set(evpos, 1, gsl_vector_get(delr, 1) - minev + (2 * j + 1) * minev / ny);
						gsl_vector_set(evpos, 2, gsl_vector_get(delr, 2) - minev + (2 * k + 1) * minev / nz);
						
						argexp = 0;
						
						for(r = 0; r < 3; r++)
						{
							for(s = 0; s < 3; s++)
							{
								argexp += gsl_vector_get(evpos, r)*gsl_vector_get(evpos, s)*gsl_matrix_get(incov12, r, s);
							}
						}
						
						gsl_vector_free(evpos);
						
						argexp *= -0.5;
						
						overprob += conj_dens12 * exp(argexp) * pow(2 * minev / nz, 2) * (maxrad - minev) / nx;
					}
				}
			}
		}
		
		return overprob;
	}
}