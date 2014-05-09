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
	int  ipos = 298;
	int npair = 2;
 	for (i = 1;i < argc;i++) {
 		if (strcmp("-i",argv[i]) == 0) {strcpy(file_name,argv[i+1]);--help_flag;}
 		if (strcmp("-h",argv[i]) == 0) {help_flag = 1;}
 		if (strcmp("-no_const",argv[i]) == 0) {constraint_flag = 0;}
 		if (strcmp("-v",argv[i]) == 0) {verbose = 1;}
		if (strcmp("-ieig",argv[i]) == 0) {strcpy(eigen_name,argv[i+1]);}
		if (strcmp("-teig",argv[i]) == 0) {strcpy(eigen_name_two,argv[i+1]);}
		if (strcmp("-pos",argv[i]) == 0) {int temp;sscanf(argv[i+1],"%d",&temp); ipos = temp;}
 		if (strcmp("-max_rmsd",argv[i]) == 0) {float temp;sscanf(argv[i+1],"%f",&temp);rmsd_cutoff = temp;}
 		if (strcmp("-size",argv[i]) == 0) {int temp;sscanf(argv[i+1],"%d",&temp); npair = temp-1;}
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
	
	//***************************************************
	//*                                                                                                     *
	//* Load eigenvector et eigenvalue									*
	//*                                                                                                     *
	//***************************************************
	
	if (verbose == 1) {printf("Loading Eigenvector One\n");}
	
	gsl_vector *eval = gsl_vector_alloc(3*atom);
	
	gsl_matrix *evec= gsl_matrix_alloc (3*atom,3*atom);
	
	load_eigen(eval,evec,eigen_name,3*atom);
	
	if (verbose == 1) {printf("Loading Eigenvector Two\n");}
	
	gsl_vector *eval_two = gsl_vector_alloc(3*atom_t);
	
	gsl_matrix *evec_two= gsl_matrix_alloc (3*atom_t,3*atom_t);
	
	load_eigen(eval_two,evec_two,eigen_name_two,3*atom_t);	
	
	gsl_matrix *k_totinv_two = gsl_matrix_alloc(atom_t*3, atom_t*3); /* Déclare et crée une matrice qui va être le pseudo inverse */
	gsl_matrix *k_totinv_two_cpy = gsl_matrix_alloc(atom_t*3, atom_t*3); /* Déclare et crée une matrice qui va être le pseudo inverse */
	gsl_matrix *done_two = gsl_matrix_alloc(atom_t*3, atom_t*3);
	gsl_matrix_set_all(k_totinv_two, 0);
	gsl_matrix_set_all(done_two, -1);
	
	int align_t[atom_t];
	//inv_mat_align(k_totinv_two,atom_t,eval_two,evec_two,0,atom_t*3,align_t);
	k_cov_inv_matrix_stem(k_totinv_two,atom_t,eval_two,evec_two,6,atom_t*3); // Génère une matrice contenant les superéléments diagonaux de la pseudo-inverse. 
	
	gsl_matrix *k_totinv = gsl_matrix_alloc(atom*3, atom*3); /* Déclare et crée une matrice qui va être le pseudo inverse */
	gsl_matrix *k_totinv_cpy = gsl_matrix_alloc(atom*3, atom*3); /* Déclare et crée une matrice qui va être le pseudo inverse */
	gsl_matrix *done = gsl_matrix_alloc(atom*3, atom*3);
	gsl_matrix_set_all(k_totinv, 0);
	gsl_matrix_set_all(done, -1);
	//inv_mat_align(k_totinv,atom,eval,evec,0,atom*3,align,done); /* Génère une matrice contenant les superéléments diagonaux de la pseudo-inverse. */
	k_cov_inv_matrix_stem(k_totinv,atom,eval,evec,0,atom*3);
	// Copy la matrice
	gsl_matrix_memcpy (k_totinv_two_cpy,k_totinv_two);
	gsl_matrix_memcpy (k_totinv_cpy,k_totinv);
	
	
	// Set the scale
	
	int ItScale = 0;
	int NbScale = 0;
	float scale[100];
	float base;
	float power;
	for(power = -5;power<1;++power) {
		for(base = 0;base<4;++base) {
			scale[NbScale] = pow(2,base)*pow(10,power);
			++NbScale;
		}
	}

	
	
	// Look at pair in contact
	printf("Pair:%d %s%d%s\n", ipos,strc_node[ ipos].res_type,strc_node[ ipos].res_number,strc_node[ ipos].chain);
	
	if (bb_sc( ipos,strc_all,all) != 1 && bb_sc( ipos+1,strc_all,all) != 0) {printf("Pair is not backbone\n");return(0);}
	int count = 0;
	int init_pair[100];
	int init_sele[100]; // Si le residus de la pair mater
	
	int align[atom]; // Aussi on va les store dans align pour l'inversement
	for(i=0;i<atom;++i) {align[i] = 1;}
	
	
	for(i =0;i<atom-1;++i) {
			if (i ==  ipos+1) {continue;}
			if (i ==  ipos) {continue;}
			//if (bb_sc(i,strc_all,all) != 0) {continue;}
			if (bb_sc(i,strc_all,all) != 1) {continue;}
			if (gsl_matrix_get(contact, ipos+1,i)+gsl_matrix_get(contact, ipos+1,i+1) < 5) {continue;}
			printf("Ori:%d I:%d Template:%f et %f Type:%d\n", ipos+1,i,gsl_matrix_get(contact, ipos+1,i),gsl_matrix_get(contact, ipos+1,i+1),bb_sc(i,strc_all,all));
			align[i] = 1;
			align[i+1] = 1;
			init_pair[count] = i;
			if (gsl_matrix_get(contact, ipos+1,i+1) < 5) {
				init_sele[count] = 0;
			} else {
				init_sele[count] = 1;
			}
			++count;
			//if (count > 4) {break;}
	}
	
	for (i =0;i<count;++i) {
		printf("I:%d %d %d\n",i,init_pair[i],init_sele[i]);
	}
	
	// Bon il faut dequoi pour builder les pair, triplet, whatever
	
	int ipair[10000][(npair+1)*2][2];for(i=0;i<10000;++i){for(j=0;j<(npair+1)*2;++j) {ipair[i][j][0] = -1;ipair[i][j][1] = -1;}}
	int ncount[npair];for(j=0;j<npair;++j) {ncount[j] = j;}
	printf("\n");
	for(i=0;i<10000;++i) {
	
		// Populate ipair
		// ipos (position a muter)
		ipair[i][0][0] = ipos;
		ipair[i][1][0] = ipos+1;
		for(j=0;j<npair;++j) {
			ipair[i][(j+1)*2][0] = init_pair[ncount[j]];
			ipair[i][(j+1)*2+1][0] = init_pair[ncount[j]]+1;
			ipair[i][(j+1)*2][1] = init_sele[ncount[j]];
			ipair[i][(j+1)*2+1][1] = init_sele[ncount[j]];
			//printf("ipair[%d][%d][0] = %d\n",i,(j+1)*2,init_pair[ncount[j]]);
			//printf("ipair[%d][%d][0] = %d\n",i,(j+1)*2+1,1+init_pair[ncount[j]]);
			//printf("%d ",init_pair[ncount[j]]);
		}		//printf("\n");
		
		// Non-redundant pair
		ncount[npair-1] += 1;
		//
		//int l;printf("\t");for(l=0;l<npair;++l) {printf("%d ",ncount[l]);}		printf("\n");
		int k;
		for (j=1;j<npair;++j) {
			if (ncount[npair-j] > count-j) {
			
			
				// Increment le niveau en dessous
				if (npair-j-1 >= 0) {
				
					ncount[npair-j-1] += 1;
					// Reset les niveaux au dessus
					for(k=0;k<npair;++k) {
						ncount[npair-j+k] = ncount[npair-j-1+k]+1;
						//printf("\t");for(l=0;l<npair;++l) {printf("%d ",ncount[l]);}		printf("\n");
					}
					
				}
			}
			
		}
		//printf("\t%d > %d et %d\n",ncount[0], count-npair,ncount[npair-1]);
		if (ncount[0] > count-npair) {break;}
	}	
	printf("My Init strc have:%d pair\n",i);
	int tot_init_pair = i;
	// Bon j'ai toute mes pair,cluster, non-redundant dans init_pair et savoir s'il sont amino acid type dependant
	

	
	
  // Screener toute les possibilite de node, et verifier qu'il on un RMSD acceptable et respecte les restrictions de natures d'acide amine
	for(j=0;j<npair;++j) {ncount[j] = 0;}
	
	//atom_t = 10;
	float dist[atom_t/2][atom_t/2];
	for(i=0;i<atom_t/2;++i) {
		for(j=0;j<atom_t/2;++j) {
			dist[i][j] = -1;
		}
	}
	int myPair;
	struct pdb_atom strc_node_cpy[atom];
	for (myPair=0;myPair<tot_init_pair;++myPair) {
		int l;
		for(i=0;i<atom_t/2;++i) {
			for(j=0;j<npair;++j) {ncount[j] = 0;}
			int k;
			while(1) {
			//for (k=0;k<pow(npair,atom_t);++k) {
		
				// Increment le ncount
		
				ncount[npair-1] += 1;
				int m;
				for (m=1;m<npair;++m) {
					//printf("ncount[%d-%d] > %d\n",npair,m,atom_t/2);
					if (ncount[npair-m] > (atom_t-1)/2) {
						ncount[npair-m] = 0;
						ncount[npair-m-1] += 1;
					}
				}
				if (ncount[0] >atom_t/2-1){break;}
			
				// Can't have same node twice or be the I node
				int redun = 0;
				for (m=0;m<npair;++m) {
				
					for (l=m+1;l<npair;++l) {
						if (ncount[m] == ncount[l]) {++redun;break;}
						if (ncount[m] == i) {++redun;break;}
						if (ncount[l] == i) {++redun;break;}
					}
					if(redun != 0) {break;}
				}
				if (redun != 0) {continue;}
			
				// Va regarder si respecte sequence
				int seq = 0;
				for(l=0;l<npair;++l) {
					if (seq != 0) {break;}
					//printf("		L:%d -> %d Init_seq:%s vs Targ_seq:%s\n",l,ncount[l]*2,strc_node[ipair[myPair][(l+1)*2][0]].res_type,strc_node_t[ncount[l]*2].res_type);
					if (ipair[myPair][(l+1)*2][1] == 1) {
						if (strncmp(strc_node[ipair[myPair][(l+1)*2][0]].res_type,strc_node_t[ncount[l]*2].res_type,3) != 0) {
							if (l != npair-1) {++ncount[l];}
							for(m=l+1;m<npair;++m) {
								ncount[m] = 0;
							}
							++seq;
						}
					
					}
			
				}
			
				if (seq != 0) {continue;}
			
			
			
				// Regarde la distance du node i, si plus grand que 15, next
			
				int dist_flag = 0;
				for(l=0;l<npair;++l) {
				
					// Si la distance n'est pas defini, on la defini
					if (dist[i][ncount[l]] == -1) {
						dist[i][ncount[l]] = pow(strc_node_t[ncount[l]*2].x_cord-strc_node_t[i*2].x_cord,2)+pow(strc_node_t[ncount[l]*2].y_cord-strc_node_t[i*2].y_cord,2)+pow(strc_node_t[ncount[l]*2].z_cord-strc_node_t[i*2].z_cord,2);
						dist[ncount[l]][i] = dist[i][ncount[l]];
					}
					if (dist[i][ncount[l]] > 15*15) {
						//printf("L:%d Dist:%f\n",l,dist[i][ncount[l]]);
						++dist_flag;
						if (l != npair-1) {++ncount[l];}
						for(m=l+1;m<npair;++m) {
								ncount[m] = 0;
							}
					}
					if (dist_flag != 0) {break;}
				}
				if (dist_flag != 0) {continue;}
			
			
			
			
				// On va regarder le RMSD
				// Reset le align
				for(j=0;j<atom;++j) {align[j] = -1;}
				for(j=0;j<atom_t;++j) {align_t[j] = 1;}
				// Le premier node doit tjrs matcher
				align[ipos] = i;
				align[ipos+1] = i+1;
				for(l=0;l<npair;++l) {
					align[ipair[myPair][(l+1)*2][0]] = ncount[l]*2;
					align[ipair[myPair][(l+1)*2+1][0]] = ncount[l]*2+1;		
					
					align_t[ncount[l]*2] = ipair[myPair][(l+1)*2][0];
					align_t[ncount[l]*2+1] = ipair[myPair][(l+1)*2+1][0];
				}
				// Look at RMSD et si seq peut exister
			
				float myRmsd = (rmsd_no(strc_node,strc_node_t,atom, align));
				
				for(ItScale=0;ItScale<NbScale;++ItScale) {
				
					if (myRmsd > 4) {continue;}
				
			
					gsl_matrix_memcpy(k_totinv_two,k_totinv_two_cpy);
					gsl_matrix_memcpy(k_totinv,k_totinv_cpy);

					// Some scaling happen
				
					gsl_matrix_scale(k_totinv_two,1.00/scale[ItScale]);
					gsl_matrix_scale(k_totinv,1.00/scale[ItScale]);

					// Need to calculate density !!!
					int score = (npair+1)*2;
					gsl_matrix *sub_covar = gsl_matrix_alloc(score*3,score*3);
					gsl_matrix_set_all(sub_covar,0);
			
					
					copy_strc( strc_node_cpy, strc_node,  atom);

					rmsd_yes_covar(strc_node_cpy,strc_node_t,atom, align,k_totinv,sub_covar);
					//print_matrix(sub_covar);
					if (verbose == 1) {printf("\tAdding Two Cross-Correlation Matrix\n");}	
					// Variable importante pour mon inversion
					gsl_matrix *incov12 = gsl_matrix_alloc(score*3,score*3);
					double conj_dens12 = 1;
					gsl_matrix *cov12 = gsl_matrix_alloc(score*3,score*3);
					gsl_matrix_set_all(cov12,0);
					// Build mix des deux cov, il faut les additioner
					int count = 0;
	
					gsl_vector *delr = gsl_vector_alloc(score*3); // Differene entre deux positions
					gsl_vector *pos = gsl_vector_alloc(score*3); // Le vecteur de différence qu'on va vouloir evaluer
	
					gsl_vector_set_all(pos,0);
					int kept_two[score];
					for(l=0;l<atom;++l) {
						if (align[l] != -1) {
							kept_two[count] = align[l];
							gsl_vector_set(delr,count*3+0,strc_node_cpy[l].x_cord-strc_node_t[align[l]].x_cord);
							gsl_vector_set(delr,count*3+1,strc_node_cpy[l].y_cord-strc_node_t[align[l]].y_cord);
							gsl_vector_set(delr,count*3+2,strc_node_cpy[l].z_cord-strc_node_t[align[l]].z_cord);
							++count;
						}
					}
				
					for(m = 0;m<score;++m) {
						for(j = 0;j<score;++j) {
							for (l = 0;l<3;++l) {for (k = 0;k<3;++k) {
							//printf("COV12:(%d,%d) = %f + %f\n",m*3+k,j*3+l,gsl_matrix_get(sub_covar    ,m*3+k,j*3+l),gsl_matrix_get(k_totinv_two,kept_two[m]*3+k,kept_two[j]*3+l));
							gsl_matrix_set(cov12,m*3+k,j*3+l,
								 gsl_matrix_get(sub_covar    ,m*3+k,j*3+l)
								+gsl_matrix_get(k_totinv_two,kept_two[m]*3+k,kept_two[j]*3+l)
							);
							}}
	
						}
	
					}
				
					gsl_vector *eval_cov = gsl_vector_alloc(3*score); 
					gsl_vector_set_all(eval_cov,0);
					gsl_matrix *evec_cov = gsl_matrix_alloc (3*score,3*score); 
					gsl_matrix_set_all(evec_cov,0);
					diagonalyse_matrix(cov12,3*score,eval_cov,evec_cov);
				
					const double PI = 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421;
					k_cov_inv_matrix_stem(incov12,score,eval_cov,evec_cov,0,3*score);
	
	
	
					conj_dens12= -0.5*score*3*log(2*PI);
					for(m=0;m<score*3;++m) {
						 //if (m < 30) {printf("I:%d -> %g Log -> %g\n",m, gsl_vector_get (eval_cov, m),log(gsl_vector_get (eval_cov, m)));}
						 if  (gsl_vector_get (eval_cov, m) < 0.000001) {continue;}
		
						 conj_dens12 -= 0.5*log(gsl_vector_get (eval_cov, m));
					}
					float dens = density_prob_n(incov12, delr, conj_dens12,pos,score*3);
				
					printf("%d%s ",strc_node_t[i*2].res_number,strc_node_t[i*2].res_type);for(l=0;l<npair;++l) {printf("%d%s ",strc_node_t[ncount[l]*2].res_number,strc_node_t[ncount[l]*2].res_type);}	
					printf("et %d%s ",strc_node[ipos].res_number,strc_node[ipos].res_type);
					for(l=0;l<npair;++l) {
						printf("%d%s%d ",strc_node[ipair[myPair][(l+1)*2][0]].res_number,strc_node[ipair[myPair][(l+1)*2][0]].res_type,ipair[myPair][(l+1)*2][1]);
					}		
					printf("%f %.4f %g\n",scale[ItScale],sqrt(myRmsd),dens);
					gsl_matrix_free(sub_covar);
					gsl_matrix_free(evec_cov);		
					gsl_matrix_free(cov12);
					gsl_matrix_free(incov12);
		
					gsl_vector_free(eval_cov);
					gsl_vector_free(pos);
					gsl_vector_free(delr);
				}
			}
		}
	}
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
	gsl_matrix_set_all (m, 0);
	int i,j,k;
	 
	for (i=0;i<nb_atom*3;++i)	{
		if (align[i/3] == -1) {continue;}
		for (j=0;j<nb_atom*3;++j) {
			if (align[j/3] == -1) {continue;}
			//if (gsl_matrix_get(done, i,j) > 0) {continue;}
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

