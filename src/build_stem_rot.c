#include "STeM.h"

void add_vecteur(gsl_vector *a,gsl_vector *b,gsl_vector *c) {
	int i;
	
	for (i=0;i<3;++i) {
		gsl_vector_set(c,i,gsl_vector_get(a,i)+gsl_vector_get(b,i));
	}
}


int main(int argc, char *argv[]) {

	int all; 								// Nombre d'atoms dans le PDB
	int atom; 								// Nombre de nodes
	float epsilon = 0.36;					// Epsilon			
	float bond_factor = 100*epsilon;		// Facteur pour poid des bond strechcing
	float angle_factor = 20*epsilon;		// Facteur pour poid des angles
	double K_phi1 = epsilon;				// Facteurs pour angles dièdres
	double K_phi3 = 0.5*epsilon;
	long seed; seed = time_seed();			// Random
 	int help_flag = 1;						// Help Flag
	int i,j,k;								// Dummy Counter
 	char file_name[100];					// Nom du PDB
 	char templaate_name[100];				// Nom du fichier du templaate
 	char eigen_name[100] = "eigen.dat";		// Nom pour output des eigen	
 	float init_templaate = 1;
 	int verbose = 0;						// Loud
 	int templaate_flag = 0;					// Flag pour savoir si tempalte ou non
	int b_factor_flag = 0;					// Output b_factor flag
	int lig = 0;							// Flag pour savoir si tient compte ligand
	int nconn;								// Nombre de connect
	float kp_factor = 1;					// acteur pour poid des angles dièdres
	int enm = 0;							//flag pour enm instead of STeM
	float tot_factor = 1;
	int energy = 0;							// Flag pour printer l'energy
	float cutoff = 18.0;
	int torsion = 0;
 	for (i = 1;i < argc;i++) {
 		if (strcmp("-i",argv[i]) == 0) {strcpy(file_name,argv[i+1]);--help_flag;}
 		if (strcmp("-o",argv[i]) == 0) {strcpy(eigen_name,argv[i+1]);}
 		if (strcmp("-h",argv[i]) == 0) {help_flag = 1;}
 		if (strcmp("-v",argv[i]) == 0) {verbose = 1;}
 		if (strcmp("-t",argv[i]) == 0) {float temp;sscanf(argv[i+1],"%f",&temp);init_templaate = temp;}  
 		if (strcmp("-lt",argv[i]) == 0) {strcpy(templaate_name,argv[i+1]);templaate_flag = -1;}
 		if (strcmp("-kr",argv[i]) == 0) {float temp;sscanf(argv[i+1],"%f",&temp);bond_factor = temp * bond_factor;}
 		if (strcmp("-kt",argv[i]) == 0) {float temp;sscanf(argv[i+1],"%f",&temp);angle_factor = temp * angle_factor;}
 		if (strcmp("-kpf",argv[i]) == 0) {float temp;sscanf(argv[i+1],"%f",&temp); kp_factor = temp;}  
 		if (strcmp("-b",argv[i]) == 0) {b_factor_flag = 1;}
 		if (strcmp("-lig",argv[i]) == 0) {lig= 1;} 
 		if (strcmp("-enm",argv[i]) == 0) {enm= 1;}  
 		if (strcmp("-ee",argv[i]) == 0) {energy= 1;}
 		if (strcmp("-angle",argv[i]) == 0) {torsion= 1;}    
 		if (strcmp("-tot",argv[i]) == 0) {float temp;sscanf(argv[i+1],"%f",&temp); tot_factor = temp;}
 		if (strcmp("-cut",argv[i]) == 0) {float temp;sscanf(argv[i+1],"%f",&temp); cutoff = temp;}  
 		
 		
 	}

 	printf("\n********************\nFile:%s\nEpsilon:%f\nKr:\t%f\nKt:\t%f\nK phi1:\t%f\nK phi3:\t%f\n********************\n",file_name,epsilon,bond_factor,angle_factor,K_phi1,K_phi3);

 	if (help_flag == 0) { } else {
 		printf("****************************\nHelp Section\n-i\tFile Input (PDB)\n-oh\tOutput Name Hessian\n-oeig\tFile Name Eigen\n-v\tVerbose\n-w\tWeight Vector\n-t\tInitial Value of templaate (negative value for random)\n\tIf Load templaate, multiply the templaate\n-lt\tLoad input templaate\n-sp\tSuper Node Mode (CA, C et N)\n-kt\tPoid de l'angle entre les nodes (1)\n-kr\tPoid de la distance entre les nodes (1)\n-lig\tTient compte duligand (sauf HOH)\n-kp\tPoid de l'angles dièdres\n-kpf\tFacteur entre \n****************************\n");
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
	
	//write_strc("structure.pdb",strc_all,all);
	
	// Pour les besoins... on limite à 800 atomes
	
	if (atom > 800) {
		printf("Too much atom, if you want to remove the limite.... vincent.frappier@usherbrooke.ca\n");
		return(1);
	}
	
	//***************************************************
 	//*													*
 	//*Build Hessian									*
 	//*													*
 	//***************************************************


    double **hessian=(double **)malloc(3*atom*sizeof(double *)); // Matrix of the Hessian 1 2 3 (bond, angle, dihedral)
    for(i=0;i<3*atom;i++) { hessian[i]=(double *)malloc(3*atom*sizeof(double));}
    for(i=0;i<3*atom;i++)for(j=0;j<(3*atom);j++){hessian[i][j]=0;}
        
	gsl_matrix *h_matrix = gsl_matrix_alloc(3*atom, 3*atom);	/*Déclare une matrice hessian matrix de 3n par 3n*/
	gsl_matrix *templaate = gsl_matrix_alloc(3*atom, 3*atom);	/*Déclare une matrice hessian matrix de 3n par 3n*/
	
	//Setting templaate
	
	if (templaate_flag == -1) {
		printf("Loading templaate... Help flag:%d\n",help_flag);
		load_matrix(templaate,templaate_name);
		//for (i=0;i<atom;++i) {printf("%10.5f ",gsl_matrix_get(templaate,0,i));}printf("\n");
		gsl_matrix_scale (templaate, init_templaate);
	} else {
		if (init_templaate >= 0) {
			gsl_matrix_set_all (templaate, init_templaate/tot_factor);
		} else {
			for (i=0;i<atom;++i) {
				for (j=0;j<atom;++j) {
					float random =ran2(&seed);
					if (abs(i-j) < 4) {random = 0;}
					gsl_matrix_set(templaate,i,j,random);
					gsl_matrix_set(templaate,j,i,random);
				}
			} // Set random number if negatif number
		}
	}
	//write_matrix("out_templaate", templaate,atom);
	//Building Matrix
		
	if (verbose == 1) {printf("Building Hessian\n");}
	
	if (verbose == 1) {printf("\tCovalent Bond Potential\n");}		
	build_1st_matrix(strc_node,hessian,atom,bond_factor/tot_factor);
	
	if (verbose == 1) {printf("\tAngle Potential\n");}	
	build_2_matrix(strc_node,hessian,atom,angle_factor/tot_factor);
	
	if (verbose == 1) {printf("\tDihedral Potential\n");}	
	build_3_matrix(strc_node, hessian,atom,K_phi1/2+K_phi3*9/2,kp_factor/tot_factor);
	
	if (verbose == 1) {printf("\tNon Local Interaction Potential\n");}	
	build_4h_matrix(strc_node,hessian,atom,epsilon,templaate);
	
	if (enm == 1) {
		 for(i=0;i<3*atom;i++)for(j=0;j<(3*atom);j++){hessian[i][j]=0;}
		 build_enm(strc_node,hessian,atom,epsilon,templaate,cutoff);
	}
	
	if (verbose == 1) {printf("\tAssigning Array\n");}	
	assignArray(h_matrix,hessian,3*atom,3*atom);
	gsl_matrix *rot_hessian = gsl_matrix_alloc(2*atom,2*atom);
	if (torsion == 1) {
		// On restrict au torsion des angle... Alors Hessian = Jacobian^T * Hessian * Hessian ... Ce qui donne (2nx3n) * (3nX3n) * (3n*2n)
		printf("Building Jacobian\n");
		// Il faut construire la Jacobian
		// On va en construire une partie (parce que connais pas toute le gars la)
		
		gsl_matrix *jacob = gsl_matrix_alloc(3*atom,2*atom); // Matrice jacobian
		gsl_matrix *tjacob = gsl_matrix_alloc(2*atom,3*atom); // Matrice jacobian transposé
		gsl_matrix_set_all(jacob,0);
		gsl_matrix_set_all(tjacob,0);
		gsl_vector *ea = gsl_vector_alloc(3); // Vector qui defini l'axe de rotation
		gsl_vector *ri_sa = gsl_vector_alloc(3); // Vector entre les deux nodes
		gsl_vector *sa = gsl_vector_alloc(3); // Vector Sa
		gsl_vector *ra = gsl_vector_alloc(3); // Vector Ra
		gsl_vector *ta = gsl_vector_alloc(3); // Ta vector
		gsl_vector *R = gsl_vector_alloc(3); // Vector qui represente le centre de masse
		gsl_vector *temp_v = gsl_vector_alloc(3); // Vecteur temporaire !
		gsl_vector *temp_v2 = gsl_vector_alloc(3); // Vecteur temporaire !
		gsl_vector *sum = gsl_vector_alloc(3); // Vecteur qui se fait rajouter toute dedans
		gsl_vector *ri = gsl_vector_alloc(3); // Vecteur qui se fait rajouter toute dedans
		gsl_vector_set_all(sum,0);
		int bump = 0;
		int now_node = 0;
		
		
		
		// Build R		
				
		gsl_vector_set_all(R,0);
		for (i=0;i<atom;++i) {
			gsl_vector_set(R,0,gsl_vector_get(R,0)+strc_node[i].x_cord);
			gsl_vector_set(R,1,gsl_vector_get(R,1)+strc_node[i].y_cord);
			gsl_vector_set(R,2,gsl_vector_get(R,2)+strc_node[i].z_cord);
		}
		printf("R vector: ");
		print_vector(R);
		// Pour chacun des angles
		
		for (i=0;i<atom*2;++i) { 
			// Pour chacun des nodes (atomes)
			for (j=0;j<atom;++j) {
				printf("Angle:%d Node:%d\n",i,j);
				// Probleme d'indexation de i
				now_node = i/2;
				
				// Reset le vecteurs a addiotionner
				gsl_vector_set_all(sum,0);
				
				// Si bump = 0... Il s'agit de N-CA, si bump = 1, il s'agit de CA-C
				if (i % 2 == 0) {bump = 0;} else {bump = 1;}
								
				// On veut implementer Jia = ea * (ri-sa)xia + A x ri + ta
				
				// Première partie (ea * ri-sa)
				
				if (bump == 0) {
					assign_vector(strc_all,i," CA  ",all,strc_all,i," N   ",all,ea);
				} else {
					assign_vector(strc_all,i," CA  ",all,strc_all,i," C   ",all,ea);
				}
				
				// Build ri-sa
				
				gsl_vector_set(ri_sa,0,strc_node[j].x_cord-strc_node[now_node].x_cord);
				gsl_vector_set(ri_sa,1,strc_node[j].y_cord-strc_node[now_node].y_cord);
				gsl_vector_set(ri_sa,1,strc_node[j].z_cord-strc_node[now_node].z_cord);
				
				// Facteur Xia
				
				if (bump == 0 && j > now_node) {gsl_vector_scale(ri_sa,0);}
				if (bump == 1 && j < now_node) {gsl_vector_scale(ri_sa,0);}
							
				// Add la première partie
				dot_product_v(ea,ri_sa,sum);
				
				// Build Ta => ta = -ea*R - (Ma/M)ea*(Ra-sa)
				
				// -ea*R = temp
				gsl_vector_set_all(temp_v,0);
				add_vecteur(temp_v,ea,temp_v);
				

				gsl_vector_scale(temp_v,-1);
				dot_product_v(temp_v,R,temp_v);
				
				// (Ma/M) ea = temp_v2
				
				gsl_vector_set_all(temp_v2,0);
				add_vecteur(temp_v2,ea,temp_v2);
				float factor = -1;
				if (bump == 0) {
					factor = float(now_node)/float(atom);
				} else {
					factor = float(atom-(now_node+1))/float(atom);
				}
				gsl_vector_scale(temp_v2,factor);
				// (Ra - Sa)
				
				gsl_vector_set(sa,0,-strc_node[now_node].x_cord);
				gsl_vector_set(sa,1,-strc_node[now_node].y_cord);
				gsl_vector_set(sa,2,-strc_node[now_node].z_cord);
				
				// Ra la somme de tous les vecteurs d'après, on storer Ra-sa dans sa
				gsl_vector_set_all(ra,0);
				for (k=0;k<atom;++k) {
					if (bump == 0 && k >= now_node) {continue;}
					if (bump == 1 && k <= now_node) {continue;}
					gsl_vector_set(ra,0,strc_node[k].x_cord+gsl_vector_get(ra,0));
					gsl_vector_set(ra,1,strc_node[k].y_cord+gsl_vector_get(ra,1));
					gsl_vector_set(ra,2,strc_node[k].z_cord+gsl_vector_get(ra,2));
				
				}

				add_vecteur(ra,sa,sa);

				// (Ma/M)ea*(Ra-sa) = temp_v2*sa qui va être storer dans temp_v2
				
				dot_product_v(temp_v2,sa,temp_v2);
				
				// Ensuite on veut -ea*R - temp_v2 (qui est l'autre partie de l'équation)... 

				gsl_vector_scale(temp_v2,-1.0);
				add_vecteur(temp_v2,temp_v,ta);
				
				// On ajoute ta a sum
				
				add_vecteur(ta,sum,sum);
				
				// On veut builder Aa*ri
				
				// Ri est le vecteur qui correspond à la position du node à comparer
				
				gsl_vector_set(ri,0,strc_node[j].x_cord);
				gsl_vector_set(ri,1,strc_node[j].y_cord);
				gsl_vector_set(ri,2,strc_node[j].z_cord);
				
				// Aa est plus tricky, il s'agit de: Ia*ea+Ma(R-Ra)*(ea*(Ra-sa))
				
				//Ra-sa était déja calculé dans sa et ea aussi, alors on va sroter ea*(ra-sa) dans temp_v
				
				dot_product_v(ea,sa,temp_v);
				
				// R-Ra va etre storer dans temp_v2
				
				gsl_vector_set_all(temp_v2,0);
				add_vecteur(ra,temp_v2,temp_v2);
				gsl_vector_scale(temp_v2,-1.0);
				add_vecteur(temp_v2,R,temp_v2);
				
				// Ma est déja en quelques sortes storer dans factor (Ma/M)
				
				gsl_vector_scale(temp_v2,factor*float(atom));
				
				// On a alors, Ia*ea+temp_v2*temp_v
				
				dot_product_v(temp_v2,temp_v,temp_v);
				
				// Ia*ea+temp_v => ce qui reste
				
				// Je sais pas trop c est quoi Ia, alors, on laisse faire pour maintenant !
				
				add_vecteur(ea,temp_v,temp_v);
				
				// Aa = temp_v, il faut faire le dot product de Aa *ri
				
				dot_product_v(temp_v,ri,temp_v);
				
				// Rajoute à sum
				
				add_vecteur(sum,temp_v,sum);
				
				// Set Vecteur
				for (k=0;k<3;++k) {		
					gsl_matrix_set(jacob,j*3+k,i,gsl_vector_get(sum,k));
					gsl_matrix_set(tjacob,i,j*3+k,gsl_vector_get(sum,k));
				}
				
				
				
			}		
		}
		write_matrix("jacob.dat", jacob,3*atom,2*atom);
			
		// Transform l'hessian
		//gsl_matrix *rot_hessian = gsl_matrix_alloc(2*atom,2*atom);
		gsl_matrix *temp_hessian = gsl_matrix_alloc(2*atom,3*atom);
		multiplie_matrix(tjacob,2*atom,3*atom,h_matrix,3*atom,3*atom,temp_hessian);
		multiplie_matrix(temp_hessian,2*atom,3*atom,jacob,3*atom,2*atom,rot_hessian);
		for(i=0;i<2*atom;++i) {
			for(j=0;j<2*atom;++j) {
				if (gsl_matrix_get(rot_hessian,i,j) < 0.00001) {
					gsl_matrix_set(rot_hessian,i,j,0.0);
				}
			}
		}
		write_matrix("hessian.dat",rot_hessian,2*atom,2*atom);
		
	}
	
	//***************************************************
	//*													*
	//*Diagonalyse the matrix							*
	//*													*
	//***************************************************
	
	
	if (verbose == 1) {printf("Diagonalizing Hessian\n");}
	gsl_vector *eval = gsl_vector_alloc((3-torsion)*atom); /*Déclare un vector qui va contenir les eigenvalue */	
	gsl_matrix *evec= gsl_matrix_alloc ((3-torsion)*atom,(3-torsion)*atom); /*Déclare une matrice qui va contenir les eigenvector correspondant à eval*/
	if (torsion == 0) {
		diagonalyse_matrix(h_matrix,(3-torsion)*atom,eval,evec); /*Diagonalyse la matrice cartesienne*/
	} else {
		diagonalyse_matrix(rot_hessian,(3-torsion)*atom,eval,evec); /*Diagonalyse la matrice cartesienne*/
	
	}
	write_eigen(eigen_name,evec,eval,(3-torsion)*atom);
	if (energy == 1) {printf("Energy:%10.5f\n",calc_energy(atom,eval));}
	 
	//***************************************************
	//*													*
	//*Inverse matrix									*
	//*													*
	//***************************************************
	 if (b_factor_flag == 1) {
		 if (verbose == 1) {printf("Inversing Matrix\n");}
		gsl_matrix *k_inverse = gsl_matrix_alloc(atom, atom); /*Déclare et crée une matrice qui va être le pseudo inverse*/
		k_inverse_matrix_stem(k_inverse,atom,eval,evec);
	
		//write_matrix("b_factor.dat", k_inverse,atom,atom);

		printf("Correlation:%f\n",correlate(k_inverse,strc_node, atom));
		/*for (i=0;i<all;++i) {
			strc_all[i].b_factor = gsl_matrix_get(k_inverse,strc_all[i].node,strc_all[i].node);
		}*/
		//write_strc("b_factor.pdb",strc_all,all);
		//if (super_node != 3) {write_strc_b("b_factor.pdb",strc_all,all,k_inverse,super_node);}
		gsl_matrix_free(k_inverse);
		
	}
	
	// Freeing
	
	gsl_matrix_free(h_matrix);
	gsl_matrix_free(templaate);
	gsl_vector_free(eval);
	gsl_matrix_free(evec);
	
	for(i=0;i<3*atom;i++){free(hessian[i]);}
	free(hessian);
	for(i=0;i<nconn;i++) {free(connect_h[i]);}
	free(connect_h);
}


