#include "STeM.h"

double covar(float *a,float *b,int atom) {
	
	double moya = 0.0;
	double moyb = 0.0;
	double tot = 0.0;
	int i;
	int count = 0;
	// Moyenne de a et b
	
	for(i=0;i<atom;++i) {
		if (a[i] < -9999.0 || b[i] < -9999.0) {continue;}
		++count;
		moya += a[i];
		moyb += b[i];
		//printf("	%f::%f\n",a[i],b[i]);
	}
	moya /= count;
	moyb /= count;
	//printf("Moy:%f %f\n",moya,moyb);
	
	for(i=0;i<atom;++i) {
		if (a[i] < -9999.0 || b[i] < -9999.0) {continue;}
		tot += ( a[i]-moya)*(b[i]-moyb);
		//printf("	Tot:%f += (%f-%f)*(%f-%f) = %f\n",tot,a[i],moya,b[i],moyb,( a[i]-moya)*(b[i]-moyb));
	}
	
	tot /= (count-1);
	return(tot);
}

int main(int argc, char *argv[]) {
	int i,j,k;
	//int prot = 0;
	
	int nbrAT[argc-1];
	int nconn[argc-1];
	int nbrCA[argc-1];
	//int* align[argc-1];
	struct pdb_atom* allpdb[argc-1];
	struct pdb_atom* allall[argc-1];
	
	for (i = 1;i < argc;i++) {
		printf("I:%d/%d %s ",i,argc,argv[i]);
		// Build Strc principale, only CA
		
		nbrAT[i-1] = count_atom(argv[i]);
		nconn[i-1] = count_connect(argv[i]);
		
		int **connect_h=(int **)malloc(nconn[i-1]*sizeof(int *)); 
    	for(k=0;k<nconn[i-1];k++) { connect_h[k]=(int *)malloc(6*sizeof(int));}
		
		assign_connect(argv[i],connect_h);
		
		struct pdb_atom* strc_all;
		
		strc_all = (pdb_atom*)malloc(nbrAT[i-1]*sizeof(pdb_atom));
		if (strc_all == NULL) {return(1);}
		
		nbrCA[i-1] = build_all_strc(argv[i],strc_all);
		check_lig(strc_all,connect_h,nconn[i-1],nbrAT[i-1]);
		
		struct pdb_atom* strc_node;
		strc_node = (pdb_atom*)malloc(nbrCA[i-1]*sizeof(pdb_atom));
		if (strc_node == NULL) {return(1);}
		
		nbrCA[i-1] = build_cord_CA(strc_all, strc_node,nbrAT[i-1],0,connect_h,nconn[i-1]);
		
		printf("%d %d\n",nbrCA[i-1],nbrAT[i-1]);
		allpdb[i-1] = strc_node;
		
		allall[i-1] = strc_all;
		
		
		
	}

	// Build total align
	printf("ALIGNING\n");
	int align[argc-2][1000];
	for(i=1;i<argc-1;++i) {		
		int align_temp[nbrCA[i]];
		int score = node_align_low(allpdb[0],nbrCA[0],allpdb[i],nbrCA[i],align[i-1]);
		node_align_low(allpdb[i],nbrCA[i],allpdb[0],nbrCA[0],align_temp);
 		printf("I:%d %s RMSD:%8.5f Score: %d/%d\n",i,argv[i+1],sqrt(rmsd_yes(allpdb[i],allpdb[0],nbrCA[i], align_temp,allall[i],nbrAT[i])),score,nbrCA[0]);
 		/*char temp[50];
 		sprintf(temp,"strc_%d.pdb",i);
		write_strc(temp, allall[i],nbrAT[i]);
		write_strc("strc_0.pdb", allall[0],nbrAT[0]);*/
		
	}
	
	// Build Array the coord
	int master_align[nbrCA[0]];
	int index = 0;
	float array[nbrCA[0]*3][argc-1];
	for(i=0;i<nbrCA[0];++i) {
		//printf("I:%d %3d%s",i,allpdb[0][i].res_number,allpdb[0][i].res_type);
		k=0;
		for (j=0;j<argc-2;++j) {
			//printf("::%3d%s ",allpdb[j+1][align[j][i]].res_number,allpdb[j+1][align[j][i]].res_type);
			if (align[j][i] == -1) {++k;}
		}	
		//printf("K:%d\n",k);
		// Hard coding the caca pour prendre data IVET
		/*if ( (allpdb[0][i].res_number > 4    && allpdb[0][i].res_number < 32 ) ||
		     (allpdb[0][i].res_number  > 35  && allpdb[0][i].res_number < 117) ||
		     (allpdb[0][i].res_number  > 120 && allpdb[0][i].res_number < 169) ||
		     (allpdb[0][i].res_number  > 184 && allpdb[0][i].res_number < 353)) {
		     	k = 0;
		     } else {
		     	k = 100;
		     }*/
		
		master_align[i] = k;
		
		

		     
		if (k < (argc-1) * 0.1) {
			for(j=0;j<argc-1;++j) {
				if (j == 0) {
					array[index+0][j] = allpdb[0][i].x_cord;
					array[index+1][j] = allpdb[0][i].y_cord;
					array[index+2][j] = allpdb[0][i].z_cord;
				} else {
					if (align[j-1][i] == -1) {
						
						array[index+0][j] = -9999.9;
						array[index+1][j] = -9999.9;
						array[index+2][j] = -9999.9;
					} else {
						array[index+0][j] = allpdb[j][align[j-1][i]].x_cord;
						array[index+1][j] = allpdb[j][align[j-1][i]].y_cord;
						array[index+2][j] = allpdb[j][align[j-1][i]].z_cord;
					}
				}
			//	printf("ALign:%d:: %d ==  align[%d][%d] ---- array[%d][%d] = %f\n",i,align[j-1][i],j-1,i,index,j,array[index+0][j]);
			}
			
			index += 3;
		}
		
	}
	printf("Build Covariance matrix with index of:%d\n",index);
	
	gsl_matrix *cova = gsl_matrix_alloc(index,index);
	for(i=0;i<index;++i) {
		for (j=0;j<index;++j) {
			// Storer dans covariance matrix au point i,j  la covariance de array[i] qui contient autant d élément que de prot
			gsl_matrix_set(cova,i,j,covar(array[i],array[j],argc-1));
			//printf("I:%d J:%d Covar:%f\n",i,j,covar(array[i],array[j],argc-1));
						
		}
	}
	//write_matrix("covariance.mtx", cova,index,index);
	printf("Diagonalyse Matrix\n");
	
	//Diagonalyse the matrix
	
	gsl_vector *eval_temp = gsl_vector_alloc(index); /*Déclare un vector qui va contenir les eigenvalue */	
	gsl_matrix *evec_temp = gsl_matrix_alloc (index,index); /*Déclare une matrice qui va contenir les eigenvector correspondant à eval*/
	diagonalyse_matrix(cova,index,eval_temp,evec_temp); /*Diagonalyse la matrice*/
	gsl_eigen_symmv_sort (eval_temp, evec_temp,GSL_EIGEN_SORT_ABS_DESC);
	// Rebuild the eigenvector for strc one... will put 0 to Node that were not match
	
	gsl_vector *eval = gsl_vector_alloc(3*nbrCA[0]); /*Déclare un vector qui va contenir les eigenvalue */	
	gsl_matrix *evec = gsl_matrix_alloc (3*nbrCA[0],3*nbrCA[0]);
	gsl_matrix_set_all(evec,0);
	gsl_vector_set_all(eval,0);
	// J veux copien un eigen genre voir
	printf("Rebuild Eigenvector\n");
	
	float tot = 0.0000000000000000;
	
	for (i=0;i<index;++i) {
		tot += gsl_vector_get(eval_temp,i);
	
	}
	
	for (i=0;i<10;++i) {
		printf("I:%d %5.3f \n",i,gsl_vector_get(eval_temp,i)/tot*100);
	
	}
	
	for(j=0;j<index;++j) {
		if (j+6 > nbrCA[0]*3-2) {break;}
		k=0;
		gsl_vector_set(eval,j+6,gsl_vector_get(eval_temp,j));
		for(i=0;i<nbrCA[0];++i) {
			//printf("I:%d Master Align:%d < %f ?\n",i,master_align[i] , (argc-1) * 0.1);
			if (master_align[i] < (argc-1) * 0.1) {
				gsl_matrix_set(evec,i*3+0,j+6,gsl_matrix_get(evec_temp,k*3+0,j));
				gsl_matrix_set(evec,i*3+1,j+6,gsl_matrix_get(evec_temp,k*3+1,j));
				gsl_matrix_set(evec,i*3+2,j+6,gsl_matrix_get(evec_temp,k*3+2,j));
				++k;
			}
		}
	}
	write_eigen("pca_eigen.dat",evec,eval,3*nbrCA[0]);
	return(0);
}
