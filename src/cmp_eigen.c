#include "STeM.h"

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

float corr_eigen(gsl_matrix *evec,gsl_matrix *pca,int atom,int a,int b) {
	int i;
	double moy_x = 0.00000,moy_y = 0.00000;
 	double std_x = 0.00000,std_y=0.0000000;
 	double prod_tot=0.0000000000;
 	int count = 0;
	for (i=0;i<atom;++i) {
		//printf("I:%d	%f	%f\n",i,gsl_matrix_get(evec,i,a),gsl_matrix_get(pca,i,b));
		//printf("I:%d %8.5f::%8.5f\n",i,gsl_matrix_get(evec,i,b),gsl_matrix_get(pca,i,b));
		if (gsl_matrix_get(pca,i,b) == 0) {continue;}
		++count;
		moy_x += gsl_matrix_get(evec,i,a);
		moy_y += gsl_matrix_get(pca ,i,b);
	}
	moy_x /= count;
	moy_y /= count;
	
	for (i=0;i<atom;++i) {
		if (gsl_matrix_get(pca,i,b) == 0) {continue;}
		std_x += (moy_x - gsl_matrix_get(evec,i,a))*(moy_x - gsl_matrix_get(evec,i,a));
		std_y += (moy_y - gsl_matrix_get(pca ,i,b))*(moy_y - gsl_matrix_get(pca ,i,b));
		prod_tot += (gsl_matrix_get(evec,i,a)-moy_x)*(gsl_matrix_get(pca,i,b)-moy_y);
	}
	
	
	
	return(prod_tot/sqrt(std_y*std_x));

}

float over_eigen(gsl_matrix *evec,gsl_matrix *pca,int atom,int aa,int bb) {
	int i;
	float a = 0.00000,b = 0.00000,c = 0.0000;
 	float x = 0.000000,z=0.00000;
 	float la=0.0000000000,lb=0.000000000;
 	// Caculate lenght (because of the vector reject)
 	
 	for(i=0;i<atom;++i) {
 		if (gsl_matrix_get(pca,i,aa) == 0) {continue;}
 		//printf("I COMPARE: %f et %f\n",gsl_matrix_get(evec,i,aa),gsl_matrix_get(pca,i,bb));
 		la += gsl_matrix_get(evec,i,aa)*gsl_matrix_get(evec,i,aa);
 		lb += gsl_matrix_get(pca ,i,bb)*gsl_matrix_get(pca ,i,bb);
 		
 	
 	}
 	la = sqrt(la);
 	lb = sqrt(lb);
 	//printf("La:%f Lb:%f\n",la,lb);
 	for(i=0;i<atom;++i) {
 		
 		
	 		x = gsl_matrix_get(evec,i,aa);
	 		z = gsl_matrix_get(pca ,i,bb);
	 		if (z == 0) {continue;}
			a += x*z;
	 	
 	}
 	
 	
 	return(a/(la*lb));

}

int main(int argc, char *argv[]) {
	int i,j;
	char eigen_name[100] = "eigen.dat";
	char pca_name[100] = "pca_eigen.dat";
	int atom[2];
	
	for (i = 1;i < argc;i++) {
 		if (strcmp("-ieig",argv[i]) == 0) {strcpy(eigen_name,argv[i+1]);}
		if (strcmp("-pca",argv[i]) == 0)  {strcpy(pca_name,argv[i+1]);}
	}
	printf("My Eigen:%s    My Pca:%s\n",eigen_name,pca_name);
	
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
	float corr;
	printf("J:NMA I:PCA\n");
	//float overlap;
	for(i=6;i<17;++i) {
		for(j=6;j<17;++j) {
			corr = corr_eigen(evec,evecpca,atom[0],j,i);
			//overlap = over_eigen(evec,evecpca,atom[0],j,i);
			if (corr*corr > 0.1) {printf("J:%3d I:%3d %8.5f :: %8.5f\n",j,i,corr*corr,corr_eigen(evec,evecpca,atom[0],j,i));}
			
	
		}
	}
	return(1);
}
