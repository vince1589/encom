#include "STeM.h"
#include <gsl/gsl_linalg.h>

void Extend(int *old, int ne, int ce, int **connex_mat, int *cliques, int c, int max, int tresh, int **clique_mat, int mat_size, int *cliquenum);

int clique_fusion(int **clique_mat, int mat_size, int *cliquenum);

int get_connex(int nod1, int nod2, struct pdb_atom *init, int natom, float max);

int main(int argc, char *argv[])
{
	int all; /*Number of atoms in the initial PDB*/
	int atom; /*Number of initial CAs*/
	
	int help_flag = 1;
	char file_name[500];
	char eigen_name[500] = "eigen.dat";
	char correl_name[500] = "correl.dat";
	char out_name[500];
	char out_correl[500];
	int verbose = 0;
	double beta = 0.000005;
	int n_tresh = 4;
	float s_tresh = 1.645;
	float h_tresh = 0.33;
	float anti = 1.0;
	float cutoff = 5.0;
	int bresn[800];
	char bresc[800];
	int nbind = 0;
	int aresn[800];
	int aresc[800];
	int nallo = 0;
	
	int test = 1;
	
	int i,j,k,l;
	
	int lig = 0;
	int nconn;
	int print_flag = 0;
	int print_correl = 0;
	int hard_tresh = 0;
	int mat_flag = 0;
	int eig_flag = 0;
	
	for (i = 1; i < argc; i++)
	{
		if (strcmp("-i",argv[i]) == 0) {strcpy(file_name,argv[i+1]);--help_flag;}
		
		if (strcmp("-h",argv[i]) == 0) {help_flag = 1;}
		if (strcmp("-a",argv[i]) == 0) {anti = -1.0;}
		if (strcmp("-v",argv[i]) == 0) {verbose = 1;}
		if (strcmp("-lig",argv[i]) == 0) {lig= 1;}
		if (strcmp("-ieig",argv[i]) == 0) {strcpy(eigen_name,argv[i+1]); eig_flag = 1;}
		if (strcmp("-imat",argv[i]) == 0) {strcpy(correl_name,argv[i+1]); mat_flag = 1;}
		if (strcmp("-o",argv[i]) == 0) {strcpy(out_name,argv[i+1]); print_flag = 1;}
		if (strcmp("-om",argv[i]) == 0) {strcpy(out_correl,argv[i+1]); print_correl = 1;}
		if (strcmp("-cut",argv[i]) == 0) {float temp;sscanf(argv[i+1],"%f",&temp); cutoff = temp;}
		
		if(strcmp("-bres", argv[i]) == 0)
		{
			j = 0;
			
			int tmp1;
			char tmp2;
			
			while(sscanf(argv[i+j+1], "%i_%c", &tmp1, &tmp2) != EOF)
			{
				bresn[j] = tmp1;
				
				bresc[j] = tmp2;
				
				j++;
				
				nbind++;
			}
		}
		
		if(strcmp("-ares", argv[i]) == 0)
		{
			j = 0;
			
			int tmp1;
			char tmp2;
			
			while(sscanf(argv[i+j+1], "%i_%c", &tmp1, &tmp2) != EOF)
			{
				aresn[j] = tmp1;
				
				aresc[j] = tmp2;
				
				j++;
				
				nallo++;
			}
		}
		
		if (strcmp("-tr",argv[i]) == 0) {int temp;sscanf(argv[i+1],"%i",&temp); n_tresh = temp;}
		if (strcmp("-str",argv[i]) == 0) {float temp;sscanf(argv[i+1],"%f",&temp); s_tresh = temp;}
		if (strcmp("-htr",argv[i]) == 0) {float temp;sscanf(argv[i+1],"%f",&temp); h_tresh = temp; hard_tresh++;}
		
		if (strcmp("-b",argv[i]) == 0) {float temp;sscanf(argv[i+1],"%f",&temp);beta = temp;}
	}
	
	if (help_flag == 1)
	{
		printf("****************************\nHelp Section\n-i\tFile Input (PDB)\n-v\tVerbose\n-w\tWeight Vector\n-t\tInitial Value of template (negative value for random)\n\tIf Load Template, multiply the template\n-lt\tLoad input template\n-sp\tSuper Node Mode (CA, N, C)\n-kt\tPoid de l'angle entre les nodes (1)\n-kr\tPoid de la distance entre les nodes (1)\n-f\tFile to fit\n****************************\n");
		return(0);
	}
	
	if(test == 1)
	{
		printf("Bind :\n");
		
		for(i == 0; i < nbind; i++)
		{
			printf("%1i_%1c\n", bresn[i], bresc[i]);
		}
		
		printf("Allo :\n");
		
		for(i == 0; i < nallo; i++)
		{
			printf("%1i_%1c\n", aresn[i], aresc[i]);
		}
		
		return 0;
	}
	
	//***************************************************
	//*													*
	//*Builds a structure contaning information on the initial pdb structure
	//*													*
	//***************************************************
	
	all = count_atom(file_name);
	
	nconn = count_connect(file_name);
	
	if (verbose == 1) {printf("Connect:%d\n",nconn);}
	
	if (verbose == 1) {printf("Assigning Structure\n\tAll Atom\n");}
	
	// Array with all connects
	
	int **connect_h=(int **)malloc(nconn*sizeof(int *)); 
	
	for(i=0;i<nconn;i++) { connect_h[i]=(int *)malloc(6*sizeof(int));}
	
	assign_connect(file_name,connect_h);
	
	// Assigns all the atoms
	
	struct pdb_atom strc_all[all];
	
	atom = build_all_strc(file_name,strc_all); // Retourne le nombre de Node
	
	if (atom > 800) {printf("Too much nodes .... To fix, ask vincent.frappier@usherbrooke.ca\n");return(1);}
	
	if (verbose == 1) {printf("	Atom:%d\n",all);}
	
	check_lig(strc_all,connect_h,nconn,all);
	
	// Assigns all Nodes
	
	if (verbose == 1) {printf("	CA Structure\n");}
	
	if (verbose == 1) {printf("	Node:%d\n",atom);}
	
	struct pdb_atom strc_node[atom];
	
	atom = build_cord_CA(strc_all, strc_node,all,lig,connect_h,nconn);
	
	if (verbose == 1) {printf("	Assign Node:%d\n",atom);}
	
	// Free Connect
	
	//for(i=0;i<nconn;i++) {printf("I:%d\n",i);free(connect_h[i]);}
	//free(connect_h);
	
	printf("Check 1\n");
	
	gsl_matrix *hess = gsl_matrix_alloc(3*atom,3*atom);
	
	gsl_matrix_set_all(hess, 0);
	
	if(eig_flag == 1)
	{
		// Load hessian
	
		gsl_vector *eval = gsl_vector_alloc(3*atom);
	
		gsl_matrix *evec= gsl_matrix_alloc (3*atom,3*atom);
	
		load_eigen(eval,evec,eigen_name,3*atom);
	
		printf("Check 2\n");
	
		// Invert hessian(s)
	
		for(i = 0; i < 3*atom; i++)
		{
			if(gsl_vector_get(eval, i) > 0.00001)
			{
				for(j = 0; j < 3*atom; j++)
				{
					for(k = 0; k < 3*atom; k++)
					{
						gsl_matrix_set(hess, j, k, gsl_matrix_get(hess, j, k) + gsl_matrix_get(evec, j, i)*gsl_matrix_get(evec, k, i)/gsl_vector_get(eval, i));
					}
				}
			}
		}
		
		gsl_matrix_free(evec);
		gsl_vector_free(eval);
	}
	
	printf("Check 3\n");
	
	gsl_matrix *correl = gsl_matrix_alloc(3*atom, 3*atom);
	
	if(mat_flag == 1)
	{
		printf("Dimension : %1i\n", 3*atom);
		
		printf("Dimension : %1i\n", load_matrix(correl, correl_name));
	}
	
	if(eig_flag == 1)
	{
		for(i = 0; i < 3*atom; i++)
		{
			for(j = 0; j < 3*atom; j++)
			{
				gsl_matrix_set(correl, i, j, gsl_matrix_get(hess,i,j) / sqrt(gsl_matrix_get(hess,i,i) * gsl_matrix_get(hess,j,j)));
			}
		}
	}
	
	gsl_matrix_free(hess);
	
	if(print_flag == 1)
	{
		printf("Writing matrix...\n");
		
		write_matrix(out_name,correl,3*(atom-lig),3*(atom-lig));
	}
	
	printf("Check 4\n");
	
	// Calculate average and sd
	
	float average = 0.0;
	
	for(i = 1; i < 3*atom; i++)
	{
		for(j = 0; j < i; j++)
		{
				average += gsl_matrix_get(correl, i, j);
		}
	}
	
	average /= (3.0 * (float)atom - 1.0) * (3.0 * (float)atom - 2.0) / 2.0;
	
	float sd = 0.0;
	
	for(i = 1; i < 3*atom; i++)
	{
		for(j = 0; j < i; j++)
		{
				sd += (gsl_matrix_get(correl, i, j) - average) * (gsl_matrix_get(correl, i, j) - average);
		}
	}
	
	sd /= (3.0 * (float)atom - 1.0) * (3.0 * (float)atom - 2.0) / 2.0 - 1.0;
	
	float flevel = 1.0 / sd;
	
	sd = sqrt(sd);
	
	printf("Avg : %1.5f\nSd : %1.5f\nFreedom level : %1.5f\nEstimated molar entropy : %1.5f J/(mol.K)\n", average, sd, flevel, 8.314 * log(flevel));
	
	printf("Check 5\n");
	
	// Print correlation matrix
	
	if(print_correl == 1)
	{
		FILE *out_file;
		
		out_file = fopen(out_correl, "w");
		
		for(i = 0; i < 3*atom; i++)
		{
			
			fprintf(out_file, "\t%1.6f", gsl_matrix_get(correl, i, 0));
			
			for(j = 1; j < 3*atom; j++)
			{
				fprintf(out_file, "\t%1.6f", gsl_matrix_get(correl, i, j));
			}
			
			fprintf(out_file, "\n");
		}
	}
	
	printf("Check 6\n");
	
	// Build connex matrix
	
	int **connex=(int **)malloc(atom*sizeof(int *));
	for(i=0;i<atom;i++) {connex[i]=(int *)malloc(atom*sizeof(int));}
	for(i=0;i<atom;i++)for(j=0;j<atom;j++){connex[i][j]=0;}
	
	for(i = 0; i < 3*atom; i++)
	{
		for(j = 0; j <= i; j++)
		{
			if((hard_tresh == 0 && anti * (gsl_matrix_get(correl, i, j) - average) / sd >= s_tresh) || (hard_tresh == 1 && anti * gsl_matrix_get(correl, i, j) >= h_tresh) || i == j)
			{
				connex[i / 3][j / 3] ++;
				connex[j / 3][i / 3] ++;
			}
		}
	}
	
	for(i = 0; i < atom; i++)
	{
		for(j = 0; j <= i; j++)
		{
			if(connex[i][j] >= 3 && get_connex(i, j, strc_all, all, cutoff) == 1)
			{
				connex[i][j] = 1;
				connex[j][i] = 1;
			}
			else
			{
				connex[i][j] = 0;
				connex[j][i] = 0;
			}
		}
	}
	
	int clique[atom];
	
	int cpos = -1;
	
	int nodenums[atom];
	
	for(i = 0; i < atom; i++)
	{
		nodenums[i] = i;
	}
	
	// Build cliques matrix
	
	int **cliques=(int **)malloc(2*atom*sizeof(int *));
	for(i=0;i<2*atom;i++) {cliques[i]=(int *)malloc((atom + 1)*sizeof(int));}
	for(i=0;i<2*atom;i++)for(j=0;j<atom;j++){cliques[i][j]=-1;}
	
	for(i=0;i<atom;i++) {cliques[i][atom] = 0;}
	
	// Fill cliques matrix
	
	int cliquen = 0;
	
	Extend(nodenums, 0, atom, connex, clique, cpos, atom, n_tresh, cliques, atom, &cliquen);
	
	printf("Check 7\n");
	
	// Build networks from cliques
	
	int cliquen_old = cliquen;
	
	while(clique_fusion(cliques, atom, &cliquen) < cliquen_old) {cliquen_old = cliquen;}
	
	for(i = 0; i < cliquen; i++)
	{
		for(j = 0; j < cliques[i][atom]; j++)
		{
			printf("%1i_%1s\t", strc_node[cliques[i][j]].res_number, strc_node[cliques[i][j]].chain);
		}
		
		printf("\n\n");
	}
	
	return 1;
}

void Extend(int *old, int ne, int ce, int **connex_mat, int *cliques, int c, int max, int tresh, int **clique_mat, int mat_size, int *cliquenum)
{
	int *neww = NULL;
	int nod,fixp;
	int newne,newce,i,j,count,pos,p,s,sel,minnod;
	int loc;
	int l;
	int k;
	//printf("entered Extend\n");PAUSE;
	
	neww = (int *)malloc(ce*sizeof(int *));
	
	if(neww == NULL)
	{
		printf("unnable to allocate memory to new\n");
		exit(1);
	}
	
	minnod=ce;
	
	i=0;
	
	nod=0;
	
	while(i<ce && minnod != 0)
	{
		p=old[i];
		
		count=0;
		
		j=ne;
		
		while(j < ce && count < minnod)
		{    
			if(connex_mat[p][old[j]] == 0)
			{
				count++;
				
				pos=j;
			}
			
			j++;
		}
		
		if(count < minnod)
		{
			fixp=p;
			
			minnod=count;
			
			if(i < ne)
			{
				s=pos;
			}
			else
			{
				s=i;
				
				nod=1;
			}
		}
		i++;
	}
	
	for(nod=minnod+nod;nod>0;nod--)
	{
		p=old[s];
		
		old[s]=old[ne];
		
		old[ne]=p;
		
		sel=p;
		
		newne=0;
		
		i=0;
		
		while(i<ne)
		{
			if(connex_mat[sel][old[i]]==1) neww[newne++]=old[i];
			i++;
		}
		
		newce=newne;
		
		i=ne+1;
		
		while(i < ce)
		{
			if(connex_mat[sel][old[i]] == 1) neww[newce++]=old[i];
			i++;
		}
		
		cliques[++c] = sel;
		
		if(newce == 0)
		{
			if(c >= tresh && *cliquenum < 2*mat_size)
			{
				for(int loca = 0; loca <= c; loca++)
				{
					clique_mat[*cliquenum][loca] = cliques[loca];
					
					//printf("%1i\t", cliques[loca]);
				}
				
				clique_mat[*cliquenum][mat_size] = c + 1;
							
				//printf("\n", *cliquenum, clique_mat[*cliquenum][mat_size]);
				
				(*cliquenum) ++;
			}
		}
		else if(newne < newce)
		{
			Extend(neww, newne, newce, connex_mat, cliques, c, max, tresh, clique_mat, mat_size, cliquenum);
		}
		
		c--;
		
		ne++;
		
		if(nod > 1)
		{
			s=ne;
			
			while(connex_mat[fixp][old[s]] == 1 && s < max) s++;
		}    
	}
	
	free(neww);
	
	neww=NULL;
	
	return;
}

int clique_fusion(int **clique_mat, int mat_size, int *cliquenum)
{
	int i,j,k,l;	
	
	int match_found = 0;
	
	for(i = 0; i < *cliquenum; i++)
	{
		if(clique_mat[i][mat_size] == 0) {continue;}
		
		//printf("Size of %1i : %1i\n", i, clique_mat[i][mat_size]);
		
		for(j = i + 1; j < *cliquenum; j++)
		{
			if(clique_mat[j][mat_size] == 0) {continue;}
			
			//printf("Size of %1i : %1i\n", j, clique_mat[j][mat_size]);
			
			match_found = 0;
			
			for(k = 0; k < clique_mat[j][mat_size]; k++)
			{
				for(l = 0; l < clique_mat[i][mat_size]; l++)
				{
					if(clique_mat[i][l] == clique_mat[j][k])
					{
						match_found = 1;
						
						break;
					}
				}
				
				if(match_found == 1)
				{
					break;
				}
			}
			
			if(match_found == 1)
			{
				int matches = 0;
				
				for(k = 0; k < clique_mat[j][mat_size]; k++)
				{
					match_found = 0;
					
					for(l = 0; l < clique_mat[i][mat_size] - matches; l++)
					{
						if(clique_mat[i][l] == clique_mat[j][k])
						{
							match_found = 1;
							
							break;
						}
					}
					
					if(match_found == 0)
					{
						clique_mat[i][clique_mat[i][mat_size]] = clique_mat[j][k];
						
						clique_mat[i][mat_size] ++;
						
						matches ++;
					}
					
					clique_mat[j][k] = -1;
				}
				
				clique_mat[j][mat_size] = 0;
				
				//printf("Size of %1i and %1i : %1i\t%1i\nMat_size : %1i\n", i, j, clique_mat[i][mat_size], clique_mat[j][mat_size], 2*mat_size - 1);
			}
		}
		
		
	}
	
	for(i = 0; i < *cliquenum; i++)
	{
		if(clique_mat[i][mat_size] == 0)
		{
			j = i+1;
			
			while(clique_mat[j][mat_size] == 0 && j < *cliquenum - 1)
			{
				j++;
			}
			
			if(clique_mat[j][mat_size] == 0)
			{
				*cliquenum = i;
				
				return i;
			}
			
			//printf("Return %1i to %1i\n", j, i);
			
			for(k = 0; k < clique_mat[j][mat_size]; k++)
			{
				clique_mat[i][k] = clique_mat[j][k];
			}
			
			clique_mat[i][mat_size] = clique_mat[j][mat_size];
			
			clique_mat[j][mat_size] = 0;
		}
	}
}

int get_connex(int nod1, int nod2, struct pdb_atom *init, int natom, float max)
{
	int i,j;
	
	if(nod1 == nod2)
	{
		return 1;
	}
	
	for(i = 0; i < natom; i++)
	{
		if(init[i].node == nod1)
		{
			for(j = 0; j < natom; j++)
			{
				if(init[j].node == nod2)
				{
					if((init[i].x_cord - init[j].x_cord)*(init[i].x_cord - init[j].x_cord)
						+ (init[i].y_cord - init[j].y_cord)*(init[i].y_cord - init[j].y_cord)
						+ (init[i].z_cord - init[j].z_cord)*(init[i].z_cord - init[j].z_cord) <= max*max)
					{
						return 1;
					}
				}
			}
		}
	}
	
	return 0;
}
