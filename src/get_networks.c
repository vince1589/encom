#include "STeM.h"
#include <gsl/gsl_linalg.h>

void Extend(int *old, int ne, int ce, int **connex_mat, int *cliques, int c, int max, int tresh, struct pdb_atom *init);

int main(int argc, char *argv[])
{
	int all; /*Number of atoms in the initial PDB*/
	int atom; /*Number of initial CAs*/
	
	int help_flag = 1;
	char file_name[500];
	char eigen_name[500] = "eigen.dat";
	char out_name[500];
	char out_diags[500];
	int verbose = 0;
	double beta = 0.000005;
	int n_tresh = 4;
	float s_tresh = 3.0;
	float h_tresh = 0.33;
	float anti = 1.0;
	
	int i,j,k,l;
	
	int lig = 0;
	int nconn;
	int print_flag = 0;
	int print_diags = 0;
	int hard_tresh = 1;
	for (i = 1;i < argc;i++) {
		if (strcmp("-i",argv[i]) == 0) {strcpy(file_name,argv[i+1]);--help_flag;}
		
		if (strcmp("-h",argv[i]) == 0) {help_flag = 1;}
		if (strcmp("-a",argv[i]) == 0) {anti = -1.0;}
		if (strcmp("-v",argv[i]) == 0) {verbose = 1;}
		if (strcmp("-lig",argv[i]) == 0) {lig= 1;}
		if (strcmp("-ieig",argv[i]) == 0) {strcpy(eigen_name,argv[i+1]);}
		if (strcmp("-o",argv[i]) == 0) {strcpy(out_name,argv[i+1]); print_flag = 1;}
		if (strcmp("-od",argv[i]) == 0) {strcpy(out_diags,argv[i+1]); print_diags = 1;}
		
		if (strcmp("-tr",argv[i]) == 0) {int temp;sscanf(argv[i+1],"%i",&temp); n_tresh = temp;}
		if (strcmp("-str",argv[i]) == 0) {float temp;sscanf(argv[i+1],"%f",&temp); s_tresh = temp; hard_tresh--;}
		if (strcmp("-htr",argv[i]) == 0) {float temp;sscanf(argv[i+1],"%f",&temp); h_tresh = temp;}
		
		if (strcmp("-b",argv[i]) == 0) {float temp;sscanf(argv[i+1],"%f",&temp);beta = temp;}
	}
	
	if (help_flag == 1)
	{
		printf("****************************\nHelp Section\n-i\tFile Input (PDB)\n-v\tVerbose\n-w\tWeight Vector\n-t\tInitial Value of template (negative value for random)\n\tIf Load Template, multiply the template\n-lt\tLoad input template\n-sp\tSuper Node Mode (CA, N, C)\n-kt\tPoid de l'angle entre les nodes (1)\n-kr\tPoid de la distance entre les nodes (1)\n-f\tFile to fit\n****************************\n");
		return(0);
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
	
	// Load hessian
	
	gsl_vector *eval = gsl_vector_alloc(3*atom);
	
	gsl_matrix *evec= gsl_matrix_alloc (3*atom,3*atom);
	
	load_eigen(eval,evec,eigen_name,3*atom);
	
	printf("Check 2\n");
	
	// Invert hessian(s)
	
	gsl_matrix *hess = gsl_matrix_alloc(3*atom,3*atom);
	
	gsl_matrix_set_all(hess, 0);
	
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
	
	printf("Check 3\n");
	
	gsl_matrix *correl = gsl_matrix_alloc(3*atom, 3*atom);
	
	for(i = 0; i < 3*atom; i++)
	{
		for(j = 0; j < 3*atom; j++)
		{
			gsl_matrix_set(correl, i, j, gsl_matrix_get(hess,i,j) / sqrt(gsl_matrix_get(hess,i,i) * gsl_matrix_get(hess,j,j)));
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
	
	float average = 0;
	
	for(i = 0; i < 3*atom; i++)
	{
		for(j = 0; j < 3*atom; j++)
		{
			if(i != j)
			{
				average += gsl_matrix_get(correl, i, j);
			}
		}
	}
	
	average /= 3.0 * (float)atom * (3.0 * (float)atom - 1.0);
	
	float sd = 0;
	
	for(i = 0; i < 3*atom; i++)
	{
		for(j = 0; j < 3*atom; j++)
		{
			if(i != j)
			{
				sd += (gsl_matrix_get(correl, i, j) - average) * (gsl_matrix_get(correl, i, j) - average);
			}
		}
	}
	
	sd /= 3.0 * (float)atom * (3.0 * (float)atom - 1.0);
	sd = sqrt(sd);
	
	printf("Avg : %1.5f\tSd : %1.5f\n", average, sd);
	
	printf("Check 5\n");
	
	if(print_diags == 1)
	{
		FILE *out_file;
		
		out_file = fopen(out_diags, "w");
		
		gsl_matrix *covij = gsl_matrix_alloc(3,3);
		
		for(i = 0; i < atom; i++)
		{
			for(j = i; j < atom; j++)
			{
				gsl_matrix_set_all(covij, 0);
				
				fprintf(out_file, "Node %1i vs node %1i\n\n", i, j);
				
				for(k = 0; k < 3; k++)
				{
					for(l = 0; l < 3; l++)
					{
						gsl_matrix_set(covij, k, l, gsl_matrix_get(correl, 3*i + k, 3*j + l));
						
						fprintf(out_file, "%1.6f\t", gsl_matrix_get(correl, 3*i + k, 3*j + l));
					}
					
					fprintf(out_file, "\n");
				}
				
				fprintf(out_file, "\n");
				
				gsl_vector *eval2 = gsl_vector_alloc(3);
				
				gsl_matrix *evec2= gsl_matrix_alloc (3,3);
				
				diagonalyse_matrix(covij, 3, eval2, evec2);
				
				for(k = 0; k < 3; k++)
				{
					for(l = 0; l < 3; l++)
					{
						fprintf(out_file, "%1.6f\t", gsl_matrix_get(evec2, k, l));
					}
					
					fprintf(out_file, "\n");
				}
				
				fprintf(out_file, "\n");
				
				for(k = 0; k < 3; k++)
				{
					fprintf(out_file, "%1.6f\t", gsl_vector_get(eval2, k));
				}
				
				fprintf(out_file, "\n\n");
				
				gsl_matrix_free(evec2);
				
				gsl_vector_free(eval2);
			}
		}
	}
	
	printf("Check 6\n");
	
	// Build connex matrix
	
	int **connex=(int **)malloc(3*atom*sizeof(int *)); // Matrix of the Hessian 1 2 3 (bond, angle, dihedral)
	for(i=0;i<3*atom;i++) {connex[i]=(int *)malloc(3*atom*sizeof(int));}
	for(i=0;i<3*atom;i++)for(j=0;j<3*atom;j++){connex[i][j]=0;}
	
	for(i = 0; i < atom; i++)
	{
		for(j = 0; j < atom; j++)
		{
			connex[i][j] = 0;
		}
	}
	
	for(i = 0; i < 3*atom; i++)
	{
		for(j = 0; j < 3*atom; j++)
		{
			if((hard_tresh == 0 && anti * (gsl_matrix_get(correl, i, j) - average) / sd >= s_tresh) || (hard_tresh == 1 && anti * gsl_matrix_get(correl, i, j) >= h_tresh) || i == j)
			{
				connex[i / 3][j / 3] ++;
			}
		}
	}
	
	for(i = 0; i < atom; i++)
	{
		for(j = 0; j < atom; j++)
		{
			if(connex[i][j] >= 3)
			{
				connex[i][j] = 1;
			}
			else
			{
				connex[i][j] = 0;
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
	
	Extend(nodenums, 0, atom, connex, clique, cpos, atom, n_tresh, strc_node);
}

void Extend(int *old, int ne, int ce, int **connex_mat, int *cliques, int c, int max, int tresh, struct pdb_atom *init)
{
	int *neww = NULL;
	int nod,fixp;
	int newne,newce,i,j,count,pos,p,s,sel,minnod;
	int loc;
	int l;
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
			if(c >= tresh)
			{
				for(int loca = 0; loca <= c; loca++)
				{
					printf("%1i_%s ", init[cliques[loca]].res_number, init[cliques[loca]].chain);
				}
				
				printf("\n");
			}
		}
		else if(newne < newce)
		{
			Extend(neww, newne, newce, connex_mat, cliques, c, max, tresh, init);
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