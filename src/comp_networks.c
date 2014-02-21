#include "STeM.h"

void Extend(int *old, int ne, int ce, int **connex_mat, int *cliques, int c, int max, int tresh, struct pdb_atom *init, struct pdb_atom *init_t, int *linker, int *aln);

int main(int argc, char *argv[])
{
	int all; /*Number of atoms in the initial PDB*/
	int atom; /*Number of initial CAs*/
	int all_t; /*Number of atoms in the target PDB*/
	int atom_t; /*Number of target CAs*/
	
	int help_flag = 1;
	int has_target = 0;
	char file_name[500];
	char file_eigen[500];
	char check_name[500];
	char check_eigen[500];
	int verbose = 0;
	int n_tresh = 4;
	float s_tresh = 3.0;
	float h_tresh = 0.1;
	int hard_tresh = 0;
	int lig = 0;
	int lig_t = 0;
	
	int i,j,k,l;
	
	int nconn;
	int print_flag = 0;
	float ligalign = 0; // Flag/valeur pour aligner seulement les résidus dans un cutoff du ligand, 0, one le fait pas... > 0... le cutoff
	for (i = 1;i < argc;i++) {
		if (strcmp("-i",argv[i]) == 0) {strcpy(file_name,argv[i+1]);--help_flag;}
		if (strcmp("-ieig",argv[i]) == 0) {strcpy(file_eigen,argv[i+1]);}
		
		if (strcmp("-h",argv[i]) == 0) {help_flag = 1;}
		if (strcmp("-v",argv[i]) == 0) {verbose = 1;}

		if (strcmp("-t",argv[i]) == 0) {strcpy(check_name,argv[i+1]);help_flag = 0;has_target = 1;}
		if (strcmp("-teig",argv[i]) == 0) {strcpy(check_eigen,argv[i+1]);}
		
		if (strcmp("-lig",argv[i]) == 0) {lig = 1;}
		if (strcmp("-ligt",argv[i]) == 0) {lig_t = 1;}
		
		if (strcmp("-ligc",argv[i]) == 0) {float temp;sscanf(argv[i+1],"%f",&temp);ligalign = temp;}
		
		if (strcmp("-tr",argv[i]) == 0) {int temp;sscanf(argv[i+1],"%i",&temp);n_tresh = temp;}
		if (strcmp("-str",argv[i]) == 0) {float temp;sscanf(argv[i+1],"%f",&temp);s_tresh = temp;}
		if (strcmp("-htr",argv[i]) == 0) {float temp;sscanf(argv[i+1],"%f",&temp);s_tresh = temp; hard_tresh = 0;}
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
	
	//***************************************************
	//*													*
	//*Builds a structure contaning information on the target pdb structure
	//*													*
	//***************************************************
	
	nconn = 0;
	
	all_t = count_atom(check_name);
	
	nconn = count_connect(check_name);
	
	if (verbose == 1) {printf("Connect:%d\n",nconn);}
	
	if (verbose == 1) {printf("Assigning Structure\n\tAll Atom\n");}
	
	// Array with all connects
	
	int **connect_t=(int **)malloc(nconn*sizeof(int *));
	
	for(i=0;i<nconn;i++) { connect_t[i]=(int *)malloc(6*sizeof(int));}
	
	assign_connect(check_name,connect_t);
	
	// Assigns all the atoms
	
	struct pdb_atom strc_all_t[all_t];
	
	atom_t = build_all_strc(check_name,strc_all_t); // Retourne le nombre de Node
	
	if (atom_t > 800) {printf("Too much node.... To fix, ask vincent.frappier@usherbrooke.ca\n");return(1);}
	
	if (verbose == 1) {printf("	Atom:%d\n",all_t);}
	
	check_lig(strc_all_t,connect_t,nconn,all_t);
	
	// Assigns all Nodes
	
	if (verbose == 1) {printf("	CA Structure\n");}
	
	if (verbose == 1) {printf("	Node:%d\n",atom_t);}
	
	struct pdb_atom strc_node_t[atom_t];
	
	atom_t = build_cord_CA(strc_all_t, strc_node_t,all_t,lig_t,connect_t,nconn);
	
	if (verbose == 1) {printf("	Assign Node:%d\n",atom_t);}
	
	//***************************************************
	//*													*
	//*Aligns both structures
	//*													*
	//***************************************************
	
	int align[atom];
	
	int score = node_align(strc_node,atom,strc_node_t,atom_t,align);
	
	printf("RMSD:%8.5f Score: %d/%d\n",sqrt(rmsd_no(strc_node,strc_node_t,atom, align)),score,atom);
	
	if ((float)score/(float)atom < 0.8)
	{
		printf("Low Score... Will try an homemade alignement !!!\n");
		
		score = node_align_low(strc_node,atom,strc_node_t,atom_t,align);
		
		printf("RMSD:%8.5f Score: %d/%d\n",sqrt(rmsd_no(strc_node,strc_node_t,atom, align)),score,atom);
	}
	
	if (ligalign > 0)
	{
		score = node_align_lig(strc_node,atom,strc_node_t,atom_t,align,strc_all,all,strc_all_t,all_t,ligalign);
		
		printf("RMSD:%8.5f Score: %d/%d\n",sqrt(rmsd_no(strc_node,strc_node_t,atom, align)),score,atom);
	}
	
	printf("RMSD:%8.5f Score: %d/%d\n",sqrt(rmsd_yes(strc_node,strc_node_t,atom, align,strc_all,all)),score,atom);
	
	
	printf("Check 1\n");
	
	// Build hessians
	
	gsl_vector *eval = gsl_vector_alloc(3*atom);
	
	gsl_matrix *evec = gsl_matrix_alloc (3*atom,3*atom);
	
	load_eigen(eval,evec,file_eigen,3*atom);
	
	gsl_vector *eval_t = gsl_vector_alloc(3*atom_t);
	
	gsl_matrix *evec_t = gsl_matrix_alloc (3*atom_t,3*atom_t);
	
	load_eigen(eval_t,evec_t,check_eigen,3*atom_t);
	
	printf("Check 1.1\n");
	
	// Invert hessians
	
	gsl_matrix *ihess = gsl_matrix_alloc(3*atom,3*atom);
	
	gsl_matrix_set_all(ihess, 0);
	
	gsl_matrix *ihess_t = gsl_matrix_alloc(3*atom_t,3*atom_t);
	
	gsl_matrix_set_all(ihess_t, 0);
	
	for(i = 0; i < 3*atom; i++)
	{
		if(gsl_vector_get(eval, i) > 0.00001)
		{
			for(j = 0; j < 3*atom; j++)
			{
				for(k = 0; k < 3*atom; k++)
				{
					gsl_matrix_set(ihess, j, k, gsl_matrix_get(ihess, j, k) + gsl_matrix_get(evec, j, i)*gsl_matrix_get(evec, k, i)/gsl_vector_get(eval, i));
				}
			}
		}
	}
	
	for(i = 0; i < 3*atom_t; i++)
	{
		if(gsl_vector_get(eval_t, i) > 0.00001)
		{
			for(j = 0; j < 3*atom_t; j++)
			{
				for(k = 0; k < 3*atom_t; k++)
				{
					gsl_matrix_set(ihess_t, j, k, gsl_matrix_get(ihess_t, j, k) + gsl_matrix_get(evec_t, j, i)*gsl_matrix_get(evec_t, k, i)/gsl_vector_get(eval_t, i));
				}
			}
		}
	}
	
	printf("Check 2\n");
	
	// Build mini-hessians
	
	gsl_matrix *mini_hess_i = gsl_matrix_alloc(3*score, 3*score);
	gsl_matrix_set_all(mini_hess_i, 10);
	gsl_matrix *mini_hess_t = gsl_matrix_alloc(3*score, 3*score);
	gsl_matrix_set_all(mini_hess_t, 10);
	
	int sup_line = 0;
	
	int sup_to_node[score];
	
	for(i = 0; i < atom; i++)
	{
		if(align[i] == -1) {continue;}
		
		sup_to_node[sup_line] = i;
		
		int sup_col = 0;
		
		for(j = 0; j < atom; j++)
		{
			if(align[j] == -1) {continue;}
			
			for(k = 0; k < 3; k++)
			{
				for(l = 0; l < 3; l++)
				{
					gsl_matrix_set(mini_hess_i, 3*sup_line + k, 3*sup_col + l, gsl_matrix_get(ihess, 3*i + k, 3*j + l));
					gsl_matrix_set(mini_hess_t, 3*sup_line + k, 3*sup_col + l, gsl_matrix_get(ihess_t, 3*align[i] + k, 3*align[j] + l));
				}
			}
			
			sup_col ++;
		}
		
		sup_line++;
	}
	
	gsl_matrix_free(ihess);
	gsl_matrix_free(ihess_t);
	
	printf("Check 3\n");
	
	//Build comp_matrix
	
	gsl_matrix *comp_matrix = gsl_matrix_alloc(3*score, 3*score);
	
	for(i = 0; i < 3*score; i++)
	{
		for(j = 0; j < 3*score; j++)
		{
			gsl_matrix_set(comp_matrix, i, j, gsl_matrix_get(mini_hess_t, i, j) - gsl_matrix_get(mini_hess_i, i, j));
		}
	}
	
	gsl_matrix_free(mini_hess_i);
	gsl_matrix_free(mini_hess_t);
	
	// Calculate average and sd
	
	float average = 0;
	
	for(i = 0; i < 3*score; i++)
	{
		for(j = 0; j < 3*score; j++)
		{
			if(i != j)
			{
				average += gsl_matrix_get(comp_matrix, i, j);
			}
		}
	}
	
	average /= 3.0 * (float)atom * (3.0 * (float)atom - 1.0);
	
	float sd = 0;
	
	for(i = 0; i < 3*score; i++)
	{
		for(j = 0; j < 3*score; j++)
		{
			if(i != j)
			{
				sd += (gsl_matrix_get(comp_matrix, i, j) - average) * (gsl_matrix_get(comp_matrix, i, j) - average);
			}
		}
	}
	
	sd /= 3.0 * (float)atom * (3.0 * (float)atom - 1.0);
	sd = sqrt(sd);
	
	printf("Avg : %1.5f\tSd : %1.5f\n", average, sd);
	
	printf("Check 4\n");
	
	// Build connex matrix
	
	int **connex=(int **)malloc(score*sizeof(int *)); // Matrix of the Hessian 1 2 3 (bond, angle, dihedral)
	for(i=0;i<score;i++) {connex[i]=(int *)malloc(score*sizeof(int));}
	for(i=0;i<score;i++)for(j=0;j<score;j++){connex[i][j]=0;}
	
	for(i = 0; i < score; i++)
	{
		for(j = 0; j < score; j++)
		{
			connex[i][j] = 0;
		}
	}
	
	if(hard_tresh == 0)
	{
		for(i = 0; i < 3*score; i++)
		{
			for(j = 0; j < 3*score; j++)
			{
				if((gsl_matrix_get(comp_matrix, i, j) - average) / sd >= s_tresh || (gsl_matrix_get(comp_matrix, i, j) - average) / sd <= -1.0 * s_tresh || i == j)
				{
					connex[i / 3][j / 3] = 1;
				}
			}
		}
	}
	else
	{
		for(i = 0; i < 3*score; i++)
		{
			for(j = 0; j < 3*score; j++)
			{
				if(gsl_matrix_get(comp_matrix, i, j) >= h_tresh || gsl_matrix_get(comp_matrix, i, j) <= -1.0 * h_tresh || i == j)
				{
					connex[i / 3][j / 3] = 1;
				}
			}
		}
	}
	
	int clique[score];
	
	int cpos = 0;
	
	int nodenums[score];
	
	for(i = 0; i < score; i++)
	{
		nodenums[i] = i;
	}
	
	Extend(nodenums, 0, score, connex, clique, cpos, score, n_tresh, strc_node, strc_node_t, sup_to_node, align);
}

void Extend(int *old, int ne, int ce, int **connex_mat, int *cliques, int c, int max, int tresh, struct pdb_atom *init, struct pdb_atom *init_t, int *linker, int *aln)
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
			if(c >= 2)
			{
				for(int loca = 0; loca <= c; loca++)
				{
					printf("(%1i,%1i) ", init[linker[cliques[loca]]].res_number, init_t[aln[linker[cliques[loca]]]].res_number);
				}
				
				printf("\n");
			}
		}
		else if(newne < newce)
		{
			Extend(neww, newne, newce, connex_mat, cliques, c, max, tresh, init, init_t, linker, aln);
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