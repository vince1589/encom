#include "STeM.h"
#include <gsl/gsl_linalg.h>

void Extend(int *old, int ne, int ce, int **connex_mat, int *cliques, int c, int max, int tresh, int **clique_mat, int mat_size, int *cliquenum);
int clique_fusion(int **clique_mat, int mat_size, int *cliquenum);
int get_connex(int nod1, int nod2, struct pdb_atom *init, int natom, float max);
int md_print(struct pdb_atom *init, int size_all, int **conx, int size_node, char *out);
float lc(int node, int **conx_mat, int *network, int size);

int main(int argc, char *argv[])
{
	int i,j,k,l;
	int all; /* Number of atoms in the initial PDB */
	int atom; /* Number of initial CAs */
	int nconn;
	
	char file_name[500];
	
	char eigen_name[500] = "eigen.dat";
	char correl_name[500] = "correl.dat";
	
	char out_cd[500] = "cor_vs_dist.dat"; // Correlation vs distance out_name
	char out_correl[500] = "correl.dat"; // Correlation matrix out_name
	char out_lc[500] = "lc.dat"; // Local closeness out_name
	char out_net[500] = "networks.dat"; // Networks out_name
	
	/*
	The whole thing should be tested with selection of both correlated and anti-correlated residues in the same graph.
	Tresholds would then have to be increased (1.645 will be 10% most correlated or anti-correlated instead of 5% most correlated).
	*/
	
	int n_tresh = 4; // Gives the clique size treshold. SHOULD be tested with n = 3 and up. Has not been tested yet (except on 4).
	int z_tresh = 3; // Gives the number of z-matches required for connexion. SHOULD be tested with n = 1..9. Should also be tested with the trace, just in case (but won't be so physicaly relevant).
	float s_tresh = 1.645; // Gives the Z treshold (optimal to this date : 2.241 or 1.125% most correlated)
	float h_tresh = 0.33; // Gives a hard treshold (based on r instead of Z, not in use)
	float anti = 1.0; // If changed to -1.0, tells program to search for anitcorrelated nodes instead
	float cutoff = 5.0; // Gives a distance cutoff to the connex completion step (current optimum is 6.0, but should be tested from 5.0 to 6.5 with 0.1 increments)
	
	int print_correl = 0; // Tells the program to print the correlation matrix
	int print_cd = 0; // Tells the program to print correlation vs distance of residues
	int print_lc = 0; // Tells the program to print local closeness
	int print_net = 0; // Tells the program to print networks
	
	int mat_flag = 0; // Tells program to load correlation matrix
	int eig_flag = 0; // Tells program to load hessian
	
	int help_flag = 1;
	int verbose = 0;
	int lig = 0;
	int hard_tresh = 0; // Tells the program to use the hard treshold
	int cde_flag = 0; // Tells the program to end after printing correlation vs distance of residues
	int cor_flag = 0; // Tells the program to end after calculating the correlation matrix
	int full_flag = 0; // Tells the program to use full networks instead of clique networks
	
	int *bresn; // Binding site residue number array (to be freed)
	bresn = (int *) malloc(10000*sizeof(int));
	
	// Binding site residue chain array (to be freed)
	
	char **bresc;
	bresc = (char **) malloc(10000*sizeof(char *));
	
	if(bresc == NULL) {printf("malloc FAIL\n"); return 0;}
	
	for(i = 0; i < 10000; i++)
	{
		bresc[i] = (char *) malloc(6*sizeof(char));
		
		if(bresc[i] == NULL) {printf("malloc FAIL\n"); return 0;}	
	}
	
	// End of bresc declaration
	
	int nbind = 0;
	
	int *aresn; // Allosteric site residue number array (to be freed)
	aresn = (int *) malloc(10000*sizeof(int));
	
	// Allosteric site residue chain array (to be freed)
	
	char **aresc;
	aresc = (char **) malloc(10000*sizeof(char *));
	
	if(aresc == NULL) {printf("malloc FAIL\n"); return 0;}
	
	for(i = 0; i < 10000; i++)
	{
		aresc[i] = (char *) malloc(6*sizeof(char));
		
		if(aresc[i] == NULL) {printf("malloc FAIL\n"); return 0;}	
	}
	
	// End of aresc declaration
	
	int nallo = 0;
	
	for (i = 1; i < argc; i++)
	{
		if (strcmp("-i",argv[i]) == 0) {strcpy(file_name,argv[i+1]);--help_flag;}
		
		if (strcmp("-h",argv[i]) == 0) {help_flag = 1;}
		if (strcmp("-a",argv[i]) == 0) {anti = -1.0;}
		if (strcmp("-v",argv[i]) == 0) {verbose = 1;}
		if (strcmp("-lig",argv[i]) == 0) {lig= 1;}
		if (strcmp("-cde",argv[i]) == 0) {cde_flag = 1;}
		if (strcmp("-cor",argv[i]) == 0) {cor_flag = 1;}
		if (strcmp("-full",argv[i]) == 0) {full_flag = 1;}
		
		if (strcmp("-ieig",argv[i]) == 0) {strcpy(eigen_name,argv[i+1]); eig_flag = 1;}
		if (strcmp("-imat",argv[i]) == 0) {strcpy(correl_name,argv[i+1]); mat_flag = 1;}
		
		if (strcmp("-om",argv[i]) == 0) {strcpy(out_correl,argv[i+1]); print_correl = 1;}
		if (strcmp("-olc",argv[i]) == 0) {strcpy(out_lc,argv[i+1]); print_lc = 1;}
		if (strcmp("-ont",argv[i]) == 0) {strcpy(out_net,argv[i+1]); print_net = 1;}
		if (strcmp("-cvd",argv[i]) == 0) {strcpy(out_cd,argv[i+1]); print_cd = 1;}

		if (strcmp("-cut",argv[i]) == 0) {float temp;sscanf(argv[i+1],"%f",&temp); cutoff = temp;}
		if (strcmp("-tr",argv[i]) == 0) {int temp;sscanf(argv[i+1],"%i",&temp); n_tresh = temp;}
		if (strcmp("-ztr",argv[i]) == 0) {int temp;sscanf(argv[i+1],"%i",&temp); z_tresh = temp;}
		if (strcmp("-str",argv[i]) == 0) {float temp;sscanf(argv[i+1],"%f",&temp); s_tresh = temp;}
		if (strcmp("-htr",argv[i]) == 0) {float temp;sscanf(argv[i+1],"%f",&temp); h_tresh = temp; hard_tresh++;}
		
		if(strcmp("-bres", argv[i]) == 0) // Load the residues of binding site
		{
			j = 1;
			
			int tmp1;
			char tmp2;
			
			printf("Scanning bind\n");
			
			while(sscanf(argv[i+j], "%i_%c", &tmp1, &tmp2) > 0)
			{
				//printf("Match! %1i\t%1c\n", tmp1, tmp2);
				
				bresn[j - 1] = tmp1;
				
				bresc[j - 1][0] = tmp2;
				
				j++;
				
				nbind++;
				
				if(i + j == argc) {break;}
			}
		}
		
		if(strcmp("-ares", argv[i]) == 0) // Load the residues of allosteric site
		{
			j = 1;
			
			int tmp1;
			char tmp2;
			
			printf("Scanning allo\n");
			
			while(sscanf(argv[i+j], "%i_%c", &tmp1, &tmp2) > 0)
			{
				//printf("Match! %1i\t%1c\n", tmp1, tmp2);
				
				aresn[j - 1] = tmp1;
				
				aresc[j - 1][0] = tmp2;
				
				j++;
				
				nallo++;
				
				if(i + j == argc) {break;}
			}
		}
	}
	
	if (help_flag == 1)
	{
		printf("****************************\nHelp Section\n-i\tFile Input (PDB)\n-v\tVerbose\n-w\tWeight Vector\n-t\tInitial Value of template (negative value for random)\n\tIf Load Template, multiply the template\n-lt\tLoad input template\n-sp\tSuper Node Mode (CA, N, C)\n-kt\tPoid de l'angle entre les nodes (1)\n-kr\tPoid de la distance entre les nodes (1)\n-f\tFile to fit\n****************************\n");
		return(0);
	}
	
	//***************************************************
	//*													*
	//*	Build a structure contaning information on the initial pdb structure
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
	
	if (atom > 5000) {printf("Too many nodes! To fix, ask vincent.frappier@usherbrooke.ca\n");return(1);}
	
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
	
	printf("Check 0\n");
	
	// Scan strc_node for binding and allosteric nodes
	
	int bnodes[nbind];
	int anodes[nallo];
	
	for(i = 0; i < nbind; i++)
	{
		for(j = 0; j < atom; j++)
		{
			if(strc_node[j].res_number == bresn[i] && strcmp(bresc[i], strc_node[j].chain) == 0)
			{
				bnodes[i] = strc_node[j].node;
			}
		}
	}
	
	for(i = 0; i < nallo; i++)
	{
		for(j = 0; j < atom; j++)
		{
			if(strc_node[j].res_number == aresn[i] && strcmp(aresc[i], strc_node[j].chain) == 0)
			{
				anodes[i] = strc_node[j].node;
			}
		}
	}
	
	free(bresn);
	free(aresn);
	
	for(i = 0; i < 10000; i++)
	{
		free(bresc[i]);
		bresc[i] = NULL;
		
		free(aresc[i]);
		aresc[i] = NULL;
	}
	
	free(bresc);
	free(aresc);
	
	bresn = NULL;
	bresc = NULL;
	aresn = NULL;
	aresc = NULL;
	
	/*
	printf("Bind nodes :\n");
	
	for(i = 0; i < nbind; i++)
	{
		printf("%1i\n", bnodes[i]);
	}
	
	printf("\nAllo nodes :\n");
	
	for(i = 0; i < nallo; i++)
	{
		printf("%1i\n", anodes[i]);
	}
	*/
	
	printf("\nCheck 1\n");
	
	//***************************************************
	//*													*
	//*	Build connectivity matrix
	//*													*
	//***************************************************
	
	// Either load hessian or correlation matrix
	
	gsl_matrix *hess = gsl_matrix_alloc(3*atom,3*atom);
	
	gsl_matrix_set_all(hess, 0);
	
	if(eig_flag == 1) // Load hessian if eigen flag is 1
	{
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
	
	if(mat_flag == 1) // Load correlation matrix if matrix flag is 1
	{
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
	
	printf("Check 4\n");
	
	// Calculate average and sd of correlations
	
	float average = 0.0;
	
	for(i = 1; i < 3*atom; i++)
	{
		for(j = 0; j < i; j++)
		{
				average += gsl_matrix_get(correl, i, j);
		}
	}
	
	average /= (3.0 * (float)atom - 1.0) * 3.0 * (float)atom / 2.0;
	
	float sd = 0.0;
	
	for(i = 1; i < 3*atom; i++)
	{
		for(j = 0; j < i; j++)
		{
				sd += (gsl_matrix_get(correl, i, j) - average) * (gsl_matrix_get(correl, i, j) - average);
		}
	}
	
	sd /= (3.0 * (float)atom - 1.0) * 3.0 * (float)atom / 2.0 - 1.0;
	
	float flevel = 1.0 / sd;
	
	sd = sqrt(sd);
	
	printf("Avg : %1.5f\nSd : %1.5f\nFreedom level : %1.5f\nEstimated molar entropy : %1.5f J/(mol.K)\n", average, sd, flevel, 8.314 * log(flevel));
	
	printf("Check 5\n");
	
	// Print correlation matrix (if asked)
	
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
		
	if(cor_flag == 1)
	{
		return 0;
	}
	
	printf("Check 6\n");
	
	// Build connex matrix
	
	int **connex=(int **)malloc(atom*sizeof(int *));
	
	if(connex == NULL) {printf("malloc FAIL\n"); return 0;}
	
	for(i=0;i<atom;i++) {connex[i]=(int *)malloc(atom*sizeof(int)); if(connex[i] == NULL) {printf("malloc FAIL\n"); return 0;}}
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
	
	// Print correlation vs distance (if asked)
	
	if(print_cd == 1)
	{
		if(md_print(strc_all, all, connex, atom, out_cd) == 0) {printf("md_print failure\n"); return 0;};
	}
	
	gsl_matrix_free(correl);
	
	if(cde_flag == 1)
	{
		return 0;
	}
	
	// Complete connex matrix
	
	for(i = 0; i < atom; i++)
	{
		for(j = 0; j <= i; j++)
		{
			if(connex[i][j] >= z_tresh && get_connex(i, j, strc_all, all, cutoff) == 1)
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
	
	//***************************************************
	//*													*
	//*	Build network matrix
	//*													*
	//***************************************************
	
	int **cliques=(int **)malloc(2*atom*sizeof(int *));
	
	if(cliques == NULL) {printf("malloc FAIL\n"); return 0;}
	
	for(i=0;i<2*atom;i++) {cliques[i]=(int *)malloc((atom + 1)*sizeof(int)); if(cliques[i] == NULL) {printf("malloc FAIL\n"); return 0;}}
	for(i=0;i<2*atom;i++)for(j=0;j<atom;j++){cliques[i][j]=-1;}
	
	for(i=0;i<atom;i++) {cliques[i][atom] = 0;}
	
	int cliquen = 0;
	
	if(full_flag == 0)
	{
		// Fill cliques matrix
	
		Extend(nodenums, 0, atom, connex, clique, cpos, atom, n_tresh, cliques, atom, &cliquen);
	
		printf("Check 7\n");
	
		// Build networks from cliques
	
		int cliquen_old = cliquen;
	
		while(clique_fusion(cliques, atom, &cliquen) < cliquen_old) {printf("%1i\n", cliquen_old); cliquen_old = cliquen;}
	}
	else
	{
		printf("In progress!\n");
		
		return 0;
	}
	
	//***************************************************
	//*													*
	//*	Calculate local closeness
	//*													*
	//***************************************************
	
	// Print networks for other uses
	
	if(print_net == 1)
	{
		FILE *out_name_net;
	
		out_name_net = fopen(out_net, "w");
		
		printf("Writing in %s\n", out_net);
		
		for(i = 0; i < cliquen; i++)
		{
			fprintf(out_name_net, "Clique %1i:\n", i + 1);
		
			for(j = 0; j < cliques[i][atom]; j++)
			{
				fprintf(out_name_net, "%1i_%1s\t", strc_node[cliques[i][j]].res_number, strc_node[cliques[i][j]].chain);
			}
		
			fprintf(out_name_net,"\n\n");
		}
	}
	else
	{
		for(i = 0; i < cliquen; i++)
		{
			printf("Clique %1i:\n", i + 1);
		
			for(j = 0; j < cliques[i][atom]; j++)
			{
				printf("%1i_%1s\t", strc_node[cliques[i][j]].res_number, strc_node[cliques[i][j]].chain);
			}
		
			printf("\n\n");
		}
	}
	
	
	// Calculate local closeness for each combination of nodes
	
	if(print_lc == 1)
	{
		FILE *out_name_lc;
	
		out_name_lc = fopen(out_lc, "w");
	
		if(out_name_lc == NULL) {printf("Print_lc failure!\n"); return 0;}
	
		fprintf(out_name_lc, "Binding site :\n");
	
		for(i = 0; i < nbind; i++)
		{
			fprintf(out_name_lc, "%1i_%1s\t", strc_node[bnodes[i]].res_number, strc_node[bnodes[i]].chain);
		}
	
		fprintf(out_name_lc, "\n\nAllosteric site :\n");
	
		for(i = 0; i < nallo; i++)
		{
			fprintf(out_name_lc,  "%1i_%1s\t", strc_node[anodes[i]].res_number, strc_node[anodes[i]].chain);
		}
	
		fprintf(out_name_lc, "\n\n");
	
		float max_lc = 0.0;
		
		int bnode_max = 0;
		
		int anode_max = 0;
	
		for(i = 0; i < cliquen; i++)
		{
			for(j = 0; j < nbind; j++)
			{
				for(k = 0; k < nallo; k++)
				{
					float new_lc = lc(bnodes[j], connex, cliques[i], atom)*lc(anodes[k], connex, cliques[i], atom);
				
					if(new_lc > max_lc) {max_lc = new_lc; bnode_max = bnodes[j]; anode_max = anodes[k];}
				}
			}
		
			fprintf(out_name_lc, "Local closeness product in network %1i : %1.10f (Residues %1i_%1s and %1i_%1s)\n", i+1, max_lc, strc_node[bnode_max].res_number,
				strc_node[bnode_max].chain, strc_node[anode_max].res_number, strc_node[anode_max].chain);
		
			max_lc = 0.0;
		}
	}
	else
	{
		float max_lc = 0.0;
	
		for(i = 0; i < cliquen; i++)
		{
			for(j = 0; j < nbind; j++)
			{
				for(k = 0; k < nallo; k++)
				{
					float new_lc = lc(bnodes[j], connex, cliques[i], atom)*lc(anodes[k], connex, cliques[i], atom);
				
					if(new_lc > max_lc) {max_lc = new_lc;}
				}
			}
		
			printf("Local closeness product in network %1i : %1.10f\n", i+1, max_lc);
		
			max_lc = 0.0;
		}
	}
	
	return 0;
}

// XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
// XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
// XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

//***************************************************
//*													*
//*	Extend function
//*													*
//***************************************************

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

//***************************************************
//*													*
//*	clique_fusion function
//*													*
//***************************************************

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
		else
		{
			printf("Size of %1i : %1i\n", i + 1, clique_mat[i][mat_size]);
		}
	}
	
	return *cliquenum;
}

//***************************************************
//*													*
//*	get_connex function
//*													*
//***************************************************

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

//***************************************************
//*													*
//*	md_print function
//*													*
//***************************************************

int md_print(struct pdb_atom *init, int size_all, int **conx, int size_node, char *out)
{
	int i,j,k,l;
	
	int **match;
	
	match = (int **) malloc(size_node*sizeof(int *));
	
	if(match == NULL) {printf("malloc FAIL\n"); return 0;}
	
	for(i = 0; i < size_node; i++)
	{
		match[i] = (int *) malloc(size_node*sizeof(int));
		
		if(match[i] == NULL) {printf("malloc FAIL\n"); return 0;}	
	}
	
	for(i = 0; i < size_node; i++)
	{
		for(j = 0; j < size_node; j++)
		{
			match[i][j] = -1;
		}
	}
	
	float **dist;
	
	dist = (float **) malloc(size_node*sizeof(float *));
	
	if(dist == NULL) {printf("malloc FAIL\n"); return 0;}
	
	for(i = 0; i < size_node; i++)
	{
		dist[i] = (float *) malloc(size_node*sizeof(float));
		
		if(dist[i] == NULL) {printf("malloc FAIL\n"); return 0;}
	}
	
	for(i = 0; i < size_node; i++)
	{
		for(j = 0; j < size_node; j++)
		{
			dist[i][j] = 1000000.0;
		}
	}
	
	for(i = 0; i < size_all; i++)
	{
		for(j = i+1; j < size_all; j++)
		{
			if(init[i].node != init[j].node)
			{
				float sqdist = (init[i].x_cord - init[j].x_cord)*(init[i].x_cord - init[j].x_cord) + (init[i].y_cord - init[j].y_cord)*(init[i].y_cord - init[j].y_cord)
						+ (init[i].z_cord - init[j].z_cord)*(init[i].z_cord - init[j].z_cord);
				
				if(sqdist < dist[init[i].node][init[j].node])
				{
					dist[init[i].node][init[j].node] = sqdist;
					dist[init[j].node][init[i].node] = sqdist;
				}
				
				if(match[init[i].node][init[j].node] == -1 && conx[init[i].node][init[j].node] >= 3)
				{
					match[init[i].node][init[j].node] = 1;
					match[init[j].node][init[i].node] = 1;
				}
				else if(match[init[i].node][init[j].node] == -1)
				{
					match[init[i].node][init[j].node] = 0;
					match[init[j].node][init[i].node] = 0;
				}
			}
		}
	}
	
	printf("MD Check 5\n");
	
	FILE *out_md;
	
	out_md = fopen(out, "w");
	
	if(out_md == NULL) {return 0;}
	
	printf("MD Check 6\n");
	
	for(i = 0; i < size_node; i++)
	{
		for(j = i + 1; j < size_node; j++)
		{
			for(k = 0; k < size_all; k++)
			{
				if(init[k].node == i)
				{
					for(l = 0; l < size_all; l++)
					{
						if(init[l].node == j)
						{
							fprintf(out_md, "%1i_%s\t%1i_%s\t%1.6f\t%1i\n", init[k].res_number, init[k].chain, init[l].res_number, init[l].chain, sqrt(dist[i][j]), match[i][j]);
							
							break;
						}
					}
					
					break;
				}
			}
		}
	}
	
	return 1;
}

//***************************************************
//*													*
//*	local closeness function
//*													*
//***************************************************

float lc(int node, int **conx_mat, int *network, int size)
{
	int taken[network[size]];
	
	int i, j, k;
	
	int tkn_size = 1;
	
	int expl_size = 1;
	
	int level = 1;
	
	float lc = 0.0;
	
	for(i = 1; i < network[size]; i++)
	{
		taken[i] = -1;
	}
	
	taken[0] = node;
	
	int explore[network[size]];
	
	for(i = 1; i < network[size]; i++)
	{
		explore[i] = -1;
	}
	
	explore[0] = node;
	
	while(taken[network[size] - 1] == -1)
	{
		int explore_new[network[size]];
		
		int expl_newsize = 0;
		
		for(i = 0; i < expl_size; i++)
		{
			for(j = 0; j < network[size]; j++)
			{
				if(conx_mat[explore[i]][network[j]] == 1)
				{
					int istaken = 0;
					
					for(k = 0; k < tkn_size; k++)
					{
						if(taken[k] == network[j])
						{
							istaken = 1;
							
							break;
						}
					}
					
					if(istaken == 0)
					{
						explore_new[expl_newsize] = network[j];
						
						lc += 1.0 / ((float)level * (float)level);
						
						expl_newsize ++;
						
						taken[tkn_size] = network[j];
						
						tkn_size ++;
					}
				}
			}
		}
		
		if(tkn_size == 1)
		{
			break;
		}
		
		level ++;
		
		expl_size = expl_newsize;
		
		for(i = 0; i < expl_newsize; i++)
		{
			explore[i] = explore_new[i];
		}
		
		for(i = expl_size; i < network[size]; i++)
		{
			explore[i] = -1;
		}
	}
	
	return lc;
}
