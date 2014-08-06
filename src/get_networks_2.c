#include "STeM.h"
#include <gsl/gsl_linalg.h>
#include <cstdlib>
#include "math.h"

double D_B(int node1, int node2, struct pdb_atom *strc_node, gsl_matrix *cov);

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
	
	int print_correl = 0; // Tells the program to print the correlation matrix
	int print_cd = 0; // Tells the program to print correlation vs distance of residues
	
	int mat_flag = 0; // Tells program to load correlation matrix
	int eig_flag = 0; // Tells program to load hessian
	
	int help_flag = 1;
	int verbose = 0;
	int lig = 0;
	int cde_flag = 0; // Tells the program to end after printing correlation vs distance of residues
	int cor_flag = 0; // Tells the program to end after calculating the correlation matrix
	
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
	
	printf("Check 3\n");
	
	gsl_matrix *correl = gsl_matrix_alloc(atom, atom);
	
	for(i = 0; i < atom; i++)
	{
		for(j = 0; j < atom; j++)
		{
			gsl_matrix_set(correl, i, j, (abs(gsl_matrix_get(hess,3*i,3*j)) + abs(gsl_matrix_get(hess,3*i + 1,3*j + 1)) + abs(gsl_matrix_get(hess,3*i + 2,3*j + 2))) / sqrt((abs(gsl_matrix_get(hess,3*i,3*i)) + abs(gsl_matrix_get(hess,3*i + 1,3*i + 1)) + abs(gsl_matrix_get(hess,3*i + 2,3*i + 2))) * (abs(gsl_matrix_get(hess,3*j,3*j)) + abs(gsl_matrix_get(hess,3*j + 1,3*j + 1)) + abs(gsl_matrix_get(hess,3*j + 2,3*j + 2)))) / D_B(i, j, strc_node, hess));
		}
	}
	
	gsl_matrix_free(hess);
	
	
	
	return 0;
}

// XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
// XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
// XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

//***************************************************
//*													*
//*	D_B function (Bhattacharyya distance)
//*													*
//***************************************************

double D_B(int node1, int node2, struct pdb_atom *strc_node, gsl_matrix *cov)
{
	int p,q;
	
	gsl_matrix *ani1 = gsl_matrix_alloc(3,3);
	
	gsl_matrix *ani2 = gsl_matrix_alloc(3,3);
	
	gsl_matrix *ani12 = gsl_matrix_alloc(3,3);
	
	for(p = 0; p < 3; p++)
	{
		for(q = 0; q < 3; q++)
		{
			gsl_matrix_set(ani1, p, q, gsl_matrix_get(cov, 3*node1 + p, 3*node1 + q));
			gsl_matrix_set(ani2, p, q, gsl_matrix_get(cov, 3*node2 + p, 3*node2 + q));
			
			gsl_matrix_set(ani12, p, q, (gsl_matrix_get(ani1, p, q) + gsl_matrix_get(ani1, p, q)) / 2.0);
		}
	}
	
	double bd = gsl_matrix_get(ani12,0,0)*(strc_node[node2].x_cord - strc_node[node1].x_cord)*(strc_node[node2].x_cord - strc_node[node1].x_cord);
	bd += gsl_matrix_get(ani12,1,1)*(strc_node[node2].y_cord - strc_node[node1].y_cord)*(strc_node[node2].y_cord - strc_node[node1].y_cord);
	bd += gsl_matrix_get(ani12,2,2)*(strc_node[node2].z_cord - strc_node[node1].z_cord)*(strc_node[node2].z_cord - strc_node[node1].z_cord);
	bd += 2*gsl_matrix_get(ani12,0,1)*(strc_node[node2].x_cord - strc_node[node1].x_cord)*(strc_node[node2].y_cord - strc_node[node1].y_cord);
	bd += 2*gsl_matrix_get(ani12,0,2)*(strc_node[node2].x_cord - strc_node[node1].x_cord)*(strc_node[node2].z_cord - strc_node[node1].z_cord);
	bd += 2*gsl_matrix_get(ani12,1,2)*(strc_node[node2].y_cord - strc_node[node1].y_cord)*(strc_node[node2].z_cord - strc_node[node1].z_cord);
	
	gsl_linalg_chloesky_decomp(ani1);
	gsl_linalg_chloesky_decomp(ani2);
	gsl_linalg_chloesky_decomp(ani12);
	
	bd += 0.5*log(gsl_matrix_get(ani12,0,0)*gsl_matrix_get(ani12,0,0)*gsl_matrix_get(ani12,1,1)*gsl_matrix_get(ani12,1,1)*gsl_matrix_get(ani12,2,2)*gsl_matrix_get(ani12,2,2));
	bd -= 0.5*log(gsl_matrix_get(ani1,0,0)*gsl_matrix_get(ani1,1,1)*gsl_matrix_get(ani1,2,2)*gsl_matrix_get(ani2,0,0)*gsl_matrix_get(ani2,1,1)*gsl_matrix_get(ani2,2,2));
	
	return bd;
}
