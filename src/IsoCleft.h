#include "general.h"
#include <gsl/gsl_linalg.h>
#include "Random_Gen.c"

// TO COMPILE:
// In Mac OS X with GSL 1.13 in /usr/local/bin
// gcc -I/usr/local/include -c IsoCleft.c
// gcc IsoCleft.o -L/usr/local/lib -lgsl -lgslcblas -lm -o IsoCleft
// In Linux
// gcc -I/sw/arch/include -c IsoCleft.c
// gcc IsoCleft.o -L/sw/arch/lib -lgsl -lgslcblas -lm -o IsoCleft

//#define ASSIGN_ATOM_TYPE atom_types_ETY
//#define NUM_OF_ATOM_TYPES 5
//#define ASSIGN_ATOM_TYPE atom_types_LPC
//#define NUM_OF_ATOM_TYPES 9
int num_of_atom_types;

#define MAX_COMMAND_LINE_LEN 500

/*Define Boolean type */
typedef	enum { FALSE, TRUE }	bool;

typedef struct AtomStructure sAtom;
typedef sAtom *psAtom;

struct AtomStructure{
  char   line[82];   // PDB line
  char   nline[82];   // PDB line
  char   name[5];    // atom name
  float  radius;     // atomic radius
  float  coor[3];    // coordinates
  int    inum;
  int    pnum;
  int    type;
  char   rnam[4];
  int    rnum;
  psAtom next_ca,prev_ca;
  psAtom next,prev;
};

typedef struct AssocGraphVertex sAGVx;
typedef sAGVx *psAGVx;

struct AssocGraphVertex{
  psAtom a;
  psAtom b;
  int    inum;
  psAGVx next;
};

/*
  typedef struct CalphaStruct sCalpha;
  typedef sCalpha *psCalpha;
  struct CalphaStruct{
  psAtom   atm;
  psCalpha next,prev;
  };
*/

psAtom head_ca_a = NULL;
psAtom head_ca_b = NULL;
int    ca_nodes = 0;
int    ca_edges = 0;
int    nca_a = 0;
int    nca_b = 0;

char   command[MAX_COMMAND_LINE_LEN];
psAtom atoms_a = NULL;
psAtom atoms_b = NULL;
psAtom all_a   = NULL;
psAtom all_b   = NULL;
psAtom het_a   = NULL;
psAtom het_b   = NULL;
psAtom *rat_a  = NULL;
psAtom *rat_b  = NULL;
int    *num_a  = NULL;
int    *num_b  = NULL;

/*
  psCalpha ca_a  = NULL;
  psCalpha ca_b  = NULL;
  int      nca_a = 0;
  int      nca_b = 0;
*/

int het_a_isnull;
int het_b_isnull;

//psAGVx ca_agv = NULL;

psAGVx agv     = NULL;
int    get_maximal_clique;
int    max_clique_size;

//float R[3][3];                              // Rotation Matrix
//float cen_a[3],cen_b[3];
char cleftfile_a[100];
char cleftfile_b[100];
char hetfile_a[100];
char hetfile_b[100];
char pdbfile_a[100];
char pdbfile_b[100];
char outbase[100];
int  superpose_pdb=0;
int  flag_all,flag_het;

float Delta_Dist_Threshold,CA_Delta_Dist_Threshold,Threshold_Neighbourhood;
int   Skip_Neighbourhood_Test;
int   bk_mode;

typedef struct tCliqueStruct tsClique;
typedef tsClique *tClique;

struct tCliqueStruct{
  int        Vnum;
  int        *Vertices;
  float      cen_a[3];
  float      cen_b[3];
  gsl_matrix *mat_r;
  double     DetR;
  float      rmsd;
  float      Similarity;
  int        inum;
  int        reject;
  tClique    next,prev;
};



int     *connected = NULL;
int     c,numNodes,numEdges;
int     *compsub   = NULL;
tClique Cliques    = NULL;
int     Clique_threshold=3;
int     clique_inum=-1;
int     got_one = 0;
int     get_all = 1;
int     ca_only;
int     all_atom_only;
int     jtt[20][20];
int     JTT_rank_threshold;
int     output_superimposedPDBs;

int  bk();
void Extend(int old[],int ne,int ce);
void AddNewClique(int n, int list[]);
void Random_bin_graph(float fraction);

psAtom ReadCleftFile(char file[], char chn);
int    (*atom_type_function)(char rnam[], char anam[]);
int    atom_types_ETY(char rnam[], char anam[]);
int    atom_types_LPC(char rnam[], char anam[]);
float  dist3d(float a[], float b[]);
float  assign_radius(char atm[]);
void   Calculate_Rotation_GSL(tClique c);
void   Print_Superposed_Clique(tClique c,int clique_counter);
void   Output_Clique_Stats(tClique c,int clique_counter);
void   read_commandline(int argc, char *argv[]);
double SupSVD(gsl_matrix *mat_u);
double gsl_matrix_Det3D(gsl_matrix *M);
int    Filter_Cliques(int total);
void   sort_Cliques();
void   swap_Cliques(tClique b, tClique d);
void   non_redundant_Cliques();
int Accept_Correspondence(char rnam_a[],char rnam_b[]);
tClique Get_Top_Clique();
int nam2n(char rnam[]);
void Threshold_JTT_Matrix(int threshold);
