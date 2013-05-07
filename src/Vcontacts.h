
/**********************************************************************
 * PROGRAM Vcontacts.h
 * calculates the atom-atom contacts within a protein, along with the
 * solvent accessible surface.
 *
 * ====================================================================

 * NOTES ON METHODOLOGY:
 *
 * The method calculates the planes of contact the center atom makes 
 * will each atom in the contact list, then determines the points
 * of intersection of the planes and intersections with a sphere. The
 * sphere has a radius equal to the sum of the van der Waals radius of the 
 * center atom plus the radius of a solvent atom (water). The set of 
 * intersections of the planes defines a polygon surrounding the center 
 * atom. This polygon is projected onto the surface of the sphere. The
 * area of each projection is calculated as the sum of spherical triangles
 * and arc segments, where an arc segment is the difference between a 
 * angluar segment of a spherical cap and the corresponding spherical triangle.
 *
 * The points of intersection are calculated using a convex hull algorithm.
 * The algorithm operates with efficiency O(Nk), where N is the number of 
 * input atoms and k is the number of edges on the contact polyhedron.
 *
 * The convex hull algorithm effectively describes the atom contacts only if
 * no 'engulfing' atoms are present, ie., atoms which cover more than 50% of
 * the atom contact surface. This may occur between closely bonded atoms of
 * different radii, such as some C=O bonds. A correction factor accounts for 
 * this.
 **********************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <math.h>
#include "STeM.h"

#define MAX_CONT 100
#define MAX_PATH 200
#define CELLSIZE 7.0
#define PAUSE fgets(stp,sizeof(stp),stdin) /* shortcut for pausing */

// ------------------------- structure definitions ----------------------
struct AtomCalcSAS_struct {
  int   atomnum;      // record number from PDB
  float coor[3];      // xyz coordinates of atom
  char  atomname[5];  // name of atom, CA, etc
  char  res[4];       // residue name from PDB
  int   resnum;
  int   inum;         // internal number in FA structurex
  char  chn;          // chain id '-' for blank ADDED RJN 29.0.4.2008
  int   type;         // atom type (interaction)
  float radius;       // radius of atom
  int   boxnum;       // the box number atom is assigned to.
  char  source;       // source of atom, protein, ligand, 
  float SAS;          // solvent exposed surface area
  float vol;          // atom volume
  char  done;         // flag if atom contacts have already been calculated
};
typedef struct AtomCalcSAS_struct atomsas;

struct AtomIndex_struct {
  int   nument;      // number of entries in box 
  int   first;       // location of first entry in PDBlist
};
typedef struct AtomIndex_struct atomindex;

struct ContactList_struct {
  int     index;   // index to PDB atom
  double  area;    // contact area, square angstroms
  double  dist;    // distance to atom zero
  char    flag;    // to keep or not. 'X' == omit.
};
typedef struct ContactList_struct contactlist;
  
struct Plane_struct {
  double Ai[4];      // parameters A,B,C,D of contact plane Ax+By+Cz+D=0 
  double dist;       // distance from plane to origin
  int    index;      // index to which record in PDB or ligand array.
  double area;       // contact area in square angstroms
  char   flag;       // 'X' if no contact, 'E' if an engulfing atom.
};
typedef struct Plane_struct plane;

struct Vertex_struct {
  double xi[3];      // x,y,z coordinates (x1,x2,x3)
  double dist;       // distance to origin
  int    plane[3];   // identification of intersecting planes. -1 = sphere.
};
typedef struct Vertex_struct vertex;

struct ptIndex_struct {
  int  numpts;       // number of points defining face
  int  pt[MAX_CONT]; // index to polyhedron points
};
typedef struct ptIndex_struct ptindex;

struct EdgeVector_struct {
  double V[3];       // vector for edge
  int startpt;       // initial vertex
  int endpt;         // final vertex
  int plane[2];      // planes defining edge (using single atoms for now)
  int startplane;    // third plane at start point
  int endplane;      // third plane at end point
  char arc;          // flag for arc point calculations
};
typedef struct EdgeVector_struct edgevector;

struct ca_struct {
  int prev;     // previous contact location in ca_index
  int atom;     // PDBarray number of current contact (NOT PDB record number)
  float area;        // contact area
  float dist;        // distance between atoms
};
typedef struct ca_struct ca_struct;

struct VC_Global_struct{
  
  int  first;    // reference CF calculations (can recalculate)

  // ----------------- Global variables -----------------
  atomsas   *Calc;      // pointer to PDB array (dynamically allocated)
  atomindex *box;      // index to PDB atoms within cubic grid
  int       *Calclist;             // list of atoms ordered by box number
  contactlist *contlist;
  
  ptindex    ptorder[1100];  // for ordering vertices around each face // RJN 081008 changed 100 to 1100 and 200 2200 
  vertex     centerpt[1100]; // center points for each contact face    // RJN 081008
  vertex     poly[2200];     // polyhedron vertices                    // RJN 081008
  plane      cont[1100];     // atom and contact plane information     // RJN 081008
  edgevector vedge[2200];
  ca_struct  *ca_rec;        // array - contact area records
  int   *ca_index;              // array - index to first ca_recs for each atom.
  int   numcarec;          
  int   ca_recsize;
  int        dim;                   // dimension in units CELLSIZE
  char       ch;                    // for debugging 
  char       planedef;              // = X, R, or B
  char       showbonded;             // = Y or N.
  char       normalize;              // = Y or N. normalize areas to area of sphere.
  int        *seed;                  // seed vertices for new polyhedra
  char       radfilename[MAX_PATH];  // custom name for radii.dat file including path ADDED RJN 28.04.2008
  // --------------------- function prototypes ----------------------
};
typedef struct VC_Global_struct VC_Global;

float   Vcontacts(struct pdb_atom *strc,resid*,VC_Global*,int);

float   ic2cf(VC_Global*,atom*,resid*,rot*,gridpoint*,float**,int,float*);
cfstr   vcfunction(VC_Global*,atom*,resid*,int); 
float   xs2cf(VC_Global*,atom*,resid*,int,int*);
float   pb2cf(VC_Global*,atom*,resid*,rot*,int,int*, int);

double  spherical_arc(vertex, vertex, vertex, float);
double  cosPQR(double *, double *, double *);
char    test_point(double *, plane *, int, float, int, int, int);
char    order_faces(int, vertex *, vertex *, float, int, int, plane *, ptindex *);
void    project_points(VC_Global*,vertex *, vertex *, float, int, int, plane *, int);
int     voronoi_poly2(VC_Global*,int, plane *, float, int, contactlist *);
int     add_vertex(vertex *, int, double *, int, int, int);
void    add_vedge(edgevector *, int, plane *, int, int, int, vertex *, int);
int     solve_3x3(double *, double *, double *, double *);
int     solve_2xS(plane, plane, float, double *, double *);
float   calc_region(VC_Global*,atom*,int,int);
void    calc_areas(VC_Global*,vertex *, vertex *, float, int, int, plane *, ptindex *, int);
void    index_protein(atom*,resid*,VC_Global*,int,int);
void    save_areas(VC_Global*,plane *, contactlist *, int, int);
int     get_contlist4(VC_Global*,atom*,int, contactlist *, int, float, int, int);
void    save_seeds(VC_Global*,plane *, vertex *, int, int);
void    get_firstvert(VC_Global*,plane *, int *, int *, int *, int, int);


