#include <stdio.h>
#include <string.h>
#include <gsl/gsl_matrix.h>
#include <math.h>
#include <gsl/gsl_eigen.h>
#include <time.h>
#include <gsl/gsl_linalg.h>

struct pdb_atom {
	 	int atom_number; /*Numéro d'atomes*/
	 	float x_cord; /*Coordoné en x*/
	 	float y_cord; /*Coordoné en y*/
	 	float z_cord; /*Coordoné en z*/
	 	char atom_prot_type[6]; /*Si carbone alpha, oxygène, etc.*/
	 	unsigned int res_number; /*Numéro de résidus*/
	 	unsigned int atom_type; /*Atom (=1) ou hétéroatome (=2) ou ligand (=3) ou ADN/ARN (=4) ou ADN/ARN HET (=5)*/
	 	char res_type[6]; /*Type de resisdus*/
	 	float b_factor; /*b_factor*/
	 	char chain[6];
	 	unsigned int type;
	 	int node;
	 	int node_c[6];
	 	float mass;
};


// STeM_lib_boinc.c
//	STeM boinc
int boinc_load_matrix(gsl_matrix *newstrc,char filename[100]);
void boinc_read_vcon(char filename[100],int atom, struct pdb_atom *strc, gsl_matrix *m);
int boinc_count_connect(char filename[100]);
int boinc_assign_connect(char filename[100],int **con);
int boinc_count_atom(char filename[100]);
int boinc_assign_cord_all_atom(char filename[100],struct pdb_atom *structure,int node);
void boinc_apply_eigen(struct pdb_atom *strc,int atom,gsl_matrix *m,int mode,float amplitude);
float boinc_overlap(int atom, int mode,gsl_matrix *m,gsl_vector *d);
float boinc_rmsd_no(struct pdb_atom *init,struct pdb_atom *targ,int atom, int *align);
//void fit_vince(struct pdb_atom *init,struct pdb_atom *targ,int atom,int all,struct pdb_atom *all_init,struct pdb_atom *all_targ,gsl_matrix *evec, int *align, int nb_mode, int mode,gsl_vector *eval);
float fit_svd(struct pdb_atom *init,struct pdb_atom *targ,int atom,int all,int atom_t,int all_t,struct pdb_atom *all_init,struct pdb_atom *all_targ,gsl_matrix *evec, int *align, int nb_mode, int mode,gsl_vector *eval);
void boinc_assign_atom_type(struct pdb_atom *strc,int atom);
void boinc_all_interaction(struct pdb_atom *strc,int atom,int res_n, gsl_matrix *ma, int node,int flag,gsl_matrix *mvcon,gsl_matrix *mint);

// STeM_lib_fit.c
//	STeM fit
int node_align(struct pdb_atom *strc,int atom,struct pdb_atom *strc_t,int atom_t, int *align);
int node_align_low(struct pdb_atom *strc,int atom,struct pdb_atom *strc_t,int atom_t, int *align);
int node_align_lig(struct pdb_atom *strc,int atom,struct pdb_atom *strc_t,int atom_t, int *align,struct pdb_atom *strc_all,int all,struct pdb_atom *strc_all_t,int atom_all,float cutoff);
void rotate_all(gsl_matrix *rota,struct pdb_atom *all_init,int all);
void write_movie(FILE *out_file, struct pdb_atom *newstrc,int nb_atom,int model);
double gsl_matrix_Det3D(gsl_matrix *M);
void multiplie_matrix(gsl_matrix *a,int a_one,int a_two,gsl_matrix *b,int b_one,int b_two,gsl_matrix *c);
float rmsd_no(struct pdb_atom *init,struct pdb_atom *targ,int atom, int *align);
float rmsd_yes(struct pdb_atom *init,struct pdb_atom *targ,int atom, int *align,struct pdb_atom *init_all,int all);
void apply_eigen(struct pdb_atom *strc,int atom,gsl_matrix *m,int mode,float amplitude);
void center_yes(struct pdb_atom *init,struct pdb_atom *targ,int atom,int atom_t, int *align);
float vector_lenght(gsl_vector *d,int atom);
float overlap(int atom, int mode,gsl_matrix *m,gsl_vector *d, int *align);
void fit(struct pdb_atom *init,struct pdb_atom *targ,int atom,int all,struct pdb_atom *all_init,struct pdb_atom *all_targ,gsl_matrix *eval, int *align);
void fit_math(struct pdb_atom *init,struct pdb_atom *targ,int atom,int all,struct pdb_atom *all_init,struct pdb_atom *all_targ,gsl_matrix *eval, int *align);
void fit_vince(struct pdb_atom *init,struct pdb_atom *targ,int atom,int all,struct pdb_atom *all_init,struct pdb_atom *all_targ,gsl_matrix *eval, int *align, int nb_mode, int mode);
void nrg_rmsd(struct pdb_atom *init,int atom,gsl_matrix *evec, int *align, int nb_mode, int mode,gsl_vector *eval);
void fit_mc_torsion(struct pdb_atom *init,struct pdb_atom *targ,int atom,int all,int atom_t,int all_t,struct pdb_atom *all_init,struct pdb_atom *all_targ,gsl_matrix *evec, int *align, int nb_mode, int mode,gsl_vector *eval);

// STeM_lib_grid_motion.c
//	STeM grid_motion
void write_matrix_pdb(char filename[100], gsl_matrix *m,int nb_atom,int nb_atom_1);
void load_eigen_motion(gsl_vector *v,gsl_matrix *m,char filename[100],int atom, int mode);
void print_image_motion(struct pdb_atom *old,gsl_matrix *m,int mode,float max,float min,int atom,char out_name[100], int lig_f);
void print_grid_motion(struct pdb_atom *old,gsl_matrix *m,int mode,int atom,char out_name[100],gsl_matrix *g,int line_l,int nbr,int lig);
int test_grid(double **grid,int nbr_mode,int mode,gsl_matrix *matrix,struct pdb_atom *strc,int atom,float max_dist, int point);
void find_max_ampli(struct pdb_atom *strc,gsl_matrix *m,int mode,int atom,float max,float *max_value,float *min_value);
void find_max_ampli_two(struct pdb_atom *strc,gsl_matrix *m,int mode,int atom,float max,double *max_value,double *min_value,struct pdb_atom *ori);
int test_point(struct pdb_atom *strc, double *actual, float limit, int atom,int mode, int nb_mode, gsl_matrix *matrix);
int test_point_two(struct pdb_atom *strc, double *actual, float limit, int atom,int mode, int nb_mode, gsl_matrix *matrix);
int build_grid(double **grid,double *resolution,int nbr_mode,int atom,gsl_matrix *matrix,struct pdb_atom *strc,float limit,int mode, int now, double *actual);
void print_image_torsion(struct pdb_atom *old,gsl_matrix *evec,int mode,float max,float min,int atom,char out_name[100],int nbnode,struct pdb_atom *strc_node);

// STeM_lib_hessian.c
//	STeM hessian
void build_2_matrix(struct pdb_atom *strc, double **hes,int nb_atom,double K_theta);
void build_3_matrix(struct pdb_atom *strc, double **hes,int nb_atom,double K_phi,float factor);
void build_1st_matrix(struct pdb_atom *strc, double **hes,int nb_atom,double K_r);
void build_4h_matrix(struct pdb_atom *strc, double **hes,int nb_atom,double epsi,gsl_matrix *tmp);
void build_enm(struct pdb_atom *strc, double **hes,int nb_atom,double epsi,gsl_matrix *tmp,float cutoff);
void diagonalyse_matrix (gsl_matrix *m,int nb_atom, gsl_vector *evl,gsl_matrix *evc);
void mass_weight_hessian(gsl_matrix *m,int atom,struct pdb_atom *strc_all);
void adjust_weight_evec(gsl_matrix *m,int atom,struct pdb_atom *strc);

// STeM_lib_other.c
//	STeM other
void load_eigen(gsl_vector *v,gsl_matrix *m,char filename[100],int atom);
void write_eigen(char filename[100], gsl_matrix *m,gsl_vector *v,int nb_atom);
void write_matrix(char filename[100], gsl_matrix *m,int nb_atom,int nb_atom_2);
float ran2(long *idum);
long time_seed();
int load_matrix(gsl_matrix *newstrc,char filename[100]);
void load_eigen_grid(gsl_vector *v,gsl_matrix *m,char filename[100],int atom, int mode);
void write_eigen_mat(char filename[100], gsl_matrix *m,int nb_atom,int nb_atom_1,int mode, int nbr_mode);
void write_grid_mat(char filename[100], gsl_matrix *m,int nb_atom,int nb_atom_1);
double dot_p(double a[3], double b[3]);
void write_strc(char filename[100], struct pdb_atom *newstrc,int nb_atom);
void write_strc_b(char filename[100], struct pdb_atom *old,int nb_atom,gsl_matrix *m, int node);
void assignArray(gsl_matrix *m,double **a,int count,int count_1);
void correlate_sp(gsl_matrix *m,struct pdb_atom *strc, int atom);
int count_atom_CA_n(struct pdb_atom *all,int atom, int node, int ligand);
float correlate(gsl_matrix *m,struct pdb_atom *strc, int atom);
void k_inverse_matrix_stem(gsl_matrix *m,int nb_atom, gsl_vector *evl,gsl_matrix *evc);
void buid_pre(gsl_matrix *m,gsl_matrix *evc,int nb_atom);
void k_inverse_matrix_stem_temp(gsl_matrix *m,int nb_atom, gsl_vector *evl,gsl_matrix *evc);
float calc_energy(int atom,gsl_vector *eval);


// STeM_lib_rot.c
//	STeM rot
double vector_lenght_v(gsl_vector *v);
void print_vector(gsl_vector *v);
double vector_scalar_v(gsl_vector *a,gsl_vector *b);
double vector_angle_v(gsl_vector *a,gsl_vector *b);
void rotate_matrix(gsl_matrix *rota,double angle,gsl_vector *v);
void center_node(int all,struct pdb_atom *all_init,int node);
void dot_product_v(gsl_vector *a,gsl_vector *b,gsl_vector *c);
void assign_vector(struct pdb_atom *fstr,int fnode,char fstring[5],int fatom,struct pdb_atom *sstr,int snode,char sstring[5],int satom,gsl_vector *v);
void rotate_phi(gsl_matrix *rota,struct pdb_atom *all_init,int all,int node);
void rotate_psy(gsl_matrix *rota,struct pdb_atom *all_init,int all,int node);
void print_matrix(gsl_matrix *mat);
void rotate_bb(struct pdb_atom *init,struct pdb_atom *targ,int atom,int all,struct pdb_atom *all_init,struct pdb_atom *all_targ,int *align);


// STeM_lib_STRC.c
//	STeM STRC
void copy_strc(struct pdb_atom *target, struct pdb_atom *initial, int atom);
int covalent_bond(int i, int j,int **con, int ncon);
int count_connect(char filename[100]);
int assign_connect(char filename[100],int **con);
int check_covalent_CA(struct pdb_atom *CA,struct pdb_atom *strc,int atom,int all,int a,int b);
void check_lig(struct pdb_atom *strc,int **con,int ncon, int atom);
int count_atom(char filename[100]);
int build_cord_CA(struct pdb_atom *all, struct pdb_atom *CA, int atom, int ligand,int **con, int ncon);
int assign_cord_all_atom(char filename[100],struct pdb_atom *structure);
int build_all_strc(char filename[100],struct pdb_atom *structure);

// STeM_lib_templaate.c
//	STeM templaate
void read_vcon(char filename[100],int atom, struct pdb_atom *strc, gsl_matrix *m);
void print_templaate(struct pdb_atom *newstrc,int atom,gsl_matrix *m,char filename[100],float min,float max);
double templaate_average(gsl_matrix *m,int atom);
void all_interaction(struct pdb_atom *strc,int atom,int res_n, gsl_matrix *ma, int flag,gsl_matrix *mvcon,gsl_matrix *mint,struct pdb_atom *ca_strc);
void assign_atom_type(struct pdb_atom *strc,int atom);
void assign_lig_type(struct pdb_atom *strc, int atom, char inp[500]);

// STeM_lib_vcon.c
//	Vcon pimper !

int vcon_file(char FILENAME[60],struct pdb_atom *strc, gsl_matrix *m,int allatom);
int vcon_file_dom(struct pdb_atom *strc, gsl_matrix *m, int allatom);
