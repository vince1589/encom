 #include "STeM.h"
 
 double vector_lenght_v(gsl_vector *v) {
	return(sqrt(gsl_vector_get(v,0)*gsl_vector_get(v,0)+
		 gsl_vector_get(v,1)*gsl_vector_get(v,1)+
		 gsl_vector_get(v,2)*gsl_vector_get(v,2)));
}
 
 void print_vector(gsl_vector *v) {
	int i;
	for(i=0;i<3;++i) {
		printf("%10.5f ",gsl_vector_get(v,i));
	}
	printf("\n");
} 
 
 double vector_scalar_v(gsl_vector *a,gsl_vector *b) {
	double value = 0;
	int i;
	for(i=0;i<3;++i) {
		value += gsl_vector_get(a,i)*gsl_vector_get(b,i);
	}
	return(value);
}
 
 double vector_angle_v(gsl_vector *a,gsl_vector *b) {
	double pi = 3.14159265358979323846264338327950288419716939937510;
	//printf("********************************\nFUNCTION D'ANGLE\nPRINTING VECTOR:\n");
	
	//print_vector(a);
	//print_vector(b);
	//printf("SCARLAR:%f	Lenght:%f %f\n",vector_scalar_v(a,b),vector_lenght_v(a),vector_lenght_v(b));
	double value = (vector_scalar_v(a,b)/vector_lenght_v(a)/vector_lenght_v(b));
	//printf("VALUE:%f\n",value);
	if ((value < 1.00000001) && (value > 0.99999999)) {return(0);}
	if ((value < pi*2 + 0.000001) && (value > pi*2 - 0.000001)) {return(1);}
	value = acos(value);
	//printf("VALUE AFTER:%g\n",value);
	return(value);
}
 void rotate_matrix(gsl_matrix *rota,double angle,gsl_vector *v) {
	if (vector_lenght_v(v) == 0) {return;}
	gsl_vector_scale(v,1/vector_lenght_v(v));
	double xv = gsl_vector_get(v,0);
	double yv = gsl_vector_get(v,1);
	double zv = gsl_vector_get(v,2);

	gsl_matrix_set(rota,0,0,cos(angle)+xv*xv*(1-cos(angle)));
	gsl_matrix_set(rota,0,1,yv*xv*(1-cos(angle))+zv*sin(angle));
	gsl_matrix_set(rota,0,2,zv*xv*(1-cos(angle))-yv*sin(angle));
	gsl_matrix_set(rota,1,0,yv*xv*(1-cos(angle))-zv*sin(angle));
	gsl_matrix_set(rota,1,1,cos(angle)+yv*yv*(1-cos(angle)));
	gsl_matrix_set(rota,1,2,zv*yv*(1-cos(angle))+xv*sin(angle));
	gsl_matrix_set(rota,2,0,zv*xv*(1-cos(angle))+yv*sin(angle));
	gsl_matrix_set(rota,2,1,zv*yv*(1-cos(angle))-xv*sin(angle));
	gsl_matrix_set(rota,2,2,cos(angle)+zv*zv*(1-cos(angle)));
}

void center_node(int all,struct pdb_atom *all_init,int node) {
	int i = node,j,k;
	float xi,yi,zi,xt,yt,zt; 
	// Temp value
	for (i=0;i<all;++i) {
		if ((all_init[i].node == node) && (strcmp(all_init[i].atom_prot_type," CA  ") == 0)) {
			xi = all_init[i].x_cord;yi = all_init[i].y_cord; zi = all_init[i].z_cord;
			break;
		}
	}
	// On soustrait à tous les atomes
 	for(j=0;j<all;++j) {
 		all_init[j].x_cord -= xi;  all_init[j].y_cord -= yi;  all_init[j].z_cord -= zi;
 	}

}

void dot_product_v(gsl_vector *a,gsl_vector *b,gsl_vector *c) {
//	printf("********************************\nFUNCTION DOT PRODUCT\nPRINTING VECTOR:\n");
//	print_vector(a);
//	print_vector(b);
	gsl_vector *d = gsl_vector_alloc(3);
	gsl_vector_set(d,0,gsl_vector_get(a,1)*gsl_vector_get(b,2)-gsl_vector_get(a,2)*gsl_vector_get(b,1));
	gsl_vector_set(d,1,gsl_vector_get(a,2)*gsl_vector_get(b,0)-gsl_vector_get(a,0)*gsl_vector_get(b,2));
	gsl_vector_set(d,2,gsl_vector_get(a,0)*gsl_vector_get(b,1)-gsl_vector_get(a,1)*gsl_vector_get(b,0));
	
	gsl_vector_set(c,0,gsl_vector_get(d,0));
	gsl_vector_set(c,1,gsl_vector_get(d,1));
	gsl_vector_set(c,2,gsl_vector_get(d,2));
	
//	printf("PRODUCT\n");
//	print_vector(c);
}

void assign_vector(struct pdb_atom *fstr,int fnode,char fstring[5],int fatom,struct pdb_atom *sstr,int snode,char sstring[5],int satom,gsl_vector *v) {
	int i,j;
	for(i=0;i<fatom;++i) {
	//	printf("I:%d	Type:%s	Node:%d.... Targ:%d\n",i,fstr[i].atom_prot_type,fstr[i].node,fnode);
		if ((fstr[i].node == fnode) && (strcmp(fstr[i].atom_prot_type,fstring) == 0)) {
			for(j=0;j<satom;++j) {
				//printf("	J:%d	Type:%s	Node:%d...Targ:%d\n",j,sstr[j].atom_prot_type,sstr[j].node,snode);
		 		if ((sstr[j].node == snode) && (strcmp(sstr[j].atom_prot_type,sstring) == 0)) {
		 			gsl_vector_set(v,0,fstr[i].x_cord - sstr[j].x_cord);
		 			gsl_vector_set(v,1,fstr[i].y_cord - sstr[j].y_cord);
		 			gsl_vector_set(v,2,fstr[i].z_cord - sstr[j].z_cord);

		 			break;
		 		}
		 	}
		  break;
		 }
 	}
}

 void rotate_phi(gsl_matrix *rota,struct pdb_atom *all_init,int all,int node) {
	int j;
	double temp_x, temp_y, temp_z;
	for(j=0;j<all;++j) {
		//printf("%d I ROTATE NODE %d et String:%s...",j,all_init[j].node,all_init[j].atom_prot_type );
		if (node < all_init[j].node || (node == all_init[j].node && strcmp(all_init[j].atom_prot_type," O   ") == 0) ) {
			//printf("YES");
 			temp_x = all_init[j].x_cord;
 			temp_y = all_init[j].y_cord;
 			temp_z = all_init[j].z_cord;
 			all_init[j].x_cord = gsl_matrix_get(rota,0,0) * temp_x + gsl_matrix_get(rota,1,0) * temp_y + gsl_matrix_get(rota,2,0) * temp_z;
 			all_init[j].y_cord = gsl_matrix_get(rota,0,1) * temp_x + gsl_matrix_get(rota,1,1) * temp_y + gsl_matrix_get(rota,2,1) * temp_z;
 			all_init[j].z_cord = gsl_matrix_get(rota,0,2) * temp_x + gsl_matrix_get(rota,1,2) * temp_y + gsl_matrix_get(rota,2,2) * temp_z;
 			
 		}
 		//printf("\n");
 	}
}

void rotate_psy(gsl_matrix *rota,struct pdb_atom *all_init,int all,int node) {
	int j;
	double temp_x, temp_y, temp_z;
	for(j=0;j<all;++j) {
		//printf("J:%d Node:%d	Targ:%d\n",j,all_init[j].node,node);
		if (node > all_init[j].node ) {
			//printf("I ROTATE %s et node:%d\n",all_init[j].atom_prot_type,all_init[j].node);
 			temp_x = all_init[j].x_cord;
 			temp_y = all_init[j].y_cord;
 			temp_z = all_init[j].z_cord;
 			all_init[j].x_cord = gsl_matrix_get(rota,0,0) * temp_x + gsl_matrix_get(rota,1,0) * temp_y + gsl_matrix_get(rota,2,0) * temp_z;
 			all_init[j].y_cord = gsl_matrix_get(rota,0,1) * temp_x + gsl_matrix_get(rota,1,1) * temp_y + gsl_matrix_get(rota,2,1) * temp_z;
 			all_init[j].z_cord = gsl_matrix_get(rota,0,2) * temp_x + gsl_matrix_get(rota,1,2) * temp_y + gsl_matrix_get(rota,2,2) * temp_z;
 		//	printf("Node:%d	%f	%f	%f\n",all_init[j].node,all_init[j].x_cord,all_init[j].y_cord,all_init[j].z_cord );
 		}
 	}
}

void print_matrix(gsl_matrix *mat) {
	int i,j;
	for(i=0;i<mat->size1;++i) {
		for(j=0;j<mat->size2;++j) {
			printf("%10.5f ",gsl_matrix_get(mat,i,j));
		}
		printf("\n");
	}
}
 
 void rotate_bb(struct pdb_atom *init,struct pdb_atom *targ,int atom,int all,struct pdb_atom *all_init,struct pdb_atom *all_targ,int *align) {
 	
 	int mc_total = 10000;
 	int i,l,j,rand_node;
 	float old,newstrc;
 	float random;
 	float psy_phi;
 	long seed; 
 	gsl_vector *vect = gsl_vector_alloc(3);
 	gsl_matrix *rota = gsl_matrix_alloc(3,3);
 	seed = time_seed();
 	/*FILE *file;
 	file = fopen("rota.pdb","w");*/
 	
 	struct pdb_atom store_init[all];
 	struct pdb_atom temp_init[atom];
 	
 	for(i=0;i<atom;++i) {
 		temp_init[i].x_cord = init[i].x_cord;
		temp_init[i].y_cord = init[i].y_cord;
		temp_init[i].z_cord = init[i].z_cord;
		temp_init[i].res_number = init[i].res_number;
		temp_init[i].atom_number = init[i].atom_number;
		temp_init[i].node = init[i].node;
		temp_init[i].atom_type = init[i].atom_type;
		strcpy(temp_init[i].atom_prot_type,init[i].atom_prot_type);
		strcpy(temp_init[i].res_type,init[i].res_type);
		strcpy(temp_init[i].chain,init[i].chain);
 	}
 	
 	
 	for(i=0;i<all;++i) {
 		store_init[i].x_cord = all_init[i].x_cord;
		store_init[i].y_cord = all_init[i].y_cord;
		store_init[i].z_cord = all_init[i].z_cord;
		store_init[i].res_number = all_init[i].res_number;
		store_init[i].atom_number = all_init[i].atom_number;
		store_init[i].node = all_init[i].node;
		store_init[i].atom_type = all_init[i].atom_type;
		strcpy(store_init[i].atom_prot_type,all_init[i].atom_prot_type);
		strcpy(store_init[i].res_type,all_init[i].res_type);
		strcpy(store_init[i].chain,all_init[i].chain);
		//printf("I COPY:%s %s %d %f %f %f\n",store_init[i].atom_prot_type,store_init[i].res_type,store_init[i].res_number,store_init[i].x_cord,store_init[i].y_cord,store_init[i].z_cord);
 	}
 	
 
 	int m;
 	int model = 1;

 	//write_movie(file,store_init,all,model);
 	old = rmsd_no(temp_init,targ,atom,align);
 	printf("Init:%f	",old);
	random = -0.001;
 	for(m=0;m<mc_total;++m) {
	 	for(i=0;i<all;++i) {
		 	store_init[i].x_cord = all_init[i].x_cord;
			store_init[i].y_cord = all_init[i].y_cord;
			store_init[i].z_cord = all_init[i].z_cord;
	 	}
	 	
	 	for(l=0;l<atom;++l) {
			temp_init[l].x_cord = init[l].x_cord;
			temp_init[l].y_cord = init[l].y_cord;
			temp_init[l].z_cord = init[l].z_cord;
		}
	 	rand_node = ran2(&seed)*atom;
	 	random = (ran2(&seed)-0.5)/10;
	 	psy_phi = ran2(&seed);
	 	
	 	center_node(all,store_init,rand_node);
	 	center_node(atom,temp_init,rand_node);
	 	
	 	if (psy_phi < 0.5) {	 	
		 	assign_vector(all_init,rand_node," CA  ",all,all_init,rand_node," C   ",all,vect);
			rotate_matrix(rota,random,vect);
		 	rotate_phi(rota,store_init,all,rand_node);
		 	rotate_phi(rota,temp_init,atom,rand_node);
		} else {	 	
		 	assign_vector(all_init,rand_node," CA  ",all,all_init,rand_node," N   ",all,vect);
			rotate_matrix(rota,random,vect);
		 	rotate_psy(rota,store_init,all,rand_node);
		 	rotate_psy(rota,temp_init,atom,rand_node);
		 }
	 	
	 	newstrc = rmsd_no(temp_init,targ,atom,align);
	 	if (newstrc < old) {
	 		center_node(all,store_init,50);
	 		//printf("M:%d Node:%4d Rand:%8.5f	%f\n",m,rand_node,random,newstrc);
	 		++model;
	 		//write_movie(file,store_init,all,model);
	 		old = newstrc;
		 	for(i=0;i<all;++i) {
			 	all_init[i].x_cord = store_init[i].x_cord;
				all_init[i].y_cord = store_init[i].y_cord;
				all_init[i].z_cord = store_init[i].z_cord;
		 	}
		 	
		 	for(l=0;l<atom;++l) {
				init[l].x_cord = temp_init[l].x_cord;
				init[l].y_cord = temp_init[l].y_cord;
				init[l].z_cord = temp_init[l].z_cord;
			}	 		
	 	}
	 	
	 }
	 printf("Final:%f	",old);
 }
 
 
 
void center_strc(int all,struct pdb_atom *strc,gsl_vector *trans) {
	float x =0.0;
	float y= 0.0;
	float z= 0.0;
	float tot = 0.0;
	int i;
	
	for(i=0;i<all;++i) {
		x += strc[i].x_cord;
		y += strc[i].y_cord;
		z += strc[i].z_cord;
		tot += 1.0;
	}
	x /= tot;
	y /= tot;
	z /= tot;
	
	for(i=0;i<all;++i) {
		strc[i].x_cord -= x;
		strc[i].y_cord -= y;
		strc[i].z_cord -= z;
	}
	
	gsl_vector_set(trans,0,x);
	gsl_vector_set(trans,1,y);
	gsl_vector_set(trans,2,z);
	
}
