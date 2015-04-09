#include "STeM.h"
void build_2_matrix(struct pdb_atom *strc, double **hes,int nb_atom,double K_theta) {
 int i,j,k;
 double lpl,lql; //Distance de points
 double dot; // Dot product des vecteurs
 double G;
 double dGdXi,dGdYi,dGdZi;
 double dGdXj,dGdYj,dGdZj;
 double dGdXk,dGdYk,dGdZk;
 

 	for (i=1;i < nb_atom+1;++i) {
		for(j=i+1;j<nb_atom+1;++j) {
			for (k=j+1;k<nb_atom+1;++k) {
		 		
		 	int connect_flag = 0;
		 		int t;
		 		for (t=0;t<6;++t) {
		 			if (strc[i-1].node_c[t] == strc[j-1].node) {++connect_flag;}
		 			if (strc[j-1].node_c[t] == strc[k-1].node) {++connect_flag;}
		 			
		 		}
		 		
		 		if (connect_flag == 2) {} else {continue;}
		 		
		 		dot = (strc[i-1].x_cord-strc[j-1].x_cord)*(strc[k-1].x_cord-strc[j-1].x_cord)
		 			 +(strc[i-1].y_cord-strc[j-1].y_cord)*(strc[k-1].y_cord-strc[j-1].y_cord)
		 			 +(strc[i-1].z_cord-strc[j-1].z_cord)*(strc[k-1].z_cord-strc[j-1].z_cord);
		 		
		 		lpl = sqrt((strc[i-1].x_cord-strc[j-1].x_cord)*(strc[i-1].x_cord-strc[j-1].x_cord)+
		 	 		       (strc[i-1].y_cord-strc[j-1].y_cord)*(strc[i-1].y_cord-strc[j-1].y_cord)+
		 	 		       (strc[i-1].z_cord-strc[j-1].z_cord)*(strc[i-1].z_cord-strc[j-1].z_cord));//Distance I et J
		 	 		
		 	 	lql = sqrt((strc[k-1].x_cord-strc[j-1].x_cord)*(strc[k-1].x_cord-strc[j-1].x_cord)+
		 	 		       (strc[k-1].y_cord-strc[j-1].y_cord)*(strc[k-1].y_cord-strc[j-1].y_cord)+
		 	 		       (strc[k-1].z_cord-strc[j-1].z_cord)*(strc[k-1].z_cord-strc[j-1].z_cord)); // Distance J et K
		 	 		       
		 	 	G=dot/(lpl*lql);

		 	 	if ( (1-G*G) == 0 ) {printf("Valeur de G invalide (angle de 180 degrées)\nIgnore terme de l'acide aminée %d %d %d\n",i,j,k);continue;} 
		 	 	 	 	
		 	 	
		 	 	
		 	 	dGdXi=((strc[k-1].x_cord-strc[j-1].x_cord)*lpl*lql-dot*(lql/lpl)*(strc[i-1].x_cord-strc[j-1].x_cord))/((lpl*lql)*(lpl*lql));
		 	 	dGdYi=((strc[k-1].y_cord-strc[j-1].y_cord)*lpl*lql-dot*(lql/lpl)*(strc[i-1].y_cord-strc[j-1].y_cord))/((lpl*lql)*(lpl*lql));
		 	 	dGdZi=((strc[k-1].z_cord-strc[j-1].z_cord)*lpl*lql-dot*(lql/lpl)*(strc[i-1].z_cord-strc[j-1].z_cord))/((lpl*lql)*(lpl*lql));
		 	 	
		 	 	
		 	 	dGdXj=((2*strc[j-1].x_cord-strc[i-1].x_cord-strc[k-1].x_cord)*lpl*lql-dot*(lql/lpl)*(strc[j-1].x_cord-strc[i-1].x_cord)-dot*(lpl/lql)*(strc[j-1].x_cord-strc[k-1].x_cord))/((lpl*lql)*(lpl*lql));
		 	 	dGdYj=((2*strc[j-1].y_cord-strc[i-1].y_cord-strc[k-1].y_cord)*lpl*lql-dot*(lql/lpl)*(strc[j-1].y_cord-strc[i-1].y_cord)-dot*(lpl/lql)*(strc[j-1].y_cord-strc[k-1].y_cord))/((lpl*lql)*(lpl*lql));
		 	 	dGdZj=((2*strc[j-1].z_cord-strc[i-1].z_cord-strc[k-1].z_cord)*lpl*lql-dot*(lql/lpl)*(strc[j-1].z_cord-strc[i-1].z_cord)-dot*(lpl/lql)*(strc[j-1].z_cord-strc[k-1].z_cord))/((lpl*lql)*(lpl*lql));
		 		
		 		
		 		 dGdXk=((strc[i-1].x_cord-strc[j-1].x_cord)*lpl*lql-dot*(lpl/lql)*(strc[k-1].x_cord-strc[j-1].x_cord))/((lpl*lql)*(lpl*lql));
		 		 dGdYk=((strc[i-1].y_cord-strc[j-1].y_cord)*lpl*lql-dot*(lpl/lql)*(strc[k-1].y_cord-strc[j-1].y_cord))/((lpl*lql)*(lpl*lql));
		 		 dGdZk=((strc[i-1].z_cord-strc[j-1].z_cord)*lpl*lql-dot*(lpl/lql)*(strc[k-1].z_cord-strc[j-1].z_cord))/((lpl*lql)*(lpl*lql));
		 		 
				hes[3*i-3][3*j-3]        += 2*K_theta/(1-G*G)*dGdXi*dGdXj;
				hes[3*i-2][3*j-2]        += 2*K_theta/(1-G*G)*dGdYi*dGdYj;
				hes[3*i-1][3*j-1]        += 2*K_theta/(1-G*G)*dGdZi*dGdZj;
				hes[3*i-3][3*j-2]        += 2*K_theta/(1-G*G)*dGdXi*dGdYj;
				hes[3*i-3][3*j-1]        += 2*K_theta/(1-G*G)*dGdXi*dGdZj;
				hes[3*i-2][3*j-3]        += 2*K_theta/(1-G*G)*dGdYi*dGdXj;
				hes[3*i-2][3*j-1]        += 2*K_theta/(1-G*G)*dGdYi*dGdZj;
				hes[3*i-1][3*j-3]        += 2*K_theta/(1-G*G)*dGdZi*dGdXj;
				hes[3*i-1][3*j-2]        += 2*K_theta/(1-G*G)*dGdZi*dGdYj;
				hes[3*j-3][3*i-3]        += 2*K_theta/(1-G*G)*dGdXj*dGdXi;
				hes[3*j-2][3*i-2]        += 2*K_theta/(1-G*G)*dGdYj*dGdYi;
				hes[3*j-1][3*i-1]        += 2*K_theta/(1-G*G)*dGdZj*dGdZi;
				hes[3*j-3][3*i-2]        += 2*K_theta/(1-G*G)*dGdXj*dGdYi;
				hes[3*j-3][3*i-1]        += 2*K_theta/(1-G*G)*dGdXj*dGdZi;
				hes[3*j-2][3*i-3]        += 2*K_theta/(1-G*G)*dGdYj*dGdXi;
				hes[3*j-2][3*i-1]        += 2*K_theta/(1-G*G)*dGdYj*dGdZi;
				hes[3*j-1][3*i-3]        += 2*K_theta/(1-G*G)*dGdZj*dGdXi;
				hes[3*j-1][3*i-2]        += 2*K_theta/(1-G*G)*dGdZj*dGdYi;
				hes[3*j-3][3*k-3]        += 2*K_theta/(1-G*G)*dGdXj*dGdXk;
				hes[3*j-2][3*k-2]        += 2*K_theta/(1-G*G)*dGdYj*dGdYk;
				hes[3*j-1][3*k-1]        += 2*K_theta/(1-G*G)*dGdZj*dGdZk;
				hes[3*j-3][3*k-2]        += 2*K_theta/(1-G*G)*dGdXj*dGdYk;
				hes[3*j-3][3*k-1]        += 2*K_theta/(1-G*G)*dGdXj*dGdZk;
				hes[3*j-2][3*k-3]        += 2*K_theta/(1-G*G)*dGdYj*dGdXk;
				hes[3*j-2][3*k-1]        += 2*K_theta/(1-G*G)*dGdYj*dGdZk;
				hes[3*j-1][3*k-3]        += 2*K_theta/(1-G*G)*dGdZj*dGdXk;
				hes[3*j-1][3*k-2]        += 2*K_theta/(1-G*G)*dGdZj*dGdYk;
				hes[3*k-3][3*j-3]        += 2*K_theta/(1-G*G)*dGdXk*dGdXj;
				hes[3*k-2][3*j-2]        += 2*K_theta/(1-G*G)*dGdYk*dGdYj;
				hes[3*k-1][3*j-1]        += 2*K_theta/(1-G*G)*dGdZk*dGdZj;
				hes[3*k-3][3*j-2]        += 2*K_theta/(1-G*G)*dGdXk*dGdYj;
				hes[3*k-3][3*j-1]        += 2*K_theta/(1-G*G)*dGdXk*dGdZj;
				hes[3*k-2][3*j-3]        += 2*K_theta/(1-G*G)*dGdYk*dGdXj;
				hes[3*k-2][3*j-1]        += 2*K_theta/(1-G*G)*dGdYk*dGdZj;
				hes[3*k-1][3*j-3]        += 2*K_theta/(1-G*G)*dGdZk*dGdXj;
				hes[3*k-1][3*j-2]        += 2*K_theta/(1-G*G)*dGdZk*dGdYj;
				hes[3*i-3][3*k-3]        += 2*K_theta/(1-G*G)*dGdXi*dGdXk;
				hes[3*i-2][3*k-2]        += 2*K_theta/(1-G*G)*dGdYi*dGdYk;
				hes[3*i-1][3*k-1]        += 2*K_theta/(1-G*G)*dGdZi*dGdZk;
				hes[3*i-3][3*k-2]        += 2*K_theta/(1-G*G)*dGdXi*dGdYk;
				hes[3*i-3][3*k-1]        += 2*K_theta/(1-G*G)*dGdXi*dGdZk;
				hes[3*i-2][3*k-3]        += 2*K_theta/(1-G*G)*dGdYi*dGdXk;
				hes[3*i-2][3*k-1]        += 2*K_theta/(1-G*G)*dGdYi*dGdZk;
				hes[3*i-1][3*k-3]        += 2*K_theta/(1-G*G)*dGdZi*dGdXk;
				hes[3*i-1][3*k-2]        += 2*K_theta/(1-G*G)*dGdZi*dGdYk;
				hes[3*k-3][3*i-3]        += 2*K_theta/(1-G*G)*dGdXk*dGdXi;
				hes[3*k-2][3*i-2]        += 2*K_theta/(1-G*G)*dGdYk*dGdYi;
				hes[3*k-1][3*i-1]        += 2*K_theta/(1-G*G)*dGdZk*dGdZi;
				hes[3*k-3][3*i-2]        += 2*K_theta/(1-G*G)*dGdXk*dGdYi;
				hes[3*k-3][3*i-1]        += 2*K_theta/(1-G*G)*dGdXk*dGdZi;
				hes[3*k-2][3*i-3]        += 2*K_theta/(1-G*G)*dGdYk*dGdXi;
				hes[3*k-2][3*i-1]        += 2*K_theta/(1-G*G)*dGdYk*dGdZi;
				hes[3*k-1][3*i-3]        += 2*K_theta/(1-G*G)*dGdZk*dGdXi;
				hes[3*k-1][3*i-2]        += 2*K_theta/(1-G*G)*dGdZk*dGdYi;
				hes[3*i-3][3*i-3]        += 2*K_theta/(1-G*G)*dGdXi*dGdXi;
				hes[3*i-2][3*i-2]        += 2*K_theta/(1-G*G)*dGdYi*dGdYi;
				hes[3*i-1][3*i-1]        += 2*K_theta/(1-G*G)*dGdZi*dGdZi;                        
				hes[3*i-3][3*i-2]        += 2*K_theta/(1-G*G)*dGdXi*dGdYi;
				hes[3*i-3][3*i-1]        += 2*K_theta/(1-G*G)*dGdXi*dGdZi;
				hes[3*i-2][3*i-3]        += 2*K_theta/(1-G*G)*dGdYi*dGdXi;
				hes[3*i-2][3*i-1]        += 2*K_theta/(1-G*G)*dGdYi*dGdZi;
				hes[3*i-1][3*i-3]        += 2*K_theta/(1-G*G)*dGdZi*dGdXi;
				hes[3*i-1][3*i-2]        += 2*K_theta/(1-G*G)*dGdZi*dGdYi;                        
				hes[3*j-3][3*j-3]        += 2*K_theta/(1-G*G)*dGdXj*dGdXj;
				hes[3*j-2][3*j-2]        += 2*K_theta/(1-G*G)*dGdYj*dGdYj;
				hes[3*j-1][3*j-1]        += 2*K_theta/(1-G*G)*dGdZj*dGdZj;                        
				hes[3*j-3][3*j-2]        += 2*K_theta/(1-G*G)*dGdXj*dGdYj;
				hes[3*j-3][3*j-1]        += 2*K_theta/(1-G*G)*dGdXj*dGdZj;
				hes[3*j-2][3*j-3]        += 2*K_theta/(1-G*G)*dGdYj*dGdXj;
				hes[3*j-2][3*j-1]        += 2*K_theta/(1-G*G)*dGdYj*dGdZj;
				hes[3*j-1][3*j-3]        += 2*K_theta/(1-G*G)*dGdZj*dGdXj;
				hes[3*j-1][3*j-2]        += 2*K_theta/(1-G*G)*dGdZj*dGdYj;
				hes[3*k-3][3*k-3]        += 2*K_theta/(1-G*G)*dGdXk*dGdXk;
				hes[3*k-2][3*k-2]        += 2*K_theta/(1-G*G)*dGdYk*dGdYk;
				hes[3*k-1][3*k-1]        += 2*K_theta/(1-G*G)*dGdZk*dGdZk;                        
				hes[3*k-3][3*k-2]        += 2*K_theta/(1-G*G)*dGdXk*dGdYk;
				hes[3*k-3][3*k-1]        += 2*K_theta/(1-G*G)*dGdXk*dGdZk;
				hes[3*k-2][3*k-3]        += 2*K_theta/(1-G*G)*dGdYk*dGdXk;
				hes[3*k-2][3*k-1]        += 2*K_theta/(1-G*G)*dGdYk*dGdZk;
				hes[3*k-1][3*k-3]        += 2*K_theta/(1-G*G)*dGdZk*dGdXk;
				hes[3*k-1][3*k-2]        += 2*K_theta/(1-G*G)*dGdZk*dGdYk;
		 		
			}
		}
 	} 
 }
 
 
void build_3_matrix(struct pdb_atom *strc, double **hes,int nb_atom,double K_phi, float factor) {
 
	int i,j,k,l;
	double Xi,Yi,Zi,Xj,Yj,Zj,Xk,Yk,Zk,Xl,Yl,Zl;
	double v1[3],v2[3];
	double G,lv1l,lv2l;

	float K_phi_temp = K_phi * factor;
	int connect_flag = 0;	
	int t;
	
	// Recherche 4 atoms reliées qui forme des angles diedres
	
	
	for (i=1;i < nb_atom+1;++i) {
		//printf("I:%d\n",i-1);
		for(j=i+1;j<nb_atom+1;++j) {
			connect_flag = 0;
			for (t=0;t<6;++t) {
				if (strc[i-1].node_c[t] == strc[j-1].node) {++connect_flag;}
			}
			if (connect_flag == 0) {continue;}
			//printf("\tJ:%d\n",j-1);
			for (k=j+1;k<nb_atom+1;++k) {
				connect_flag = 0;
				for (t=0;t<6;++t) {
					if (strc[j-1].node_c[t] == strc[k-1].node) {++connect_flag;}
				}
				if (connect_flag == 0) {continue;}
				for(l=k+1;l<nb_atom+1;++l) {
					//printf("%d %d %d %d\n",i-1,j-1,k-1,l-1);
					connect_flag = 0;
			 		
			 		for (t=0;t<6;++t) {
			 			if (strc[i-1].node_c[t] == strc[j-1].node) {++connect_flag;}
			 			if (strc[j-1].node_c[t] == strc[k-1].node) {++connect_flag;}
			 			if (strc[k-1].node_c[t] == strc[l-1].node) {++connect_flag;}
			 		}
			 		
			 		if (connect_flag == 3) {} else {continue;}
					
					//printf("	%d %d %d %d\n",i-1,j-1,k-1,l-1);
					Xi=strc[i-1].x_cord;
					Yi=strc[i-1].y_cord;
					Zi=strc[i-1].z_cord;
					Xj=strc[j-1].x_cord;
					Yj=strc[j-1].y_cord;
					Zj=strc[j-1].z_cord;
					Xk=strc[k-1].x_cord;
					Yk=strc[k-1].y_cord;
					Zk=strc[k-1].z_cord;
					Xl=strc[l-1].x_cord;
					Yl=strc[l-1].y_cord;
					Zl=strc[l-1].z_cord;
		
					//printf("%d\ti:%s\tj:%s\tk:%s\tl:%s\n",m,strc[i-1].atom_prot_type,strc[j-1].atom_prot_type,strc[k-1].atom_prot_type,strc[l-1].atom_prot_type);
					//if ((strncmp(strc[j-1].atom_prot_type," C ",3) == 0) || (strncmp(strc[k-1].atom_prot_type," N ",3) == 0) ) {K_phi_temp = K_phi *2;} else { K_phi_temp = K_phi/factor;}
					 //printf("K_phi = %f\n",K_phi_temp);          

					v1[0] = (Yj-Yi)*(Zk-Zj)-(Zj-Zi)*(Yk-Yj);
					v1[1] = (Zj-Zi)*(Xk-Xj)-(Xj-Xi)*(Zk-Zj);
					v1[2] = (Xj-Xi)*(Yk-Yj)-(Yj-Yi)*(Xk-Xj);
				
					v2[0] = (Yk-Yj)*(Zl-Zk)-(Zk-Zj)*(Yl-Yk);
					v2[1] = (Zk-Zj)*(Xl-Xk)-(Xk-Xj)*(Zl-Zk);
					v2[2] = (Xk-Xj)*(Yl-Yk)-(Yk-Yj)*(Xl-Xk);

					lv1l=sqrt(v1[0]*v1[0]+v1[1]*v1[1]+v1[2]*v1[2]);
					lv2l=sqrt(v2[0]*v2[0]+v2[1]*v2[1]+v2[2]*v2[2]);
				
					G=dot_p(v1,v2)/(lv1l*lv2l);
		
					double dv1dXi[3] ={0,Zk-Zj,Yj-Yk};
					double dv1dYi[3] ={Zj-Zk,0,Xk-Xj};
					double dv1dZi[3] ={Yk-Yj,Xj-Xk,0};
					double dv1dXj[3] ={0,Zi-Zk,Yk-Yi};
					double dv1dYj[3] ={Zk-Zi,0,Xi-Xk};
					double dv1dZj[3] ={Yi-Yk,Xk-Xi,0};       
					double dv1dXk[3] ={0,Zj-Zi,Yi-Yj};
					double dv1dYk[3] ={Zi-Zj,0,Xj-Xi};
					double dv1dZk[3] ={Yj-Yi,Xi-Xj,0};
					double dv1dXl[3] ={0,0,0};
					double dv1dYl[3] ={0,0,0};
					double dv1dZl[3] ={0,0,0};
					double dv2dXi[3] ={0,0,0};
					double dv2dYi[3] ={0,0,0};
					double dv2dZi[3] ={0,0,0};        
					double dv2dXj[3] ={0,Zl-Zk,Yk-Yl};
					double dv2dYj[3] ={Zk-Zl,0,Xl-Xk};
					double dv2dZj[3] ={Yl-Yk,Xk-Xl,0};
					double dv2dXk[3] ={0,Zj-Zl,Yl-Yj};
					double dv2dYk[3] ={Zl-Zj,0,Xj-Xl};
					double dv2dZk[3] ={Yj-Yl,Xl-Xj,0};
					double dv2dXl[3] ={0,Zk-Zj,Yj-Yk};
					double dv2dYl[3] ={Zj-Zk,0,Xk-Xj};
					double dv2dZl[3] ={Yk-Yj,Xj-Xk,0};
		
					double K1=(Yj-Yi)*(Zk-Zj)-(Yk-Yj)*(Zj-Zi);
					double K2=(Xk-Xj)*(Zj-Zi)-(Xj-Xi)*(Zk-Zj);
					double K3=(Xj-Xi)*(Yk-Yj)-(Xk-Xj)*(Yj-Yi);
				
				
					double dlv1ldXi=(2*K2*(Zk-Zj)+2*K3*(Yj-Yk))/(2*sqrt(K1*K1+K2*K2+K3*K3));
					double dlv1ldYi=(2*K1*(Zj-Zk)+2*K3*(Xk-Xj))/(2*sqrt(K1*K1+K2*K2+K3*K3));
					double dlv1ldZi=(2*K1*(Yk-Yj)+2*K2*(Xj-Xk))/(2*sqrt(K1*K1+K2*K2+K3*K3));               
				   
					double dlv1ldXj=(2*K2*(Zi-Zk)+2*K3*(Yk-Yi))/(2*sqrt(K1*K1+K2*K2+K3*K3));
					double dlv1ldYj=(2*K1*(Zk-Zi)+2*K3*(Xi-Xk))/(2*sqrt(K1*K1+K2*K2+K3*K3));
					double dlv1ldZj=(2*K1*(Yi-Yk)+2*K2*(Xk-Xi))/(2*sqrt(K1*K1+K2*K2+K3*K3));        
				   
					double dlv1ldXk=(2*K2*(Zj-Zi)+2*K3*(Yi-Yj))/(2*sqrt(K1*K1+K2*K2+K3*K3));     
					double dlv1ldYk=(2*K1*(Zi-Zj)+2*K3*(Xj-Xi))/(2*sqrt(K1*K1+K2*K2+K3*K3));   
					double dlv1ldZk=(2*K1*(Yj-Yi)+2*K2*(Xi-Xj))/(2*sqrt(K1*K1+K2*K2+K3*K3));
						  
					double dlv1ldXl=0; 
					double dlv1ldYl=0;    
					double dlv1ldZl=0;
				
					double dlv2ldXi=0;
					double dlv2ldYi=0;
					double dlv2ldZi=0;

					double L1=(Yk-Yj)*(Zl-Zk)-(Yl-Yk)*(Zk-Zj);
					double L2=(Xl-Xk)*(Zk-Zj)-(Xk-Xj)*(Zl-Zk);
					double L3=(Xk-Xj)*(Yl-Yk)-(Xl-Xk)*(Yk-Yj);
				
				
					double dlv2ldXj=(2*L2*(Zl-Zk)+2*L3*(Yk-Yl))/(2*sqrt(L1*L1+L2*L2+L3*L3));
					double dlv2ldYj=(2*L1*(Zk-Zl)+2*L3*(Xl-Xk))/(2*sqrt(L1*L1+L2*L2+L3*L3));
					double dlv2ldZj=(2*L1*(Yl-Yk)+2*L2*(Xk-Xl))/(2*sqrt(L1*L1+L2*L2+L3*L3));

					double dlv2ldXk=(2*L2*(Zj-Zl)+2*L3*(Yl-Yj))/(2*sqrt(L1*L1+L2*L2+L3*L3));
					double dlv2ldYk=(2*L1*(Zl-Zj)+2*L3*(Xj-Xl))/(2*sqrt(L1*L1+L2*L2+L3*L3));
					double dlv2ldZk=(2*L1*(Yj-Yl)+2*L2*(Xl-Xj))/(2*sqrt(L1*L1+L2*L2+L3*L3));

					double dlv2ldXl=(2*L2*(Zk-Zj)+2*L3*(Yj-Yk))/(2*sqrt(L1*L1+L2*L2+L3*L3));
					double dlv2ldYl=(2*L1*(Zj-Zk)+2*L3*(Xk-Xj))/(2*sqrt(L1*L1+L2*L2+L3*L3));
					double dlv2ldZl=(2*L1*(Yk-Yj)+2*L2*(Xj-Xk))/(2*sqrt(L1*L1+L2*L2+L3*L3));
		
					double dGdXi=((dot_p(dv1dXi,v2)+dot_p(dv2dXi,v1))*lv1l*lv2l-dot_p(v1,v2)*(dlv1ldXi*lv2l+dlv2ldXi*lv1l))/((lv1l*lv2l)*(lv1l*lv2l));
					double dGdYi=((dot_p(dv1dYi,v2)+dot_p(dv2dYi,v1))*lv1l*lv2l-dot_p(v1,v2)*(dlv1ldYi*lv2l+dlv2ldYi*lv1l))/((lv1l*lv2l)*(lv1l*lv2l));
					double dGdZi=((dot_p(dv1dZi,v2)+dot_p(dv2dZi,v1))*lv1l*lv2l-dot_p(v1,v2)*(dlv1ldZi*lv2l+dlv2ldZi*lv1l))/((lv1l*lv2l)*(lv1l*lv2l));

					double dGdXj=((dot_p(dv1dXj,v2)+dot_p(dv2dXj,v1))*lv1l*lv2l-dot_p(v1,v2)*(dlv1ldXj*lv2l+dlv2ldXj*lv1l))/((lv1l*lv2l)*(lv1l*lv2l));
					double dGdYj=((dot_p(dv1dYj,v2)+dot_p(dv2dYj,v1))*lv1l*lv2l-dot_p(v1,v2)*(dlv1ldYj*lv2l+dlv2ldYj*lv1l))/((lv1l*lv2l)*(lv1l*lv2l));
					double dGdZj=((dot_p(dv1dZj,v2)+dot_p(dv2dZj,v1))*lv1l*lv2l-dot_p(v1,v2)*(dlv1ldZj*lv2l+dlv2ldZj*lv1l))/((lv1l*lv2l)*(lv1l*lv2l));

					double dGdXk=((dot_p(dv1dXk,v2)+dot_p(dv2dXk,v1))*lv1l*lv2l-dot_p(v1,v2)*(dlv1ldXk*lv2l+dlv2ldXk*lv1l))/((lv1l*lv2l)*(lv1l*lv2l));
					double dGdYk=((dot_p(dv1dYk,v2)+dot_p(dv2dYk,v1))*lv1l*lv2l-dot_p(v1,v2)*(dlv1ldYk*lv2l+dlv2ldYk*lv1l))/((lv1l*lv2l)*(lv1l*lv2l));
					double dGdZk=((dot_p(dv1dZk,v2)+dot_p(dv2dZk,v1))*lv1l*lv2l-dot_p(v1,v2)*(dlv1ldZk*lv2l+dlv2ldZk*lv1l))/((lv1l*lv2l)*(lv1l*lv2l));

					double dGdXl=((dot_p(dv1dXl,v2)+dot_p(dv2dXl,v1))*lv1l*lv2l-dot_p(v1,v2)*(dlv1ldXl*lv2l+dlv2ldXl*lv1l))/((lv1l*lv2l)*(lv1l*lv2l));
					double dGdYl=((dot_p(dv1dYl,v2)+dot_p(dv2dYl,v1))*lv1l*lv2l-dot_p(v1,v2)*(dlv1ldYl*lv2l+dlv2ldYl*lv1l))/((lv1l*lv2l)*(lv1l*lv2l));
					double dGdZl=((dot_p(dv1dZl,v2)+dot_p(dv2dZl,v1))*lv1l*lv2l-dot_p(v1,v2)*(dlv1ldZl*lv2l+dlv2ldZl*lv1l))/((lv1l*lv2l)*(lv1l*lv2l));
		
					hes[3*i-3][3*j-3]     +=(2*K_phi_temp)/(1-G*G)*dGdXi*dGdXj;
					hes[3*i-2][3*j-2]     +=(2*K_phi_temp)/(1-G*G)*dGdYi*dGdYj;
					hes[3*i-1][3*j-1]     +=(2*K_phi_temp)/(1-G*G)*dGdZi*dGdZj;
					hes[3*i-3][3*j-2]     +=(2*K_phi_temp)/(1-G*G)*dGdXi*dGdYj;
					hes[3*i-3][3*j-1]     +=(2*K_phi_temp)/(1-G*G)*dGdXi*dGdZj;
					hes[3*i-2][3*j-3]     +=(2*K_phi_temp)/(1-G*G)*dGdYi*dGdXj;
					hes[3*i-2][3*j-1]     +=(2*K_phi_temp)/(1-G*G)*dGdYi*dGdZj;
					hes[3*i-1][3*j-3]     +=(2*K_phi_temp)/(1-G*G)*dGdZi*dGdXj;
					hes[3*i-1][3*j-2]     +=(2*K_phi_temp)/(1-G*G)*dGdZi*dGdYj;
					hes[3*j-3][3*i-3]     +=(2*K_phi_temp)/(1-G*G)*dGdXj*dGdXi;
					hes[3*j-2][3*i-2]     +=(2*K_phi_temp)/(1-G*G)*dGdYj*dGdYi;
					hes[3*j-1][3*i-1]     +=(2*K_phi_temp)/(1-G*G)*dGdZj*dGdZi;
					hes[3*j-3][3*i-2]     +=(2*K_phi_temp)/(1-G*G)*dGdXj*dGdYi;
					hes[3*j-3][3*i-1]     +=(2*K_phi_temp)/(1-G*G)*dGdXj*dGdZi;
					hes[3*j-2][3*i-3]     +=(2*K_phi_temp)/(1-G*G)*dGdYj*dGdXi;
					hes[3*j-2][3*i-1]     +=(2*K_phi_temp)/(1-G*G)*dGdYj*dGdZi;
					hes[3*j-1][3*i-3]     +=(2*K_phi_temp)/(1-G*G)*dGdZj*dGdXi;
					hes[3*j-1][3*i-2]     +=(2*K_phi_temp)/(1-G*G)*dGdZj*dGdYi;
					hes[3*i-3][3*l-3]     +=(2*K_phi_temp)/(1-G*G)*dGdXi*dGdXl;
					hes[3*i-2][3*l-2]     +=(2*K_phi_temp)/(1-G*G)*dGdYi*dGdYl;
					hes[3*i-1][3*l-1]     +=(2*K_phi_temp)/(1-G*G)*dGdZi*dGdZl;
					hes[3*i-3][3*l-2]     +=(2*K_phi_temp)/(1-G*G)*dGdXi*dGdYl;
					hes[3*i-3][3*l-1]     +=(2*K_phi_temp)/(1-G*G)*dGdXi*dGdZl;
					hes[3*i-2][3*l-3]     +=(2*K_phi_temp)/(1-G*G)*dGdYi*dGdXl;
					hes[3*i-2][3*l-1]     +=(2*K_phi_temp)/(1-G*G)*dGdYi*dGdZl;
					hes[3*i-1][3*l-3]     +=(2*K_phi_temp)/(1-G*G)*dGdZi*dGdXl;
					hes[3*i-1][3*l-2]     +=(2*K_phi_temp)/(1-G*G)*dGdZi*dGdYl;
					hes[3*l-3][3*i-3]     +=(2*K_phi_temp)/(1-G*G)*dGdXl*dGdXi;
					hes[3*l-2][3*i-2]     +=(2*K_phi_temp)/(1-G*G)*dGdYl*dGdYi;
					hes[3*l-1][3*i-1]     +=(2*K_phi_temp)/(1-G*G)*dGdZl*dGdZi;
					hes[3*l-3][3*i-2]     +=(2*K_phi_temp)/(1-G*G)*dGdXl*dGdYi;
					hes[3*l-3][3*i-1]     +=(2*K_phi_temp)/(1-G*G)*dGdXl*dGdZi;
					hes[3*l-2][3*i-3]     +=(2*K_phi_temp)/(1-G*G)*dGdYl*dGdXi;
					hes[3*l-2][3*i-1]     +=(2*K_phi_temp)/(1-G*G)*dGdYl*dGdZi;
					hes[3*l-1][3*i-3]     +=(2*K_phi_temp)/(1-G*G)*dGdZl*dGdXi;
					hes[3*l-1][3*i-2]     +=(2*K_phi_temp)/(1-G*G)*dGdZl*dGdYi;
					hes[3*k-3][3*j-3]     +=(2*K_phi_temp)/(1-G*G)*dGdXk*dGdXj;
					hes[3*k-2][3*j-2]     +=(2*K_phi_temp)/(1-G*G)*dGdYk*dGdYj;
					hes[3*k-1][3*j-1]     +=(2*K_phi_temp)/(1-G*G)*dGdZk*dGdZj;
					hes[3*k-3][3*j-2]     +=(2*K_phi_temp)/(1-G*G)*dGdXk*dGdYj;
					hes[3*k-3][3*j-1]     +=(2*K_phi_temp)/(1-G*G)*dGdXk*dGdZj;
					hes[3*k-2][3*j-3]     +=(2*K_phi_temp)/(1-G*G)*dGdYk*dGdXj;
					hes[3*k-2][3*j-1]     +=(2*K_phi_temp)/(1-G*G)*dGdYk*dGdZj;
					hes[3*k-1][3*j-3]     +=(2*K_phi_temp)/(1-G*G)*dGdZk*dGdXj;
					hes[3*k-1][3*j-2]     +=(2*K_phi_temp)/(1-G*G)*dGdZk*dGdYj;
					hes[3*j-3][3*k-3]     +=(2*K_phi_temp)/(1-G*G)*dGdXj*dGdXk;
					hes[3*j-2][3*k-2]     +=(2*K_phi_temp)/(1-G*G)*dGdYj*dGdYk;
					hes[3*j-1][3*k-1]     +=(2*K_phi_temp)/(1-G*G)*dGdZj*dGdZk;
					hes[3*j-3][3*k-2]     +=(2*K_phi_temp)/(1-G*G)*dGdXj*dGdYk;
					hes[3*j-3][3*k-1]     +=(2*K_phi_temp)/(1-G*G)*dGdXj*dGdZk;
					hes[3*j-2][3*k-3]     +=(2*K_phi_temp)/(1-G*G)*dGdYj*dGdXk;
					hes[3*j-2][3*k-1]     +=(2*K_phi_temp)/(1-G*G)*dGdYj*dGdZk;
					hes[3*j-1][3*k-3]     +=(2*K_phi_temp)/(1-G*G)*dGdZj*dGdXk;
					hes[3*j-1][3*k-2]     +=(2*K_phi_temp)/(1-G*G)*dGdZj*dGdYk;
					hes[3*i-3][3*k-3]     +=(2*K_phi_temp)/(1-G*G)*dGdXi*dGdXk;
					hes[3*i-2][3*k-2]     +=(2*K_phi_temp)/(1-G*G)*dGdYi*dGdYk;
					hes[3*i-1][3*k-1]     +=(2*K_phi_temp)/(1-G*G)*dGdZi*dGdZk;
					hes[3*i-3][3*k-2]     +=(2*K_phi_temp)/(1-G*G)*dGdXi*dGdYk;
					hes[3*i-3][3*k-1]     +=(2*K_phi_temp)/(1-G*G)*dGdXi*dGdZk;
					hes[3*i-2][3*k-3]     +=(2*K_phi_temp)/(1-G*G)*dGdYi*dGdXk;
					hes[3*i-2][3*k-1]     +=(2*K_phi_temp)/(1-G*G)*dGdYi*dGdZk;
					hes[3*i-1][3*k-3]     +=(2*K_phi_temp)/(1-G*G)*dGdZi*dGdXk;
					hes[3*i-1][3*k-2]     +=(2*K_phi_temp)/(1-G*G)*dGdZi*dGdYk;
					hes[3*k-3][3*i-3]     +=(2*K_phi_temp)/(1-G*G)*dGdXk*dGdXi;
					hes[3*k-2][3*i-2]     +=(2*K_phi_temp)/(1-G*G)*dGdYk*dGdYi;
					hes[3*k-1][3*i-1]     +=(2*K_phi_temp)/(1-G*G)*dGdZk*dGdZi;
					hes[3*k-3][3*i-2]     +=(2*K_phi_temp)/(1-G*G)*dGdXk*dGdYi;
					hes[3*k-3][3*i-1]     +=(2*K_phi_temp)/(1-G*G)*dGdXk*dGdZi;
					hes[3*k-2][3*i-3]     +=(2*K_phi_temp)/(1-G*G)*dGdYk*dGdXi;
					hes[3*k-2][3*i-1]     +=(2*K_phi_temp)/(1-G*G)*dGdYk*dGdZi;
					hes[3*k-1][3*i-3]     +=(2*K_phi_temp)/(1-G*G)*dGdZk*dGdXi;
					hes[3*k-1][3*i-2]     +=(2*K_phi_temp)/(1-G*G)*dGdZk*dGdYi;
					hes[3*l-3][3*j-3]     +=(2*K_phi_temp)/(1-G*G)*dGdXl*dGdXj;
					hes[3*l-2][3*j-2]     +=(2*K_phi_temp)/(1-G*G)*dGdYl*dGdYj;
					hes[3*l-1][3*j-1]     +=(2*K_phi_temp)/(1-G*G)*dGdZl*dGdZj;
					hes[3*l-3][3*j-2]     +=(2*K_phi_temp)/(1-G*G)*dGdXl*dGdYj;
					hes[3*l-3][3*j-1]     +=(2*K_phi_temp)/(1-G*G)*dGdXl*dGdZj;
					hes[3*l-2][3*j-3]     +=(2*K_phi_temp)/(1-G*G)*dGdYl*dGdXj;
					hes[3*l-2][3*j-1]     +=(2*K_phi_temp)/(1-G*G)*dGdYl*dGdZj;
					hes[3*l-1][3*j-3]     +=(2*K_phi_temp)/(1-G*G)*dGdZl*dGdXj;
					hes[3*l-1][3*j-2]     +=(2*K_phi_temp)/(1-G*G)*dGdZl*dGdYj;
					hes[3*j-3][3*l-3]     +=(2*K_phi_temp)/(1-G*G)*dGdXj*dGdXl;
					hes[3*j-2][3*l-2]     +=(2*K_phi_temp)/(1-G*G)*dGdYj*dGdYl;
					hes[3*j-1][3*l-1]     +=(2*K_phi_temp)/(1-G*G)*dGdZj*dGdZl;
					hes[3*j-3][3*l-2]     +=(2*K_phi_temp)/(1-G*G)*dGdXj*dGdYl;
					hes[3*j-3][3*l-1]     +=(2*K_phi_temp)/(1-G*G)*dGdXj*dGdZl;
					hes[3*j-2][3*l-3]     +=(2*K_phi_temp)/(1-G*G)*dGdYj*dGdXl;
					hes[3*j-2][3*l-1]     +=(2*K_phi_temp)/(1-G*G)*dGdYj*dGdZl;
					hes[3*j-1][3*l-3]     +=(2*K_phi_temp)/(1-G*G)*dGdZj*dGdXl;
					hes[3*j-1][3*l-2]     +=(2*K_phi_temp)/(1-G*G)*dGdZj*dGdYl;
					hes[3*l-3][3*k-3]     +=(2*K_phi_temp)/(1-G*G)*dGdXl*dGdXk;
					hes[3*l-2][3*k-2]     +=(2*K_phi_temp)/(1-G*G)*dGdYl*dGdYk;
					hes[3*l-1][3*k-1]     +=(2*K_phi_temp)/(1-G*G)*dGdZl*dGdZk;
					hes[3*l-3][3*k-2]     +=(2*K_phi_temp)/(1-G*G)*dGdXl*dGdYk;
					hes[3*l-3][3*k-1]     +=(2*K_phi_temp)/(1-G*G)*dGdXl*dGdZk;
					hes[3*l-2][3*k-3]     +=(2*K_phi_temp)/(1-G*G)*dGdYl*dGdXk;
					hes[3*l-2][3*k-1]     +=(2*K_phi_temp)/(1-G*G)*dGdYl*dGdZk;
					hes[3*l-1][3*k-3]     +=(2*K_phi_temp)/(1-G*G)*dGdZl*dGdXk;
					hes[3*l-1][3*k-2]     +=(2*K_phi_temp)/(1-G*G)*dGdZl*dGdYk;
					hes[3*k-3][3*l-3]     +=(2*K_phi_temp)/(1-G*G)*dGdXk*dGdXl;
					hes[3*k-2][3*l-2]     +=(2*K_phi_temp)/(1-G*G)*dGdYk*dGdYl;
					hes[3*k-1][3*l-1]     +=(2*K_phi_temp)/(1-G*G)*dGdZk*dGdZl;
					hes[3*k-3][3*l-2]     +=(2*K_phi_temp)/(1-G*G)*dGdXk*dGdYl;
					hes[3*k-3][3*l-1]     +=(2*K_phi_temp)/(1-G*G)*dGdXk*dGdZl;
					hes[3*k-2][3*l-3]     +=(2*K_phi_temp)/(1-G*G)*dGdYk*dGdXl;
					hes[3*k-2][3*l-1]     +=(2*K_phi_temp)/(1-G*G)*dGdYk*dGdZl;
					hes[3*k-1][3*l-3]     +=(2*K_phi_temp)/(1-G*G)*dGdZk*dGdXl;
					hes[3*k-1][3*l-2]     +=(2*K_phi_temp)/(1-G*G)*dGdZk*dGdYl;
					hes[3*i-3][3*i-3]     +=2*K_phi_temp/(1-G*G)*dGdXi*dGdXi;
					hes[3*i-2][3*i-2]     +=2*K_phi_temp/(1-G*G)*dGdYi*dGdYi;
					hes[3*i-1][3*i-1]     +=2*K_phi_temp/(1-G*G)*dGdZi*dGdZi;
					hes[3*i-3][3*i-2]     +=2*K_phi_temp/(1-G*G)*dGdXi*dGdYi;
					hes[3*i-3][3*i-1]     +=2*K_phi_temp/(1-G*G)*dGdXi*dGdZi;
					hes[3*i-2][3*i-3]     +=2*K_phi_temp/(1-G*G)*dGdYi*dGdXi;
					hes[3*i-2][3*i-1]     +=2*K_phi_temp/(1-G*G)*dGdYi*dGdZi;
					hes[3*i-1][3*i-3]     +=2*K_phi_temp/(1-G*G)*dGdZi*dGdXi;
					hes[3*i-1][3*i-2]     +=2*K_phi_temp/(1-G*G)*dGdZi*dGdYi;
					hes[3*j-3][3*j-3]     +=2*K_phi_temp/(1-G*G)*dGdXj*dGdXj;
					hes[3*j-2][3*j-2]     +=2*K_phi_temp/(1-G*G)*dGdYj*dGdYj;
					hes[3*j-1][3*j-1]     +=2*K_phi_temp/(1-G*G)*dGdZj*dGdZj;
					hes[3*j-3][3*j-2]     +=2*K_phi_temp/(1-G*G)*dGdXj*dGdYj;
					hes[3*j-3][3*j-1]     +=2*K_phi_temp/(1-G*G)*dGdXj*dGdZj;
					hes[3*j-2][3*j-3]     +=2*K_phi_temp/(1-G*G)*dGdYj*dGdXj;
					hes[3*j-2][3*j-1]     +=2*K_phi_temp/(1-G*G)*dGdYj*dGdZj;
					hes[3*j-1][3*j-3]     +=2*K_phi_temp/(1-G*G)*dGdZj*dGdXj;
					hes[3*j-1][3*j-2]     +=2*K_phi_temp/(1-G*G)*dGdZj*dGdYj;
					hes[3*k-3][3*k-3]     +=2*K_phi_temp/(1-G*G)*dGdXk*dGdXk;
					hes[3*k-2][3*k-2]     +=2*K_phi_temp/(1-G*G)*dGdYk*dGdYk;
					hes[3*k-1][3*k-1]     +=2*K_phi_temp/(1-G*G)*dGdZk*dGdZk;
					hes[3*k-3][3*k-2]     +=2*K_phi_temp/(1-G*G)*dGdXk*dGdYk;
					hes[3*k-3][3*k-1]     +=2*K_phi_temp/(1-G*G)*dGdXk*dGdZk;
					hes[3*k-2][3*k-3]     +=2*K_phi_temp/(1-G*G)*dGdYk*dGdXk;
					hes[3*k-2][3*k-1]     +=2*K_phi_temp/(1-G*G)*dGdYk*dGdZk;
					hes[3*k-1][3*k-3]     +=2*K_phi_temp/(1-G*G)*dGdZk*dGdXk;
					hes[3*k-1][3*k-2]     +=2*K_phi_temp/(1-G*G)*dGdZk*dGdYk;
					hes[3*l-3][3*l-3]     +=2*K_phi_temp/(1-G*G)*dGdXl*dGdXl;
					hes[3*l-2][3*l-2]     +=2*K_phi_temp/(1-G*G)*dGdYl*dGdYl;
					hes[3*l-1][3*l-1]     +=2*K_phi_temp/(1-G*G)*dGdZl*dGdZl;
					hes[3*l-3][3*l-2]     +=2*K_phi_temp/(1-G*G)*dGdXl*dGdYl;
					hes[3*l-3][3*l-1]     +=2*K_phi_temp/(1-G*G)*dGdXl*dGdZl;
					hes[3*l-2][3*l-3]     +=2*K_phi_temp/(1-G*G)*dGdYl*dGdXl;
					hes[3*l-2][3*l-1]     +=2*K_phi_temp/(1-G*G)*dGdYl*dGdZl;
					hes[3*l-1][3*l-3]     +=2*K_phi_temp/(1-G*G)*dGdZl*dGdXl;
					hes[3*l-1][3*l-2]     +=2*K_phi_temp/(1-G*G)*dGdZl*dGdYl;
				}
			}
		}
	}
}


 void build_1st_matrix(struct pdb_atom *strc, double **hes,int nb_atom,double K_r) {

	int i,j;
 	double bx,by,bz,dist,distijsqr;

	for (i=1;i < nb_atom+1;++i) {
		for(j=i+1;j<nb_atom+1;++j) {
 		
	 		int connect_flag = 0;
	 		int t;
	 		for (t=0;t<6;++t) {
	 			if (strc[i-1].node_c[t] == strc[j-1].node) {connect_flag = 1;}
	 		}
	 		
	 		if (connect_flag == 0) {
	 			if (strc[i-1].atom_type == 3 && strc[j-1].atom_type == 3) {} else {
	 				if (j-i == 1) {
	 					printf("		Covalent Break:%s%4d et %s%4d\n",strc[i-1].res_type,strc[i-1].res_number,strc[j-1].res_type,strc[j-1].res_number);
	 				}
	 				continue;
	 			}
	 		}
	 		//printf("I:%d J:%d\n",i,j);

	 		
			bx=strc[i-1].x_cord - strc[j-1].x_cord;
   		bz=strc[i-1].z_cord - strc[j-1].z_cord;
   		by=strc[i-1].y_cord - strc[j-1].y_cord;
   		dist = sqrt(bx*bx+by*by+bz*bz);
   		distijsqr = dist*dist;
	   		

	   		// diagonals of off-diagonal super elements (1st term)
	   		
	   		hes[3*i-3][3*j-3]        += -2*K_r*bx*bx/distijsqr;
			hes[3*i-2][3*j-2]        += -2*K_r*by*by/distijsqr;
			hes[3*i-1][3*j-1]        += -2*K_r*bz*bz/distijsqr;


		
			// off-diagonals of off-diagonal super elements (1st term)
		
			hes[3*i-3][3*j-2]        += -2*K_r*bx*by/distijsqr;
			hes[3*i-3][3*j-1]        += -2*K_r*bx*bz/distijsqr;
			hes[3*i-2][3*j-3]        += -2*K_r*by*bx/distijsqr;
			hes[3*i-2][3*j-1]        += -2*K_r*by*bz/distijsqr;
			hes[3*i-1][3*j-3]        += -2*K_r*bz*bx/distijsqr;
			hes[3*i-1][3*j-2]        += -2*K_r*bz*by/distijsqr;  
		
			//diagonals of off-diagonal super elements (1st term)
		
			hes[3*j-3][3*i-3]        += -2*K_r*bx*bx/distijsqr;
			hes[3*j-2][3*i-2]        += -2*K_r*by*by/distijsqr;
			hes[3*j-1][3*i-1]        += -2*K_r*bz*bz/distijsqr;
		
			//off-diagonals of off-diagonal super elements (1st term)
		
			hes[3*j-3][3*i-2]        += -2*K_r*bx*by/distijsqr;
			hes[3*j-3][3*i-1]        += -2*K_r*bx*bz/distijsqr;
			hes[3*j-2][3*i-3]        += -2*K_r*by*bx/distijsqr;
			hes[3*j-2][3*i-1]        += -2*K_r*by*bz/distijsqr;
			hes[3*j-1][3*i-3]        += -2*K_r*bz*bx/distijsqr;
			hes[3*j-1][3*i-2]        += -2*K_r*bz*by/distijsqr;
		
			// update the diagonals of diagonal super elements
		
			hes[3*i-3][3*i-3]        += 2*K_r*bx*bx/distijsqr;
			hes[3*i-2][3*i-2]        += 2*K_r*by*by/distijsqr;
			hes[3*i-1][3*i-1]        += 2*K_r*bz*bz/distijsqr;
		
			// update the off-diagonals of diagonal super elements
		
			hes[3*i-3][3*i-2]        += 2*K_r*bx*by/distijsqr;
			hes[3*i-3][3*i-1]        += 2*K_r*bx*bz/distijsqr;
			hes[3*i-2][3*i-3]        += 2*K_r*by*bx/distijsqr;
			hes[3*i-2][3*i-1]        += 2*K_r*by*bz/distijsqr;
			hes[3*i-1][3*i-3]        += 2*K_r*bz*bx/distijsqr;
			hes[3*i-1][3*i-2]        += 2*K_r*bz*by/distijsqr;
		
			// update the diagonals of diagonal super elements
		
			hes[3*j-3][3*j-3]        += 2*K_r*bx*bx/distijsqr;
			hes[3*j-2][3*j-2]        += 2*K_r*by*by/distijsqr;
			hes[3*j-1][3*j-1]        += 2*K_r*bz*bz/distijsqr;
		
			// update the off-diagonals of diagonal super elements
		
			hes[3*j-3][3*j-2]        += 2*K_r*bx*by/distijsqr;
			hes[3*j-3][3*j-1]        += 2*K_r*bx*bz/distijsqr;
			hes[3*j-2][3*j-3]        += 2*K_r*by*bx/distijsqr;
			hes[3*j-2][3*j-1]        += 2*K_r*by*bz/distijsqr;
			hes[3*j-1][3*j-3]        += 2*K_r*bz*bx/distijsqr;
			hes[3*j-1][3*j-2]        += 2*K_r*bz*by/distijsqr;                          
		}
	}
}
 
void build_4h_matrix(struct pdb_atom *strc, double **hes,int nb_atom,double epsi,gsl_matrix *tmp) {
 	int i,j;
 	int l,m,n;
 	double bx,by,bz,dist,distijsqr;
 	double temp;
 	int temp_1;
 	int temp_2;
 	int temp_3;
 	for (i=1;i<nb_atom+1;++i) {
 		for (j=1;j<nb_atom+1;++j) {
 			
			if (i==j) {continue;}
			
			int connect_flag = 0;
			if(i>nb_atom - 20) {
				
			}
			for (l=0;l<6;++l) {
				temp_1 = strc[i-1].node_c[l];
				if (temp_1 > nb_atom-1) {break;}
				//printf("Node:%d Connect:%d Node:%d\n",strc[i-1].node,temp_1,strc[j-1].node);
				if (temp_1 == -1) {continue;}
				if (temp_1 == strc[j-1].node) {connect_flag = 1;break;}
				for(m=0;m<6;++m) {
					temp_2 = strc[temp_1].node_c[m];
					if (temp_2 > nb_atom-1) {break;}					
					//printf("	Node:%d Connect:%d Connect:%d Node:%d\n",strc[i-1].node,temp_1,temp_2,strc[j-1].node);
					if (temp_2 == -1) {continue;}
					if (temp_2 == strc[j-1].node) {connect_flag = 1;break;}
					for(n=0;n<6;++n) {
						temp_3 = strc[temp_2].node_c[n];
						if (temp_3 > nb_atom-1) {break;}
						//printf("	Node:%d Connect:%d Connect:%d Node:%d\n",strc[i-1].node,temp_1,temp_2,strc[j-1].node);
						if (temp_3 == -1) {continue;}
					if (temp_3 == strc[j-1].node) {connect_flag = 1;break;}
				}
				}
				if (connect_flag == 1) {break;}			
					
			}
			

			if (connect_flag == 1) {continue;}

 			if ((strc[i-1].atom_type == 3) || (strc[j-1].atom_type == 3) ) {} else {
 			 			
	 			if (abs(strc[i-1].res_number - strc[j-1].res_number) < 4) {
	 				if (strcmp(strc[i-1].chain, strc[j-1].chain) != 0) {
	 					////printf("Chain Break between Chain %s and %s\n",strc[i-1].chain,strc[j-1].chain);
	 				} else {
	 				continue;
	 				}
	 			}
 			}

 			//printf("I'M DOING %d and %d index\t",i-1,j-1);	
 			bx=strc[i-1].x_cord - strc[j-1].x_cord;
       		bz=strc[i-1].z_cord - strc[j-1].z_cord;
       		by=strc[i-1].y_cord - strc[j-1].y_cord;
       		dist = sqrt(bx*bx+by*by+bz*bz);
       		distijsqr = dist*dist*dist*dist;
       		temp = gsl_matrix_get(tmp,i-1,j-1);
       		//printf("%d	%d	%10.5f\n",i,j,temp);
       		      		
     		// diagonals of off-diagonal super elements (1st term)
       		
			hes[3*i-3][3*j-3]        += -120*epsi*bx*bx/distijsqr*temp;
			hes[3*i-2][3*j-2]        += -120*epsi*by*by/distijsqr*temp;
			hes[3*i-1][3*j-1]        += -120*epsi*bz*bz/distijsqr*temp;
		
				
       		// off-diagonals of off-diagonal super elements (1st term)
       		
			hes[3*i-3][3*j-2]        += -120*epsi*bx*by/distijsqr*temp;
			hes[3*i-3][3*j-1]        += -120*epsi*bx*bz/distijsqr*temp;
			hes[3*i-2][3*j-3]        += -120*epsi*by*bx/distijsqr*temp;
			hes[3*i-2][3*j-1]        += -120*epsi*by*bz/distijsqr*temp;
			hes[3*i-1][3*j-3]        += -120*epsi*bx*bz/distijsqr*temp;
			hes[3*i-1][3*j-2]        += -120*epsi*by*bz/distijsqr*temp; 
     
			//update the diagonals of diagonal super elements
			
			hes[3*i-3][3*i-3]        += 120*epsi*bx*bx/distijsqr*temp;
			hes[3*i-2][3*i-2]        += 120*epsi*by*by/distijsqr*temp;
			hes[3*i-1][3*i-1]        += 120*epsi*bz*bz/distijsqr*temp; 
          
            //update the off-diagonals of diagonal super elements
            
			hes[3*i-3][3*i-2]        += 120*epsi*bx*by/distijsqr*temp;
			hes[3*i-3][3*i-1]        += 120*epsi*bx*bz/distijsqr*temp;
			hes[3*i-2][3*i-3]        += 120*epsi*by*bx/distijsqr*temp;
			hes[3*i-2][3*i-1]        += 120*epsi*by*bz/distijsqr*temp;
			hes[3*i-1][3*i-3]        += 120*epsi*bz*bx/distijsqr*temp;
			hes[3*i-1][3*i-2]        += 120*epsi*bz*by/distijsqr*temp;  
			
			//printf("OK\n");
		}
 	} 
 }
 
 void build_enm(struct pdb_atom *strc, double **hes,int nb_atom,gsl_matrix *tmp,float cutoff,int pfanm) {
 	int i,j;
 	double bx,by,bz,dist,distijsqr;
 	double temp;
	//printf("IN BUILD ENM\n");
	double epsi = 1.0;
 	for (i=1;i<nb_atom+1;++i) {
 		for (j=1;j<nb_atom+1;++j) {
 			
			if (i==j) {continue;}
			
			

 			//printf("I'M DOING %d and %d index\t",i-1,j-1);	
 			bx=strc[i-1].x_cord - strc[j-1].x_cord;
       		bz=strc[i-1].z_cord - strc[j-1].z_cord;
       		by=strc[i-1].y_cord - strc[j-1].y_cord;
       		dist = (bx*bx+by*by+bz*bz);
       		if (dist > cutoff*cutoff) {continue;}
       		distijsqr = dist;
       		if (pfanm == 1) {
       			distijsqr = dist*dist;
       		}
       		temp = gsl_matrix_get(tmp,i-1,j-1);
       		epsi = 1.0/120.0;
       		//printf("%d	%d	%10.5f %10.5f\n",i,j,temp,dist);
       		//temp = 1;  		
       		//distijsqr = 1;
     		// diagonals of off-diagonal super elements (1st term)
       		
			hes[3*i-3][3*j-3]        += -120*epsi*bx*bx/distijsqr*temp;
			hes[3*i-2][3*j-2]        += -120*epsi*by*by/distijsqr*temp;
			hes[3*i-1][3*j-1]        += -120*epsi*bz*bz/distijsqr*temp;
			//printf("Value:%10.5f * %10.5f = %10.5f\n",-120.0,epsi,-120*epsi);
				
       		// off-diagonals of off-diagonal super elements (1st term)
       		
			hes[3*i-3][3*j-2]        += -120*epsi*bx*by/distijsqr*temp;
			hes[3*i-3][3*j-1]        += -120*epsi*bx*bz/distijsqr*temp;
			hes[3*i-2][3*j-3]        += -120*epsi*by*bx/distijsqr*temp;
			hes[3*i-2][3*j-1]        += -120*epsi*by*bz/distijsqr*temp;
			hes[3*i-1][3*j-3]        += -120*epsi*bx*bz/distijsqr*temp;
			hes[3*i-1][3*j-2]        += -120*epsi*by*bz/distijsqr*temp; 
     		// printf("2 hessian[%d][%d] = %f\n", 3*i-1,3*j-2,hes[3*i-1][3*j-2]);
			//update the diagonals of diagonal super elements
			//
			hes[3*i-3][3*i-3]        += 120*epsi*bx*bx/distijsqr*temp;
			hes[3*i-2][3*i-2]        += 120*epsi*by*by/distijsqr*temp;
			hes[3*i-1][3*i-1]        += 120*epsi*bz*bz/distijsqr*temp; 
           
            //update the off-diagonals of diagonal super elements
            
			hes[3*i-3][3*i-2]        += 120*epsi*bx*by/distijsqr*temp;
			hes[3*i-3][3*i-1]        += 120*epsi*bx*bz/distijsqr*temp;
			hes[3*i-2][3*i-3]        += 120*epsi*by*bx/distijsqr*temp;
			hes[3*i-2][3*i-1]        += 120*epsi*by*bz/distijsqr*temp;
			hes[3*i-1][3*i-3]        += 120*epsi*bz*bx/distijsqr*temp;
			hes[3*i-1][3*i-2]        += 120*epsi*bz*by/distijsqr*temp;  
			//printf("1 hessian[%d][%d] = %f\n", 3*i-1,3*j-2,hes[3*i-1][3*j-2]);
			//printf("OK\n");
		}
 	} 
 }
 
 
 
  void diagonalyse_matrix (gsl_matrix *m,int nb_atom, gsl_vector *evl,gsl_matrix *evc) {
	
 	gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc (nb_atom);
 	gsl_eigen_symmv (m, evl, evc, w);
 
 	gsl_eigen_symmv_free (w);
 	gsl_eigen_symmv_sort (evl, evc,GSL_EIGEN_SORT_ABS_ASC);
 	
 }
 
 void mass_weight_hessian(gsl_matrix *m,int atom,struct pdb_atom *strc) {
 	// Fonction qui normalize la hessian avec les masses des résidus
  	int i,j;

 	for (i=0;i<atom*3;++i) {
 		//printf("I:%d Index:%d Mass:%f\n",i, i/3 ,mass[i/3]);
 		for (j=0;j<atom*3;++j) {
 			//printf("I:%4d J:%4d Before:%8f ",i,j,gsl_matrix_get(m,i,j));
 			//printf("Mass:%8f %8f ",mass[i/3],mass[j/3]);
 			gsl_matrix_set(m,i,j,gsl_matrix_get(m,i,j)/sqrt(strc[i/3].mass*strc[j/3].mass));
 			//printf("After:%8f\n",gsl_matrix_get(m,i,j));
 		}
 	} 
 }


void adjust_weight_evec(gsl_matrix *m,int atom,struct pdb_atom *strc) {
	int i,j;
	for (i=0;i<atom*3;++i) {
 		for (j=0;j<atom*3;++j) {
 			gsl_matrix_set(m,i,j,gsl_matrix_get(m,i,j)/sqrt(strc[i/3].mass));
 		}
 	} 	


}
