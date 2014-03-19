#include "IsoCleft.h"

/***********************************************************************/
/*        1         2         3         4         5         6         7*/
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/*        1         2         3         4         5         6         7*/
/***********************************************************************/
int main(int argc, char *argv[]){
  psAtom a,b;
  int i,j,k,l,type,inum,m;
  psAGVx v,w;
  float f,Delta_Dist;
  tClique c,top_clique;
  int *egraph;
  int *ca_egraph;
  int ca_node_count;
  int clique_counter=0;
  char filename[100];
  FILE *outfile_ptr;
  gsl_matrix *mat_r  = gsl_matrix_alloc(3,3);
  double DetR;
  int total_cliques;
  int out_size_flag;
  float Similarity=0.0;
  float ca_counter=0.0;
  float coor[3];
  char pdbout[100];
  float rmsd;
  int n;
  int get_all_ori;
  int noat;

  //printf("There are %d ca in a and %d in b",nca_a,nca_b);
  //PAUSE;

  read_commandline(argc, argv);
  //printf("Input Files: %s %s\n",cleftfile_a,cleftfile_b);

  Threshold_JTT_Matrix(JTT_rank_threshold);

  //printf("ok so far\n");PAUSE;

  sprintf(filename,"%s.isd",outbase);
  outfile_ptr = fopen(filename, "w");
  if(outfile_ptr == NULL){
    fprintf(stderr, "ERROR: unable to open output file %s for write\n",filename);
    exit(8);
  }
  fprintf(outfile_ptr,"Output of IsoClefts\n");
  fprintf(outfile_ptr,"Input Files: %s %s\n",cleftfile_a,cleftfile_b);
  fclose(outfile_ptr);
  
  rat_a=malloc(num_of_atom_types*sizeof(psAtom));
  rat_b=malloc(num_of_atom_types*sizeof(psAtom));
  if(rat_a == NULL || rat_b == NULL){
    printf("unable to allocate memory for rat_a or rat_b\n");
    exit(8);
  }

  for(i=0;i<num_of_atom_types;i++) rat_a[i]=NULL;
  for(i=0;i<num_of_atom_types;i++) rat_b[i]=NULL;

  //printf("cleftfiles: <%s> <%s>\n",cleftfile_a,cleftfile_b);PAUSE;
  atoms_a=ReadCleftFile(cleftfile_a,'A');
  atoms_b=ReadCleftFile(cleftfile_b,'B');
  if(atoms_a == NULL || atoms_b == NULL){
    printf("atoms_a OR atoms_b is NULL\n");
    return(8);
  }

  

  if(flag_het == 1){
    het_a=ReadCleftFile(hetfile_a,'A');
    het_b=ReadCleftFile(hetfile_b,'B');

    het_a_isnull=0;
    het_b_isnull=0;
    if(het_a == NULL){
      het_a_isnull=1;
      printf("het_a is NULL\n");
    }
    if(het_b == NULL){
      het_b_isnull=1;
      printf("het_b is NULL\n");
    }

  }

  if(superpose_pdb == 1){
    //printf("pdbfiles: <%s> <%s>\n",pdbfile_a,pdbfile_b);PAUSE;
    all_a=ReadCleftFile(pdbfile_a,'*');
    all_b=ReadCleftFile(pdbfile_b,'*');
    if(all_a == NULL || all_b == NULL){
      printf("all_a OR all_b is NULL\n");
      return(8);
    }
  }
  
  
  //--------------- new section on C-alpha pre-processing-------------
  if(all_atom_only == 1){Skip_Neighbourhood_Test = 1;}  
  if(all_atom_only == 0){
    Skip_Neighbourhood_Test = 0;  
    a=atoms_a;
    do{
      if(strcmp(a->name," CA ") == 0){
	ca_counter += 1.0;
	if(head_ca_a == NULL){
	  head_ca_a = a;
	  a->next_ca = a;
	}else{
	  a->next_ca = head_ca_a->next_ca;
	  head_ca_a->next_ca = a;
	}
	
	b=atoms_b;
	do{
	  if(strcmp(b->name," CA ") == 0 && 
	     Accept_Correspondence(a->rnam,b->rnam) == 1){
	    ca_nodes++;
	  }
	  b=b->next;
	}while(b != atoms_b);  
	
      }
      a=a->next;
    }while(a != atoms_a);  
    
    b=atoms_b;
    do{
      //printf("residues %s %d\n",b->rnam,b->rnum);
      if(strcmp(b->name," CA ") == 0){
	ca_counter += 1.0;
	//nca_b++;
	if(head_ca_b == NULL){
	  head_ca_b = b;
	  b->next_ca = b;
	}else{
	  b->next_ca = head_ca_b->next_ca;
	  head_ca_b->next_ca = b;
	}
      }
      b=b->next;
    }while(b != atoms_b);  
    
    
    agv=malloc(ca_nodes*sizeof(sAGVx));
    if(agv == NULL){
      printf("unable to allocate memory for agv\n");
      exit(8);
    }
    
    ca_node_count=0;
    a=head_ca_a;
    do{
      b=head_ca_b;
      do{
	//printf("residues %s %d and %s %d ",a->rnam,a->rnum,b->rnam,b->rnum);
	if(Accept_Correspondence(a->rnam,b->rnam) == 1){
	  //printf("accepted");
	  agv[ca_node_count].a=a;
	  agv[ca_node_count].b=b;
	  ca_node_count++;
	}
	//printf("\n");
	//PAUSE;
	b=b->next_ca;
      }while(b != head_ca_b);
      a=a->next_ca;
    }while(a != head_ca_a);
    
    printf("there are %d nodes in the C-alpha graph\n",ca_nodes);
    //PAUSE;
    
    ca_edges=0;
    ca_egraph = malloc((ca_nodes*ca_nodes)*sizeof(int));
    if(ca_egraph == NULL){
      printf("unnable to allocate memory to ca_egraph\n");
      exit(1);
    }
    
    
    for(i=0;i<ca_nodes*ca_nodes;i++) ca_egraph[i]=0;
    
    for(i=0;i<ca_nodes;i++){
      ca_egraph[i*ca_nodes+i]=1;
      for(j=0;j<i;j++){
	if(agv[i].a == agv[j].a || agv[i].b == agv[j].b) continue;
	k=i*ca_nodes + j;
	l=j*ca_nodes + i;
	Delta_Dist=fabs(dist3d(agv[i].a->coor,agv[j].a->coor) -
			dist3d(agv[i].b->coor,agv[j].b->coor));
	//printf("Deslta_Dist=%f Thres=%f\n",Delta_Dist,CA_Delta_Dist_Threshold);
	if(Delta_Dist <= CA_Delta_Dist_Threshold){
	  //printf("added edge\n");
	  ca_egraph[k]=1;
	  ca_egraph[l]=1;
	  ca_edges++;
	}
	//PAUSE;
      }
    }
    
    f=ca_nodes*(ca_nodes-1)/2;
    printf("ca_edges: %d (%g)\n",ca_edges,(float)ca_edges/f);
    
    connected=ca_egraph;
    numEdges=ca_edges;
    numNodes=ca_nodes;
    
    /*
      for(i=0;i<numNodes;i++){
      for(j=0;j<numNodes;j++) 
      if(connected[i*numNodes+j] == 1){
      printf("%d a=(%d %d) b=(%d %d)\n",connected[i*numNodes+j],
      ca_agv[i].a->inum,ca_agv[j].a->inum,
      ca_agv[i].b->inum,ca_agv[j].b->inum
      ); 
      printf("\n");
      PAUSE;
      }
      }
    */
    
    // get_all = 0 -> simplified BK
    // get_all = 1 -> exact BK
    
    max_clique_size = 0;
    get_all = 0;                 
    if(bk_mode < 3) get_all = 1; // bk_mode 0,1,2 -> exact
    //printf("get_all = %d\n",get_all);
    total_cliques=Bron_Kerbosch();
    if(total_cliques == 0 && bk_mode > 5){
      printf("1nd stage with 0 cliques, S/E mode\n");
      get_all = 1;
      got_one = 0;
      total_cliques=Bron_Kerbosch();
    }
    printf("C-alpha Bron Kerbosch [bk_mode = %d] finished: %d Total Cliques\n",bk_mode,total_cliques);
    
    // skip rest of Ca procedure of no clique was found...
    if(total_cliques != 0){
      top_clique=Get_Top_Clique();
      
      // Tanimoto Coeficient of Similarity (C/(A+B-C))
      printf("Similarity = %5.3f Vnum=%d\n",ca_counter,top_clique->Vnum);
      top_clique->Similarity = top_clique->Vnum/(ca_counter - top_clique->Vnum);
      
      
      printf("Top Clique: size=%d rmsd=%5.3f \n",top_clique->Vnum, top_clique->rmsd);
      //PAUSE;
      
      printf("Maximal Clique Size=%d\n",max_clique_size);
      
      /*
	c=Cliques;
	if(c != NULL){
	do{
	printf("%d %g\n",c->Vnum,c->DetR);
	PAUSE;
	c=c->next;
	}while(c != Cliques);
	}
	PAUSE;
      */
      
      a=atoms_a;
      do{
	//printf("%s\n",a->nline);
	for(i=0;i<3;i++){
	  coor[i]=top_clique->cen_b[i];
	  for(j=0;j<3;j++) coor[i] += 
			     (a->coor[j]-top_clique->cen_a[j])*gsl_matrix_get(top_clique->mat_r,i,j);
	  //printf("i=%d old=%8.3f new=%8.3f ",i,a->coor[i],coor);
	}
	for(i=0;i<3;i++) a->coor[i]=coor[i];
	//printf("\n");
	a=a->next;
      }while(a != atoms_a);
      
      
      if(superpose_pdb == 1){
	a=all_a;
	do{
	  //printf("%s\n",a->nline);
	  for(i=0;i<3;i++){
	    coor[i]=top_clique->cen_b[i];
	    for(j=0;j<3;j++) coor[i] += 
			       (a->coor[j]-top_clique->cen_a[j])*gsl_matrix_get(top_clique->mat_r,i,j);
	    //printf("i=%d old=%8.3f new=%8.3f ",i,a->coor[i],coor);
	  }
	  for(i=0;i<3;i++) a->coor[i]=coor[i];
	  //printf("\n");
	  a=a->next;
	}while(a != all_a);
      }
      /*
	for(i=0;i<3;i++){
	top_clique->cen_a[i] = 0.0;
	top_clique->cen_b[i] = 0.0;
	for(l=0;l<top_clique->Vnum;l++){
	top_clique->cen_a[i] += agv[top_clique->Vertices[l]].a->coor[i];
	top_clique->cen_b[i] += agv[top_clique->Vertices[l]].b->coor[i];
	}
	top_clique->cen_a[i] /= top_clique->Vnum;
	top_clique->cen_b[i] /= top_clique->Vnum;
	}
      */
      
      Output_Clique_Stats(top_clique,0);
      
      // if no Ca clique was found, reset distance threshold to do all-atom nontheless.
      }else{
      Skip_Neighbourhood_Test = 1;
    }
    free(agv);
    Cliques = NULL;
    printf("end of new C-alpha section\n");//PAUSE;
    if(ca_only == 1) return;
  }
  //---------------------------------------------------------------------------
  num_a=malloc(num_of_atom_types*sizeof(int));
  num_b=malloc(num_of_atom_types*sizeof(int));
  if(num_a == NULL || num_b == NULL){
    printf("unable to allocate memory for num_a or num_b\n");
    exit(8);
  }
  for(i=0;i<num_of_atom_types;i++) num_a[i]=0;
  for(i=0;i<num_of_atom_types;i++) num_b[i]=0;

  a=atoms_b;
  do{
    if(a->type != 7){
      Similarity += 1.0;
    }
    a=a->next;
  }while(a != atoms_b);
  a=atoms_a;
  do{
    if(a->type != 7){
      Similarity += 1.0; 
    }
    a=a->next;
  }while(a != atoms_a);

  // ATENTION AT THIS POINT THE -1 MEANS WE ARE IGNORINNG ALL TYPE 7 CA atoms.
  noat=num_of_atom_types;
  if(num_of_atom_types == 9){noat=7;}
  //printf("%d\n",noat);PAUSE;

  for(type=0;type<noat;type++){
    a=atoms_a;
    do{
      if(a->type == type){
	num_a[type]++;
	if(rat_a[type] == NULL){
	  a->prev=a;
	}else{
	  a->prev=rat_a[type];
	}
	rat_a[type]=a;
      }
      a=a->next;
    }while(a != atoms_a);
    
    if(rat_a[type] != NULL){
      a=rat_a[type];
      do{
	if(a->prev == a) a->prev=rat_a[type];
	a=a->prev;
      }while(a != rat_a[type]);
    }
    
    a=atoms_b;
    do{
      if(a->type == type){
	num_b[type]++;
	if(rat_b[type] == NULL){
	  a->prev=a;
	}else{
	  a->prev=rat_b[type];
	}
	rat_b[type]=a;
      }
      a=a->next;
    }while(a != atoms_b);
    
    if(rat_b[type] != NULL){
      a=rat_b[type];
      do{
	if(a->prev == a) a->prev=rat_b[type];
	a=a->prev;
      }while(a != rat_b[type]);
    }
  }
  
  
  
  numNodes=0;
  for(type=0;type<num_of_atom_types;type++){
      if(rat_a[type] == NULL || rat_b[type] == NULL) continue;
      a=rat_a[type];
      do{
	b=rat_b[type];
	do{
	  if(Skip_Neighbourhood_Test == 1){
	    numNodes++;
	  }else{
	    if(dist3d(a->coor,b->coor) < Threshold_Neighbourhood) numNodes++;
	  }
	  b=b->prev;
	}while(b != rat_b[type]);	
	a=a->prev;
      }while(a != rat_a[type]);
  }
  
  printf("numNodes=%d\n",numNodes);    
  
  agv=malloc(numNodes*sizeof(sAGVx));
  if(agv == NULL){
    printf("unable to allocate memory for agv\n");
    exit(8);
  }
  
  inum=0;
    for(type=0;type<num_of_atom_types;type++){
      if(rat_a[type] == NULL || rat_b[type] == NULL) continue;
      a=rat_a[type];
      do{
	b=rat_b[type];
	do{
	  if(Skip_Neighbourhood_Test == 1){
	    agv[inum].a=a;
	    agv[inum].b=b;
	    inum++;
	  }else{
	    if(dist3d(a->coor,b->coor) < Threshold_Neighbourhood){
	      agv[inum].a=a;
	      agv[inum].b=b;
	      inum++;
	    }
	  }
	  b=b->prev;
	}while(b != rat_b[type]);	
	a=a->prev;
      }while(a != rat_a[type]);
    }
    
    
    /*
      numNodes=0;
      for(type=0;type<NUM_OF_ATOM_TYPES;type++){
      numNodes += num_a[type]*num_b[type];
      //printf("num_a[%d]=%d num_b[%d]=%d\n",type,num_a[type],type,num_b[type]);
      }
      
      printf("numNodes=%d\n",numNodes);
      
      
      agv=malloc(numNodes*sizeof(sAGVx));
      if(agv == NULL){
      printf("unable to allocate memory for agv\n");
      exit(8);
    }
    
    inum=0;
    for(type=0;type<NUM_OF_ATOM_TYPES;type++){
    if(rat_a[type] == NULL || rat_b[type] == NULL) continue;
    a=rat_a[type];
    do{
    b=rat_b[type];
    do{
    agv[inum].a=a;
    agv[inum].b=b;
    inum++;
    b=b->prev;
    }while(b != rat_b[type]);	
    a=a->prev;
    }while(a != rat_a[type]);
    }
    */
    
    /*
      for(i=0;i<numNodes;i++){
      printf("i=%d a=%d (%d) b=%d (%d)\n",i,agv[i].a->inum,agv[i].a->type,
      agv[i].b->inum,agv[i].b->type);
      PAUSE;
      }
    */
    
    numEdges=0;
    egraph = malloc((numNodes*numNodes)*sizeof(int));
    if(egraph == NULL){
      printf("unnable to allocate memory to connected\n");
      exit(1);
    }
    
    
    for(i=0;i<numNodes*numNodes;i++) egraph[i]=0;
    
    for(i=0;i<numNodes;i++){
      egraph[i*numNodes+i]=1;
      for(j=0;j<i;j++){
	if(agv[i].a == agv[j].a || agv[i].b == agv[j].b) continue;
	k=i*numNodes + j;
	l=j*numNodes + i;
	Delta_Dist=fabs(dist3d(agv[i].a->coor,agv[j].a->coor) -
			dist3d(agv[i].b->coor,agv[j].b->coor));
	if(Delta_Dist <= Delta_Dist_Threshold){
	  egraph[k]=1;
	  egraph[l]=1;
	  numEdges++;
	}
      }
    }
    
    f=numNodes*(numNodes-1)/2;
    printf("numEdges: %d (%g)\n",numEdges,(float)numEdges/f);
    
    
    /*
      for(i=0;i<numNodes;i++){
      for(j=0;j<numNodes;j++) 
      if(egraph[i*numNodes+j] == 1){
      printf("%d a=(%d %d) b=(%d %d)\n",egraph[i*numNodes+j],
      agv[i].a->inum,agv[j].a->inum,
      agv[i].b->inum,agv[j].b->inum
      ); 
      //printf("\n");
      PAUSE;
      }
      }
    */
    
    connected=egraph;

    got_one = 0;
    max_clique_size = 0;
    clique_inum=-1;

    get_all = 0;
    if(bk_mode == 0 || bk_mode == 3 || bk_mode == 6) get_all = 1; // Exact 2nd stage
    printf("get_all = %d\n",get_all);
    total_cliques=Bron_Kerbosch();
    if(total_cliques == 0 && (bk_mode == 2 || bk_mode == 5 || bk_mode == 8)){
      printf("2nd stage with 0 cliques, S/E mode\n");
      get_all = 1;
      got_one = 0;
      printf("get_all = %d\n",get_all);
      total_cliques=Bron_Kerbosch();
    }
    printf("Cleft Bron Kerbosch [bk_mode = %d] finished: %d Total Cliques\n",bk_mode,total_cliques);
    
    
    
    //c=Cliques;
    //if(c != NULL){
    //do{
	//printf("%d %g\n",c->Vnum,c->DetR);
      //PAUSE;
    //c=c->next;
    //}while(c != Cliques);
    //}else{
    //printf("no cliques...\n");
    //}
    //PAUSE;
    
    
    /*
      total_cliques=Filter_Cliques(total_cliques);
      printf("After Filtering: %d Cliques left\n",total_cliques);
      
      sprintf(filename,"%s.isd",outbase);
      outfile_ptr = fopen(filename, "w");
      if(outfile_ptr == NULL){
      fprintf(stderr, "ERROR: unable to open output file %s for write\n",filename);
      exit(8);
      }
      fprintf(outfile_ptr,"Output of IsoClefts\n");
      fprintf(outfile_ptr,"Input Files: %s %s\n",cleftfile_a,cleftfile_b);
      fclose(outfile_ptr);
      
      
      //max_clique_size=filter_cliques();
      printf("Maximal Clique Size=%d\n",max_clique_size);
      
      
      //c=Cliques;
      //if(c != NULL){
      //do{
      //printf("%d %g\n",c->Vnum,c->DetR);
      //PAUSE;
      //c=c->next;
      //}while(c != Cliques);
      //}
      //PAUSE;

      c=Cliques;
      if(c != NULL){
      do{
      if(c->reject == 0){
      out_size_flag=1;
      if(get_maximal_clique == 1 && c->Vnum != max_clique_size) out_size_flag=0;
      if(out_size_flag == 1){
      //printf("inum=%d counter=%d\n",c->inum,clique_counter);
      Output_Clique_Stats(c,clique_counter);
      clique_counter++;
      }
      }
      c=c->next;
      }while(c != Cliques);
      }
    */
    
    
    //max_clique_size=filter_cliques();
    if(total_cliques != 0){
      printf("Maximal Clique Size=%d\n",max_clique_size);
      
      
      top_clique=Get_Top_Clique();
      
      // Tanimoto Coeficient of Similarity (C/(A+B-C))
      printf("Similarity = %5.3f Vnum=%d\n",Similarity,top_clique->Vnum);
      Similarity = top_clique->Vnum/(Similarity - top_clique->Vnum);
      printf("Similarity = %5.3f\n",Similarity);
      top_clique->Similarity = Similarity;
      
      Output_Clique_Stats(top_clique,1);
      //Print_Superposed_Clique(top_clique,1);
    }else{
      printf("No Cliques Found!\n");
    }

    return;
}
/***********************************************************************/
/*        1         2         3         4         5         6         7*/
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/*        1         2         3         4         5         6         7*/
/***********************************************************************/
void Threshold_JTT_Matrix(int threshold){
  int i,j;
  int jtt_ranks[20][20] =
    {
      {0,8,10,7,14,10,6,4,15,12,9,11,13,15,5,2,1,17,16,3},
      {6,0,10,14,9,2,11,3,4,12,8,1,13,16,7,5,7,11,15,12},
      {7,11,0,2,15,9,8,6,5,10,13,3,15,16,14,1,4,17,12,13},
      {4,10,2,0,13,8,1,3,6,11,11,9,12,13,10,5,7,14,9,8},
      {4,3,10,13,0,14,15,4,9,11,8,14,12,6,12,1,7,9,2,5},
      {7,3,8,10,17,0,1,11,4,15,6,2,13,17,5,8,9,16,14,12},
      {5,10,7,1,15,2,0,4,12,12,11,3,13,14,11,8,9,15,14,6},
      {2,4,6,5,11,12,3,0,15,16,14,9,16,17,10,1,8,13,16,7},
      {10,3,4,6,13,1,12,11,0,14,7,9,15,12,5,8,9,16,2,15},
      {6,10,8,13,16,15,12,15,15,0,2,9,4,5,14,7,3,17,11,1},
      {8,9,15,16,15,7,14,15,11,1,0,12,4,3,5,6,10,13,12,2},
      {6,1,3,9,15,2,3,7,11,11,10,0,8,14,11,5,4,15,13,12},
      {5,7,11,13,15,9,12,12,14,2,1,6,0,10,13,8,4,16,15,3},
      {7,13,12,14,6,14,13,12,9,4,1,14,8,0,9,3,11,10,2,5},
      {2,6,13,12,16,5,11,8,7,14,3,10,15,13,0,1,4,17,16,9},
      {1,7,3,10,8,13,16,4,17,14,6,11,18,9,5,0,2,19,15,12},
      {1,8,4,11,13,10,11,9,12,3,9,6,7,15,6,2,0,16,14,5},
      {9,1,12,13,5,8,10,3,13,11,2,10,10,6,12,4,9,0,5,7},
      {12,11,5,6,4,11,15,13,2,8,7,14,15,1,14,3,9,11,0,10},
      {2,12,14,10,11,13,7,6,17,1,3,13,4,9,11,8,5,16,15,0}
    };

  //printf("inside JTT threshold=%d\n",threshold);PAUSE;
  for(i=0;i<20;i++){
    for(j=0;j<20;j++){
      jtt[i][j] =0;
      //printf("jtt_ranks[%d][%d]=%d",i,j,jtt_ranks[i][j]);PAUSE;
      if(jtt_ranks[i][j] < threshold) jtt[i][j] = 1;
    }
  }
  
  //printf("inside JTT threshold=%d\n",threshold);PAUSE;
  return;
}
/***********************************************************************/
/*        1         2         3         4         5         6         7*/
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/*        1         2         3         4         5         6         7*/
/***********************************************************************/
int nam2n(char rnam[]){
  
  if(strcmp(rnam,"ALA") == 0) return 0;
  if(strcmp(rnam,"ARG") == 0) return 1;
  if(strcmp(rnam,"ASN") == 0) return 2;
  if(strcmp(rnam,"ASP") == 0) return 3;
  if(strcmp(rnam,"CYS") == 0) return 4;
  if(strcmp(rnam,"GLN") == 0) return 5;
  if(strcmp(rnam,"GLU") == 0) return 6;
  if(strcmp(rnam,"GLY") == 0) return 7;
  if(strcmp(rnam,"HIS") == 0) return 8;
  if(strcmp(rnam,"ILE") == 0) return 9;
  if(strcmp(rnam,"LEU") == 0) return 10;
  if(strcmp(rnam,"LYS") == 0) return 11;
  if(strcmp(rnam,"MET") == 0) return 12;
  if(strcmp(rnam,"PHE") == 0) return 13;
  if(strcmp(rnam,"PRO") == 0) return 14;
  if(strcmp(rnam,"SER") == 0) return 15;
  if(strcmp(rnam,"THR") == 0) return 16;
  if(strcmp(rnam,"TRP") == 0) return 17;
  if(strcmp(rnam,"TYR") == 0) return 18;
  if(strcmp(rnam,"VAL") == 0) return 19;

  return -1;

}
/***********************************************************************/
/*        1         2         3         4         5         6         7*/
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/*        1         2         3         4         5         6         7*/
/***********************************************************************/
int Accept_Correspondence(char rnam_a[],char rnam_b[]){
  int accept = 0;

  //if(strcmp(rnam_a,rnam_b) == 0) accept = 1;

  if(jtt[nam2n(rnam_a)][nam2n(rnam_b)] == 1 || 
     jtt[nam2n(rnam_b)][nam2n(rnam_a)] == 1) accept = 1;

  return accept;
}
/***********************************************************************/
/*        1         2         3         4         5         6         7*/
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/*        1         2         3         4         5         6         7*/
/***********************************************************************/
void Output_Clique_Stats(tClique c,int clique_counter){
  int i,j,m;
  float coor;
  psAtom a,b;
  char filename[100],pdbout[100],hetout[100];
  FILE *outfile_ptr;


  sprintf(filename,"%s.isd",outbase);
  //sprintf(pdbout,"%s_%d.pdb",outbase,clique_counter);
  sprintf(hetout,"%s_%d_het.pdb",outbase,clique_counter);

  //if(output_superimposedPDBs == 1){
  //printf("in Output_Cliqe_Stats will print to %s\n",pdbout);
  //PAUSE;

  printf("CLIQUE %3d SIZE %3d RMSD %5.3f\n",clique_counter,c->Vnum,c->rmsd);

  outfile_ptr = fopen(filename, "a");
  if(outfile_ptr == NULL){
    fprintf(stderr, "ERROR: unable to open output file %s for write\n",filename);
    exit(8);
  }
  
  fprintf(outfile_ptr,"-- CLIQUE %3d inum %d\n",clique_counter,c->inum);
  fprintf(outfile_ptr,"REMARK  COMMAND %s\n",command);
  if(flag_het == 1) fprintf(outfile_ptr,"HFILE %s\n",hetout);
  fprintf(outfile_ptr,"REMARK  Centre A: ");
  for(m=0;m<3;m++) fprintf(outfile_ptr,"%8.3f",c->cen_a[m]);
  fprintf(outfile_ptr,"\n");
  fprintf(outfile_ptr,"REMARK  Centre B: ");
  for(m=0;m<3;m++) fprintf(outfile_ptr,"%8.3f",c->cen_b[m]);
  fprintf(outfile_ptr,"\n");
  fprintf(outfile_ptr,"REMARK  Rotation Matrix:\n");
  for (i=0;i<3;i++) {
    fprintf(outfile_ptr,"REMARK  ");
    for (j=0;j<3;j++) fprintf(outfile_ptr,"%8.3f",gsl_matrix_get(c->mat_r,i,j));
    fprintf(outfile_ptr,"\n");
  }
  fprintf(outfile_ptr,"REMARK  Rotation Matrix Determinant: %g\n",c->DetR);
  fprintf(outfile_ptr,"REMARK  Clique %3d SIZE %3d RMSD %5.3f Similarity %5.3f\n",
	  clique_counter,c->Vnum,c->rmsd,c->Similarity);
  for(i=0;i<c->Vnum;i++)
    fprintf(outfile_ptr,"REMARK  NODE %7d [%5d,%5d]\n",c->Vertices[i],
	    agv[c->Vertices[i]].a->pnum,
	    agv[c->Vertices[i]].b->pnum);
  fclose(outfile_ptr);
  
  return;
}
/***********************************************************************/
/*        1         2         3         4         5         6         7*/
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/*        1         2         3         4         5         6         7*/
/***********************************************************************/
void Print_Superposed_Clique(tClique c,int clique_counter){
  int i,j,m,rnum,jc;
  float coor;
  psAtom a,b;
  char filename[100],pdbout[100],icmscript[100],achn,anam[4],orstr[20];
  char pdbname_a[100],pdbname_b[100];
  FILE *outfile_ptr;

  sprintf(filename,"%s_%d.pdb",outbase,clique_counter);

  printf("File: %s \n",filename);

  outfile_ptr = fopen(filename, "w");
  if(outfile_ptr == NULL){
    fprintf(stderr, "ERROR: unable to open output file %s for write\n",filename);
    exit(8);
  }
  
  fprintf(outfile_ptr,"REMARK  Output of IsoClefts\n");
  fprintf(outfile_ptr,"REMARK  COMMAND %s\n",command);
  fprintf(outfile_ptr,"REMARK  Input Files: %s %s\n",cleftfile_a,cleftfile_b);

  fprintf(outfile_ptr,"REMARK  Centre A: ");
  for(m=0;m<3;m++) fprintf(outfile_ptr,"%8.3f",c->cen_a[m]);
  fprintf(outfile_ptr,"\n");
  fprintf(outfile_ptr,"REMARK  Centre B: ");
  for(m=0;m<3;m++) fprintf(outfile_ptr,"%8.3f",c->cen_b[m]);
  fprintf(outfile_ptr,"\n");

  fprintf(outfile_ptr,"REMARK  Rotation Matrix:\n");
  for (i=0; i<3; i++) {
    fprintf(outfile_ptr,"REMARK  ");
    for (j=0; j<3; j++) fprintf(outfile_ptr,"%8.3f",gsl_matrix_get(c->mat_r,i,j));
    fprintf(outfile_ptr,"\n");
  }
  fprintf(outfile_ptr,"REMARK  Rotation Matrix Determinant: %g\n",c->DetR);
  fprintf(outfile_ptr,"REMARK  Clique %3d SIZE %3d RMSD %5.3f Similarity %5.3f\n",
	  clique_counter,c->Vnum,c->rmsd,c->Similarity);
  //fprintf(outfile_ptr,"REMARK  Clique (size=%d):\n",c->Vnum);
  for(i=0;i<c->Vnum;i++)
    fprintf(outfile_ptr,"REMARK  NODE %7d [%5d,%5d]\n",c->Vertices[i],
	    agv[c->Vertices[i]].a->pnum,
	    agv[c->Vertices[i]].b->pnum);
  a=atoms_a;
  do{
    fprintf(outfile_ptr,"%s",a->nline);
    for(i=0;i<3;i++){
      coor=c->cen_b[i];
      for(j=0;j<3;j++) coor += (a->coor[j]-c->cen_a[j])*gsl_matrix_get(c->mat_r,i,j);
      fprintf(outfile_ptr,"%8.3f",coor);
    }
    fprintf(outfile_ptr,"\n");
    a=a->next;
  }while(a != atoms_a);

  b=atoms_b;
  do{
    fprintf(outfile_ptr,"%s",b->nline);
    for(m=0;m<3;m++) fprintf(outfile_ptr,"%8.3f",b->coor[m]);
    fprintf(outfile_ptr,"\n");
    b=b->next;
  }while(b != atoms_b);

  fclose(outfile_ptr);

  if(flag_het == 1){
    sprintf(pdbout,"%s_%d_het.pdb",outbase,clique_counter);
    //printf("<%s>\n",pdbout);PAUSE;

    outfile_ptr = fopen(pdbout, "w");
    if(outfile_ptr == NULL){
      fprintf(stderr, "ERROR: unable to open output file %s for write\n",pdbout);
      exit(8);
    }
    
    fprintf(outfile_ptr,"MODEL     1\n");
    b=het_b;
    do{
      fprintf(outfile_ptr,"%s",b->nline);
      for(m=0;m<3;m++) fprintf(outfile_ptr,"%8.3f",b->coor[m]);
      fprintf(outfile_ptr,"\n");
      b=b->next;
    }while(b != het_b);
    fprintf(outfile_ptr,"ENDMDL\n");
    fprintf(outfile_ptr,"MODEL     2\n");
    a=het_a;
    do{
      fprintf(outfile_ptr,"%s",a->nline);
      for(i=0;i<3;i++){
	coor=c->cen_b[i];
	for(j=0;j<3;j++) coor += (a->coor[j]-c->cen_a[j])*gsl_matrix_get(c->mat_r,i,j);
	fprintf(outfile_ptr,"%8.3f",coor);
      }
      fprintf(outfile_ptr,"\n");
      a=a->next;
    }while(a != het_a);
    fprintf(outfile_ptr,"ENDMDL\n");

    fclose(outfile_ptr);
  }

  //printf("superpose_pdb=%d\n",superpose_pdb);PAUSE;
  if(superpose_pdb == 1){
    strncpy(pdbname_a,pdbfile_a,strlen(pdbfile_a)-4);
    pdbname_a[strlen(pdbfile_a)-4]='\0';
    strncpy(pdbname_b,pdbfile_b,strlen(pdbfile_b)-4);
    pdbname_b[strlen(pdbfile_b)-4]='\0';
    //printf("pdbname_a=<%s> pdbname_b=<%s>\n",pdbname_a,pdbname_b);PAUSE;

    sprintf(pdbout,"%s_B.pdb",outbase);
    //printf("<%s>\n",pdbout);PAUSE;

    outfile_ptr = fopen(pdbout, "w");
    if(outfile_ptr == NULL){
      fprintf(stderr, "ERROR: unable to open output file %s for write\n",pdbout);
      exit(8);
    }
    
    b=all_b;
    do{
      //printf("%s",b->nline);
      //PAUSE;
      fprintf(outfile_ptr,"%s",b->nline);
      for(m=0;m<3;m++) fprintf(outfile_ptr,"%8.3f",b->coor[m]);
      fprintf(outfile_ptr,"\n");
      b=b->prev;
    }while(b != all_b);
    fclose(outfile_ptr);

    sprintf(pdbout,"%s_A.pdb",outbase);
    //printf("<%s>\n",pdbout);PAUSE;
    outfile_ptr = fopen(pdbout, "w");
    if(outfile_ptr == NULL){
      fprintf(stderr, "ERROR: unable to open output file %s for write\n",pdbout);
      exit(8);
    }
    a=all_a;
    do{
      fprintf(outfile_ptr,"%s",a->nline);
      for(i=0;i<3;i++){
	coor=c->cen_b[i];
	for(j=0;j<3;j++) coor += (a->coor[j]-c->cen_a[j])*gsl_matrix_get(c->mat_r,i,j);
	fprintf(outfile_ptr,"%8.3f",coor);
      }
      fprintf(outfile_ptr,"\n");
      a=a->prev;
    }while(a != all_a);
    fclose(outfile_ptr);

    // icm script to visualize clique

    sprintf(icmscript,"%s.icm",outbase);
    outfile_ptr = fopen(icmscript, "w");
    if(outfile_ptr == NULL){
      fprintf(stderr, "ERROR: unable to open output file %s for write\n",icmscript);
      exit(8);
    }
    fprintf(outfile_ptr,"read pdb \"%s_A.pdb\" name=\"%s\"\n",outbase,pdbname_a);
    fprintf(outfile_ptr,"read pdb \"%s_B.pdb\" name=\"%s\"\n",outbase,pdbname_b);
    fprintf(outfile_ptr,"assign sstructure a_*.*\n");

    for(i=0;i<c->Vnum;i++){
      achn=tolower(agv[c->Vertices[i]].a->line[21]);
      if(isspace(achn)) achn = 'm';
      rnum=agv[c->Vertices[i]].a->rnum;
      jc=0;
      for(j=0;j<4;j++){
	if(!isspace(agv[c->Vertices[i]].a->name[j])){
	  anam[jc]=tolower(agv[c->Vertices[i]].a->name[j]);
	  jc++;
	}
      }
      anam[jc]='\0';
      sprintf(orstr,"CLQ_A | ");
      if(i==0) sprintf(orstr,"");
      fprintf(outfile_ptr,"CLQ_A = %sa_1.%c/%d/%s\n",orstr,achn,rnum,anam);
    }
    fprintf(outfile_ptr,"CLQ_A_r = Res(CLQ_A)\n");    
    for(i=0;i<c->Vnum;i++){
      achn=tolower(agv[c->Vertices[i]].b->line[21]);
      if(isspace(achn)) achn = 'm';
      rnum=agv[c->Vertices[i]].b->rnum;
      jc=0;
      for(j=0;j<4;j++){
	if(!isspace(agv[c->Vertices[i]].b->name[j])){
	  anam[jc]=tolower(agv[c->Vertices[i]].b->name[j]);
	  jc++;
	}
      }
      anam[jc]='\0';
      sprintf(orstr,"CLQ_B | ");
      if(i==0) sprintf(orstr,"");
      fprintf(outfile_ptr,"CLQ_B = %sa_2.%c/%d/%s\n",orstr,achn,rnum,anam);
    }
    fprintf(outfile_ptr,"CLQ_B_r = Res(CLQ_B)\n");    
    
    fprintf(outfile_ptr,"ds Name(a_1.) blue, -0.9 0.9\n");
    fprintf(outfile_ptr,"ds Name(a_2.) red, -0.9 0.8\n");
    fprintf(outfile_ptr,"ds ribbon a_*.*\n");    
    fprintf(outfile_ptr,"ds wire CLQ_A_r 35\n");    
    fprintf(outfile_ptr,"ds wire CLQ_B_r 230\n");    
    fprintf(outfile_ptr,"ds ball CLQ_A blue\n");    
    fprintf(outfile_ptr,"ds ball CLQ_B red\n");    
    fprintf(outfile_ptr,"color ribbon a_1.* 35\n");
    fprintf(outfile_ptr,"color ribbon a_2.* 230\n");
    fprintf(outfile_ptr,"ds residue label CLQ_A_r 35\n");
    fprintf(outfile_ptr,"ds residue label CLQ_B_r 230\n");
    fprintf(outfile_ptr,"ds atom label CLQ_A blue\n");
    fprintf(outfile_ptr,"ds atom label CLQ_B red\n");
    fprintf(outfile_ptr,"center CLQ_A\n");
    fclose(outfile_ptr);
  }



  return;
}
/***********************************************************************/
/*        1         2         3         4         5         6         7*/
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/*        1         2         3         4         5         6         7*/
/***********************************************************************/
void Calculate_Rotation_GSL(tClique c){
  int i,j,l,k,m;
  double val;

  for(i=0;i<3;i++){
    c->cen_a[i] = 0.0;
    c->cen_b[i] = 0.0;
    for(l=0;l<c->Vnum;l++){
      c->cen_a[i] += agv[c->Vertices[l]].a->coor[i];
      c->cen_b[i] += agv[c->Vertices[l]].b->coor[i];
    }
    c->cen_a[i] /= c->Vnum;
    c->cen_b[i] /= c->Vnum;
  }

  for(i=0;i<3;i++){
    for(j=0;j<3;j++){
      val = 0.0;
      for(l=0;l<c->Vnum;l++){
	val += (agv[c->Vertices[l]].a->coor[i] - c->cen_a[i])*
	  (agv[c->Vertices[l]].b->coor[j] - c->cen_b[j]);
      }
      gsl_matrix_set(c->mat_r,i,j,val);
    }
  }


  c->DetR=SupSVD(c->mat_r);
  
  return ;
}
/***********************************************************************/
/*        1         2         3         4         5         6         7*/
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/*        1         2         3         4         5         6         7*/
/***********************************************************************/
double SupSVD(gsl_matrix *mat_u){
  int i,j,k,m;
  double val;
  double eps = 0.0001;

  //gsl_matrix *mat_u  = gsl_matrix_alloc(3,3); 
  gsl_matrix *mat_ut = gsl_matrix_alloc(3,3); 
  gsl_matrix *mat_v  = gsl_matrix_alloc(3,3); 
  //gsl_matrix *mat_vt = gsl_matrix_alloc(3,3); 
  //gsl_matrix *mat_s  = gsl_matrix_alloc(3,3); 
  //gsl_matrix *mat_us = gsl_matrix_alloc(3,3); 

  gsl_vector *vec_s  = gsl_vector_alloc(3);
  gsl_vector *vec_w  = gsl_vector_alloc(3);

  /*
    for(i=0;i<3;i++)
    for(j=0;j<3;j++){
    m=3*i+j;
    gsl_matrix_set(mat_u,i,j,*(mat+m));
    }
  */

  /*
    printf("\n");
    for(i=0;i<3;i++){
    for(j=0;j<3;j++)
    printf("A(%d,%d)=%g ",i,j,gsl_matrix_get(mat_u,i,j));
    printf("\n");
    }
  */

  gsl_linalg_SV_decomp(mat_u,mat_v,vec_s,vec_w);
  /*
    This function factorizes the M-by-N matrix A (mat_u) into the singular 
    value  decomposition A = USV^T. On output the matrix A is replaced by U. 
    The diagonal elements of the singular value matrix S are stored in the 
    vector S (vec_s). The singular values are non-negative and form a 
    nonincreasing sequence from S_1 to S_N. The matrix V (mat_v) contains the 
    elements of V in untransposed form. To form the product USV^T it is 
    necessary to take the transpose of V. A workspace of length N is required 
    in work (vec_w). This routine uses the Golub-Reinsch SVD algorithm. 
  */
 
  /*
    printf("\n");
    for(i=0;i<3;i++){
    for(j=0;j<3;j++)
    printf("u[%d,%d]=%g ",i,j,gsl_matrix_get(mat_u,i,j));
    printf("\n");
    }
    
    
    for(i=0;i<3;i++){
    for(j=0;j<3;j++){
    gsl_matrix_set(mat_s,i,j,0.0);
    if(i == j) gsl_matrix_set(mat_s,i,j,gsl_vector_get(vec_s,i));
    }
    }
  */

  gsl_matrix_transpose_memcpy(mat_ut,mat_u);
  //gsl_matrix_transpose_memcpy(mat_vt,mat_v);

  /*
    for(i=0;i<3;i++){
    for(j=0;j<3;j++){
    m=3*i+j;
    val=0.0;
    for(k=0;k<3;k++){
    val += gsl_matrix_get(mat_u,i,k)*gsl_matrix_get(mat_ut,k,j);
    }
    if(val*val < eps) val=0.0;
    printf("uut[%d,%d]=%g ",i,j,val);
    }
    printf("\n");
    }
    for(i=0;i<3;i++){
    for(j=0;j<3;j++){
    m=3*i+j;
    val=0.0;
    for(k=0;k<3;k++){
    val += gsl_matrix_get(mat_v,i,k)*gsl_matrix_get(mat_vt,k,j);
    }
    if(val*val < eps) val=0.0;
    printf("vvt[%d,%d]=%g ",i,j,val);
    }
    printf("\n");
    }
    
    for(i=0;i<3;i++){
    for(j=0;j<3;j++){
    m=3*i+j;
    val=0.0;
    for(k=0;k<3;k++){
    val += gsl_matrix_get(mat_u,i,k)*gsl_matrix_get(mat_s,k,j);
    }
    //if(val*val < eps) val=0.0;
    gsl_matrix_set(mat_us,i,j,val);
    }
    }

    for(i=0;i<3;i++){
    for(j=0;j<3;j++){
    m=3*i+j;
    val=0.0;
    for(k=0;k<3;k++){
    val += gsl_matrix_get(mat_us,i,k)*gsl_matrix_get(mat_vt,k,j);
    }
    if(val*val < eps) val=0.0;
    printf("usvt[%d,%d]=%g ",i,j,val);
    }
    printf("\n");
    }
    
    
    printf("gnu_Det3D(u)=%g\n",gsl_matrix_Det3D(mat_u));
  */

  for(i=0;i<3;i++)
    for(j=0;j<3;j++){
      m=3*i+j;
      val=0.0;
      for(k=0;k<3;k++){
	val += gsl_matrix_get(mat_v,i,k)*gsl_matrix_get(mat_ut,k,j);
      }
      if(val*val < eps) val = 0.0;
      gsl_matrix_set(mat_u,i,j,val);      
    }
  
  //printf("gnu_Det3D(r)=%g\n",gsl_matrix_Det3D(mat_u));
  //PAUSE;

  gsl_matrix_free(mat_ut);
  gsl_matrix_free(mat_v);
  gsl_vector_free(vec_s);
  gsl_vector_free(vec_w);
 
  return gsl_matrix_Det3D(mat_u);

}
/***********************************************************************/
/*        1         2         3         4         5         6         7*/
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/*        1         2         3         4         5         6         7*/
/***********************************************************************/
double gsl_matrix_Det3D(gsl_matrix *M){
  int i,j;
  double det;

  //  guide to indexes: 0=x, 1=y, 2=z
  //  00 01 02
  //  10 11 12
  //  20 21 22

  det = gsl_matrix_get(M,0,0)*gsl_matrix_get(M,1,1)*gsl_matrix_get(M,2,2)
    + gsl_matrix_get(M,0,1)*gsl_matrix_get(M,1,2)*gsl_matrix_get(M,2,0)
    + gsl_matrix_get(M,0,2)*gsl_matrix_get(M,2,1)*gsl_matrix_get(M,1,0)
    - gsl_matrix_get(M,0,2)*gsl_matrix_get(M,1,1)*gsl_matrix_get(M,2,0)
    - gsl_matrix_get(M,1,2)*gsl_matrix_get(M,2,1)*gsl_matrix_get(M,0,0)
    - gsl_matrix_get(M,2,2)*gsl_matrix_get(M,1,0)*gsl_matrix_get(M,0,1);

  return det;
}
/***********************************************************************/
/*        1         2         3         4         5         6         7*/
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/*        1         2         3         4         5         6         7*/
/***********************************************************************/
psAtom ReadCleftFile(char file[], char chn){
  FILE  *infile_ptr;        // pointer to input file
  int i,j;
  char  buffer[100];         // a line from the INPUT file
  char  coor_char[10];      // string used to read the coordinates
  char  field[7];
  psAtom a;
  psAtom atoms=NULL;
  int inum=0;

  infile_ptr = fopen(file, "r");

  if(infile_ptr == NULL){
    fprintf(stderr, "ERROR: unable to open input cleft file: %s\n",file);
    exit(8);
  }
  while (fgets(buffer, sizeof(buffer),infile_ptr)){
    for(i=0;i<6;i++) field[i]=buffer[i];
    field[6]='\0';
    if(strcmp(field,"ATOM  ") != 0 && strcmp(field,"HETATM") != 0) continue;
    if(buffer[13] == 'H') continue; // IGNORING HYDROGEN ATOMS

    NEW( a, sAtom );
    a->inum=inum;

    for(i=0;i<81;i++) a->line[i]=buffer[i];
    a->line[80]='\0';
    //printf("<%s>\n",a->line);

    for(i=0;i<30;i++) a->nline[i]=buffer[i];
    if(chn  != '*') a->nline[21]=chn;
    //a->nline[29]='\0';
    //printf("<%s>",a->nline);PAUSE;

    for(i=0;i<5;i++) coor_char[i]=buffer[i+6];
    coor_char[5]='\0';
    sscanf(coor_char,"%d",&a->pnum);

    for(i=0;i<=3;i++){a->name[i]=buffer[i+12];}
    a->name[4]='\0';

    a->radius=assign_radius(a->name);
    
    for(i=0;i<3;i++) a->rnam[i]=buffer[i+17];
    a->rnam[4]='\0';

    //a->type=ASSIGN_ATOM_TYPE(a->rnam,a->name);
    a->type=atom_type_function(a->rnam,a->name);
    //printf(" ->%d\n",a->type);

    for(i=0;i<4;i++) coor_char[i]=buffer[i+22];
    coor_char[4]='\0';
    sscanf(coor_char,"%d",&a->rnum);
    //printf("rnum = %d (%s)\n",a->rnum,coor_char);
    //PAUSE;

    for (j=0;j<=2;j++){
      for(i=0;i<8;i++){coor_char[i]=buffer[30+j*8+i];}
      coor_char[8]='\0';
      sscanf(coor_char,"%f",&a->coor[j]);
    }

    a->next_ca=NULL;

    if(atoms){
      a->next=atoms->next;
      if(chn  == '*'){
	a->prev=atoms;
	atoms->next->prev = a;
      }
      atoms->next = a;
    }else{
      atoms = a;
      a->next = atoms;
      if(chn  == '*') a->prev = atoms;
    }

    inum++;
  }

  return atoms;
}
/***********************************************************************/
/*        1         2         3         4         5         6         7*/
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/*        1         2         3         4         5         6         7*/
/***********************************************************************/
float dist3d(float a[], float b[]){
  float dist=0.0;
  int i;

  for(i=0;i<3;i++){
    dist += (a[i]-b[i])*(a[i]-b[i]);
  }
  dist = sqrt(dist);

  return dist;
}
/***********************************************************************/
/*        1         2         3         4         5         6         7*/
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/*        1         2         3         4         5         6         7*/
/***********************************************************************/
void read_commandline(int argc, char *argv[]){
  int   i,j;
  int   die_flag=0;
  char  usage[1000],helper_str[200];
  int   manc=1;
  char* manv[]={"-i"}; // obligatory argument handles
  int   man_ok=0;
  int   incr;
  int   len_command=0;
  char atom_type_function_char[4];

  // assignment of default values to optional parameters
  strcpy(outbase,"ICO");
  flag_het=0;
  get_maximal_clique=1;
  JTT_rank_threshold=5;
  CA_Delta_Dist_Threshold=3.50;
  Threshold_Neighbourhood=4.00;
  Delta_Dist_Threshold=4.00;
  bk_mode = 0;
  all_atom_only = 0;
  ca_only = 0;
  output_superimposedPDBs = 0;

  strcpy(atom_type_function_char,"LPC\0");

  strcpy(usage,"Usage:\n      IsoClefts [obligatory-arguments list] ");
  strcpy(usage,"[optional-arguments list]\n\n");
  strcat(usage,"Obligatory Arguments:\n");
  strcat(usage,"-i <cleft_a> <cleft_b>: Cleft files to be compared\n");
  strcat(usage,"Optional Arguments:\n");
  strcat(usage,"-h\t\t\t: Include Hetero Groups in full output\n");
  sprintf(helper_str,"-oO\t\t\t: Output base filename[%s]\n",outbase);
  strcat(usage,helper_str);

  sprintf(helper_str,"-r\t\t\t: JTT rank Threshold, between 1->20 [%d]\n",JTT_rank_threshold);
  strcat(usage,helper_str);

  sprintf(helper_str,"-c\t\t\t: C-aplpha Node Distance Threshold [%4.2f A]\n",CA_Delta_Dist_Threshold);
  strcat(usage,helper_str);

  sprintf(helper_str,"-n\t\t\t: Neighborhood sphere radius [%4.2f A]\n",Threshold_Neighbourhood);
  strcat(usage,helper_str);

  sprintf(helper_str,"-d\t\t\t: Node Distance Threshold [%4.2f A]\n",Delta_Dist_Threshold);
  strcat(usage,helper_str);

  strcat(usage,"-p <1st pdb> <2nd pdb>\t: PDB filenames for clique based superimposition\n");
  sprintf(helper_str,"-s \t\t\t: Level of BK search simplification [%d]\n",bk_mode);
  strcat(usage,helper_str);
  strcat(usage,"-f \t\t\t: Ca atom (first stage) only\n");
  strcat(usage,"-a \t\t\t: All atom (second stage) only\n");
  strcat(usage,"-t \t\t\t: Chemical element atom types (N,C,O,S,P) [LPC types]\n");
  
  for(i=0;i<argc;i++){
    len_command += strlen(argv[i]);
    //printf("argv[%d]=%s len=%d len_command=%d \n",i,argv[i],strlen(argv[i]),len_command);
    //PAUSE;
  }

  strcat(command,argv[0]);
  if(len_command < MAX_COMMAND_LINE_LEN){
    for(i=1;i<argc;i++){
      strcat(command," ");
      strcat(command,argv[i]);
    }
    strcat(command,"\0");
  }else{
    printf("len_command %d > than max %d\n",len_command,MAX_COMMAND_LINE_LEN);
  }

  printf("command: %s\n",command);
  // check existence of mandatory arguments
  for(j=0;j<manc;j++){
    for(i=1;i<argc-1;i++){
      if(strcmp(argv[i],manv[j]) == 0){man_ok++;}
    } 
  }
  if(man_ok < manc){die_flag=1;}

  // die or let live...
  if(die_flag != 0){
    printf("INPUT ERROR\n%s\n",usage);
    exit(8);
  }

  i=1;
  // copy argv values to the respective global arguments
  while(i < argc){
    //printf("so far ok, reading arg %d\n",i);
  
    if(strcmp(argv[i],"-i")==0){
      strcpy(cleftfile_a,argv[i+1]);
      strcpy(cleftfile_b,argv[i+2]);
      //printf("Input Files: %s %s\n",cleftfile_a,cleftfile_b);PAUSE;
      incr=3;      
    }

    if(strcmp(argv[i],"-o")==0 || strcmp(argv[i],"-O")==0){
      strcpy(outbase,argv[i+1]);
      if(strcmp(argv[i],"-O")==0){output_superimposedPDBs = 1;}
      incr=2;      
    }

    if(strcmp(argv[i],"-h")==0){
      flag_het=1;
      strcpy(hetfile_a,argv[i+1]);
      strcpy(hetfile_b,argv[i+2]);
      incr=3;      
    }

    if(strcmp(argv[i],"-r")==0){
      sscanf(argv[i+1],"%d",&JTT_rank_threshold);
      incr=2;
    }

    if(strcmp(argv[i],"-c")==0){
      sscanf(argv[i+1],"%f",&CA_Delta_Dist_Threshold);
      incr=2;
    }

    if(strcmp(argv[i],"-n")==0){
      sscanf(argv[i+1],"%f",&Threshold_Neighbourhood);
      incr=2;
    }

    if(strcmp(argv[i],"-d")==0){
      sscanf(argv[i+1],"%f",&Delta_Dist_Threshold);
      incr=2;
    }

    if(strcmp(argv[i],"-p")==0){
      superpose_pdb=1;
      strcpy(pdbfile_a,argv[i+1]);
      strcpy(pdbfile_b,argv[i+2]);
      incr=3;  
    }

    if(strcmp(argv[i],"-s")==0){
      sscanf(argv[i+1],"%d",&bk_mode);
      incr=2;  
    }

    if(strcmp(argv[i],"-a")==0){
      all_atom_only = 1;
      incr=1;  
    }

    if(strcmp(argv[i],"-f")==0){
      ca_only = 1;
      incr=1;  
    }

    if(strcmp(argv[i],"-t")==0){
      strcpy(atom_type_function_char,"ETY");
      incr=1;
    }

    i = i + incr;
  }

  if(strcmp(atom_type_function_char,"LPC")==0){
    //printf("using LPC\n");
    atom_type_function = atom_types_LPC;
    num_of_atom_types=9;
  }
  if(strcmp(atom_type_function_char,"ETY")==0){
    //printf("using ETY\n");
    atom_type_function = atom_types_ETY;
    num_of_atom_types=6;
  }

  return;
}

/***********************************************************************/
/*        1         2         3         4         5         6         7*/
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/*        1         2         3         4         5         6         7*/
/***********************************************************************/
float assign_radius(char atm[]){
  float r;
  char  c;

  c=atm[1];
  
  switch(c){
  case 'N':
    r=1.70;
    break;
  case 'O':
    r=1.50;
    break;
  case 'C':
    r=1.90;
    break;
  case 'S':
    r=1.90;
    break;
  case 'P':
    r=1.23;
    break;
  case 'L':
    r=1.75;
    break;
  case 'F':
    r=1.70;
    break;
  case 'R':
    r=1.85;
    break;
  case 'I':
    r=2.00;
    break;
  case 'E':
    r=1.50;
    break;
  case 'A':
    r=1.50;
    break;
  default :
    r=1.50;
    break;
  }

  return (r);
}
/***********************************************************************/
/*        1         2         3         4         5         6         7*/
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/*        1         2         3         4         5         6         7*/
/***********************************************************************/
int atom_types_ETY(char rnam[], char anam[]){
  //printf("in ETY, %s <%s>",rnam,anam);

  if(anam[0] == 'N' || anam[1] == 'N') return 1; // does not work for Na
  if(anam[1] == 'C' || anam[1] == 'C') return 2; // does not work for Ca
  if(anam[0] == 'O' || anam[1] == 'O') return 3;
  if(anam[0] == 'S' || anam[1] == 'S') return 4;
  if(anam[0] == 'P' || anam[1] == 'P') return 5;
  return 0;
}
/***********************************************************************/
/*        1         2         3         4         5         6         7*/
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/*        1         2         3         4         5         6         7*/
/***********************************************************************/
int atom_types_LPC(char rnam[], char anam[]){
  //printf("in LPC, %s <%s>",rnam,anam);

  if(strcmp(anam," N  ") == 0 && strcmp(rnam,"ALA") == 0) return 3;
  if(strcmp(anam," CA ") == 0 && strcmp(rnam,"ALA") == 0) return 7;
  if(strcmp(anam," C  ") == 0 && strcmp(rnam,"ALA") == 0) return 6;
  if(strcmp(anam," O  ") == 0 && strcmp(rnam,"ALA") == 0) return 2;
  if(strcmp(anam," CB ") == 0 && strcmp(rnam,"ALA") == 0) return 4;

  if(strcmp(anam," N  ") == 0 && strcmp(rnam,"ARG") == 0) return 3;
  if(strcmp(anam," CA ") == 0 && strcmp(rnam,"ARG") == 0) return 7;
  if(strcmp(anam," C  ") == 0 && strcmp(rnam,"ARG") == 0) return 6;
  if(strcmp(anam," O  ") == 0 && strcmp(rnam,"ARG") == 0) return 2;
  if(strcmp(anam," CB ") == 0 && strcmp(rnam,"ARG") == 0) return 4;
  if(strcmp(anam," CG ") == 0 && strcmp(rnam,"ARG") == 0) return 4;
  if(strcmp(anam," CD ") == 0 && strcmp(rnam,"ARG") == 0) return 7;
  if(strcmp(anam," NE ") == 0 && strcmp(rnam,"ARG") == 0) return 3;
  if(strcmp(anam," CZ ") == 0 && strcmp(rnam,"ARG") == 0) return 6;
  if(strcmp(anam," NH1") == 0 && strcmp(rnam,"ARG") == 0) return 3;
  if(strcmp(anam," NH2") == 0 && strcmp(rnam,"ARG") == 0) return 3;

  if(strcmp(anam," N  ") == 0 && strcmp(rnam,"ASN") == 0) return 3;
  if(strcmp(anam," CA ") == 0 && strcmp(rnam,"ASN") == 0) return 7;
  if(strcmp(anam," C  ") == 0 && strcmp(rnam,"ASN") == 0) return 6;
  if(strcmp(anam," O  ") == 0 && strcmp(rnam,"ASN") == 0) return 2;
  if(strcmp(anam," CB ") == 0 && strcmp(rnam,"ASN") == 0) return 4;
  if(strcmp(anam," CG ") == 0 && strcmp(rnam,"ASN") == 0) return 6;
  if(strcmp(anam," OD1") == 0 && strcmp(rnam,"ASN") == 0) return 2;
  if(strcmp(anam," ND2") == 0 && strcmp(rnam,"ASN") == 0) return 3;

  if(strcmp(anam," N  ") == 0 && strcmp(rnam,"ASP") == 0) return 3;
  if(strcmp(anam," CA ") == 0 && strcmp(rnam,"ASP") == 0) return 7;
  if(strcmp(anam," C  ") == 0 && strcmp(rnam,"ASP") == 0) return 6;
  if(strcmp(anam," O  ") == 0 && strcmp(rnam,"ASP") == 0) return 2;
  if(strcmp(anam," CB ") == 0 && strcmp(rnam,"ASP") == 0) return 4;
  if(strcmp(anam," CG ") == 0 && strcmp(rnam,"ASP") == 0) return 6;
  if(strcmp(anam," OD1") == 0 && strcmp(rnam,"ASP") == 0) return 2;
  if(strcmp(anam," OD2") == 0 && strcmp(rnam,"ASP") == 0) return 2;

  if(strcmp(anam," N  ") == 0 && strcmp(rnam,"CYS") == 0) return 3;
  if(strcmp(anam," CA ") == 0 && strcmp(rnam,"CYS") == 0) return 7;
  if(strcmp(anam," C  ") == 0 && strcmp(rnam,"CYS") == 0) return 6;
  if(strcmp(anam," O  ") == 0 && strcmp(rnam,"CYS") == 0) return 2;
  if(strcmp(anam," CB ") == 0 && strcmp(rnam,"CYS") == 0) return 4;
  if(strcmp(anam," SG ") == 0 && strcmp(rnam,"CYS") == 0) return 6;

  if(strcmp(anam," N  ") == 0 && strcmp(rnam,"GLN") == 0) return 3;
  if(strcmp(anam," CA ") == 0 && strcmp(rnam,"GLN") == 0) return 7;
  if(strcmp(anam," C  ") == 0 && strcmp(rnam,"GLN") == 0) return 6;
  if(strcmp(anam," O  ") == 0 && strcmp(rnam,"GLN") == 0) return 2;
  if(strcmp(anam," CB ") == 0 && strcmp(rnam,"GLN") == 0) return 4;
  if(strcmp(anam," CG ") == 0 && strcmp(rnam,"GLN") == 0) return 4;
  if(strcmp(anam," CD ") == 0 && strcmp(rnam,"GLN") == 0) return 6;
  if(strcmp(anam," OE1") == 0 && strcmp(rnam,"GLN") == 0) return 2;
  if(strcmp(anam," NE2") == 0 && strcmp(rnam,"GLN") == 0) return 3;

  if(strcmp(anam," N  ") == 0 && strcmp(rnam,"GLU") == 0) return 3;
  if(strcmp(anam," CA ") == 0 && strcmp(rnam,"GLU") == 0) return 7;
  if(strcmp(anam," C  ") == 0 && strcmp(rnam,"GLU") == 0) return 6;
  if(strcmp(anam," O  ") == 0 && strcmp(rnam,"GLU") == 0) return 2;
  if(strcmp(anam," CB ") == 0 && strcmp(rnam,"GLU") == 0) return 4;
  if(strcmp(anam," CG ") == 0 && strcmp(rnam,"GLU") == 0) return 4;
  if(strcmp(anam," CD ") == 0 && strcmp(rnam,"GLU") == 0) return 6;
  if(strcmp(anam," OE1") == 0 && strcmp(rnam,"GLU") == 0) return 2;
  if(strcmp(anam," OE2") == 0 && strcmp(rnam,"GLU") == 0) return 2;

  if(strcmp(anam," N  ") == 0 && strcmp(rnam,"GLY") == 0) return 3;
  if(strcmp(anam," CA ") == 0 && strcmp(rnam,"GLY") == 0) return 6;
  if(strcmp(anam," C  ") == 0 && strcmp(rnam,"GLY") == 0) return 6;
  if(strcmp(anam," O  ") == 0 && strcmp(rnam,"GLY") == 0) return 2;

  if(strcmp(anam," N  ") == 0 && strcmp(rnam,"HIS") == 0) return 3;
  if(strcmp(anam," CA ") == 0 && strcmp(rnam,"HIS") == 0) return 7;
  if(strcmp(anam," C  ") == 0 && strcmp(rnam,"HIS") == 0) return 6;
  if(strcmp(anam," O  ") == 0 && strcmp(rnam,"HIS") == 0) return 2;
  if(strcmp(anam," CB ") == 0 && strcmp(rnam,"HIS") == 0) return 4;
  if(strcmp(anam," CG ") == 0 && strcmp(rnam,"HIS") == 0) return 5;
  if(strcmp(anam," ND1") == 0 && strcmp(rnam,"HIS") == 0) return 1;
  if(strcmp(anam," CD2") == 0 && strcmp(rnam,"HIS") == 0) return 5;
  if(strcmp(anam," CE1") == 0 && strcmp(rnam,"HIS") == 0) return 5;
  if(strcmp(anam," NE2") == 0 && strcmp(rnam,"HIS") == 0) return 1;

  if(strcmp(anam," N  ") == 0 && strcmp(rnam,"ILE") == 0) return 3;
  if(strcmp(anam," CA ") == 0 && strcmp(rnam,"ILE") == 0) return 7;
  if(strcmp(anam," C  ") == 0 && strcmp(rnam,"ILE") == 0) return 6;
  if(strcmp(anam," O  ") == 0 && strcmp(rnam,"ILE") == 0) return 2;
  if(strcmp(anam," CB ") == 0 && strcmp(rnam,"ILE") == 0) return 4;
  if(strcmp(anam," CG1") == 0 && strcmp(rnam,"ILE") == 0) return 4;
  if(strcmp(anam," CG2") == 0 && strcmp(rnam,"ILE") == 0) return 4;
  if(strcmp(anam," CD1") == 0 && strcmp(rnam,"ILE") == 0) return 4;

  if(strcmp(anam," N  ") == 0 && strcmp(rnam,"LEU") == 0) return 3;
  if(strcmp(anam," CA ") == 0 && strcmp(rnam,"LEU") == 0) return 7;
  if(strcmp(anam," C  ") == 0 && strcmp(rnam,"LEU") == 0) return 6;
  if(strcmp(anam," O  ") == 0 && strcmp(rnam,"LEU") == 0) return 2;
  if(strcmp(anam," CB ") == 0 && strcmp(rnam,"LEU") == 0) return 4;
  if(strcmp(anam," CG ") == 0 && strcmp(rnam,"LEU") == 0) return 4;
  if(strcmp(anam," CD1") == 0 && strcmp(rnam,"LEU") == 0) return 4;
  if(strcmp(anam," CD2") == 0 && strcmp(rnam,"LEU") == 0) return 4;

  if(strcmp(anam," N  ") == 0 && strcmp(rnam,"LYS") == 0) return 3;
  if(strcmp(anam," CA ") == 0 && strcmp(rnam,"LYS") == 0) return 7;
  if(strcmp(anam," C  ") == 0 && strcmp(rnam,"LYS") == 0) return 6;
  if(strcmp(anam," O  ") == 0 && strcmp(rnam,"LYS") == 0) return 2;
  if(strcmp(anam," CB ") == 0 && strcmp(rnam,"LYS") == 0) return 4;
  if(strcmp(anam," CG ") == 0 && strcmp(rnam,"LYS") == 0) return 4;
  if(strcmp(anam," CD ") == 0 && strcmp(rnam,"LYS") == 0) return 4;
  if(strcmp(anam," CE ") == 0 && strcmp(rnam,"LYS") == 0) return 7;
  if(strcmp(anam," NZ ") == 0 && strcmp(rnam,"LYS") == 0) return 3;

  if(strcmp(anam," N  ") == 0 && strcmp(rnam,"MET") == 0) return 3;
  if(strcmp(anam," CA ") == 0 && strcmp(rnam,"MET") == 0) return 7;
  if(strcmp(anam," C  ") == 0 && strcmp(rnam,"MET") == 0) return 6;
  if(strcmp(anam," O  ") == 0 && strcmp(rnam,"MET") == 0) return 2;
  if(strcmp(anam," CB ") == 0 && strcmp(rnam,"MET") == 0) return 4;
  if(strcmp(anam," CG ") == 0 && strcmp(rnam,"MET") == 0) return 4;
  if(strcmp(anam," SD ") == 0 && strcmp(rnam,"MET") == 0) return 6;
  if(strcmp(anam," CE ") == 0 && strcmp(rnam,"MET") == 0) return 4;

  if(strcmp(anam," N  ") == 0 && strcmp(rnam,"PHE") == 0) return 3;
  if(strcmp(anam," CA ") == 0 && strcmp(rnam,"PHE") == 0) return 7;
  if(strcmp(anam," C  ") == 0 && strcmp(rnam,"PHE") == 0) return 6;
  if(strcmp(anam," O  ") == 0 && strcmp(rnam,"PHE") == 0) return 2;
  if(strcmp(anam," CB ") == 0 && strcmp(rnam,"PHE") == 0) return 4;
  if(strcmp(anam," CG ") == 0 && strcmp(rnam,"PHE") == 0) return 5;
  if(strcmp(anam," CD1") == 0 && strcmp(rnam,"PHE") == 0) return 5;
  if(strcmp(anam," CD2") == 0 && strcmp(rnam,"PHE") == 0) return 5;
  if(strcmp(anam," CE1") == 0 && strcmp(rnam,"PHE") == 0) return 5;
  if(strcmp(anam," CE2") == 0 && strcmp(rnam,"PHE") == 0) return 5;
  if(strcmp(anam," CZ ") == 0 && strcmp(rnam,"PHE") == 0) return 5;

  if(strcmp(anam," N  ") == 0 && strcmp(rnam,"PRO") == 0) return 6;
  if(strcmp(anam," CA ") == 0 && strcmp(rnam,"PRO") == 0) return 4;
  if(strcmp(anam," C  ") == 0 && strcmp(rnam,"PRO") == 0) return 6;
  if(strcmp(anam," O  ") == 0 && strcmp(rnam,"PRO") == 0) return 2;
  if(strcmp(anam," CB ") == 0 && strcmp(rnam,"PRO") == 0) return 4;
  if(strcmp(anam," CG ") == 0 && strcmp(rnam,"PRO") == 0) return 4;
  if(strcmp(anam," CD ") == 0 && strcmp(rnam,"PRO") == 0) return 4;

  if(strcmp(anam," N  ") == 0 && strcmp(rnam,"SER") == 0) return 3;
  if(strcmp(anam," CA ") == 0 && strcmp(rnam,"SER") == 0) return 7;
  if(strcmp(anam," C  ") == 0 && strcmp(rnam,"SER") == 0) return 6;
  if(strcmp(anam," O  ") == 0 && strcmp(rnam,"SER") == 0) return 2;
  if(strcmp(anam," CB ") == 0 && strcmp(rnam,"SER") == 0) return 6;
  if(strcmp(anam," OG ") == 0 && strcmp(rnam,"SER") == 0) return 1;

  if(strcmp(anam," N  ") == 0 && strcmp(rnam,"THR") == 0) return 3;
  if(strcmp(anam," CA ") == 0 && strcmp(rnam,"THR") == 0) return 7;
  if(strcmp(anam," C  ") == 0 && strcmp(rnam,"THR") == 0) return 6;
  if(strcmp(anam," O  ") == 0 && strcmp(rnam,"THR") == 0) return 2;
  if(strcmp(anam," CB ") == 0 && strcmp(rnam,"THR") == 0) return 6;
  if(strcmp(anam," OG1") == 0 && strcmp(rnam,"THR") == 0) return 1;
  if(strcmp(anam," CG2") == 0 && strcmp(rnam,"THR") == 0) return 4;

  if(strcmp(anam," N  ") == 0 && strcmp(rnam,"TRP") == 0) return 3;
  if(strcmp(anam," CA ") == 0 && strcmp(rnam,"TRP") == 0) return 7;
  if(strcmp(anam," C  ") == 0 && strcmp(rnam,"TRP") == 0) return 6;
  if(strcmp(anam," O  ") == 0 && strcmp(rnam,"TRP") == 0) return 2;
  if(strcmp(anam," CB ") == 0 && strcmp(rnam,"TRP") == 0) return 4;
  if(strcmp(anam," CG ") == 0 && strcmp(rnam,"TRP") == 0) return 5;
  if(strcmp(anam," CD1") == 0 && strcmp(rnam,"TRP") == 0) return 5;
  if(strcmp(anam," CD2") == 0 && strcmp(rnam,"TRP") == 0) return 5;
  if(strcmp(anam," NE1") == 0 && strcmp(rnam,"TRP") == 0) return 3;
  if(strcmp(anam," CE2") == 0 && strcmp(rnam,"TRP") == 0) return 5;
  if(strcmp(anam," CE3") == 0 && strcmp(rnam,"TRP") == 0) return 5;
  if(strcmp(anam," CZ2") == 0 && strcmp(rnam,"TRP") == 0) return 5;
  if(strcmp(anam," CZ3") == 0 && strcmp(rnam,"TRP") == 0) return 5;
  if(strcmp(anam," CH2") == 0 && strcmp(rnam,"TRP") == 0) return 5;

  if(strcmp(anam," N  ") == 0 && strcmp(rnam,"TYR") == 0) return 3;
  if(strcmp(anam," CA ") == 0 && strcmp(rnam,"TYR") == 0) return 7;
  if(strcmp(anam," C  ") == 0 && strcmp(rnam,"TYR") == 0) return 6;
  if(strcmp(anam," O  ") == 0 && strcmp(rnam,"TYR") == 0) return 2;
  if(strcmp(anam," CB ") == 0 && strcmp(rnam,"TYR") == 0) return 4;
  if(strcmp(anam," CG ") == 0 && strcmp(rnam,"TYR") == 0) return 5;
  if(strcmp(anam," CD1") == 0 && strcmp(rnam,"TYR") == 0) return 5;
  if(strcmp(anam," CD2") == 0 && strcmp(rnam,"TYR") == 0) return 5;
  if(strcmp(anam," CE1") == 0 && strcmp(rnam,"TYR") == 0) return 5;
  if(strcmp(anam," CE2") == 0 && strcmp(rnam,"TYR") == 0) return 5;
  if(strcmp(anam," CZ ") == 0 && strcmp(rnam,"TYR") == 0) return 5;
  if(strcmp(anam," OH ") == 0 && strcmp(rnam,"TYR") == 0) return 1;

  if(strcmp(anam," N  ") == 0 && strcmp(rnam,"VAL") == 0) return 3;
  if(strcmp(anam," CA ") == 0 && strcmp(rnam,"VAL") == 0) return 7;
  if(strcmp(anam," C  ") == 0 && strcmp(rnam,"VAL") == 0) return 6;
  if(strcmp(anam," O  ") == 0 && strcmp(rnam,"VAL") == 0) return 2;
  if(strcmp(anam," CB ") == 0 && strcmp(rnam,"VAL") == 0) return 4;
  if(strcmp(anam," CG1") == 0 && strcmp(rnam,"VAL") == 0) return 4;
  if(strcmp(anam," CG2") == 0 && strcmp(rnam,"VAL") == 0) return 4;

  return 0;
}
/***********************************************************************/
/*        1         2         3         4         5         6         7*/
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/*        1         2         3         4         5         6         7*/
/***********************************************************************/
void Random_bin_graph(float fraction){
  int i,j,k;
  int l[]={0,1};

  seed=time_seed();
  printf("seed=%ld\n",seed);

  if(fraction > 0.5){
    fraction = 1.0 - fraction;
    l[0]=1;
    l[1]=0;
  }

  numEdges=fraction*numNodes*(numNodes-1)/2;

  printf("numEdges=%d\n",numEdges);

  for(i=0;i<numNodes;i++){
    for(j=0;j<numNodes;j++){
      k=l[0];
      if(i == j) k=1;
      connected[i*numNodes+j]=k;
    }
  }


  while(numEdges > 0){
    i=ran2(&seed)*numNodes;
    j=ran2(&seed)*numNodes;
    if(i == j) continue;
    if(connected[i*numNodes+j] == l[0]){
      connected[i*numNodes+j]=l[1];
      connected[j*numNodes+i]=l[1];
      numEdges--;
    }
  }

  return;
}

/***********************************************************************/
/*        1         2         3         4         5         6         7*/
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/*        1         2         3         4         5         6         7*/
/***********************************************************************/
int Bron_Kerbosch(){
  int c,i,j;
  int *all = NULL;

  //printf("numNodes=%d\n",numNodes);
  //PAUSE;
  all = malloc(numNodes*sizeof(int));
  if(all == NULL){
    printf("unnable to allocate memory to all\n");
    exit(1);
  }
  compsub = malloc(numNodes*sizeof(int));
  if(compsub == NULL){
    printf("unnable to allocate memory to compsub\n");
    exit(1);
  }
  
  
  for(c=0;c<numNodes;c++) all[c]=c;
  c=0;
  Extend(all,0,numNodes);
  
  return (clique_inum+1);
}
/***********************************************************************/
/*        1         2         3         4         5         6         7*/
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/*        1         2         3         4         5         6         7*/
/***********************************************************************/
 void Extend(int old[],int ne,int ce){
   int *new = NULL;
   int nod,fixp;
   int newne,newce,i,j,count,pos,p,s,sel,minnod;
   int loc;
   int l;
   //printf("entered Extend\n");PAUSE;

  new = malloc(ce*sizeof(int));
  if(new == NULL){
    printf("unnable to allocate memory to new\n");
    exit(1);
  }

  minnod=ce;
  i=0;
  nod=0;
  while(i<ce && minnod != 0){
    p=old[i];
    count=0;
    j=ne;
    while(j < ce && count < minnod){    
      if(connected[p*numNodes+old[j]] == 0){
	count++;
	pos=j;
      }
      j++;
    }
    if(count < minnod){
      fixp=p;
      minnod=count;
      if(i < ne){
	s=pos;
      }else{
	s=i;
	nod=1;
      }
    }
    i++;
  }

  for(nod=minnod+nod;nod>0;nod--){
    p=old[s];
    old[s]=old[ne];
    old[ne]=p;
    sel=p;

    newne=0;
    i=0;
    while(i<ne){
      if(connected[sel*numNodes+old[i]]==1) new[newne++]=old[i];
      i++;
    }

    newce=newne;
    i=ne+1;
    while(i < ce){
      if(connected[sel*numNodes+old[i]] == 1) new[newce++]=old[i];
      i++;
    }
    compsub[c++]=sel;
    if(newce == 0){
      if(c >= Clique_threshold && got_one == 0){
	//printf("adding new clique\n");PAUSE;
	AddNewClique(c, compsub);
	if(get_all == 0) got_one = 1;
      }
    }else if(newne < newce){
      if(got_one == 0) Extend(new,newne,newce);
    }
    c--;
    ne++;
    if(nod > 1){
      s=ne;
      while(connected[fixp*numNodes+old[s]] == 1 && s < numNodes) s++;
    }    
  }

  free(new);
  new=NULL;
  return;
}
/***********************************************************************/
/*        1         2         3         4         5         6         7*/
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/*        1         2         3         4         5         6         7*/
/***********************************************************************/
void AddNewClique(int n, int list[]){
  int l,i,j;
  tClique c;
  int *clist;
  gsl_matrix *mat_r  = NULL;
  double DetR;

  if(n < max_clique_size)  return;

  mat_r  = gsl_matrix_alloc(3,3);
  NEW(c, tsClique);
  clist=malloc(n*sizeof(int));
  if(clist == NULL){
    printf("unable to create clist: OUT OF MEMORY!\n");
    exit(1);
  }
  c->mat_r=mat_r;
  c->Vnum=n;
  memcpy(clist,list,n*sizeof(int));
  c->Vertices=clist;

  Calculate_Rotation_GSL(c);
  //printf("inside add new clque");PAUSE;

  //printf("DetR=%g\n",c->DetR);

  if(fabs(c->DetR - 1.0) > 0.001){
    //printf("Not adding as a clique\n");
    FREE(c);
    gsl_matrix_free(mat_r);
    FREE(clist);
    //PAUSE;
    return;
  }


  //for(i=0;i<c->Vnum;i++){printf("%d ",c->Vertices[i]);}
  //PAUSE;

  /*
    else{
    d=Cliques;
    do{
    min_incommon=0;
    //for(i=0;i<
    d=d->next;
    }while(d != Cliques);
    }
  */

  //printf("Adding as a clique\n");
  //PAUSE;
  
  clique_inum++;

  c->inum = clique_inum;
  c->reject = 0;

  if(c->Vnum > max_clique_size){
    max_clique_size=c->Vnum;
    //printf("max_clique_size=%d\n",max_clique_size);
  }

  ADD(Cliques,c);

  return;
}
/***********************************************************************/
/*        1         2         3         4         5         6         7*/
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/*        1         2         3         4         5         6         7*/
/***********************************************************************/
tClique Get_Top_Clique(){
  tClique c,m,f;
  int i,j,n;
  float coor;
  float min_rmsd=1000.0;

  if(Cliques == NULL) return NULL;

  // Calculate RMSD values for all cliques with max_clique_size 
  c=Cliques;
  do{
    c->rmsd=999.0; // standard value for those with less than max_clique_size nodes
    if(c->Vnum == max_clique_size){
      c->rmsd=0.0; // standard value for those with less than max_clique_size nodes
      for(n=0;n<c->Vnum;n++){	
	for(i=0;i<3;i++){
	  coor=c->cen_b[i];
	  for(j=0;j<3;j++){
	    coor += (agv[c->Vertices[n]].a->coor[j] - c->cen_a[j])*gsl_matrix_get(c->mat_r,i,j);
	  }

	  c->rmsd +=(agv[c->Vertices[n]].b->coor[i] - coor)*(agv[c->Vertices[n]].b->coor[i] - coor);
	}
      }
      c->rmsd = sqrt(c->rmsd/c->Vnum);
      //printf("RMSD=%g\n",c->rmsd);
      //PAUSE;
    }
    c=c->next;
  }while(c != Cliques);


  // find Clique with Vnum = max_clique_size & lowest RMSD
  m=Cliques;
  f=Cliques->next;
  do{
    if(f->Vnum == max_clique_size && f->rmsd < min_rmsd){
      //printf("RMSD m=%5.3f (%d) f=%5.3f(%d)",m->rmsd,m->Vnum,f->rmsd,f->Vnum);PAUSE;
      m=f;
      min_rmsd = f->rmsd;
    }
    f=f->next;
  }while(f != Cliques);
  //Cliques=m;

  //printf("RMSD=%5.3f min=%5.3f Vnum=%d",m->rmsd,min_rmsd,m->Vnum);PAUSE;

  return m;
}
/***********************************************************************/
/*        1         2         3         4         5         6         7*/
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/*        1         2         3         4         5         6         7*/
/***********************************************************************/
int Filter_Cliques(int total){
  tClique c;
  int i,j,n;
  float coor;
  
  total=0;
  if(Cliques == NULL) return total;

 if(Cliques->next != Cliques){
    sort_Cliques();
    non_redundant_Cliques();
  }

  // Calculate RMSD values for accepted cliques
  c=Cliques;
  do{
    if(c->reject == 0){
      total++;
      //printf("inum=%d Vnum=%d reject=%d\n",c->inum,c->Vnum,c->reject);
      c->rmsd=0.0;
      for(n=0;n<c->Vnum;n++){
        for(i=0;i<3;i++){
          coor=c->cen_b[i];
          for(j=0;j<3;j++) coor +=
            (agv[c->Vertices[n]].a->coor[j] - c->cen_a[j])*gsl_matrix_get(c->mat_r,i,j);
          c->rmsd += (agv[c->Vertices[n]].b->coor[i] - coor)*
            (agv[c->Vertices[n]].b->coor[i] - coor);
        }
      }
      c->rmsd = sqrt(c->rmsd/c->Vnum);
      //printf("RMSD=%g\n",c->rmsd);
      //PAUSE;
    }
    c=c->next;
  }while(c != Cliques);

  return total;
}
/***********************************************************************/
/*        1         2         3         4         5         6         7*/
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/*        1         2         3         4         5         6         7*/
/***********************************************************************/
void non_redundant_Cliques(){
  tClique c,d,e;
  int min_incommon_thres=3;
  int min_incommon;
  int found_node;
  int i,j;

  c=Cliques;
  do{
    if(c->reject == 0){
      /*
	printf("Accepted %d (size %d reject=%d):\n",c->inum,c->Vnum,c->reject);
	printf("Next %d (size %d reject=%d):\n",c->next->inum,c->next->Vnum,c->next->reject);
	printf("Prev %d (size %d reject=%d):\n",c->prev->inum,c->prev->Vnum,c->prev->reject);
	for(i=0;i<c->Vnum;i++) printf("%d ",c->Vertices[i]);
	printf("\n");
	PAUSE;
      */
      d=c->next;
      do{
	if(d->reject == 0){
	  //printf("Candidate %d:\n",d->inum);
	  //for(i=0;i<d->Vnum;i++) printf("%d ",d->Vertices[i]);
	  //printf("\n");
	  //PAUSE;
	  min_incommon=0;
	  i=0;
	  while(i<c->Vnum && min_incommon < min_incommon_thres){    
	    found_node=0;
	    j=0;
	    while(found_node == 0 && j<d->Vnum){
	      if(d->Vertices[j] == c->Vertices[i]){
		found_node=1;
		min_incommon++;
	      }
	      j++;
	    }
	    i++;
	  }
	  
	  if(min_incommon >= min_incommon_thres){
	    /*
	      printf("Accepted inum=%d size=%d:\n",c->inum,c->Vnum);
	      for(i=0;i<c->Vnum;i++) printf("%d ",c->Vertices[i]);
	      printf("\n");
	      printf("Rejected inum=%d size=%d:\n",d->inum,d->Vnum);
	      for(i=0;i<d->Vnum;i++) printf("%d ",d->Vertices[i]);
	      printf("\n")
	      };
	    */
	    //printf(" -> rejected\n");
	    d->reject=1;
	    /*
	      e=d;
	      d->prev->next=e->next;
	      d->next->prev=e->prev;
	      d=e->next;
	      free(e->Vertices);
	      gsl_matrix_free(e->mat_r);
	      free(e);
	      //e=NULL;
	      */
	  }
	  //else{
	  //printf(" -> accepted\n");
	  //}
	}
	d=d->next;
      }while(d != c);
      /*
	printf("Accepted %d (size %d reject=%d):\n",c->inum,c->Vnum,c->reject);
	printf("Next %d (size %d reject=%d):\n",c->next->inum,c->next->Vnum,c->next->reject);
	printf("Prev %d (size %d reject=%d):\n",c->prev->inum,c->prev->Vnum,c->prev->reject);
	PAUSE;
      */
    }
    c=c->next;
  }while(c != Cliques);




  return;
}
/***********************************************************************/
/*        1         2         3         4         5         6         7*/
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/*        1         2         3         4         5         6         7*/
/***********************************************************************/
void sort_Cliques(){
  tClique f,m,x;

  // Shift start of circular list to top cleft
  m=Cliques;
  f=Cliques->next;
  do{
    if(f->Vnum > m->Vnum) m=f;
    f=f->next;
  }while(f != Cliques);
  Cliques=m;


  // sort rest by each time finding themaixum in remaining of list and
  // swapping the maximum with current point if maximum > current

  f=Cliques->next;
  do{
    x=f;
    m=f->next;
    do{
      if(m->Vnum > x->Vnum) x=m;
      m=m->next;
    }while(m != Cliques);

    if(x != f){
      swap_Cliques(f,x);
      f=x;
    }
    f=f->next;
  }while(f != Cliques->prev);

  return;
}
/***********************************************************************/
/*        1         2         3         4         5         6         7*/
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/*        1         2         3         4         5         6         7*/
/***********************************************************************/
void swap_Cliques(tClique b, tClique d){
  tClique dprev,dnext,f;
  bool  consec=FALSE;
  
  f=d;
  do{
    if(f == b){
      b=d;
      d=f;
      break;
    } 
    f=f->next;
  }while(f != Cliques);


  if(b->next == d || b->prev == d){consec=TRUE;}

  dprev=d->prev;
  dnext=d->next;

  d->next->prev = d->prev;
  d->prev->next = d->next;

  if(consec==FALSE){
    b->next->prev = b->prev;
    b->prev->next = b->next;
    
    d->prev=b->prev;
    d->next=b->next;
    b->prev->next=d;
    b->next->prev=d;
    
    dprev->next=b;
    dnext->prev=b;
    
    b->next=dnext;
    b->prev=dprev;
  }else{
    b->prev->next=d;
    b->prev=d;
    d->prev=b->prev;
    d->next=b;
  }

  return;
}
