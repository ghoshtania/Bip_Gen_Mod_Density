#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"body.h"
#include"help.h"
#include <time.h>
#include <string.h>  // For memcpy
extern double cp;
void printGraph(struct graph *G);
int findxiny(int x, int *y)
{
    // Check if the pointer y is not NULL
    if (y == NULL) {
        printf("Error: y is NULL!\n");
        return -1;  // or handle the error appropriately
    }

    // Check if the first element of y (y[0]) is valid
    if (y[0] <= 0) {
        printf("Error: Invalid size in y[0]. Size should be positive.\n");
        return -1;  // or handle the error appropriately
    }

    int j;
    // Loop through y[1] to y[y[0]]
    for (j = 1; j <= y[0]; j++) {
        if (y[j] == x) {
            return j;  // Return index if x is found
        }
    }
    
    // Return -1 if x is not found
    return -1;
}

double power(double a, double b) {
    // Edge case: 0 raised to a negative power → return 0 instead of infinity
    if (a == 0.0 && b < 0.0) {
        return 0.0;
    }
    // Normal case
    return pow(a, b);
}

int outpart(struct part ans,double n, int N, char *fname)
{
    int i;
    printf("%d\t%lf\n",ans.com,ans.Q/(n));
	FILE *fo;
	fo=fopen(fname,"w");
    	for (i=0;i<N;i++)
    		fprintf(fo,"%d\n",ans.pa[i]);
	fclose(fo);
    return 0;
}

int updateG(struct graph *G, int x, int y) {
    int i, j,k;
    int *newlist_blue, *newlist_red; //These are new adjacency lists for blue and red nodes, respectively, after merging.
    double *newl2_blue, *newl2_red;//These hold the corresponding weights for connections in newlist_blue and newlist_red.
    int new_size_blue, new_size_red;//These determine the size of the new adjacency and weight lists by adding the sizes of nodes x and y.

    // Determine the sizes for the new adjacency and weight lists
    new_size_blue = G[0].nl_blue[x][0] + G[0].nl_blue[y][0];
    new_size_red = G[0].nl_red[x][0] + G[0].nl_red[y][0];

    newlist_blue = (int*)malloc(sizeof(int) * (new_size_blue + 1));
    newlist_red = (int*)malloc(sizeof(int) * (new_size_red + 1));
    newl2_blue = (double*)malloc(sizeof(double) * (new_size_blue + 1));
    newl2_red = (double*)malloc(sizeof(double) * (new_size_red + 1));

    // Initialize new lists
    newlist_blue[0] = 0;
    newlist_red[0] = 0;
    int direct_connection_redX = 0;//if direct conn. from red node of X cluster
    int direct_connection_blueX = 0;//if direct conn from blue node of X cluster
    k=0;
    newl2_red[0]=G->wl_red[x][0]+G->wl_red[y][0];
    newl2_blue[0]=G->wl_blue[x][0]+G->wl_blue[y][0];
   
//from x group
    //for red node of x group
    if (G[0].nl_red[x][0]!=0)
    {
        for (i=1;i<=G[0].nl_red[x][0];i++)//for all blue neighbors of red node
        {
            if (G[0].nl_red[x][i]==y) 
            {	
                direct_connection_redX = 1 ;k=-1;	newl2_red[0]+=G[0].wl_red[x][i]; newl2_blue[0]+=G[0].wl_red[x][i];  
            }
            else
            {
                newlist_red[i+k]=G[0].nl_red[x][i];
                newl2_red[i+k]=G[0].wl_red[x][i];
            }
        }
    }
    
    k=0;
    //for blue node of x group
    if (G[0].nl_blue[x][0]!=0)
    {
        for (i=1;i<=G[0].nl_blue[x][0];i++)//for all red neighbors of red node
        {
            if (G[0].nl_blue[x][i]==y) 
            { 
                direct_connection_blueX = 1 ;	k=-1;	newl2_blue[0]+=G[0].wl_blue[x][i]; newl2_red[0]+=G[0].wl_blue[x][i]; 
            }
            else{
                newlist_blue[i+k]=G[0].nl_blue[x][i];
                newl2_blue[i+k]=G[0].wl_blue[x][i];
                }
        }
    }
    // only update newlist_blue[0] if there was a direct connection
    if (direct_connection_redX) 
    {
        newlist_red[0] = G[0].nl_red[x][0] - 1;
    }
    else {newlist_red[0] = G[0].nl_red[x][0]; newl2_red[0]=newl2_blue[0];}
    if (direct_connection_blueX) 
    {
        newlist_blue[0] = G[0].nl_blue[x][0] - 1;
    }
    else { newlist_blue[0] = G[0].nl_blue[x][0]; newl2_blue[0] = newl2_red[0];}

//from the y cluster
    //for red node of y cluster
    if (G[0].nl_red[y][0]!=0)
    {
        for (i=1;i<=G[0].nl_red[y][0];i++)
        if (G[0].nl_red[y][i]!=x)
        {
        
            int sg_blue,k1_blue,k2_blue;
            sg_blue=G[0].nl_red[y][i];
            k1_blue=findxiny(x,G[0].nl_blue[sg_blue]);// the index of x in the adjacency list of s
            k2_blue=findxiny(y,G[0].nl_blue[sg_blue]);//the index of y in the adjacency list of sg
            if (k1_blue>0)
            {
                G[0].wl_blue[sg_blue][k1_blue]+=G[0].wl_blue[sg_blue][k2_blue];
                int test2;
                test2=findxiny(sg_blue,newlist_red);
                if (test2>0)
                newl2_red[test2]=G[0].wl_blue[sg_blue][k1_blue];
                else
                {
                    //printf("Error!\n");
                    exit(0);
                }
                if (k2_blue!=G[0].nl_blue[sg_blue][0])
                {
                    G[0].nl_blue[sg_blue][k2_blue]=G[0].nl_blue[sg_blue][G[0].nl_blue[sg_blue][0]];
                    G[0].wl_blue[sg_blue][k2_blue]=G[0].wl_blue[sg_blue][G[0].nl_blue[sg_blue][0]];
                }
                G[0].nl_blue[sg_blue][0]--;
            }
            else
            {
                if (k2_blue>0)
                G[0].nl_blue[sg_blue][k2_blue]=x;
                else
                {
                    //printf("Error1!\n");
                    exit (0);
                }
                newlist_red[0]+=1;
                newlist_red[newlist_red[0]]=sg_blue;
                newl2_red[newlist_red[0]]=G[0].wl_blue[sg_blue][k2_blue];
            } 


        }
    }
    //for blue node of y cluster
    if (G[0].nl_blue[y][0]!=0)
    {
        for (i=1;i<=G[0].nl_blue[y][0];i++)
        if (G[0].nl_blue[y][i]!=x)
        {
            int sg_red,k1_red,k2_red;
            sg_red=G[0].nl_blue[y][i];
            k1_red=findxiny(x,G[0].nl_red[sg_red]);// the index of x in the adjacency list of sg
            k2_red=findxiny(y,G[0].nl_red[sg_red]);//the index of y in the adjacency list of sg
            if (k1_red>0)
            {
                G[0].wl_red[sg_red][k1_red]+=G[0].wl_red[sg_red][k2_red];
                int test_red;
                test_red=findxiny(sg_red,newlist_blue);
                if (test_red>0)
                newl2_blue[test_red]=G[0].wl_red[sg_red][k1_red];
                else
                {
                    //printf("Error2!\n");
                    exit(0);
                }
                if (k2_red!=G[0].nl_red[sg_red][0])
                {
                    G[0].nl_red[sg_red][k2_red]=G[0].nl_red[sg_red][G[0].nl_red[sg_red][0]];
                    G[0].wl_red[sg_red][k2_red]=G[0].wl_red[sg_red][G[0].nl_red[sg_red][0]];
                }
                G[0].nl_red[sg_red][0]--;
            }
            else
            {
                if (k2_red>0)
                G[0].nl_red[sg_red][k2_red]=x;
                else
                {
                    //printf("Error4!\n");
                    exit (0);
                }
                newlist_blue[0]+=1;
                newlist_blue[newlist_blue[0]]=sg_red;
                newl2_blue[newlist_blue[0]]=G[0].wl_red[sg_red][k2_red];
            } 
        }
    }
    free(G[0].nl_red[x]);	free(G[0].wl_red[x]);
    free(G[0].nl_blue[x]);	free(G[0].wl_blue[x]);
    G->nl_red[x]=newlist_red;
    G->wl_red[x]=newl2_red;
    G->nl_blue[x]=newlist_blue;
    G->wl_blue[x]=newl2_blue;
    free(G[0].nl_red[y]);	free(G[0].wl_red[y]);
    G->nl_red[y]=NULL;	G->wl_red[y]=NULL;

    free(G[0].nl_blue[y]);	free(G[0].wl_blue[y]);
    G->nl_blue[y]=NULL;	G->wl_blue[y]=NULL;
    printf ("");
    if (y!=G[0].com-1)
    {
        if (G[0].nl_red[G[0].com-1][0]!=0)
            {   for (i=1;i<=G[0].nl_red[G[0].com-1][0];i++)
                {
                int last_blue,k2a;
                last_blue=G[0].nl_red[G[0].com-1][i];//it is a blue node
                k2a=findxiny(G[0].com-1,G[0].nl_blue[last_blue]);
                G[0].nl_blue[last_blue][k2a]=y;
                }
            }
        
        if (G[0].nl_blue[G[0].com-1][0]!=0)
            {
                for (i=1;i<=G[0].nl_blue[G[0].com-1][0];i++)
                {
                int last_red,k2b;
                if (G[0].nl_blue == NULL || G[0].nl_blue[G[0].com - 1] == NULL) {
                    printf("ERROR: nl_blue is not allocated properly!\n");
                }
                last_red=G[0].nl_blue[G[0].com-1][i];//it is a red node
                fflush(stdout);
                k2b=findxiny(G[0].com-1,G[0].nl_red[last_red]);
                G[0].nl_red[last_red][k2b]=y;
                }
            }
    G->nl_blue[y]=G->nl_blue[G->com-1];
    G->wl_blue[y]=G->wl_blue[G->com-1];
    G->nl_red[y]=G->nl_red[G->com-1];
    G->wl_red[y]=G->wl_red[G->com-1];
    G->nl_blue[G->com-1]=NULL;
    G->wl_blue[G->com-1]=NULL;
    G->nl_red[G->com-1]=NULL;
    G->wl_red[G->com-1]=NULL;
    
    }
    
    G[0].dl_red[x]+=G[0].dl_red[y];
    G[0].dl_blue[x]+=G[0].dl_blue[y];
    if (y!=G[0].com-1)
    {
        G[0].dl_red[y]=G[0].dl_red[G[0].com-1];
        G[0].dl_blue[y]=G[0].dl_blue[G[0].com-1];
    }
    G->ml_red[x][0]+=G->ml_red[y][0];
    G->ml_blue[x][0]+=G->ml_blue[y][0];

    if (y!=G->com-1)
    {
	G->ml_red[y][0]=G->ml_red[G->com-1][0];
    G->ml_blue[y][0]=G->ml_blue[G->com-1][0];
    }
    
    return 0;
}

void extrafG(struct part ans,struct graph G,int *s)
{
    int i,j;

    for (i=0;i<G.com;i++)
    {
        if (G.ml_red[i]!=0)
            {
                for (j=1;j<=G.ml_red[i][0];j++)
                    {   if (G.ml_red[i][j]<=G.N&&G.ml_red[i][j]>=0)
                            ans.pa[G.ml_red[i][j]]=s[i]+1;
                        else
                        {
                            printf("wrong\n");
                            exit(1);
                        }
                    }
            }
        if (G.ml_blue[i]!=0)
            {
                for (j=1;j<=G.ml_blue[i][0];j++)
                    {   if (G.ml_blue[i][j]<=G.N&&G.ml_blue[i][j]>=0)
                        ans.pa[G.ml_blue[i][j]]=s[i]+1;
                        else
                        {
                            printf("extrafG wrong when memberadsadasdaddasd\n");
                            exit(1);
                        }

                    }
                
            }
        else
            {
            printf("extrafG wrong because null member listasdasdas\n");
            exit(1);
            }

    }

}


void gcopy(struct graph G,struct graph *Ga)
{
   	Ga->com=G.com; Ga->n=G.n; Ga->N_red = G.N_red; Ga->N_blue = G.N_blue;
    int i,j;

    Ga->ml_red=(int**)malloc(sizeof(int*)*G.com);
    Ga->nl_red=(int**)malloc(sizeof(int*)*G.com);
    Ga->wl_red=(double**)malloc(sizeof(double*)*G.com);
    Ga->dl_red=(double*)malloc(sizeof(double)*G.com);

    Ga->ml_blue=(int**)malloc(sizeof(int*)*G.com);
    Ga->nl_blue=(int**)malloc(sizeof(int*)*G.com);
    Ga->wl_blue=(double**)malloc(sizeof(double*)*G.com);
    Ga->dl_blue=(double*)malloc(sizeof(double)*G.com);
    
    for (i=0;i<G.com;i++)
    {
        if (G.nl_red[i]!=NULL)
        {
        ///red nodes..
        Ga->ml_red[i]=(int*)malloc(sizeof(int));
        Ga[0].ml_red[i][0]=G.ml_red[i][0];

        Ga->nl_red[i]=(int*)malloc(sizeof(int)*(G.nl_red[i][0]+1));
        for (j=0;j<=G.nl_red[i][0];j++)
        Ga->nl_red[i][j]=G.nl_red[i][j];
        
        Ga[0].wl_red[i]=(double*)malloc(sizeof(double)*(G.nl_red[i][0]+1));
        for (j=0;j<=G.nl_red[i][0];j++)
        Ga[0].wl_red[i][j]=G.wl_red[i][j];
  
        Ga[0].dl_red[i]=G.dl_red[i];
        }
        
        if (G.nl_blue[i]!=NULL)
        {
        //blue nodes
        Ga->ml_blue[i]=(int*)malloc(sizeof(int));
        Ga[0].ml_blue[i][0]=G.ml_blue[i][0];

        Ga->nl_blue[i]=(int*)malloc(sizeof(int)*(G.nl_blue[i][0]+1));
        for (j=0;j<=G.nl_blue[i][0];j++)
        Ga->nl_blue[i][j]=G.nl_blue[i][j];
        
        Ga[0].wl_blue[i]=(double*)malloc(sizeof(double)*(G.nl_blue[i][0]+1));
        for (j=0;j<=G.nl_blue[i][0];j++)
        Ga[0].wl_blue[i][j]=G.wl_blue[i][j];
  
        Ga[0].dl_blue[i]=G.dl_blue[i];
        }
        
    }
}

void printGraph(struct graph *G) 
{
    printf("Graph structure:\n");
    printf("G.N_red %d: ", G->N_red);
    printf("G.N_blue %d: ", G->N_blue);
    printf("G.n: %f\n: ", G->n);
    printf("G.com %d\n: ", G->com);
    for (int i = 0; i < G->com; i++) 
    {
        {
            printf("RedNode %d: ", i);
            printf("Degree of red node: %f ", G->dl_red[i]);
            printf("Neighborlist of red node %d: ", i);
            for (int j = 0; j <= G->nl_red[i][0]; j++) 
            {
                printf("%d ", G->nl_red[i][j]);
            }
            printf("weightlist of red node %d: ", i);
            for (int j = 0; j <= G->nl_red[i][0]; j++) 
            {
                printf("%f ", G->wl_red[i][j]);
            }
    
            printf("\n");
            printf("BlueNode %d: ", i);
            printf("Degree of blue node: %f ", G->dl_blue[i]);
            printf("\n");
            printf("Neighborlist of blue node %d: ", i);
                for (int j = 0; j <= G->nl_blue[i][0]; j++) 
                {
                    printf("%d ", G->nl_blue[i][j]);
                }
            printf("\n");
            printf("weightlist of blue node %d: ", i);
            for (int j = 0; j <= G->nl_blue[i][0]; j++) 
            {
                printf("%f ", G->wl_blue[i][j]);

            }
            printf("\n");
        }
    }
}


void freeG(struct graph G)
{
    int i;
    for (i=0;i<G.com;i++)
    {

        free(G.nl_red[i]);
        free(G.wl_red[i]);
		free(G.nl_blue[i]);
        free(G.wl_blue[i]);
		free(G.ml_blue[i]);
		free(G.ml_red[i]);
    }
    free(G.dl_red);
	free(G.dl_blue);

    free(G.nl_red);
	free(G.nl_blue);
    free(G.wl_red);
	free(G.wl_blue);
	free(G.ml_red);
	free(G.ml_blue);
}


double movedQ(struct graph G, int child, int home, int school, int *s, double *deg_blue, double *deg_red, int *rednodes, int *bluenodes, double *links,double *links1)
{   int i, j, mem;    // child is the node, home is the assigned community, school is the new community
    double ans = 0;
    double L1t_r=0, L1t_b=0; //L1t_r(L1t_b)=total weight b/w red(blue) child and blue(red) nodes in home
    double L2t_r=0, L2t_b=0; //L2t_r(L2t_b)=total weight b/w red(blue) child and blue(red) nodes in school

    //red part of child
    for (i = 1; i <= G.nl_red[child][0]; i++)
    {
        mem = s[G.nl_red[child][i]];
        if (mem == home)
            L1t_r+=G.wl_red[child][i];
        else if (mem == school)
            L2t_r+=G.wl_red[child][i];
        
    }
    
    //blue part of child
    for (i = 1; i <= G.nl_blue[child][0]; i++)
    {
        mem = s[G.nl_blue[child][i]];
        if (mem == home)
            L1t_b+=G.wl_blue[child][i];
        else if (mem == school)
            L2t_b+=G.wl_blue[child][i];
        
    }
    double q1,q2,q1p,q2p,d1,d2,d1p,d2p;
    q1=links[home]-deg_blue[home]*deg_red[home]/(G.n); //normal modularity of home
	q2=links[school]-deg_blue[school]*deg_red[school]/(G.n); //normal modularity of school
    if (rednodes[home] == 0 || bluenodes[home] == 0)
		d1=0;
	else 
		d1=(links[home])/rednodes[home]/bluenodes[home];
	if (rednodes[school] == 0 || bluenodes[school] == 0)
		d2=0;
	else
		d2=(links[school])/rednodes[school]/bluenodes[school];
    double home_first, home_second;
    //taking care of red and blue part together
    if (G.wl_red[child][0]!=G.wl_blue[child][0]) printf("break1\n");
    home_first=links[home]-(L1t_r+L1t_b)-G.wl_red[child][0];
    home_second=(deg_blue[home]-G.dl_blue[child])*(deg_red[home]-G.dl_red[child]);
    if ((rednodes[home]-G.ml_red[child][0])==0 || (bluenodes[home]-G.ml_blue[child][0])==0)
        d1p=0;
    else
        d1p = (links[home]-(L1t_r+L1t_b)-G.wl_red[child][0])/(rednodes[home]-G.ml_red[child][0])/(bluenodes[home]-G.ml_blue[child][0]);


  
    q1p = home_first-home_second/G.n;
    double school_first, school_second;
    school_first=links[school]+L2t_r+L2t_b+G.wl_red[child][0];
    school_second=(deg_blue[school]+G.dl_blue[child])*(deg_red[school]+G.dl_red[child]);
    if ((rednodes[school]+G.ml_red[child][0])==0 || (bluenodes[school]+G.ml_red[child][0])==0)
        d2p=0;
    else
        d2p = (links[school]+G.wl_red[child][0]+L2t_r+L2t_b)/(rednodes[school]+G.ml_red[child][0])/(bluenodes[school]+G.ml_blue[child][0]);
    q2p = school_first-school_second/G.n;
    ans=q1p*power(d1p,cp)+q2p*power(d2p,cp)-q1*power(d1,cp)-q2*power(d2,cp);
    return  (ans);
}




struct part RG(struct graph G, int ke, unsigned int *seed)
{
    int i;
    struct part ans;
    struct graph Ga;
	
    gcopy(G,&Ga);
    ans.pa=(int*)malloc(sizeof(int)*G.N);;

    int *s,*ns,kee;
    double *deg_red,*deg_blue;
	s=(int*)malloc(sizeof(int)*G.com);
	ns=(int*)malloc(sizeof(int)*G.com);
	deg_red=(double*)malloc(sizeof(double)*G.com);
    deg_blue=(double*)malloc(sizeof(double)*G.com);
    double Q;
    Q=compQ(Ga);
    double dQmax,dQ;
    int j,k,g1,g2,lr,randg,color;
    double Q1,Q2,den1,den2,den12, first_term_dQ;
    int temp_blue,temp_red;

    for (i=0;i<G.com;i++)
	{
		s[i]=i; //Here, each node starts in its own community (s[i] = i) and ns is initially the same as s
		ns[i]=i;
	 	deg_red[i]=G.dl_red[i];   //deg[i] is initialized to the degree of the node i from the input graph G
        deg_blue[i]=G.dl_blue[i];   //deg[i] is initialized to the degree of the node i from the input graph G
    }

    ans.com=G.com;
    ans.Q=Q;
    // Community Merging Loop: G.com = initial number of communities (each node is its own community initially).
    int maxiter=G.com-1;
    
    for (i=1;i<=maxiter;i++)
        {
            dQmax=-2*G.n-1;
            lr=0; //lr is a counter for ties in modularity gain.
            if (Ga.com<ke)
                kee=Ga.com;
            else
                kee=ke;
            for (j=0;j<kee;j++)
            {
                randg=randint(Ga.com,seed);
                color=randcolor(seed);
                if (color==0 && Ga.nl_red[randg][0]!=0) 
                    {
                        for (k=1;k<=Ga.nl_red[randg][0];k++)//when blue node is picked up
                        {
                            temp_blue = Ga.nl_red[randg][k];
                            Q1=Ga.wl_red[randg][0]-Ga.dl_red[randg]*Ga.dl_blue[randg]/(double)(Ga.n);
                            Q2=Ga.wl_red[temp_blue][0]-Ga.dl_red[temp_blue]*Ga.dl_blue[temp_blue]/(double)(Ga.n);
                            if (Ga.ml_red[randg][0] == 0 || Ga.ml_blue[randg][0] == 0) 
                                {den1 = 0;} 
                            else 
                                {
                                // Calculate density using the formula: den = m_c / (n_r * n_b)
                                den1 = Ga.wl_red[randg][0] / (Ga.ml_red[randg][0] * Ga.ml_blue[randg][0]);
                                }
                            if (Ga.ml_red[temp_blue][0] == 0 || Ga.ml_blue[temp_blue][0] == 0) 
                                {den2 = 0;} 
                            else {
                                // Calculate density using the formula: den = m_c / (n_r * n_b)
                                den2 = Ga.wl_red[temp_blue][0] / (Ga.ml_red[temp_blue][0] * Ga.ml_blue[temp_blue][0]);
                                }
                            int m1;
                            if (Ga.nl_blue[randg] == NULL || Ga.nl_blue[randg][0] <= 0) {
                                m1 = -1;  // Assign -1 if the array is empty or invalid
                            } else {
                                m1=findxiny(Ga.nl_red[randg][k],Ga.nl_blue[randg]);
                            }
                            if (m1>0) 
                            {
                                den12=(Ga.wl_red[randg][0]+Ga.wl_red[temp_blue][0]+Ga.wl_red[randg][k]+Ga.wl_blue[randg][m1])/(Ga.ml_red[randg][0]+Ga.ml_red[temp_blue][0])/(Ga.ml_blue[randg][0]+Ga.ml_blue[temp_blue][0]);
                                first_term_dQ = Ga.wl_red[randg][k]+Ga.wl_blue[randg][m1]-(Ga.dl_red[randg]*Ga.dl_blue[Ga.nl_red[randg][k]]+Ga.dl_blue[randg]*Ga.dl_red[Ga.nl_red[randg][k]])/(double)(Ga.n);
                            }
                            else
                            {
                                den12=(Ga.wl_red[randg][0]+Ga.wl_red[temp_blue][0]+Ga.wl_red[randg][k])/(Ga.ml_red[randg][0]+Ga.ml_red[temp_blue][0])/(Ga.ml_blue[randg][0]+Ga.ml_blue[temp_blue][0]);
                                first_term_dQ = Ga.wl_red[randg][k]-(Ga.dl_red[randg]*Ga.dl_blue[Ga.nl_red[randg][k]]+Ga.dl_blue[randg]*Ga.dl_red[Ga.nl_red[randg][k]])/(double)(Ga.n);   
                            }
                            dQ = first_term_dQ*power(den12,cp)+Q1*(power(den12,cp)-power(den1,cp))+Q2*(power(den12,cp)-power(den2,cp));
                            if (abso(dQ-dQmax)<Inf_Sma)
                            {
                                lr++; 
                                if (rand1(seed)*lr<1)
                                {
                                    g1=randg;	g2=Ga.nl_red[randg][k];
                                }
                            }
                            else if (dQ>dQmax)
                            { 
                                lr=1;
                                dQmax=dQ;
                                g1=randg;	g2=Ga.nl_red[randg][k];
                            }
                        }  
                    } 
                else if (color==1 && Ga.nl_blue[randg][0]!=0) 
                    {
                        for (k=1;k<=Ga.nl_blue[randg][0];k++)//when blue node is picked up
                        {
                            temp_red = Ga.nl_blue[randg][k];
                            Q1=Ga.wl_blue[randg][0]-Ga.dl_blue[randg]*Ga.dl_red[randg]/(double)(Ga.n);
                            Q2=Ga.wl_blue[temp_red][0]-Ga.dl_blue[temp_red]*Ga.dl_red[temp_red]/(double)(Ga.n);
                            if (Ga.ml_blue[randg][0] == 0 || Ga.ml_red[randg][0] == 0) {
                                den1 = 0;} 
                            else {
                                // Calculate density using the formula: den = m_c / (n_r * n_b)
                                den1 = Ga.wl_blue[randg][0] / (Ga.ml_blue[randg][0] * Ga.ml_red[randg][0]);
                                }
                            if (Ga.ml_blue[temp_red][0] == 0 || Ga.ml_red[temp_red][0] == 0) {
                                den2 = 0;} 
                            else {
                                // Calculate density using the formula: den = m_c / (n_r * n_b)
                                den2 = Ga.wl_red[temp_red][0] / (Ga.ml_red[temp_red][0] * Ga.ml_blue[temp_red][0]);
                                } 
                            int m2;
                            if (Ga.nl_red[randg] == NULL || Ga.nl_red[randg][0] <= 0) {
                                m2 = -1;  // Assign -1 if the array is empty or invalid
                            } else {
                                m2 = findxiny(Ga.nl_blue[randg][k], Ga.nl_red[randg]);
                            }
                            if (m2>0) 
                            {
                                den12=(Ga.wl_blue[randg][0]+Ga.wl_blue[temp_red][0]+Ga.wl_blue[randg][k]+Ga.wl_red[randg][m2])/(Ga.ml_blue[randg][0]+Ga.ml_blue[temp_red][0])/(Ga.ml_red[randg][0]+Ga.ml_red[temp_red][0]);
                                first_term_dQ = Ga.wl_blue[randg][k]+Ga.wl_red[randg][m2]-(Ga.dl_blue[randg]*Ga.dl_red[Ga.nl_blue[randg][k]]+Ga.dl_red[randg]*Ga.dl_blue[Ga.nl_blue[randg][k]])/(double)(Ga.n);
                            }
                            else
                            {
                                den12=(Ga.wl_blue[randg][0]+Ga.wl_blue[temp_red][0]+Ga.wl_blue[randg][k])/(Ga.ml_blue[randg][0]+Ga.ml_blue[temp_red][0])/(Ga.ml_red[randg][0]+Ga.ml_red[temp_red][0]);
                                first_term_dQ = Ga.wl_blue[randg][k]-(Ga.dl_blue[randg]*Ga.dl_red[Ga.nl_blue[randg][k]]+Ga.dl_red[randg]*Ga.dl_blue[Ga.nl_blue[randg][k]])/(double)(Ga.n);
                            }
                            dQ = first_term_dQ*power(den12,cp)+Q1*(power(den12,cp)-power(den1,cp))+Q2*(power(den12,cp)-power(den2,cp));
                            if (abso(dQ-dQmax)<Inf_Sma)
                            {
                                lr++; 
                                if (rand1(seed)*lr<1)
                                {
                                    g1=randg;	g2=Ga.nl_blue[randg][k];
                                }
                            }
                        
                            else if (dQ>dQmax)
                            {  
                                lr=1;
                                dQmax=dQ;
                                g1=randg;	g2=Ga.nl_blue[randg][k];
                            }
                            
                        }
                    }
            
                else 
                    {j--; 
                    continue;
                    }
                
            }
           
            if ((color == 0 && Ga.nl_red[randg][0] != 0) || (color == 1 && Ga.nl_blue[randg][0] != 0))
                {
                    if (lr==0) break;
                    else	
                        {
                            if ((Ga.nl_red[g1][0]+Ga.nl_blue[g1][0])<(Ga.nl_red[g2][0]+Ga.nl_blue[g2][0]))
                                {
                                    g1=g1+g2;	g2=g1-g2;	g1=g1-g2;
                                }    	
                            if (Ga.com==2) break;
                            updateG(&Ga,g1,g2);

                            Ga.com--;
                    
                            for (j=0;j<G.com;j++)
                                if (s[j]==g2)
                                    s[j]=g1;
                            if (g2!=Ga.com)
                            for (j=0;j<G.com;j++)
                                if (s[j]==Ga.com)
                                    s[j]=g2;
                        }
                    Q+=dQmax;
                    if (Q>ans.Q)
                        {
                            ans.com=Ga.com;
                            ans.Q=Q;
                            for (j=0;j<G.com;j++) ns[j]=s[j];
                            for (j=0;j<Ga.com;j++) 
                            {
                                deg_red[j]=Ga.dl_red[j];
                                deg_blue[j]=Ga.dl_blue[j];
                            }
                        }
                }     
        }
    int *rednodes, *bluenodes;
    double *linksfromRed, *linksfromBlue;
    rednodes=(int*)malloc(sizeof(int)*ans.com);
    bluenodes=(int*)malloc(sizeof(int)*ans.com);
    linksfromRed=(double*)malloc(sizeof(double)*ans.com);
    linksfromBlue=(double*)malloc(sizeof(double)*ans.com);

    for (i=0;i<ans.com;i++)
        {
            rednodes[i]=0; bluenodes[i]=0;	linksfromRed[i]=0; linksfromBlue[i]=0;
        }
    
    
    for (i=0;i<G.com;i++)
        {
            rednodes[ns[i]]+=G.ml_red[i][0]; // Update red node count for the community
            
            bluenodes[ns[i]]+=G.ml_blue[i][0]; // Update red node count for the community
            for (j=1;j<=G.nl_red[i][0];j++)
                {   
                    if (ns[i]==ns[G.nl_red[i][j]])
                    {
                        linksfromRed[ns[i]]+=G.wl_red[i][j]; // Update link weight from red node for the community
                    }
                }
            for (j=1;j<=G.nl_blue[i][0];j++)
                if (ns[i]==ns[G.nl_blue[i][j]])
                    linksfromBlue[ns[i]]+=G.wl_blue[i][j]; // Update link weight for the community
        }
    for (i=0;i<G.com;i++)
    {
        linksfromRed[ns[i]]+=G.wl_red[i][0];
        linksfromBlue[ns[i]]+=G.wl_blue[i][0];
    }
    struct part ans_before_refine = ans;  
    int *ns_before = malloc(sizeof(int)*G.com);
    memcpy(ns_before, ns, sizeof(int)*G.com);
    time_t refine_start = time(NULL);
    double refine_timeout = 20.0;  
        //refine section
        int change=1; //change is a flag to indicate if any changes were made during an iteration.

    while (change)
        { 
        if (difftime(time(NULL), refine_start) > refine_timeout) 
        {
            ans = ans_before_refine;
            memcpy(ns, ns_before, sizeof(int)*G.com);
            
            free(ns_before);
            ns_before = NULL; 
            break;
        }
            change=0;
            for (i=0; i<G.com; i++)
                {   
                    g1=ns[i];
                    dQmax=-1;
                    for (j=0; j<ans.com; j++)
                        if (j!=g1)
                        {
                            dQ=movedQ(G,i,g1,j,ns,deg_blue,deg_red,rednodes,bluenodes,linksfromBlue,linksfromRed);
                            if (dQ>dQmax)
                                {
                                    dQmax=dQ;
                                    g2=j;
                                }
                        }
                    if (dQmax>Inf_Sma)
                    {
                        ans.Q+=dQmax;
                        deg_red[g1]-=G.dl_red[i];
                        deg_blue[g1]-=G.dl_blue[i];//Adjust the degree totals for the old community g1 and the new community g2.
                        deg_red[g2]+=G.dl_red[i];
                        deg_blue[g2]+=G.dl_blue[i];

                        rednodes[g1]-=G.ml_red[i][0];
                        bluenodes[g1]-=G.ml_blue[i][0];
                        rednodes[g2]+=G.ml_red[i][0];
                        bluenodes[g2]+=G.ml_blue[i][0];

                        for (j=1;j<=G.nl_red[i][0];j++)
                            if (ns[G.nl_red[i][j]]==g1)
                                {
                                    linksfromRed[g1]-=G.wl_red[i][j];
                                    linksfromBlue[g1]-=G.wl_red[i][j];
                                }
                            else if (ns[G.nl_red[i][j]]==g2)
                                {
                                    linksfromRed[g2]+=G.wl_red[i][j];
                                    linksfromBlue[g2]+=G.wl_red[i][j];
                                }
                        linksfromBlue[g1]-=G.wl_red[i][0];
                        linksfromBlue[g2]+=G.wl_red[i][0];
                        for (j=1;j<=G.nl_blue[i][0];j++)
                            if (ns[G.nl_blue[i][j]]==g1)
                                {
                                    linksfromRed[g1]-=G.wl_blue[i][j];
                                    linksfromBlue[g1]-=G.wl_blue[i][j];
                                }
                            else if (ns[G.nl_blue[i][j]]==g2)
                                {
                                    linksfromRed[g2]+=G.wl_blue[i][j];
                                    linksfromBlue[g2]+=G.wl_blue[i][j];
                                }
                        
                        linksfromRed[g1]-=G.wl_blue[i][0];
                        linksfromRed[g2]+=G.wl_blue[i][0];
                        ns[i]=g2;
                        change=1;
                        }
                    

                }

        }



free(linksfromBlue);
free(linksfromRed);
free(rednodes);
free(bluenodes);
    
extrafG(ans,G,ns);

//below is to reorder ans.partition
	j = 0;
	int *lab;
	lab=(int*)malloc(sizeof(int)*G.N);
	for (i=0;i<G.N;i++)
		lab[i]=0;
	for (i = 0;i < G.N;i++)
	{
		if (!lab[i]) //lab[i]==0 Checks if the node i hasn't been assigned a new label yet.
		{
			j++;
			lab[i]=j;
		
		for (k = i + 1;k < G.N;k++)
			if (ans.pa[k] == ans.pa[i])
				lab[k] = lab[i];
		}
	}
	for (i=0;i<G.N;i++)
		{
            ans.pa[i]=lab[i];
        }

	ans.com=j;
	free(lab);
	free(s);
	free(ns);
	free(deg_blue);
    free(deg_red);
    freeG(Ga);
    
    return ans;
}
