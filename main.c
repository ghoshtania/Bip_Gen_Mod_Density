#include <stdlib.h>
#include <time.h>
#include "help.h"
#include "rg.h"
#include <stdio.h>
#include <unistd.h>
#include <sys/wait.h>
#include <signal.h>
#define Tsmall 1e-28
#define TIME_LIMIT 5  // seconds

int *cols;
char *filename;

void inputGraph(struct graph *G);
int cleansort(struct part *ensem,int kmax,int N);
int getscore(int *score,int *edgecc, int N,int know,int *ref);
void renorm(struct graph G, int *ref);
int mergroup(struct graph *G, int x, int y);
void remolast(struct part *ensemble, int ens, int *edgecc, int N)
{
	long i;
	long j;
	for (i=0;i<N;i++)		
		cols[i]=0;
	
	for (i = 0;i<N-1;i++)
		if (!cols[i])
		{
			cols[i]=1;
			for (j = i+1;j<N;j++)
			{
				if (ensemble[0].pa[i] == ensemble[0].pa[j])
				edgecc[i*N + j]--;
				if (edgecc[i*N+j]==ens-1)
				cols[j]=1;
			}
		
		}
	//delete in ensemble
	i = 0;
	free(ensemble[0].pa);
	for (i = 0;i<ens - 1;i++)
	{
		ensemble[i] = ensemble[i + 1];
	}
	ensemble[ens - 1].pa = NULL;
}


int comps(int *p, int *q, int N)
{
    int i;
	for (i = 0;i < N;i++)
		if (p[i] != q[i])
			return 0;
    
    return 1;
}
void replaceone(struct part* ensemble,int know, struct part* inensem, int max, int *edgecc, int N)
{
	long i,j;
	for (i=0;i<N;i++)	cols[i]=0;
	for (i = 0;i<N-1;i++)
		if (!cols[i])
		{
			cols[i]=1;
			for (j = i+1;j<N;j++)
			{
				if (inensem[max].pa[i] == inensem[max].pa[j])
				edgecc[i*N + j]++;
				if (ensemble[0].pa[i]==ensemble[0].pa[j])
				edgecc[i*N+j]--;
				if (edgecc[i*N+j]==know)
				cols[j]=1;
			}
		}
	//now we need to add p' best at the ordering position
	j = -1;
	for (i = 1;i<know;i++)
		if (ensemble[i].Q>inensem[max].Q)
		{
			j = i;    break;
		}
	if (j<0)
	{
		free(ensemble[0].pa);
		for (i=1;i<know;i++)
			ensemble[i-1]=ensemble[i];
		ensemble[know-1]=inensem[max];
	}
	else
	{
		free(ensemble[0].pa);
		for (i = 1;i < j;i++)
			ensemble[i - 1] = ensemble[i];
		ensemble[j-1]=inensem[max];
	}
}

void addone(struct part *ensemble, int know, struct part *inensem, int max, int *edgecc, int N)
{
	long i, j;
	for (i=0;i<N;i++)		cols[i]=0;

	for (i = 0;i<N-1;i++)
		if (!cols[i])
		{
			cols[i]=1;
			for (j = i+1;j<N;j++)
			{
				if (inensem[max].pa[i] == inensem[max].pa[j])
				edgecc[i*N + j]++;
				if (edgecc[i*N+j]==know+1)
				cols[j]=1;
			}
		}

	j = -1;
	for (i = 1;i<know;i++)
		if (ensemble[i].Q>inensem[max].Q)
		{
			j = i;    break;
		}
	if (j<0)
	{
		ensemble[know]=inensem[max];
	}
	else
	{
		for (i = know - 1;i >= j;i--)
		{
			ensemble[i + 1] = ensemble[i];
		}
		ensemble[j]=inensem[max];	

	}

}

clock_t TIME;
void puttime(FILE *fo,char *mess)
{
	fprintf(fo,mess);
	fprintf(fo,": %lf\n",(double)(clock()-TIME)/CLOCKS_PER_SEC);
	printf(mess);
	printf(": %lf\n",(double)(clock()-TIME)/CLOCKS_PER_SEC);
	TIME=clock();
}

double cp;

int main(int argc, char *argv[])
{
    clock_t T1;
	time_t sec;
	sec=time(NULL);
    double ti;
	TIME=clock();
	
    long i,j,k;
    int size;
	FILE *fout;
	
    char fname[40];
    sprintf(fname,"results_%s",filename);
	fout=fopen(fname,"w");

    //size=omp_get_num_procs();
	size=1;
	fprintf(fout,"number of processors: %d\nclock per second: %d\n",size,CLOCKS_PER_SEC);
	printf("number of processors: %d\nclock per second: %d\n",size,CLOCKS_PER_SEC);

	unsigned int seed, *seedlist;
	seed=(unsigned int)time(NULL);
    int krg,know,copy1,copy2;
    int kmax,kp;
    krg=atoi(argv[1]); //parameter for Randomized Greedy algorithm  
    copy1=atoi(argv[2]);//number of partitions for original network
    kmax=copy1*size;
    copy2=atoi(argv[3]);//number of partitions for reduced network
    kp=copy2*size;
    cp=atof(argv[4]);//exponent for density in modularity density
    filename=(argv[5]); // input filename
	fprintf(fout,"seed=%ld\n",seed);
	fprintf(fout,"krg=%d\nkmax=%d\nkp=%d\n",krg,kmax,kp);
	printf("seed=%ld\n",seed);
	printf("krg=%d\nkmax=%d\nkp=%d\n",krg,kmax,kp);
    printf ("cp %lf\n",cp);

    seedlist=(unsigned int*)malloc(sizeof(unsigned int)*size);
    for (i=0;i<size;i++)
	seedlist[i]=seed+i;

	struct graph G;
	inputGraph(&G);
	printf("links=%lf nodes=%d\n",G.n,G.N);
    
   	puttime(fout,"reading time");
    struct part *ensemble,*iterensem;
    ensemble=(struct part*)malloc(sizeof(struct part)*kmax);
    iterensem=(struct part*)malloc(sizeof(struct part)*kp);
    printf("working on initialization\n");
    fprintf(fout,"working on initialization\n");
   
    int *edgecc, *score;
    score=(int*)malloc(sizeof(int)*G.N);
    edgecc=(int*)malloc(sizeof(int)*G.N*G.N);
    
    for (j=0;j<copy1;j++)
	{
		{
            ensemble[j*size]=RG(G,krg,seedlist);
        }
		//get set of partitions
	}

    puttime(fout,"initialization time");
    know=cleansort(ensemble,kmax,G.N);
    puttime(fout,"sort time");

    cols=(int*)malloc(sizeof(int)*G.N);
	i=0;
	while (i<G.N) {cols[i]=0; i++;}
    for (i=0;i<G.N-1;i++)
       if (!cols[i])
	{
		cols[i]=1;
        for (j=i+1;j<G.N;j++)
		{
			edgecc[i*G.N+j]=0;
			for (k=0;k<know;k++)
                {
                    if (ensemble[k].pa[i]==ensemble[k].pa[j])
                	{
                        
                        edgecc[i*G.N+j]++;
                    }
                }
            
            if (edgecc[i*G.N+j]==know)
                {
                    cols[j]=1;
                }
                
		}
	}
	
     //initialization
     puttime(fout,"count edgecc time");
    
    int iteration=0;
    int *ref, bp;
    ref=(int*)malloc(sizeof(int)*G.N);
	for (i = 0;i < G.N;i++)
		score[i] = i + 1;

    while (know>1)
    {
        iteration++;
        printf("\n");
        printf("%d:size=%d,ensem=%d,%lf,%lf\n",iteration,G.com,know,ensemble[0].Q/G.n,ensemble[know-1].Q/G.n);
        int sizeofRN;
        sizeofRN=getscore(score,edgecc,G.N,know,ref);
        renorm(G,ref);
        G.com=sizeofRN;
        puttime(fout,"reducing time");
    
        for (j=0;j<copy2;j++)
        {
        //#pragma omp parallel private (i)
        {
            iterensem[j*size]=RG(G,krg,seedlist);
    
        }
        }
        puttime(fout,"RG time");
        //pick the best
		bp = 0;	j = 1;
		for (i=1;i<kp;i++)
			if (abso(iterensem[i].Q - iterensem[bp].Q) < Tsmall)
			{
				j++;
				if (rand1(seedlist)*j < 1)
				{
					free(iterensem[bp].pa);	
					bp = i;
				}
				else
					free(iterensem[i].pa);
			}
			else if (iterensem[i].Q > iterensem[bp].Q)
			{
				j = 1;
				free(iterensem[bp].pa);
				bp=i;
			}
			else
			{
				free(iterensem[i].pa);
			}


        //update ensemble and edgecc
		i = -1;
		for (j = 0;j<know;j++)
			if (abso(iterensem[bp].Q -ensemble[j].Q )<Tsmall)
			{
				if (iterensem[bp].com == ensemble[j].com)
					if (comps(iterensem[bp].pa, ensemble[j].pa, G.N))
					{
						i = j;    break;
					}
			}
		if (i >= 0)
		{
			remolast(ensemble, know, edgecc, G.N);
			know--;
		}
		else
		{
			if (iterensem[bp].Q<ensemble[0].Q)
			{
				remolast(ensemble, know, edgecc, G.N); know--;
        
			}
			else
			{
				if (know<kmax)
				{
					addone(ensemble, know, iterensem, bp, edgecc, G.N);
					know++;
				}
				else
				{
					replaceone(ensemble, know, iterensem, bp, edgecc, G.N);
				}
			}
        }   
        puttime(fout,"update ensemble time");
    }
    fprintf(fout,"rest time=%lf\n",(double)(clock()-T1)/CLOCKS_PER_SEC);
	fprintf(fout,"real walltime=%d\n",time(NULL)-sec);
	printf("rest time=%lf\n",(double)(clock()-T1)/CLOCKS_PER_SEC);
	printf("real walltime=%d\n",time(NULL)-sec);
    //get score of final reduced network
    int sizeofRN;
    sizeofRN=getscore(score,edgecc,G.N,know,ref);
    renorm(G,ref);
    G.com=sizeofRN;
    sprintf(fname,"partition_%s",filename);    
	outpart(ensemble[0],G.n, G.N,fname);
    double Qfinal;
	Qfinal=compQ(G);
	printf("Qfinal here=%lf\n",Qfinal/G.n);
    fprintf(fout,"Qfinal=%lf\n",Qfinal/(G.n));
    fclose(fout);
	return 0;
}






int mergroup(struct graph *G, int x, int y)
{
    // printf("at first rednode neighborlist\n");
    // for (int i = 0; i <= G[0].nl_red[x][0]; i++) 
    // {
    //     printf("%d ", G[0].nl_red[x][i]);
    // }
    // printf("\n");
    int i;
    int *newlist_red, *newlist_blue;
    //newlist will count for member list and neighborlist
    double *newl2_red, *newl2_blue;
    
    //first with the members of x and y
    newlist_red=(int*)malloc(sizeof(int)*(G[0].ml_red[x][0]+G[0].ml_red[y][0]+1));
    newlist_blue=(int*)malloc(sizeof(int)*(G[0].ml_blue[x][0]+G[0].ml_blue[y][0]+1));
    for (i=1;i<=G[0].ml_red[x][0];i++)	newlist_red[i]=G[0].ml_red[x][i];
    for (i=1;i<=G[0].ml_red[y][0];i++)	newlist_red[i+G[0].ml_red[x][0]]=G[0].ml_red[y][i];
    for (i=1;i<=G[0].ml_blue[x][0];i++)	newlist_blue[i]=G[0].ml_blue[x][i];
    for (i=1;i<=G[0].ml_blue[y][0];i++)	newlist_blue[i+G[0].ml_blue[x][0]]=G[0].ml_blue[y][i];
    newlist_red[0]=G[0].ml_red[x][0]+G[0].ml_red[y][0];
    newlist_blue[0]=G[0].ml_blue[x][0]+G[0].ml_blue[y][0];
    free(G[0].ml_red[x]);
    free(G[0].ml_red[y]);
    free(G[0].ml_blue[x]);
    free(G[0].ml_blue[y]);
    G->ml_red[x]=newlist_red;
    G->ml_red[y]=NULL;
    G->ml_blue[x]=newlist_blue;
    G->ml_blue[y]=NULL;
    // printf("rednodelist\n");
    // for (int i = 0; i <= newlist_red[0]; i++) 
    // {
    //     printf("%d ", G[0].ml_red[x][i]);
    // }
    // printf("\n");
    // printf("bluenodelist\n");
    // for (int i = 0; i <= newlist_blue[0]; i++) 
    // {
    //     printf("%d ", G[0].ml_blue[x][i]);
    // }
    // printf("\n");
    //secondly, with the neighbors of x and y
    newlist_red=(int*)malloc(sizeof(int)*(G[0].nl_red[x][0]+G[0].nl_red[y][0]+1));
    newl2_red  =(double*)malloc(sizeof(double)*(G[0].nl_red[x][0]+G[0].nl_red[y][0]+1));
    newl2_red[0]=G->wl_red[x][0]+G->wl_red[y][0];
    newlist_red[0]=0;

    newlist_blue=(int*)malloc(sizeof(int)*(G[0].nl_blue[x][0]+G[0].nl_blue[y][0]+1));
    newl2_blue  =(double*)malloc(sizeof(double)*(G[0].nl_blue[x][0]+G[0].nl_blue[y][0]+1));
    newl2_blue[0]=G->wl_blue[x][0]+G->wl_blue[y][0];
    newlist_blue[0]=0;

    int direct_connection_redX = 0;//if direct conn. from red node of X cluster
    int direct_connection_blueX = 0;//if direct conn from blue node of X cluster

    //all neighbors of red node of x will be taken care of
    //for red node of x group
    if (G[0].nl_red[x][0]!=0)
    {
        //printf ("G[0].nl_red[x][0]!=0\n");
        for (i=1;i<=G[0].nl_red[x][0];i++)//for all blue neighbors of red node
        {
            if (G[0].nl_red[x][i]==y) 
            {	direct_connection_redX = 1 ;	newl2_red[0]+=G[0].wl_red[x][i]; newl2_blue[0]+=G[0].wl_red[x][i];  
                //printf ("direct connection, neighbor and newl2_red[0] and newl2_blue[0] are %d %lf %lf\n", G[0].nl_red[x][i], newl2_red[0], newl2_blue[0]);
            }
            else
            {
                newlist_red[0]++;
                newlist_red[newlist_red[0]]=G[0].nl_red[x][i];
                newl2_red[newlist_red[0]]=G[0].wl_red[x][i];
            }
        }
    }

    //for blue node of x group
    if (G[0].nl_blue[x][0]!=0)
    {
        for (i=1;i<=G[0].nl_blue[x][0];i++)//for all red neighbors of red node
        {
            if (G[0].nl_blue[x][i]==y) 
            { 
                direct_connection_blueX = 1 ; newl2_blue[0]+=G[0].wl_blue[x][i]; newl2_red[0]+=G[0].wl_blue[x][i]; 
                //printf ("direct connection, neighbor and newl2_blue[0] and newl2_red[0] are %d %lf %lf\n", G[0].nl_blue[x][i], newl2_blue[0],  newl2_red[0]);
            }
            else
            {
                newlist_blue[0]++;
                newlist_blue[newlist_blue[0]]=G[0].nl_blue[x][i];
                newl2_blue[newlist_blue[0]]=G[0].wl_blue[x][i];
            }
        }
    }

    //from the y cluster
    //for red node of y cluster
    if (G[0].nl_red[y][0]!=0)
    {
        for (i=1;i<=G[0].nl_red[y][0];i++)
        if (G[0].nl_red[y][i]!=x)
        {
            int sg_blue,k1_blue,k2_blue;
            sg_blue=G[0].nl_red[y][i];
            //printf ("not x %d\n",G[0].nl_red[y][i]);
            k1_blue=findxiny(x,G[0].nl_blue[sg_blue]);// the index of x in the adjacency list of sg
            k2_blue=findxiny(y,G[0].nl_blue[sg_blue]);//the index of y in the adjacency list of sg
            //printf("k1_blue %d k2_blue %d\n", k1_blue, k2_blue);
            if (k1_blue>0)
            {
                //printf("x ia a neighbor of sg_blue\n");
                G[0].wl_blue[sg_blue][k1_blue]+=G[0].wl_blue[sg_blue][k2_blue];
                //printf ("wl_blue[sg_blue][k1_blue] %d\n", G[0].wl_blue[sg_blue][k1_blue]);
                int test2;
                test2=findxiny(sg_blue,newlist_red);
                if (test2>0)
                newl2_red[test2]=G[0].wl_blue[sg_blue][k1_blue];
                else
                {
                    printf("Error!\n");
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
                    printf("Error1!\n");
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
            //printf ("not x %d\n",G[0].nl_blue[y][i]);
            k1_red=findxiny(x,G[0].nl_red[sg_red]);// the index of x in the adjacency list of sg
            k2_red=findxiny(y,G[0].nl_red[sg_red]);//the index of y in the adjacency list of sg
            //printf("k1_red %d k2_red %d\n", k1_red, k2_red);
            if (k1_red>0)
            {
                //printf("x ia a neighbor of sg_red\n");
                G[0].wl_red[sg_red][k1_red]+=G[0].wl_red[sg_red][k2_red];
                //printf ("wl_red[sg_red][k1_red] %d\n", G[0].wl_red[sg_red][k1_red]);
                int test;
                test=findxiny(sg_red,newlist_blue);
                if (test>0)
                newl2_blue[test]=G[0].wl_red[sg_red][k1_red];
                else
                {
                    printf("Error2!\n");
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
                    printf("Error4!\n");
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
    free(G[0].nl_blue[y]);	free(G[0].wl_blue[y]);
    G->nl_red[y]=NULL;
    G->wl_red[y]=NULL;
    G->nl_blue[y]=NULL;
    G->wl_blue[y]=NULL;
    
    G[0].dl_red[x]+=G[0].dl_red[y];
    G[0].dl_blue[x]+=G[0].dl_blue[y];

    //printf("degree red %f\n",G[0].dl_red[x]);
    //printf("degree blue %d\n",G[0].dl_blue[x]);
    // for(i=1;i<=G[0].nl_red[x][0];i++)
    //     printf("blue neighbor %d\n",G[0].nl_red[x][i]);
    // for(i=1;i<=G[0].nl_blue[x][0];i++)
    //     printf("neighbor %d\n",G[0].nl_blue[x][i]);
    return 0;
}


void movegroup(struct graph *G, int f, int t) //here ml[t] is null
{
    int i,j;

    G->ml_red[t]=G->ml_red[f];
    G->ml_red[f] = NULL; //now ml[f] is null

    G->ml_blue[t]=G->ml_blue[f];
    G->ml_blue[f] = NULL; //now ml[f] is null

    for (i=1;i<=G->nl_red[f][0];i++)
    {
        int temp_Blue;
        temp_Blue = G->nl_red[f][i];
        j=findxiny(f,G->nl_blue[temp_Blue]);
        G->nl_blue[temp_Blue][j]=t;
    }  
    for (i=1;i<=G->nl_blue[f][0];i++)
    {
        int temp_red;
        temp_red = G->nl_blue[f][i];
        j=findxiny(f,G->nl_red[temp_red]);
        G->nl_red[temp_red][j]=t;
    }
    
    G->nl_red[t]=G->nl_red[f];
    G->wl_red[t]=G->wl_red[f];
    G->dl_red[t]=G->dl_red[f];
    G->nl_blue[t]=G->nl_blue[f];
    G->wl_blue[t]=G->wl_blue[f];
    G->dl_blue[t]=G->dl_blue[f];
    G->nl_blue[f] = NULL;
    G->wl_blue[f] = NULL;
    G->nl_red[f] = NULL;
    G->wl_red[f] = NULL;

}

void renorm(struct graph G, int *ref)
{
    int i,j;
    
    //renormlize G
    for (i=0;i<G.com;i++)
    {
		//printf("i %d\n",i);
        if (ref[i]!=i)
        {
            //printf("i %d ref[i] %d\n", i, ref[i]);
			if (G.ml_red[ref[i]] == NULL || G.ml_blue[ref[i]] == NULL)
				{
                    // printf("blue neighbors\n");
                    // for (j=0; j<G.nl_red[i][0]+1; j++) printf("%d ",G.nl_red[i][j]);
                    // printf("\n");
                    // printf("red neighbors\n");
                    //for (j=0; j<G.nl_blue[i][0]+1; j++) printf("%d ",G.nl_blue[i][j]);
                    //printf("\n");
                    movegroup(&G,i,ref[i]);
                    // printf("moving %d and %d\n", i ,ref[i]);
                    // printf("blue neighbors\n");
                    //for (j=0; j<G.nl_red[ref[i]][0]+1; j++) printf("%d ",G.nl_red[ref[i]][j]);
                    // printf("\n");
                    // printf("red neighbors\n");
                    // for (j=0; j<G.nl_blue[ref[i]][0]+1; j++) printf("%d ",G.nl_blue[ref[i]][j]);
                    // printf("\n");
                }
			else
                {
					mergroup(&G,ref[i],i);
					//printf("merging %d to %d \n",i, ref[i]);
				}
        }
    }
}


int getscore(int *score,int *edgecc, int N,int know,int *ref)
{
    long i,j,k;
    int groupsnumber;
	int *newsc;
	newsc = (int*)malloc(sizeof(int)*N);
    for (i=0;i<N;i++)
        newsc[i]=0;
    groupsnumber=0;
    for (i=0;i<N;i++)
    {
        if (!newsc[i])
        {
            ref[score[i]-1]=groupsnumber;
            groupsnumber++;
            newsc[i]=groupsnumber;

        	for (j=i+1;j<N;j++)
        	{
            if (edgecc[i*N+j]==know)
            	{
                ref[score[j] - 1] = newsc[i]-1;
        		newsc[j]=newsc[i];
                
            	}
        	}
        }
    }
	for (i = 0;i < N;i++)
		score[i] = newsc[i];
	free(newsc);
    return (groupsnumber);
}

int cleansort(struct part *ensem, int kmax,int N)
{
    int i,j;
    struct part temp;
    for (i=0;i<kmax-1;i++)
        for (j=i+1;j<kmax;j++)
        {
            if (ensem[i].Q>ensem[j].Q)
            {
                temp=ensem[i];  ensem[i]=ensem[j]; ensem[j]=temp;
            }
        }
    int *del,len;
    del=(int*)malloc(sizeof(int)*kmax);
    for (i=0;i<kmax;i++)
        del[i]=0;
	
	//Remove duplicates based on partition similarity
    for (i=0;i<kmax-1;i++)
        if (del[i]==0)
        {
            j=i+1;
            while (j<kmax&&(ensem[j].Q-ensem[i].Q)<Tsmall)
            {
                if (ensem[j].com==ensem[i].com)
                    if (comps(ensem[i].pa,ensem[j].pa,N))
                        del[j]=1;
                j++;
            }
        }
    for (i=0;i<kmax-1;i++)
        for (j=i+1;j<kmax;j++)
        {
            if (del[i]>del[j])
            {
                del[i]=0;   del[j]=1;
                temp=ensem[i];  ensem[i]=ensem[j]; ensem[j]=temp;
            }
            else if (del[i]+del[j]==0)
            {
                if (ensem[i].Q>ensem[j].Q)
                {
                    temp=ensem[i];  ensem[i]=ensem[j]; ensem[j]=temp;
                }
            }
        }
	//Find the valid length after cleaning
    len=0;
    while (len<kmax&&del[len]==0) len++;

    i=len;
    while (i<kmax)
    {
        free(ensem[i].pa);
        i++;
    }
    free(del);
    return len;
}



void inputGraph(struct graph *G)
{

	
	int i,n1,n2,N,Gn,*nnl_red,*nnl_blue;
	double w;
	FILE *fi;
	fi=fopen("info.txt","r");
	fscanf(fi,"%d%d%d",&G->N_red,&G->N_blue,&Gn);
	fclose(fi);
	G->n=0;
	G->N=G->N_red+G->N_blue;
	G->com=G->N;
	G->ml_blue=(int**)malloc(sizeof(int*)*G->N);
	G->ml_red=(int**)malloc(sizeof(int*)*G->N);
	G->nl_red=(int**)malloc(sizeof(int*)*G->N);
	G->nl_blue=(int**)malloc(sizeof(int*)*G->N);
	G->wl_blue=(double**)malloc(sizeof(double*)*G->N);
	G->wl_red=(double**)malloc(sizeof(double*)*G->N);
	G->dl_red=(double*)malloc(sizeof(double)*G->N);
	G->dl_blue=(double*)malloc(sizeof(double)*G->N);
	nnl_red=(int*)malloc(sizeof(int)*G->N);
	nnl_blue=(int*)malloc(sizeof(int)*G->N);
    // for (i = 0; i < G->N; i++) {
    //     nnl_blue[i] = 0;  // Ensures no garbage values
    // }
    // for (i = 0; i < G->N; i++) {
    //     nnl_red[i] = 0;  // Ensures no garbage values
    // }
    // Initialize the arrays with 0
    for (int i = 0; i < G->N; i++) {
        G->dl_red[i] = 0;
        G->dl_blue[i] = 0;
		nnl_red[i] = G->N;
		nnl_blue[i] = G->N;
    }

	fi=fopen("degree_red.txt","r");
	for (i=0;i<G->N_red;i++)
	{
		fscanf(fi,"%d%lf",nnl_red+i,G->dl_red+i);
		G->n+=G->dl_red[i];
	}

    
	fclose(fi);

	fi=fopen("degree_blue.txt","r");
	for (i=G->N_red;i<G->N;i++)
		{
			//int index=G->N_red+i;
			fscanf(fi,"%d %lf",nnl_blue+i,G->dl_blue+i);
            
			G->n+=G->dl_blue[i];
		}
    G->n/=2;
	printf ("\n");
	fclose(fi);
    // for (i=0;i<G->N;i++) printf("nnl_blue[%d] = %d\n", i, nnl_blue[i]); // Add print statement
    // printf("\n");

	fi=fopen("clean.txt","r");
	for (i=0;i<G->N;i++)
	{
		// G->ml[i]=(int*)malloc(sizeof(int)*2);
		// G->ml[i][0]=1;
		// G->ml[i][1]=i+1;
        //printf("nnl_red[%d] = %d\n", i, nnl_red[i]);
		G->nl_red[i]=(int*)malloc(sizeof(int)*(nnl_red[i]+1));
		G->nl_red[i][0]=0;
        //printf("nnl_blue[%d] = %d\n", i, nnl_blue[i]);
		G->nl_blue[i]=(int*)malloc(sizeof(int)*(nnl_blue[i]+1));
		G->nl_blue[i][0]=0;
		G->wl_red[i]=(double*)malloc(sizeof(double)*(nnl_red[i]+1));
		G->wl_red[i][0]=0;
		G->wl_blue[i]=(double*)malloc(sizeof(double)*(nnl_blue[i]+1));
		G->wl_blue[i][0]=0;

		// Initialize singleton clusters
    if (i < G->N_red)
    {
        G->ml_red[i] = (int *)malloc(sizeof(int) * 2);
        G->ml_red[i][0] = 1;
        G->ml_red[i][1] = i;
        G->ml_blue[i] = (int *)malloc(sizeof(int));
        G->ml_blue[i][0] = 0;
    }
    else
    {
        G->ml_blue[i] = (int *)malloc(sizeof(int) * 2);
        G->ml_blue[i][0] = 1;
        G->ml_blue[i][1] = i;
        G->ml_red[i] = (int *)malloc(sizeof(int));
        G->ml_red[i][0] = 0;
    }

	}


	for (i=0;i<Gn;i++)
	{
		fscanf(fi,"%d %d %lf",&n1,&n2, &w);
		


		if (n1 < G->N_red) { // n1 is a red node
            G->nl_red[n1][0]++;
            //printf("yes\n");
            //printf("n1 %d %d\n",n1,G->nl_red[n1][0]);
            G->nl_red[n1][G->nl_red[n1][0]] = n2 ; // Store blue node index in red adjacency list
            G->wl_red[n1][G->nl_red[n1][0]] = w;
        } else { // n1 is a blue node
            G->nl_blue[n1][0]++;
            G->nl_blue[n1][G->nl_blue[n1][0]] = n2; // Store red node index in blue adjacency list
            G->wl_blue[n1][G->nl_blue[n1][0]] = w;
        }

        if (n2 < G->N_red) { // n2 is a red node
            G->nl_red[n2][0]++;
            G->nl_red[n2][G->nl_red[n2][0]] = n1 + G->N_red; // Store blue node index in red adjacency list
            G->wl_red[n2][G->nl_red[n2][0]] = w;
        } else { // n2 is a blue node
            G->nl_blue[n2][0]++;
            G->nl_blue[n2][G->nl_blue[n2][0]] = n1; // Store red node index in blue adjacency list
            G->wl_blue[n2][G->nl_blue[n2][0]] = w;
        }
        n1--;	n2--;
    }
	fclose(fi);
	free (nnl_red);
	free (nnl_blue);
}


