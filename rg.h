#include"body.h"


int findxiny(int x, int *y);
int outpart (struct part ans,double n,int N, char *fname);
int updateG(struct graph *G, int x, int y);
double power(double a, double b);
//void extrafG(struct part ans,struct graph G,int *s);
void gcopy(struct graph G,struct graph *Ga);
void printGraph(struct graph *G);
struct part RG(struct graph G, int ke,unsigned int *seed);
//double RG(struct graph G, int ke,unsigned int *seed);