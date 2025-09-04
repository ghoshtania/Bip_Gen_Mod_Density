#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"body.h"
extern double cp;
double ran2(long *idum);


int randint(int K,unsigned int *seed)
{
	return ((int)(rand_r(seed)/(RAND_MAX*1.0+1)*K));
	
}
int randcolor(unsigned int *seed)
{
    // Choose randomly between 0 (red) and 1 (blue)
    return (int)(rand_r(seed) / (RAND_MAX * 1.0 + 1) * 2);
}

double rand1(unsigned int* seed)
{
	return (rand_r(seed)/(RAND_MAX*1.0+1));
	
}

double abso(double x)
{	if (x<0)  return (-x); return (x);	}

double compQ(struct graph G)
{
    double ans=0;
    int i;
    double den;
    double power;
    for (i=0;i<G.com;i++)
    {

        if (G.wl_red[i][0]==G.wl_blue[i][0])
            {   // Skip communities with no red and no blue nodes
                if (G.ml_red[i][0] == 0 && G.ml_blue[i][0] == 0)
                    {
                        continue;
                    }
                
                // If one type of node is missing (either n_r or n_b is zero), set density to 1
                if (G.ml_red[i][0] == 0 || G.ml_blue[i][0] == 0) 
                {
                    den = 0;
                } else 
                {
                    // Calculate density using the formula: den = m_c / (n_r * n_b)
                    den = G.wl_red[i][0] /(G.ml_red[i][0] * G.ml_blue[i][0]);
                    //printf("density1 %lf\n",den);
                }
                power=pow(den,cp);
                if (den==0) power=0;
                ans+=((G.wl_red[i][0])-G.dl_red[i]*G.dl_blue[i]/(G.n))*power;
            }
        else
            printf("Error in Mod!\n");
    }
    return (ans);
}


void lcopy(int *p, int *q, int N)
{
	int i;
	for (i=0;i<N;i++)
		p[i]=q[i];
}
