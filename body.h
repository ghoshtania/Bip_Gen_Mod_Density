#define Inf_Sma 1e-6

struct graph
{
    int N;
    int N_red;
    int N_blue;
    double n;
    int com;
    int **ml_red;
    int **ml_blue;
    int **nl_red;
    int **nl_blue;
    double **wl_red;
    double **wl_blue;
    double *dl_blue;
    double *dl_red;
};

struct part
{
    double Q;
    int com;
    int *pa;
};