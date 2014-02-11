#include <mpi.h>
#include <omp.h>
#include "body_struct.h"

typedef struct sim_param_t {
    int   npart;
    int   niter;
    int   npframe;
    float dt;
    float sig;
} sim_param_t;


static int rank;
static int nproc;
static MPI_Datatype pairtype;

/*
 * Impartim problema in subprobleme, fiecare dintre procesoare
 * va primi cel putin npart/nproc particule.
 */

void partition(int *parts, int *count, int npart)
{
    int num = npart/nproc;
    int toPartition = npart - num * nproc;
    parts[0] = 0;
    
    for (int i = 0; i < nproc; i++) {
        count[i] = num + (i < toPartition ? 1 : 0);
        parts[i + 1] = parts[i] + count[i];
    }
}

/*
 * Calculam fortele ce actioneaza asupra particulelor, avem un vector
 * de pozitii locale si un vector de pozitii globale si calculam
 * fortele intr-un vector de forte locale.
 */

void calculateForces(int n, const float *x, 
                    int istart, int iend, 
                    const float *xlocal, float *Flocal,
                    sim_param_t *params)
{
    float sig  = params->sig;
    float sig2 = sig*sig;

    for (int i = istart; i < iend; i++) {
        int ii = i - istart;
        for (int j = 0; j < n; ++j) {
            if (i != j) {
                float dx = x[2*j+0]-xlocal[2*ii+0];
                float dy = x[2*j+1]-xlocal[2*ii+1];
                
                float rad = dx * dx + dy * dy;
                float potential = 0;
                if (rad < GRAVITY * sig2) {
                	float val = sig2 / rad;
                	float val3 = val * val * val;
                	potential = 24 / rad * val3 * (1 - 2 * val3); 
                } else {
                	potential = 0;
                }
                
                Flocal[2 * ii + 0] += (potential * dx);
                Flocal[2 * ii + 1] += (potential * dy);
            }
        }
    }
}


/*
 * Pentru celelalte implementari am ales distributia particulelor
 * si vitezelor lor in acelasi timp, acum nu mai este nevoie sa facem asa.
 * Generam o distributie random a particulelor avand coordonatele
 * punctelor in spatiu si mai apoi generam velocitatile, pentru
 * particulele generate, deoarece putem respinge unele dintre 
 * cele NUM_BODY particule si nu are rost sa generam viteze si pentru acestea.
 */

int initializeParticleGeneral(int n, float *x, sim_param_t *params)
{
    const int MAX_INIT_TRIALS = 1000;

    float sig    = params->sig;
    float sig2 = sig * sig;

    for (int i = 0; i < n; ++i) {
        float r2 = 0;

        for (int trial = 0; r2 < sig2 && trial < MAX_INIT_TRIALS; trial++) {
            x[2*i+0] = (float) drand48();
            x[2*i+1] = (float) drand48();
            for (int j = 0; j < i; ++j) {
                float dx = x[2 * i + 0] - x[2 * j + 0];
                float dy = x[2 * i + 1] - x[2 * j + 1];
                r2 = dx * dx + dy * dy;
                if (r2 < sig2)
                    break;
            }
        }

        /*
         * Daca depasim numarul maxim de incercari de distribuire atunci
         * intoarcem numarul de particule distribuite pana in acest moment si asta este.
         */
        if (i > 0 && r2 < sig2) 
            return i;
    }
    return n;
}


void initializeParticlesVelocities(int n, float *v)
{
    for (int i = 0; i < n; i++) {
        double R = sqrt((-2) * log(drand48()));
        double T = 2 * M_PI * drand48();
        v[2 * i + 0] = (float) (R * cos(T));
        v[2 * i + 1] = (float) (R * sin(T));
    }
}

/*
 * Functie de simulare pentru nbody.
 *
 *
 */

void calculateVelocity1(int n, float dt, float *x, float *v, float *a)
{
    for (int i = 0; i < n; ++i, x += 2, v += 2, a += 2) {
        v[0] += a[0] * dt / 2;
        v[1] += a[1] * dt / 2;
        x[0] += v[0] * dt;
        x[1] += v[1] * dt;
    }
}


void runSimulation(int n, int nlocal,         /* numarul de body-uri*/
             int *parts, int *counts,  /* offset pentru fiecare procesor si nr de body-uri pe care le ia fiecare */
             int npframe,               /* numarul de pasi pe iteratie */
             int niter,               /* numar iteratii */
             float dt,                  /* delta t */
             float *nBodyArray,         /* pozitiile globale */
             float *nBodyArrayLocal,    /* pozitiile locale */
             float *velocityArrayLocal,    /* viteze locale */
             sim_param_t *params)       /* structura de mpi */
{
    float *alocal = (float*) malloc(2 * nlocal * sizeof(float));
    memset(alocal, 0, 2 * nlocal * sizeof(float));

    /*
     * calculam acceleratia local
     */
    calculateForces(n, nBodyArray, parts[rank], parts[rank + 1],
                   nBodyArrayLocal, alocal, params);

    for (int frame = 1; frame < niter; frame++) {
        for (int i = 0; i < npframe; ++i) {
            calculateVelocity1(nlocal, dt, nBodyArrayLocal, velocityArrayLocal, alocal);
            MPI_Allgatherv(nBodyArrayLocal, nlocal, pairtype,
                           nBodyArray, counts, parts, pairtype,
                           MPI_COMM_WORLD);
            calculateForces(n, nBodyArray, parts[rank], parts[rank + 1],
                           nBodyArrayLocal, alocal, params);
        }
    }

    free(alocal);
}


int main(int argc, char** argv)
{
    sim_param_t params;
    double tTime;
    float *nBodyArray;
    float *nBodyArrayLocal;
    float *velocityArrayLocal;
    int npart;
    int *iparts;
    int *counts;
    int nlocal;
    
    FILE * f_in;
    
    f_in = fopen("time_values.txt", "a");

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Type_vector(1, 2, 1, MPI_FLOAT, &pairtype);
    MPI_Type_commit(&pairtype);
    
   	params.npart = NUM_BODY;   /* Numarul de particule (1000)  */
    params.niter = NUM_ITER; /* Numarul de frame-uri (100) */
    params.npframe = NUM_STEP_FRAME; /* Numarul de pasi pentru fiecare frame (10) */
    params.dt = 1e-4;      /* Time step (1e-4) */
    params.sig = 1e-2;
    
    /*
     * Initializam totul in procesul master, procesul 0.
     */
    nBodyArray = (float *)malloc(2 * params.npart * sizeof(float));
    if (rank == 0) {
        npart = initializeParticleGeneral(params.npart, nBodyArray, &params);
    	/* incepem simulare dam drumul la timp */
    	tTime = omp_get_wtime();
    }
    
    /* Broadcast informatie initiala de la master */
    MPI_Bcast(&npart, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(nBodyArray, npart, pairtype, 0, MPI_COMM_WORLD);
    params.npart = npart;

    /* Decidem cine este responsabil epntru ce body */
    counts = (int *)malloc( nproc * sizeof(int));
    iparts = (int *)malloc((nproc + 1) * sizeof(int));
    partition(iparts, counts, npart);
    nlocal = counts[rank];

    /* Allocam spatiu pentru storage local si copiem datele */
    nBodyArrayLocal = (float *)malloc(2 * nlocal * sizeof(float));
    velocityArrayLocal = (float *)malloc(2 * nlocal * sizeof(float));
    memcpy(nBodyArrayLocal, nBodyArray + 2 * iparts[rank], 2 * nlocal * sizeof(float));
    initializeParticlesVelocities(nlocal, velocityArrayLocal);

    runSimulation(npart, nlocal, iparts, counts,
            params.npframe, params.niter, 
            params.dt, nBodyArray, nBodyArrayLocal, velocityArrayLocal, &params);

    free(velocityArrayLocal);
    free(nBodyArrayLocal);
    free(iparts);
    free(nBodyArray);

    MPI_Finalize();
    
     /*
     * stop the time and display it
     */
    if (rank == 0) {
    	tTime = omp_get_wtime() - tTime;
    	fprintf(f_in, "MPI run tTime for %d iterations is %f\n", NUM_ITER, tTime);
    	fclose(f_in);
    }
    
    return 0;
}
