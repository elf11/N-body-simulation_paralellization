#include <omp.h>
#include "body_struct.h"

int main(int argc, char **argv) {

	int i, j, k;
    struct body nBodyArray[NUM_BODY];
    double velocityXArray[NUM_BODY];
    double velocityYArray[NUM_BODY];
    double tTime;
    FILE * f_in;
    
    f_in = fopen("time_values.txt", "a");

    /*
     * intialize the simulation
     */
    srand(time(NULL));
    for (i = 0; i < NUM_BODY; i += 1) {
        nBodyArray[i].x = rand() % (MAX_X - MIN_X) + MIN_X;
        nBodyArray[i].y = rand() % (MAX_Y - MIN_Y) + MIN_Y;
        nBodyArray[i].vx = velocityXArray[i] = rand() % (MAX_V - MIN_V) - (MAX_V + MIN_V) / 2;
        nBodyArray[i].vy = velocityYArray[i] = rand() % (MAX_V - MIN_V) - (MAX_V + MIN_V) / 2;
        nBodyArray[i].w = rand() % (MAX_W - MIN_W) + MIN_W;
    }

    tTime = omp_get_wtime();
    
     /* start to simulate N-body */
    for(k = 0; k < NUM_ITER; k++){

		/* set the number of threads */
		omp_set_num_threads(NUM_THS);
		#pragma omp parallel private(j)
		{
			/* Calculate the position of bodies by point-to-point in each  iteration */
			#pragma omp for schedule(static)
		    for(i = 0; i < NUM_BODY; i++) {
		        for(j = 0; j < NUM_BODY; j++) {
			        if(j == i){ /* there is no need to calculate the effect from itself */
			            continue;
			        }
			        double delta_x = nBodyArray[j].x - nBodyArray[i].x;
			        double delta_y = nBodyArray[j].y - nBodyArray[i].y;
			        double dist = sqrt(pow(delta_x,2) + pow(delta_y,2));
			        if(dist == 0) { /* if two bodies have the same position, skip the force calculation */
		            	continue;
			        }
			        double force = GRAVITY * nBodyArray[j].w / (dist * dist);
			        velocityXArray[i] = velocityXArray[i] + force * delta_x / dist;
			        velocityYArray[i] = velocityYArray[i] + force * delta_y / dist;
			    }
			}
			
		    /* update the new data */
			#pragma omp for schedule(static)
		    for(i = 0; i < NUM_BODY; i++){
		        nBodyArray[i].x = nBodyArray[i].x + velocityXArray[i];
		        nBodyArray[i].y = nBodyArray[i].y + velocityYArray[i];
		        nBodyArray[i].vx = velocityXArray[i];
		        nBodyArray[i].vy = velocityYArray[i];
		    }
		}
  	}
  	
  	/*
     * stop the time and display it
     */
    tTime = omp_get_wtime() - tTime;
    fprintf(f_in, "Omp run tTime for %d iterations is %f\n", NUM_ITER, tTime);
    fclose(f_in);

	return 0;
}
