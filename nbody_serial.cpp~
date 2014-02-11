#include <omp.h>
#include "body_struct.h"

int main(int argc, char **argv) 
{
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

    for (k = 0; k < NUM_ITER; k += 1) {
        for (i = 0; i < NUM_BODY; i += 1) {
            for (j = 0; j < NUM_BODY; j += 1) {
                if (j == i) {
                    continue;
                }
                double delta_x = nBodyArray[j].x - nBodyArray[i].x;
                double delta_y = nBodyArray[j].y - nBodyArray[i].y;
                double dist = sqrt(pow(delta_x,2) + pow(delta_y,2));
                if (dist == 0) {
                    continue;
                }
                double force = GRAVITY * nBodyArray[j].w / (dist * dist);
                velocityXArray[i] = velocityXArray[i] + force * delta_x / dist;
                velocityYArray[i] = velocityYArray[i] + force * delta_y / dist;
            }
        }

        /*
         * update the data
         */
        for (i = 0; i < NUM_BODY; i += 1) {
            nBodyArray[i].x = nBodyArray[i].x + velocityXArray[i];
            nBodyArray[i].y = nBodyArray[i].y + velocityYArray[i];
            nBodyArray[i].vx = velocityXArray[i];
            nBodyArray[i].vy = velocityYArray[i];
        }

    }

    /*
     * stop the time and display it
     */
    tTime = omp_get_wtime() - tTime;
    fprintf(f_in, "Serial run tTime for %d iterations is %f\n", NUM_ITER, tTime);
    fclose(f_in);

    return 0;
}
