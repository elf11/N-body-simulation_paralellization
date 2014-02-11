#include "body_struct.h"
#include <pthread.h>

/*
 * Function to set the partition for each processor.
 */
void setStartEnd(int size,int rank,int oriStart,int oriEnd,int *start,int *count) {

	int length, block, rem, end;
	
	length = oriEnd - oriStart + 1;
	block = (int)(length / size);
	rem = length - block * size;
	
	if (length > 0) {
		if (rank < rem) {
			*start = oriStart + rank * (block + 1);
			end = *start + block;
		} else {
			*start = oriStart + rank * block + rem;
			end = *start + block - 1;
		}
	} else {
		*start = 1;
		end = 0;
	}
	
	*count = end - *start + 1;
}

/*
 * Function to calculate new information for each body.
 */
void *calculateNewInfo(void *sendInfo) {
	
	int i, j;
	struct body *old = ((struct info *)sendInfo)->NBodyData;
	int start = ((struct info *)sendInfo)->start;
	int count = ((struct info *)sendInfo)->count;
	double *velocityX = ((struct info *)sendInfo)->velocityX;
	double *velocityY = ((struct info *)sendInfo)->velocityY;
	
	/*
	 * initialize the new velocities
	 */
	 for (i = 0; i < count; i++) {
	 	velocityX[start + i] = old[start + i].vx;
	 	velocityY[start + i] = old[start + i].vy;
	 }
	 
	 for (i = 0; i < count; i++) {
	 	for (j = 0; j < NUM_BODY; j++) {
	 		if (j == (start + i)) {
	 			continue;
	 		}
	 		
	 		double delta_x = old[j].x - old[start + i].x;
            double delta_y = old[j].y - old[start + i].y;
            double dist = sqrt(pow(delta_x,2) + pow(delta_y,2));
            if(dist == 0){ /* if two bodies have the same position, skip the force calculation */
                continue;
            }
            double force = GRAVITY * old[j].w / (dist * dist);
            velocityX[start + i] = velocityX[start + i] + force * delta_x / dist;
            velocityY[start + i] = velocityY[start + i] + force * delta_y / dist;
	 		
	 	}
	 }
	
}

/*
 * Function to update the information for each body.
 */
void *updateNewInfo(void *sendInfo) {
	
	int i;
	struct body *old = ((struct info *)sendInfo)->NBodyData;
	int start = ((struct info *)sendInfo)->start;
	int count = ((struct info *)sendInfo)->count;
	double *velocityX = ((struct info *)sendInfo)->velocityX;
	double *velocityY = ((struct info *)sendInfo)->velocityY;
	
	for (i = 0; i < count; i++) {
		old[start + i].x = old[start + i].x + velocityX[start + i];
        old[start + i].y = old[start + i].y + velocityY[start + i];
        old[start + i].vx = velocityX[start + i];
        old[start + i].vy = velocityY[start + i];
	}
}

int main(int argc, char **argv) {
	
	int i, k;
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
    
    /*
     * initialize the variables for pthread
     */
     pthread_t threads[NUM_PTH];
     pthread_attr_t pattr;
     pthread_attr_init(&pattr);
     pthread_attr_setdetachstate(&pattr, PTHREAD_CREATE_JOINABLE);
     struct info sendInfo[NUM_PTH];
     
     /*
      * initialize the sending info part
      */
      for (i = 0; i < NUM_PTH; i++) {
      	sendInfo[i].NBodyData = (struct body *)nBodyArray;
      	setStartEnd(NUM_PTH, i, 0, NUM_BODY-1, &sendInfo[i].start, &sendInfo[i].count);
        sendInfo[i].velocityX = velocityXArray;
        sendInfo[i].velocityY = velocityYArray;
      }
      
      
      tTime = omp_get_wtime();
      /* start to simulate N-body */
    for(k = 0; k < NUM_ITER; k++) {
        /* Calculate the position of bodies by point-to-point in each iteration */
        for(i = 0; i < NUM_PTH; i++) {
            pthread_create(&threads[i], &pattr, calculateNewInfo, &sendInfo[i]);
        }
        
        for(i = 0; i < NUM_PTH; i++) {
            pthread_join(threads[i], NULL);
        }
        /* update the new data */
        for(i = 0; i < NUM_PTH; i++){
            pthread_create(&threads[i], &pattr, updateNewInfo, &sendInfo[i]);
        }
        for(i = 0; i < NUM_PTH; i++){
            pthread_join(threads[i], NULL);
        }
    }
    
    
    /*
     * stop the time and display it
     */
    tTime = omp_get_wtime() - tTime;
    fprintf(f_in, "PThread run tTime for %d iterations is %f\n", NUM_ITER, tTime);
    fclose(f_in);
    
    pthread_attr_destroy(&pattr);
    pthread_exit(NULL);
    return 0;
}
