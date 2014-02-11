#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <omp.h>

/*
 *  * Define the height of the window where the bodies can reside.
 *   */

#define NUM_BODY    3000
#define MAX_X       900
#define MAX_Y       600
#define MIN_X       100
#define MIN_Y       100
#define NUM_FRAMES	100
#define NUM_STEP_FRAME	1

/*
 *  * Define values for the weights of the bodies.
 *   */
#define MAX_W       1000
#define MIN_W       250

/*
 *  * Define values for the velocities of the bodies.
 *   */
#define MAX_V       20
#define MIN_V       0

#define NUM_ITER    100

/*
 * Define the number of threads in OpenMP.
 *
 */

#define NUM_THS		4

/*
 * Define number of threads for PTHREAD.
 *
 */
#define NUM_PTH		4


/*
 *  * Define gravitational constant
 *   
 */
const double GRAVITY = 6.67259;

/*
 * Define body structure.
 *
 */
struct body {
    /*
     * position = (x,y)
     */
    double x;
    double y;
    /*
     * velocities
     */
    double vx;
    double vy;
    double w;
};

/*
 * Define a structure for calculating info.
 *
 */
struct info {
	struct body *NBodyData;
	int start, count;
	double *velocityX, *velocityY;
};
