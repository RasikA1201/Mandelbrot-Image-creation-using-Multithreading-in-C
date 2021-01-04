/*
Student Name: Rasika Hedaoo
Student ID: 1001770527

*/

#include "bitmap.h"
#include <stdio.h>
#include <getopt.h>
#include <stdlib.h>

#include <math.h>
#include <errno.h>
#include <string.h>
#include <pthread.h>
#include <sys/time.h>

int iteration_to_color( int i, int max );
int iterations_at_point( double x, double y, int max );
/* void compute_image( struct bitmap *bm, double xmin, double xmax, double ymin, double ymax, int num_t, int max ); */

void show_help()
{
	printf("Use: mandel [options]\n");
	printf("Where options are:\n");
	printf("-m <max>    The maximum number of iterations per point. (default=1000)\n");
	printf("-x <coord>  X coordinate of image center point. (default=0)\n");
	printf("-y <coord>  Y coordinate of image center point. (default=0)\n");
	printf("-s <scale>  Scale of the image in Mandlebrot coordinates. (default=4)\n");
	printf("-W <pixels> Width of the image in pixels. (default=500)\n");
	printf("-H <pixels> Height of the image in pixels. (default=500)\n");
	printf("-n <pixels> Number of threads. \n");
	printf("-o <file>   Set output file. (default=mandel.bmp)\n");
	printf("-h          Show this help text.\n");
	printf("\nSome examples are:\n");
	printf("mandel -x -0.5 -y -0.5 -s 0.2\n");
	printf("mandel -x -.38 -y -.665 -s .05 -m 100\n");
	printf("mandel -x 0.286932 -y 0.014287 -s .0005 -m 1000\n\n");
}

typedef struct 
{
	struct bitmap *bm;
	double xmin;
	double xmax;
	double ymin;
	double ymax;
	int max;
	int ht_frst;
	int ht_lst;
}send_argument;

int num_t; //for implementing n
	
void *compute_image(void *ptr);

int main( int argc, char *argv[] )
{
	
	struct timeval begin;
	struct timeval end;

	gettimeofday( &begin, NULL );
  
	char c;

	// These are the default configuration values used
	// if no command line arguments are given.

	const char *outfile = "mandel.bmp";
	double xcenter = 0;
	double ycenter = 0;
	double scale = 4;
	int    image_width = 500;
	int    image_height = 500;
	int    max = 1000;
	int num_t = 1; //number of threads
	

	// For each command line argument given,
	// override the appropriate configuration value.

	while((c = getopt(argc,argv,"x:y:s:W:H:m:n:o:h"))!=-1) {
		switch(c) {
			case 'x':
				xcenter = atof(optarg);
				break;
			case 'y':
				ycenter = atof(optarg);
				break;
			case 's':
				scale = atof(optarg);
				break;
			case 'W':
				image_width = atoi(optarg);
				break;
			case 'H':
				image_height = atoi(optarg);
				break;
			case 'm':
				max = atoi(optarg);
				break;
			case 'n':
				num_t = atoi(optarg); //specify the number of threads
				break;
			case 'o':
				outfile = optarg;
				break;
			case 'h':
				show_help();
				exit(1);
				break;
		}
	}
	
	// Display the configuration of the image.
	printf("mandel: x=%lf y=%lf scale=%lf max=%d num_t=%d outfile=%s\n",xcenter,ycenter,scale,max,num_t,outfile);
	
	
	// Array of thread start here
	pthread_t *threadp = malloc(num_t *sizeof(pthread_t));
	// Array of thread ends here
	
	// Array of argument start here
	send_argument *arg = malloc(num_t * sizeof(send_argument));
	// Array of argument start here
	
	//thread math starts here
	int i;
	
	//dividing height by number of n
	int ht_new = image_height / num_t; 
	
	// Create a bitmap of the appropriate size.
	struct bitmap *bm = bitmap_create(image_width,image_height);
	
	for(i = 0; i < num_t; ++i) 
	{
		arg[i].bm = bm;
		arg[i].xmin = xcenter-scale;
		arg[i].xmax = xcenter+scale;
		arg[i].ymin = ycenter-scale;
		arg[i].ymax = ycenter+scale;
		arg[i].max = max;

		if(i == 0) {
			arg[i].ht_frst = 0;
			arg[i].ht_lst = ht_new;
		}

		else {
			arg[i].ht_frst = arg[i-1].ht_lst;
			arg[i].ht_lst= arg[i-1].ht_lst + ht_new;
		}

		pthread_create(&threadp[i], NULL, compute_image, (void*) &arg[i]);


	}
	//thread math ends here
	
	//thread join starts here
	for(i = 0; i < num_t; ++i) {
		pthread_join(threadp[i], NULL);
	}
	//thread join ends here
	
	
	

	

	// Fill it with a dark blue, for debugging
	//bitmap_reset(bm,MAKE_RGBA(0,0,255,0));

	// Compute the Mandelbrot image
	/* compute_image(bm,xcenter-scale,xcenter+scale,ycenter-scale,ycenter+scale,num_t,max);*/

	// Save the image in the stated file.
	if(!bitmap_save(bm,outfile)) {
		fprintf(stderr,"mandel: couldn't write to %s: %s\n",outfile,strerror(errno));
		return 1;
	}
	
	gettimeofday( &end, NULL );

	int time_duration = ( ( end.tv_sec - begin.tv_sec ) * 1000000 + ( end.tv_usec - begin.tv_usec ) );

	printf("Duration: %d\n", time_duration );
	return 0;
}

/*
Compute an entire Mandelbrot image, writing each point to the given bitmap.
Scale the image to the range (xmin-xmax,ymin-ymax), limiting iterations to "max"
*/

void *compute_image(void *ptr)
{
	send_argument *arg = (send_argument*) ptr;
	int i,j;
	
	
	int width = bitmap_width(arg -> bm);
	int height = bitmap_height(arg -> bm);
	int ht_new = arg -> ht_frst;
	int ht_new_last = arg -> ht_lst;
	// For every pixel in the image...

	for(j=ht_new;j<ht_new_last;j++) {

		for(i=0;i<width;i++) {

			// Determine the point in x,y space for that pixel.
			double x = arg->xmin + i*(arg->xmax - arg->xmin)/width;
			double y = arg->ymin + j*(arg->ymax - arg->ymin)/height;

			// Compute the iterations at that point.
			int iters = iterations_at_point(x,y,arg->max);

			// Set the pixel in the bitmap.
			bitmap_set(arg->bm,i,j,iters);
		}
	}
	
	return 0;
}

/*
Return the number of iterations at point x, y
in the Mandelbrot space, up to a maximum of max.
*/

int iterations_at_point( double x, double y, int max )
{
	double x0 = x;
	double y0 = y;

	int iter = 0;

	while( (x*x + y*y <= 4) && iter < max ) {

		double xt = x*x - y*y + x0;
		double yt = 2*x*y + y0;

		x = xt;
		y = yt;

		iter++;
	}

	return iteration_to_color(iter,max);
}

/*
Convert a iteration number to an RGBA color.
Here, we just scale to gray with a maximum of imax.
Modify this function to make more interesting colors.
*/

int iteration_to_color( int i, int max )
{
	int gray = 255*i/max;
	return MAKE_RGBA(gray,gray,gray,0);
}

