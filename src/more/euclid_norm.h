#ifndef _EUCLID_NORM_H_
#define _EUCLID_NORM_H_

#include <math.h>
#include <stdlib.h>
#include <ctime>

struct norm {
	double diff_norm; //euclidean norm
	double actual_norm; //relative euclidean norm
};

//overcomes numerical errors
static void calc_sum(double *array, int nof_points)
{
	int i, n;
	
	n = nof_points;
	while (n!=1) {
		int is_odd = n%2;
		n /= 2;
		for (i=0; i<n; i++) {
			array[i] = array[2*i] + array[2*i+1];
		}
		if (is_odd) {array[rand()%n] += array[2*n];}
	}
}

//euclidean norm
struct norm calc_norm(const double *actual_values, const double *est_values, int nof_points)
{
	int i;
	struct norm norm;
	double norm_assist[nof_points]; //norm_assist could be of log(nof_points) size,
	                                  //but it's done this way 'cause it's easier
	
	srand(time(NULL));
	
	/* difference */
	for (i=0; i<nof_points; i++) {
		norm_assist[i] = (actual_values[i]-est_values[i]) * (actual_values[i]-est_values[i]);
	}
	calc_sum(norm_assist, nof_points);
	norm.diff_norm = sqrt(norm_assist[0]);
	
	/* actual values */
	for (i=0; i<nof_points; i++) {
		norm_assist[i] = actual_values[i] * actual_values[i];
	}
	calc_sum(norm_assist, nof_points);
	norm.actual_norm = sqrt(norm_assist[0]);
	
	return norm;
}

#endif
