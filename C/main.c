/* PROGRAM: main.c
 * PURPOSE: Module for 3D convex hull algorithm
 * AUTHOR:  Geoffrey Card
 * DATE:    2014-10-16 - 
 * NOTES:   Time: O(n^2*log(n))    Algorithm requires sort(n) for each n, hence n*n*log(n).
 *          Space: O(n)
 * TO DO:   
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "ch3d.h"

int main(void) {
	float* p = NULL;
	unsigned int lenp = 0;
	unsigned int* l = NULL;
	unsigned int lenl = 0;
	unsigned int leni = 0;
	unsigned int i = 0;
	
	printf("\n");
	
	// seed random with time
	srand((int) time(NULL));

	lenp = 10;
	p = (float*) malloc(3*lenp*sizeof(float));
	
	GeneratePoints(p, lenp, FALSE);
	/*
	for (i = 0; i < lenp; i++) {
		//p[i*3+0] = 0;
		p[i*3+1] = 0.5*p[i*3+0];
		//p[i*3+2] = 0;
	}
	
	p[0*3+0] = -1;
	p[0*3+1] =  0;
	p[0*3+2] =  0;
	p[1*3+0] =  1;
	p[1*3+1] =  0;
	p[1*3+2] =  0;
	p[2*3+0] =  0;
	p[2*3+1] = -1;
	p[2*3+2] =  0;
	p[3*3+0] =  0;
	p[3*3+1] =  1;
	p[3*3+2] =  0;
	p[4*3+0] =  0;
	p[4*3+1] =  0;
	p[4*3+2] = -1;
	p[5*3+0] =  0;
	p[5*3+1] =  0;
	p[5*3+2] =  1;
	*/
	
	printf("%p %u", p, lenp);
	printf("\n\n\n");

	PrintPoints(p, lenp);
	printf("\n\n\n");

	ConvexHull3D(&p, &lenp, &l, &lenl, &leni, VERTICES);
	PrintResults(p, lenp, l, lenl, leni);
	printf("\n\n\n");
	free(l);
	l = NULL;
	
	ConvexHull3D(&p, &lenp, &l, &lenl, &leni, EDGES);
	PrintResults(p, lenp, l, lenl, leni);
	printf("\n\n\n");
	free(l);
	l = NULL;

	ConvexHull3D(&p, &lenp, &l, &lenl, &leni, FACETS);
	PrintResults(p, lenp, l, lenl, leni);
	printf("\n\n\n");
	free(l);
	l = NULL;
	
	free(p);
	p = NULL;
	lenp = 0;
	
	return 0;
}
