/* PROGRAM: ch3d.c
 * PURPOSE: Module for 3D convex hull algorithm
 * AUTHOR:  Geoffrey Card
 * DATE:    2014-10-13 - 
 * NOTES:   Time: O(n^2*log(n))    Algorithm requires sort(n) for each n, hence n*n*log(n).
 *          Space: O(n)
 *          Untested.
 * TO DO:   Write 2D convex hull for planar case.
 *          Test.
 */

#include "ch3d.h"

/* void CrossProduct(float* a, float* b, float* c)
 * 
 * DESCRIPTION:
 *     Cross product: c = a x b
 * 
 * PARAMETERS:
 *     float*           a        3D vector
 * 
 *     float*           b        3D vector
 * 
 *     float*           c        3D vector, receives result
 *
 * RETURNS:
 *     void
 * 
 * NOTE:
 *     Consider making this inline, due to its small size.
 */
void CrossProduct(float* a, float* b, float* c) {
	c[0] = a[1]*b[2]-a[2]*b[1];
	c[1] = a[2]*b[0]-a[0]*b[2];
	c[2] = a[0]*b[1]-a[1]*b[0];
}

/* float DotProduct(float* a, float* b)
 * 
 * DESCRIPTION:
 *     Dot product: a . b
 * 
 * PARAMETERS:
 *     float*           a        3D vector
 * 
 *     float*           b        3D vector
 *
 * RETURNS:
 *     float                     result
 * 
 * NOTE:
 *     Consider making this inline, due to its small size.
 */
float DotProduct(float* a, float* b) {
	return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

/* void UnitVector(float* u)
 * 
 * DESCRIPTION:
 *     Turns u into a unit vector.
 * 
 * PARAMETERS:
 *     float*           u        3D vector
 *
 * RETURNS:
 *     void
 * 
 * NOTE:
 *     Consider making this inline, due to its small size.
 */
void UnitVector(float* u) {
	float temp;
	
	temp = sqrt(u[0]*u[0] + u[1]*u[1] + u[2]*u[2]);
	u[0] /= temp;
	u[1] /= temp;
	u[2] /= temp;
}

/* float TripleProduct(float* a, float* b, float* c)
 * 
 * DESCRIPTION:
 *     Triple product: a . (b x c) == b . (c x a) == c . (a x b)
 * 
 * PARAMETERS:
 *     float*           a        3D vector
 * 
 *     float*           b        3D vector
 * 
 *     float*           c        3D vector
 *
 * RETURNS:
 *     float                     result
 * 
 * NOTE:
 *     Consider making this inline, due to its small size.
 */
float TripleProduct(float* a, float* b, float* c) {
	return a[0] * (b[1]*c[2]-b[2]*c[1])  +  a[1] * (b[2]*c[0]-b[0]*c[2])  +  a[2] * (b[0]*c[1]-b[1]*c[0]);
}

/* void GeneratePoints(float* p, unsigned int size, int convexset)
 * 
 * DESCRIPTION:
 *     Generates a set of size 3D points. If convexset is TRUE, the points are a convex set.
 * 
 * IMPORTANT:
 *     perform srand before this
 *     assumes p has been PRE-ALLOCATED
 * 
 * PARAMETERS:
 *     float*           p        array for storing points
 * 
 *     unsigned int     size     number of 3D points (actual size of p is 3*size)
 *
 * RETURNS:
 *     void
 */
void GeneratePoints(float* p, unsigned int size, int convexset) {
	float a1 = 0;
	float a2 = 0;
	unsigned int i = 0;
	
	if (convexset == TRUE) {
		// convex set
		for (i = 0; i < size; i++) {
			a1 = 2*PI*((float) rand()/RAND_MAX);
			a2 = asin(2*((float) rand()/RAND_MAX)-1);
			p[i*3+0] = sin(a2)*cos(a1);
			p[i*3+1] = sin(a2)*sin(a1);
			p[i*3+2] = cos(a2);
		}
	} else {
		// random points
		for (i = 0; i < size; i++) {
			p[i*3+0] = ((float) rand()/RAND_MAX);
			p[i*3+1] = ((float) rand()/RAND_MAX);
			p[i*3+2] = ((float) rand()/RAND_MAX);
		}
	}
}

/* void PrintPoints(float* p, unsigned int lenp)
 * 
 * DESCRIPTION:
 *     Prints an indexed list of all points.
 * 
 * PARAMETERS:
 *     float*           p        array of 3D points
 * 
 *     unsigned int     lenp     number of 3D points (actual size of p is 3*lenp)
 *
 * RETURNS:
 *     void
 */
void PrintPoints(float* p, unsigned int lenp) {
	unsigned int i = 0;	
	
	printf("POINTS\n");
	printf("--------+---------------------------\n");
	for (i = 0; i < lenp; i++) {
		printf("%5u   |   %6.2f   %6.2f   %6.2f\n", i, p[i*3+0], p[i*3+1], p[i*3+2]);
	}
	printf("\n");
}

/* void PrintResults(float* p, unsigned int lenp, unsigned int* l, unsigned int lenl, unsigned int leni)
 * 
 * DESCRIPTION:
 *     Prints an indexed list of all points.
 *     Prints a list of all vertex/edge/facet indices and points.
 * 
 * PARAMETERS:
 *     float*           p        array of 3D points
 * 
 *     unsigned int     lenp     number of 3D points (actual size of p is 3*lenp)
 * 
 *     float*           l        array of point indices
 * 
 *     unsigned int     lenl     number of listed objects (actual size of p is 3*lenl*leni)
 * 
 *     unsigned int     leni     number of points per list item (vertex/edge/facet)
 *                               this determines whether the list if of vertices, edges or facets
 *
 * RETURNS:
 *     void
 */
void PrintResults(float* p, unsigned int lenp, unsigned int* l, unsigned int lenl, unsigned int leni) {
	unsigned int i = 0;	
	unsigned int j = 0;
	
	// ... print points
	PrintPoints(p, lenp);
	
	if (leni == 1) {
		printf("VERTICES\n");
	} else if (leni == 2) {
		printf("EDGES\n");
	} else if (leni == 3) {
		printf("FACETS\n");
	} else {
		printf("List contains %u unknown items of length %u.", lenl, leni);
	}

	// seperator
	printf("--------+");
	for (j = 0; j < leni; j++) {
		printf("--------");
	}
	printf("---+");
	for (j = 0; j < leni; j++) {
		printf("---------");
		printf("---------");
		printf("---------");
		printf("---+");	
	}
	printf("\n");
	
	// list
	for (i = 0; i < lenl; i++) {
		// index
		printf("%5u   |", i);
		
		
		// indices
		for (j = 0; j < leni; j++) {
			printf("   %5u", l[i*leni+j]);
		}
		printf("   |");

		// points
		for (j = 0; j < leni; j++) {
			if (l[i*leni+j] < lenp) {
				printf("   %6.2f   %6.2f   %6.2f   |", p[l[i*leni+j]*3+0], p[l[i*leni+j]*3+1], p[l[i*leni+j]*3+2]);
			} else {
				break;
			}
		}
		printf("\n");
	}
	printf("\n");
}

/* void RemoveRepetitions (float** p_p, unsigned int* lenp_p)
 * 
 * DESCRIPTION:
 *     Removes repeated points (which would cause problems for the 3D convex hull algorithm).
 *     Returns new number of points.
 * 
 * IMPORTANT:
 *     I'm quite sure that this has a MEMORY LEAK
 *         but gcc give me errors when I free the memory
 * 
 * PARAMETERS:
 *     float*           p_p      pointer to array for storing points
 * 
 *     unsigned int     lenp_p   pointer to number of 3D points (actual size of p is lenp*3)
 *
 * RETURNS:
 *     void
 */
void RemoveRepetitions (float** p_p, unsigned int* lenp_p) {
	float* p = NULL;
	unsigned int lenp = 0;
	unsigned int* s1 = NULL;
	unsigned int* s2 = NULL;
	unsigned int* temp = NULL;
	unsigned int a = 0;
	unsigned int i = 0;
	unsigned int i1 = 0;
	unsigned int i2 = 0;
	unsigned int lim1 = 0;
	unsigned int lim2 = 0;
	
	// get points
	p = *p_p;
	lenp = *lenp_p;
	
	// initialize indices
	s1 = (unsigned int*) malloc(lenp*sizeof(unsigned int));
	for (i = 0; i < lenp; i++) {
		s1[i] = i;
	}
	// s2 can be garbage, will be replaced
	s2 = (unsigned int*) malloc(lenp*sizeof(unsigned int));
	
	// sort
	a = 1;
	while (a < lenp*2) {
		i = 0;
		while (i < lenp) {
			lim1 = i+1*a < lenp ? i+1*a : lenp;
			lim2 = i+2*a < lenp ? i+2*a : lenp;
			i1 = i;
			i2 = lim1;
			while (i < lim2) {
				if (i1 == lim1) {
					s2[i] = s1[i2];
					i2++;
				} else if (i2 == lim2) {
					s2[i] = s1[i1];
					i1++;
				} else if (p[s1[i1]*3+0] < p[s1[i2]*3+0]) { 
					s2[i] = s1[i1];
					i1++;
				} else if (p[s1[i1]*3+0] > p[s1[i2]*3+0]) {
					s2[i] = s1[i2];
					i2++;
				} else if (p[s1[i1]*3+1] < p[s1[i2]*3+1]) { 
					s2[i] = s1[i1];
					i1++;
				} else if (p[s1[i1]*3+1] > p[s1[i2]*3+1]) { 
					s2[i] = s1[i2];
					i2++;
				} else if (p[s1[i1]*3+2] < p[s1[i2]*3+2]) { 
					i1++;
				} else if (p[s1[i1]*3+2] > p[s1[i2]*3+2]) { 
					s2[i] = s1[i2];
					s2[i] = s1[i1];
					i2++;
				} else {
					s2[i] = s1[i1];
					i1++;
				}
				i++;
			}
		}
		a *= 2;
		// exchange addresses
		temp = s1;
		s1 = s2;
		s2 = temp;
	}

	// repurpose s2 to mark repeats
	// TRUE == repeat
	for (i = 0; i < lenp-1; i++) {
		if (p[s1[i]*3+0] == p[s1[i+1]*3+0] && p[s1[i]*3+1] == p[s1[i+1]*3+1] && p[s1[i]*3+2] == p[s1[i+1]*3+2]) {
			s2[s1[i]] = TRUE;
		} else {
			s2[s1[i]] = FALSE;
		}
	}
	s2[s1[lenp-1]] = FALSE; // this point isn't covered in the loop

	// repurpose a as last valid index
	a = lenp-1;
	for (i = lenp-1; i < lenp; i--) {
		if (s2[i] == TRUE) {
			// place point out of range
			// really, just copy out-of-range point into it's place
			p[i*3+0] = p[a*3+0];
			p[i*3+1] = p[a*3+1];
			p[i*3+2] = p[a*3+2];
			a--;
		}
	}

	free(s1);
	free(s2);
	s1 = NULL;
	s2 = NULL;
	temp = NULL;

	lenp = a+1;
	p = (float*) realloc(p, 3*lenp*sizeof(float));
	*p_p = p;
	*lenp_p = lenp;
}

/* void SortForIndices(float* l1, unsigned int* s1, unsigned int* s2, unsigned int size)
 * 
 * DESCRIPTION:
 *     Sorts s1 according to l1 (ascending). l1 remains unsorted.
 * 
 * IMPORTANT:
 *     s1 MUST CONTAIN {0, 1, ..., size-1}
 * 
 * PARAMETERS:
 *     float*           l1       unsorted array (there used to be an l2, so I kept the name)
 * 
 *     unsigned int*    s1       indices, {0, 1, ..., size-1}
 * 
 *     unsigned int*    s2       empty list (garbage OK)
 * 
 *     unsigned int     size     size of l1, s1, s2
 *
 * RETURNS:
 *     void
 * 
 * NOTES:
 *     Why is s2 pre-allocated when it's not used outside of this function? Because this function is heavily reused with the same size.
 */
void SortForIndices (float* l1, unsigned int* s1, unsigned int* s2, unsigned int size) {
	unsigned int* temps = NULL; // just a pointer, not an array
	unsigned int a = 0;
	unsigned int i = 0;
	unsigned int i1 = 0;
	unsigned int i2 = 0;
	unsigned int lim1 = 0;
	unsigned int lim2 = 0;

	a = 1;
	while (a < size*2) {
		i = 0;
		while (i < size) {
			lim1 = i+1*a < size ? i+1*a : size;
			lim2 = i+2*a < size ? i+2*a : size;
			i1 = i;
			i2 = lim1;
			while (i < lim2) {
				if (i1 == lim1) {
					s2[i] = s1[i2];
					i2++;
				} else if (i2 == lim2) {
					s2[i] = s1[i1];
					i1++;
				} else if (l1[s1[i1]] <= l1[s1[i2]]) {
					s2[i] = s1[i1];
					i1++;
				} else {
					s2[i] = s1[i2];
					i2++;
				}
				i++;
			}
		}
		a *= 2;
		temps = s1;
		s1 = s2;
		s2 = temps;
	}
}

/* unsigned int ConvexityTest (float* v, float* v0, float* v1, float* v2)
 * 
 * DESCRIPTION:
 *     Checks whether v1 is above, in, or below the plane of v0 and v2: v1 . (v0 x v2) > 0
 *     If it's in the plane, it checks whether v1 is farther away from v than the line between v0 and v2: 
 * 
 * PARAMETERS:
 *     float*           v        3D vector, vii - c (if you don't understand that, read the convex hull algorithm)
 * 
 *     float*           v0       3D vector
 * 
 *     float*           v1       3D vector
 * 
 *     float*           v2       3D vector
 *
 * RETURNS:
 *     unsigned int              return CONCAVE, PLANAR, or CONVEX
 * 
 * NOTES:
 *     This function is designed for use inside of ConvexHull3D, nowhere else.
 */
unsigned int ConvexityTest (float* v, float* v0, float* v1, float* v2) {
	float v4[3];
	float v5[3];
	float test;
	
	// check whether v_1 is above (+), in (0), or below (-) the plane of v_0 and v_2
	// test = v1 . (v0 x v2)
	test = TripleProduct(v1, v0, v2);

	if (test > EPSILON) {
		// above the plane
		return CONVEX;
	} else if (abs(test) <= EPSILON) {
		// planar
		// is it farther out than v_0 and v_2 (convex in the plane)?
		// v4 = v2 - v0
		// v5 = v4 x v
		v4[0] = v2[0]-v0[0];
		v4[1] = v2[1]-v0[1];
		v4[2] = v2[2]-v0[2];
		CrossProduct(v4, v, v5);
		// if v5 . v1 > v5 . v0 (or v5 . v2) then it is convex planar
		if (DotProduct(v5, v1) > DotProduct(v5, v1)) {
			return PLANAR;
		} else {
			return CONCAVE;
		}
	} else {
		return CONCAVE;
	}
}

/* void ConvexHull3D(float** p_p, unsigned int* lenp_p, unsigned int** l_p, unsigned int* lenl_p, unsigned int* leni_p, int list_type, int debug)
 * 
 * DESCRIPTION:
 *     Removes repeated points from p.
 *     If list_type is EDGES,  generates a convex hull, of points p, out of points.
 *     If list_type is POINTS, generates a convex hull, of points p, out of edges.
 *     If list_type is FACETS, generates a convex hull, of points p, out of triangluar facets.
 *     Else, list_type treated like FACETS.
 *     If debug is TRUE, prints messages to console to help with debugging.
 * 
 * IMPORTANT:
 *     Planar sets are not yet implemented.
 *     Not tested beyond triganles.
 * 
 * PARAMETERS:
 *     float**          p_p      pointer to array of 3D points, size 3*lenp
 * 
 *     unsigned int*    lenp_p   pointer to number of 3D points (actual size of p is 3*lenp)
 *                               it's a pointer so that it can change the original value
 * 
 *     unsigned int**   l_p      point to array that receives array of indices, size lenl*leni
 * 
 *     unsigned int*    lenl_p   pointer to the length of the list (actual size of l is leni*len_p)
 *                               it's a pointer so that it can change the original value
 * 
 *     unsigned int*    leni_p   pointer to the length of the list type
 *                               1 if it's points, 2 if it's edges, 3 if it's facets
 *                               it's a pointer so that it can change the original value
 *
 * RETURNS:
 *     void
 */
void ConvexHull3D(float** p_p, unsigned int* lenp_p, unsigned int** l_p, unsigned int* lenl_p, unsigned int* leni_p, int list_type, int debug) {
	// arrays and lengths
	float* p = NULL;
	unsigned int lenp = 0;
	unsigned int* l = NULL;
	unsigned int leni = 0;
	unsigned int lenl = 0;
	
	// queue
	int* visited = NULL;
	unsigned int* q = NULL;
	unsigned int q1 = 0;
	unsigned int q2 = 0;
	
	// points and vectors
	float c[3] = {0};
	float vii[3] = {0};
	float u0[3] = {0};
	float u1[3] = {0};
	float* v = NULL;
	
	// lists
	float* theta = NULL;
	unsigned int* s = NULL;
	unsigned int* s2 = NULL;
	unsigned int* t = NULL;
	unsigned int tcount = 0;
	float* n = NULL;
	unsigned int lenn1 = 0;
	unsigned int lenn2 = 0;
	
	// indices
	unsigned int i = 0;
	unsigned int ii = 0;
	unsigned int i9 = 0;
	unsigned int i0 = 0;
	unsigned int i1 = 0;
	unsigned int i2 = 0;
	
	// ????
	float tempf = 0;
	
	// bools for checks
	unsigned int chk1 = 0;
	unsigned int chk2 = 0;
	
	if (debug == TRUE) 	printf("DEBUG: ConvexHull3D: BEGIN\n");
	
	// get points and lengths
	p = *p_p;
	lenp = *lenp_p;
	leni = *leni_p;
	lenl = *lenl_p;
	
	if (debug == TRUE) 	printf("DEBUG: ConvexHull3D: RemoveRepetitions\n");
	// get rid of pesky repetitions get rid of pesky repetitions
	RemoveRepetitions(&p, &lenp);
	*p_p = p;
	*lenp_p = lenp;

	//
	// POINT
	//
	if (lenp == 1) {
		if (debug == TRUE) 	printf("DEBUG: ConvexHull3D: Point\n");
		lenl = 1;
		leni = 1;
		l = (unsigned int*) malloc(leni*lenl*sizeof(unsigned int));
		l[0] = 0;
		*leni_p = leni;
		*lenl_p = lenl;
		*l_p = l;
		return;
	}
	
	//
	// LINE
	//
	if (lenp == 2) {
		if (debug == TRUE) 	printf("DEBUG: ConvexHull3D: Line\n");
		if (list_type == VERTICES) {
			lenl = 2;
			leni = 1;
		} else {
			lenl = 1;
			leni = 2;
		}
		l = (unsigned int*) malloc(leni*lenl*sizeof(unsigned int));
		l[0] = 0;
		l[1] = 1;
		*leni_p = leni;
		*lenl_p = lenl;
		*l_p = l;
		return;
	}
	
	if (debug == TRUE) 	printf("DEBUG: ConvexHull3D: Non-Linearity Test\n");
	// NON-LINEARITY TEST
	// vi = pi - p0
	// for i:
	// 	if v0 x vi != 0:
	//		non-linear
	// 		break
	// repurposeing chk1, u0, u1
	chk1 = TRUE;
	u0[0] = p[1*3+0] - p[0*3+0];
	u0[1] = p[1*3+1] - p[0*3+1];
	u0[2] = p[1*3+2] - p[0*3+2];
	for (i = 2; i < lenp; i++) {
		u1[0] = p[i*3+0] - p[0*3+0];
		u1[1] = p[i*3+1] - p[0*3+1];
		u1[2] = p[i*3+2] - p[0*3+2];
		// if the cross product is NOT ZERO
		// then it's NON-LINEAR
		if (abs(u0[1]*u1[2]-u1[1]*u0[2]) > EPSILON || 
		    abs(u0[2]*u1[0]-u1[2]*u0[0]) > EPSILON || 
		    abs(u0[0]*u1[1]-u1[0]*u0[1]) > EPSILON) {
			chk1 = FALSE;
			break;
		}
	}
	
	//
	// LINEAR
	//
	if (chk1 == TRUE) {
		if (debug == TRUE) 	printf("DEBUG: ConvexHull3D: Linear\n");
		if (list_type == VERTICES) {
			lenl = 2;
			leni = 1;
		} else {
			lenl = 1;
			leni = 2;
		}
		l = (unsigned int*) malloc(leni*lenl*sizeof(unsigned int));

		// find extremes
		// repurposing q1, q2 as indices for extremes
		// repurposing c[0], c[1] as values of extremes
		q1 = 0;
		c[0] = 0;
		q2 = 1;
		c[1] = DotProduct(u0, u0);
		for (i = 2; i < lenp; i++) {
			u1[0] = p[i*3+0] - p[0*3+0];
			u1[1] = p[i*3+1] - p[0*3+1];
			u1[2] = p[i*3+2] - p[0*3+2];
			// get relative position
			tempf = DotProduct(u0, u1);
			if (tempf < c[0]) {
				// if farthest negative
				q1 = i;
				c[0] = tempf;
			} else if (tempf > c[1]) {
				// if farthest forward
				q2 = i;
				c[1] = tempf;
			}
		}
		
		l[0] = q1;
		l[1] = q2;
		*leni_p = leni;
		*lenl_p = lenl;
		*l_p = l;
		return;
	}
	
	//
	// TRIANGLE
	//
	if (lenp == 3) {
		if (debug == TRUE) 	printf("DEBUG: ConvexHull3D: Triangle\n");
		if (list_type == VERTICES) {
			lenl = 3;
			leni = 1;
			l = (unsigned int*) malloc(leni*lenl*sizeof(unsigned int));
			l[0] = 0;
			l[1] = 1;
			l[2] = 2;
		} else if (list_type == EDGES) {
			lenl = 3;
			leni = 2;
			l = (unsigned int*) malloc(leni*lenl*sizeof(unsigned int));
			l[0] = 0;
			l[1] = 1;
			l[2] = 0;
			l[3] = 2;
			l[4] = 1;
			l[5] = 2;
		} else {
			lenl = 1; // from a pure maths perspective, this should be two (one for each side), but frak it
			leni = 3;
			l = (unsigned int*) malloc(leni*lenl*sizeof(unsigned int));
			l[0] = 0;
			l[1] = 1;
			l[2] = 2;
		}
		*leni_p = leni;
		*lenl_p = lenl;
		*l_p = l;
		return;
	}
	
	if (debug == TRUE) 	printf("DEBUG: ConvexHull3D: Non-Planarity Test\n");
	// NON-PLANARITY TEST
	//	if (v0 x v1) . vi != 0:
	//		non-planar
	// 		break
	// repurposeing chk1, u0, u1, c
	chk1 = TRUE;
	u0[0] = p[1*3+0] - p[0*3+0];
	u0[1] = p[1*3+1] - p[0*3+1];
	u0[2] = p[1*3+2] - p[0*3+2];
	u1[0] = p[2*3+0] - p[0*3+0];
	u1[1] = p[2*3+1] - p[0*3+1];
	u1[2] = p[2*3+2] - p[0*3+2];
	c[0]  = u0[1]*u1[2]-u1[1]*u0[2];
	c[1]  = u0[2]*u1[0]-u1[2]*u0[0];
	c[2]  = u0[0]*u1[1]-u1[0]*u0[1];
	for (i = 3; i < lenp; i++) {
		// if |vi.v(i)| > 0
		// then it's NON-PLANAR 
		tempf = (p[i*3+0] - p[0*3+0])*c[0] + p[i*3+1] - p[0*3+1]*c[1] +  p[i*3+2] - p[0*3+2]*c[2];
		if (abs(tempf) > EPSILON) {
			chk1 = FALSE;
			break;
		}
	}
	
	//
	// PLANE
	//
	if (chk1 == TRUE) {
		if (debug == TRUE) 	printf("DEBUG: ConvexHull3D: Planar\n");
		// needs a 2D convexity test
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		
		*leni_p = leni;
		*lenl_p = lenl;
		*l_p = l;
		return;
	}
	
	//
	// TETRAHEDRON
	//
	if (lenp == 4) {
		if (debug == TRUE) 	printf("DEBUG: ConvexHull3D: Tetrahedron\n");
		if (list_type == VERTICES) {
			lenl = 4;
			leni = 1;
			l = (unsigned int*) malloc(leni*lenl*sizeof(unsigned int));
			l[0] = 0;
			l[1] = 1;
			l[2] = 2;
			l[3] = 3;
		} else if (list_type == EDGES) {
			lenl = 6;
			leni = 2;
			l = (unsigned int*) malloc(leni*lenl*sizeof(unsigned int));
			l[ 0] = 0;
			l[ 1] = 1;
			l[ 2] = 0;
			l[ 3] = 2;
			l[ 4] = 0;
			l[ 5] = 3;
			l[ 6] = 1;
			l[ 7] = 2;
			l[ 8] = 1;
			l[ 9] = 3;
			l[10] = 2;
			l[11] = 3;
		} else {
			lenl = 1;
			leni = 3;
			l = (unsigned int*) malloc(leni*lenl*sizeof(unsigned int));
			l[0] = 0;
			l[1] = 1;
			l[2] = 2;
		}
		*leni_p = leni;
		*lenl_p = lenl;
		*l_p = l;
		return;
	}
	
	//
	// CONVEX HULL
	//
	if (debug == TRUE) 	printf("DEBUG: ConvexHull3D: Convex Hull\n");
		
	// Determine maximum possible number of edges (E) and facets (F)
	// E <= V + (V-1) + (V-2) + ... + 0 = V*(V-1)/2
	// F = 2 + E - V
	// F <= 2 + V*(V-1)/2 - V
	// F <= 2+V*(V-3)/2
	// Then I did some maths and discovered that if you're building polyhedra from a minimal number of triangles:
	// E = 3*V-6
	// F = 2*V-4
	if (list_type == VERTICES) {
		leni = 1;
		lenl = lenp;
		l = (unsigned int*) malloc(leni*lenl*sizeof(unsigned int));		
	} else if (list_type == EDGES) {
		leni = 2;
		lenl = leni*(3*lenp-6);
		l = (unsigned int*) malloc(leni*lenl*sizeof(unsigned int));
	} else {
		leni = 3;
		lenl = leni*(2*lenp-4);
		l = (unsigned int*) malloc(leni*lenl*sizeof(unsigned int));
	}
	lenl = 0; // becomes incrementor for l
	

	// allocate everything
	visited = (int*) malloc(lenp*sizeof(int));
	q = (unsigned int*) malloc(lenp*sizeof(unsigned int));
	v = (float*) malloc(3*lenp*sizeof(float));
	theta = (float*) malloc(lenp*sizeof(float));
	s  = (unsigned int*) malloc(lenp*sizeof(unsigned int));
	s2 = (unsigned int*) malloc(lenp*sizeof(unsigned int));
	t  = (unsigned int*) malloc(lenp*sizeof(unsigned int));
	n = (float*) malloc(3*(2*lenp-4)*sizeof(float));

	// list of used variables
	// calculate approximate centroid and get first point
	c[0] = p[0];
	c[1] = p[1];
	c[2] = p[2];
	// find index of the lowerest point p_x_min
	ii = 0;
	for (i = 1; i < lenp; i++) {
		// visited
		visited[i] = UNVISITED;
		// center
		c[0] += p[i*3+0];
		c[1] += p[i*3+1];
		c[2] += p[i*3+2];
		// find min x
		if (p[i*3+0] < p[ii*3+0]) ii = i;
	}
	c[0] /= lenp;
	c[1] /= lenp;
	c[2] /= lenp;
	
	
	// push first point
	visited[ii] = QUEUED;
	q1 = 0;
	q2 = 0;
	q[q2++] = ii;
	
	// set list of normals, n, as empty
	lenn1 = 0;
	lenn2 = 0;
	
	// while the queue is not empty
	while (q1 < q2) {
		// pull point
		ii = q[q1++];

		if (debug == TRUE) 	printf("DEBUG: ConvexHull3D: POINT LOOP %u\n", ii);
		
		// get directional vectors
		vii[0] = p[ii*3+0] - c[0];
		vii[1] = p[ii*3+1] - c[1];
		vii[2] = p[ii*3+2] - c[2];
		if (abs(vii[0]) <= EPSILON) {
			u0[0] = 1;
			u0[1] = 0;
			u0[2] = 0;
		} else if (abs(vii[1]) <= EPSILON) {
			u0[1] = 1;
			u0[1] = 0;
			u0[2] = 0;
		} else {
			u0[0] = vii[1];
			u0[1] = -vii[0];
			u0[2] = 0;
			// this verison of UnitVector is slightly more optimized
			// don't both using the function here
			tempf = sqrt(u0[0]*u0[0] + u0[1]*u0[1]);
			u0[0] /= tempf;
			u0[1] /= tempf;
		}
		// u1 = v x u0
		CrossProduct(v, u0, u1);
		UnitVector(u1);
		
		// create a list of vectors outward from v
		// calculate their angle around v
		// create a list of test results
		for (i = 0; i < lenp; i++) {
			v[i*3+0] = p[i*3+0] - p[ii*3+0];
			v[i*3+1] = p[i*3+1] - p[ii*3+1];
			v[i*3+2] = p[i*3+2] - p[ii*3+2];
			theta[i] = atan2(
				u1[0]*v[i*3+0] + u1[1]*v[i*3+1] + u1[2]*v[i*3+2],
				u0[0]*v[i*3+0] + u0[1]*v[i*3+1] + u0[2]*v[i*3+2]);
			t[i] = CONVEX;
		}
		// create list of indices of sorted theta
		SortForIndices(theta, s, s2, lenp);
		// and mark vii as failed
		t[ii] = CONCAVE;
		
		// tcount is a failsafe
		// it counts the number of failures
		// if tcount >= lenp, there will be an infinite loop
		// therefore break, and skip to the end
		tcount = 1;
					
		if (debug == TRUE) 	{
			printf("DEBUG: ConvexHull3D: POINT LOOP %u\n", ii);
			for (i = 0; i < lenp; i++) {
				printf(" %4d", i);
			}
			printf("    | ");
			for (i = 0; i < lenp; i++) {
				printf(" %4d", s[i]);
			}
			printf("\n");
			// separator
			printf(" ");
			for (i = 0; i < lenp; i++) {
				printf("----");
			}
			printf("---+");
			for (i = 0; i < lenp; i++) {
				printf("----");
			}
			printf("\n");
		}
			
		// create a set of sequential indices
		// only count values that haven't failed the test
		//  i_???? = i9
		// 	i_prev = i0
		// 	i_curr = i1
		// 	i_next = i2
		i0 = -1;
		while (t[s[i0]] == CONCAVE) i0--;
		i9 = i0-1;
		while (t[s[i9]] == CONCAVE) i9--;
		i1 = 0;
		while (t[s[i1]] == CONCAVE) i1++;
		i2 = i1+1;
		while (t[s[i2]] == CONCAVE) i2++;
		// go through the whole list once
		// check forward and backward at each failure
		// (until success in both directions)
		while (i1 < lenp) {
			// test v_1
			t[s[i1]] = ConvexityTest(vii, &v[s[i0%lenp]*3], &v[s[i1%lenp]*3], &v[s[i2%lenp]*3]);

			if (t[s[i1%lenp]] == CONCAVE) {
				// failsafe
				tcount++;
				if (tcount >= lenp) {
					break;
					break;
				}
				// concavity detected

				// increment forward
				i1 = i2;
				i2++;
				while (t[s[i2%lenp]] == CONCAVE) i2++;

				// check forward and backward until coast is clear
				chk1 = CONCAVE;
				chk2 = CONCAVE;
				while (chk1 == CONCAVE || chk2 == CONCAVE) {
					// check backward
					chk1 = ConvexityTest(vii, &v[s[i9%lenp]*3], &v[s[i0%lenp]*3], &v[s[i1%lenp]*3]);
					if (chk1 == CONCAVE) {
						t[s[i0%lenp]] = CONCAVE;
						// failsafe
						tcount++;
						if (tcount >= lenp) {
							break;
							break;
							break;
						}
						i0 = i9;
						i9--;
						while (t[s[i9%lenp]] == CONCAVE) i9--;
					}
					// check forward
					chk2 = ConvexityTest(vii, &v[s[i0%lenp]*3], &v[s[i1%lenp]*3], &v[s[i2%lenp]*3]);
					if (chk2 == CONCAVE) {
						t[s[i1%lenp]] = CONCAVE;
						// failsafe
						tcount++;
						if (tcount >= lenp) {
							break;
							break;
							break;
						}
						i1 = i2;
						i2++;
						while (t[s[i2%lenp]] == CONCAVE) i2++;
					}
				}
				// concavity removed
			} else {
				// no concavity, continue
				i9 = i0;
				i0 = i1;
				i1 = i2;
				i2++;
				while (t[s[i2%lenp]] == CONCAVE) i2++;
			}
			if (debug == TRUE) {
				// test results, of this round
				for (i = 0; i < lenp; i++) {
					printf("%4u", t[i]);
				}
				printf("    | ");
				for (i = 0; i < lenp; i++) {
					printf("%4u", t[s[i]]);
				}
				printf("\n");
			}
		}
		
		if (debug == TRUE) {
			// test results
			for (i = 0; i < lenp; i++) {
				printf("%4u", t[i]);
			}
			printf("    | ");
			for (i = 0; i < lenp; i++) {
				printf("%4u", t[s[i]]);
			}
			printf("\n");
		}
		
		// failsafe
		tcount++;
		if (tcount >= lenp) {
			break;
		}

		// create vertices, edges, facets
		if (list_type == VERTICES) {
			// VERTICES
			// list this vertex, push the rest to the queue (they'll get added on their turn
			l[lenl*1+0] = ii;
			lenl ++;
			for (i = 0; i < lenp; i++) {
				// push all unvisited points to the queue
				if (t[i] != CONCAVE && visited[i] == UNVISITED) {
					visited[i] = QUEUED;
					q[q2++] = i;
				}
			}
		} else if (list_type == EDGES) {
			// EDGES
			// list edges
			for (i = 0; i < lenp; i++) {
				// if a point has been visited, the line already exisits
				// if a point is planar, it's line is unnecessary
				if (t[i] == CONVEX && visited[i] != VISITED) {
					l[lenl*2+0] = ii;
					l[lenl*2+1] = i;
					lenl++;
					if (visited[i] != QUEUED) {
						visited[i] = QUEUED;
						q[q2++] = i;
					}
				}
			}
		} else {
			// FACETS
			// create a triangle fan of points that pass
			i1 = 0;
			while (t[s[i1%lenp]] == CONCAVE) i1++;
			i2 = i1+1;
			while (t[s[i2%lenp]] == CONCAVE) i2++;
			while (i1 < lenp) {
				// if a triangle contains a visited point, it's already in the list
				if (visited[s[i1%lenp]] != VISITED && visited[s[i2%lenp]] != VISITED) {
					// if planar, confirm that the normal is not yet disallowed
					chk1 = TRUE;
					if (t[s[i1%lenp]] == PLANAR || t[s[i2%lenp]] == PLANAR) {
						CrossProduct(&v[s[i1%lenp]], &v[s[i1%lenp]], u0);
						UnitVector(u0);
						for (i = 0; i < lenn1; i++) {
							if (abs(u0[0]-n[i*3+0]) <= EPSILON && 
							    abs(u0[1]-n[i*3+1]) <= EPSILON && 
							    abs(u0[2]-n[i*3+2]) <= EPSILON) {
							    chk1 = FALSE;
							    break;
							}
						}
						// if planar, but normal not disallowed
						// add point, and check if normal is listed
						// if not, add it
						chk2 = FALSE;
						for (i = lenn1; i < lenn2; i++) {
							if (abs(u0[0]-n[i*3+0]) <= EPSILON && 
							    abs(u0[1]-n[i*3+1]) <= EPSILON && 
							    abs(u0[2]-n[i*3+2]) <= EPSILON) {
							    chk2 = TRUE;
							    break;
							}
						}
						if (chk2 == FALSE) {
							n[lenn2*3+0] = u0[0];
							n[lenn2*3+1] = u0[1];
							n[lenn2*3+2] = u0[2];
							lenn2++;
						}
					}
					if (chk1 == TRUE) {
						l[lenl*3+0] = ii;
						l[lenl*3+1] = s[i1%lenp];
						l[lenl*3+2] = s[i2%lenp];
						lenl++;
						if (visited[s[i1%lenp]] != QUEUED) {
							visited[s[i1%lenp]] = QUEUED;
							q[q2++] = s[i1%lenp];
						}
						if (visited[s[i2%lenp]] != QUEUED) {
							visited[s[i2%lenp]] = QUEUED;
							q[q2++] = s[i2%lenp];
						}
					}
				}
				i1 = i2;
				i2++;
				while (t[s[i2%lenp]] == CONCAVE) i2++;
			}
		}
		// make any new normals off limits
		lenn1 = lenn2;
		// mark point as off limits
		visited[ii] = VISITED;
	}
	
	if (debug == TRUE && tcount >= lenp) {
		printf("DEBUG: ConvexHull3D: FAILSAFE ENGAGED\n");
	}
	
	// free everything
	free(visited);
	free(q);
	free(v);
	free(v);
	free(theta);
	free(s);
	free(s2);
	free(t);
	free(n);
	
	// resize, in case of sub-optimal/non-convex set
	l = (unsigned int*) realloc(l, leni*lenl*sizeof(unsigned int));
	*leni_p = leni;
	*lenl_p = lenl;
	*l_p = l;

	if (debug == TRUE) 	printf("DEBUG: ConvexHull3D: END\n");
}
