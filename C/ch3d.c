/* PROGRAM: CH3D.py
 * PURPOSE: Module for 3D convex hull algorithm
 * AUTHOR:  Geoffrey Card
 * DATE:    2014-10-13
 * NOTES:   Time: O(n^2*log(n))    Algorithm requires sort(n) for each n, hence n*n*log(n).
 *          Space: O(n)
 *          Has trouble with non-triangluar facets.
 *          And points in the middle of facets.
 * TO DO:   - ConvexHull3D is incomplete. Follow the description.
 *          - Figure out what to do with ConvexityTest, now that we're testing for normals.
 *          - Make all non-trivial memory allocation happen in the main algorithm, it'll save reallocating it.
 *          - And don't forget to group initialization for loops.
 *          - And mention all this in each function's description.
 *          - Need to write reallocation parts (look for Python "del" code).
 */

#include <math.c>
#include <time.c>

#define PI 3.14159

#define TRUE  1
#define FALSE 0

#define UNVISITED 0
#define QUEUED    1
#define VISITED   2

#define VERTICES 0
#define EDGES    1
#define FACETS   2

/* void GeneratePoints(float* p, unsigned int count, int convexset)
 * 
 * DESCRIPTION:
 *     Generates a set of count points. If convexset is TRUE, the points are a convex set.
 * 
 * IMPORTANT:
 *     perform srand before this
 *     p MUST BE PRE-ALLOCATED with size 3*count
 *     p IS NOT FREED IN THIS FUNCTION
 * 
 * PARAMETERS:
 *     float*           p        array for storing points
 *     unsigned int     count    number of 3D points (actual size of p is count*3)
 *
 * RETURNS:
 *     void
 */
void GeneratePoints(float* p, unsigned int count, int convexset) {
	float a1 = 0;
	float a2 = 0;
	unsigned int i = 0;
	
	if (convexset == TRUE) {
		// convex set
		for (i = 0; i < count; i++) {
			a1 = 2*PI*((float) rand()/RAND_MAX);
			a2 = asin(2*((float) rand()/RAND_MAX)-1);
			p[i*3+0] = sin(a2)*cos(a1);
			p[i*3+1] = sin(a2)*sin(a1);
			p[i*3+2] = cos(a2);
		}
	} else {
		// random points
		for (i = 0; i < count; i++) {
			p[i*3+0] = ((float) rand()/RAND_MAX);
			p[i*3+1] = ((float) rand()/RAND_MAX);
			p[i*3+2] = ((float) rand()/RAND_MAX);
		}
	}
}

/* unsigned int RemoveRepetitions (float* p, unsigned int count)
 * 
 * DESCRIPTION:
 *     Removes repeated points (which would cause problems for the 3D convex hull algorithm).
 *     Returns new number of points.
 * 
 * IMPORTANT:
 *     p MUST BE PRE-ALLOCATED with size 3*count
 *     p IS NOT FREED IN THIS FUNCTION
 * 
 * PARAMETERS:
 *     float*           p        array for storing points
 *     unsigned int     count    number of 3D points (actual size of p is count*3)
 *
 * RETURNS:
 *     unsigned int              new value for count
 */
unsigned int RemoveRepetitions (float* p, unsigned int size) {
	unsigned int* s1 = NULL;
	unsigned int* s2 = NULL;
	unsigned int* temp = NULL;
	unsigned int a = 0;
	unsigned int i = 0;
	unsigned int i1 = 0;
	unsigned int i2 = 0;
	unsigned int lim1 = 0;
	unsigned int lim2 = 0;
	
	s1 = (unsigned int*) malloc(count*sizeof(unsigned int));
	for (i = 0; i < size; i++) {
		s1[i] = i;
	}
	
	// s2 can be garbage, will be replaced
	s2 = (unsigned int*) malloc(count*sizeof(unsigned int));
	
	// temp doesn't require data, you'll see later
	
	// merge(ish) sort
	a = 1;
	while (a < count*2) {
		i = 0
		while (i < count) {
			lim1 = i+1*a < size ? i+1*a : size;
			lim1 = i+2*a < size ? i+2*a : size;
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
					s2[i] = s1[i1];
					i1++;
				} else if (p[s1[i1]*3+2] > p[s1[i2]*3+2]) { 
					s2[i] = s1[i2];
					i2++;
				} else {
					s2[i] = s1[i1];
					i1++;
				}
				i++;
			}
		}
		a *= 2;
		// exchange addresses, cool trick
		temp = s1;
		s1 = s2;
		s2 = temp;
	}

	// repurpose s2 to mark repeats
	// TRUE == repeat, i.e. remove this
	for (i = 0; i < size-1; i++) {
		if (p[s1[i]*3+0] == p[s1[i+1]*3+0] && p[s1[i]*3+1] == p[s1[i+1]*3+1] && p[s1[i]*3+2] == p[s1[i+1]*3+2]) {
			s2[s1[i]] = TRUE;
		} else {
			s2[s1[i]] = FALSE;
		}
	}
	s2[s1[count-1]] = FALSE; // this point isn't covered in the loop

	// repurpose a as new size-1
	a = count-1;
	for (i = count-1; i < count; i--) { // i < size still works, because i is unsigned, i.e. overflow
		if (s2[i] == TRUE) {
			// place point out of range
			// really, just copy out-of-range point into it's place
			p[i*3+0] = p[a*3+0];
			p[i*3+1] = p[a*3+1];
			p[i*3+2] = p[a*3+2];
			a--;
		}
	}
	count = a+1;
	p = (float*) realloc(p, count*sizeof(float));
	free(s1);
	free(s2);
	return count;
}

/* void SortForIndices(float* l1, unsigned int* s1, unsigned int* s2, unsigned int size)
 * 
 * DESCRIPTION:
 *     Sorts s1 according to l1 (ascending). l1 remains unsorted.
 * 
 * IMPORTANT:
 *     l1, s1, s2 MUST BE PRE-ALLOCATED with size size
 *     l1, s1, s2 ARE NOT FREED IN THIS FUNCTION
 *     s1 MUST CONTAIN [0, 1, ..., size-1]
 * 
 * PARAMETERS:
 *     float*           l1       unsorted array (there used to be an l2, so I kept the name)
 *     unsigned int*    s1       indices, [0, 1, ..., size-1]
 *     unsigned int*    s2       empty list (garbage OK)
 *     unsigned int     size     size of l1, s1, s2
 *
 * RETURNS:
 *     void
 * 
 * NOTES:
 *     Why is s2 pre-allocated when it's not used outside of this function? Because this function is heavily reused with the same size.
 */
void SortForIndices (float* l1, float* l2, unsigned int* s1, unsigned int* s2, unsigned int size) {
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
		while i < size:
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
	return s1;
}

/* int ConvexityTest (float* v, float* v_0, float* v_1, float* v_2)
 * 
 * DESCRIPTION:
 *     Checks whether v_1 is above the plane of v_0 and v_2: v_1 . (v_0 x v_2) > 0
 *     If it's in the plane, it checks whether v_1 is farther away from v than the line between v_0 and v_2: 
 * 
 * IMPORTANT:
 * 
 * PARAMETERS:
 *     float*           v        3D vector, vii - c (if you don't understand that, read the convex hull algorithm)
 *     float*           v_0      3D vector
 *     float*           v_1      3D vector
 *     float*           v_2      3D vector
 *
 * RETURNS:
 *     int                       I'm reconsidering this, a float might work better
 * 
 * NOTES:
 *     need to check that float* v works for contant length array
 */
int ConvexityTest (float* v, float* v_0, float* v_1, float* v_2) {
		float v4[3];
		float test = 0;
		
		// test = (v_0 x v_2) . v_1
		test = (v_0[1]*v_2[2] - v_0[2]*v_2[1]) * v_1[0] + (v_0[2]*v_2[0] - v_0[0]*v_2[2]) * v_1[1] + (v_0[0]*v_2[1] - v_0[1]*v_2[0]) * v_1[2];
		
		// if it's in the same plane as v_0 and v_2
		// check whether it extends past them
		if (test == 0) {
			// v4 = (v_2 - v_0) x v
			// if v4 . v_1 > v_4 . v_0 (or v_4 . v_2)
			v4[0] = v_2[0]-v_0[0];
			v4[1] = v_2[1]-v_0[1];
			v4[2] = v_2[2]-v_0[2];
			v4[0] = v4[1]*v[2]-v4[2]*v[1];
			v4[1] = v4[2]*v[0]-v4[0]*v[2];
			v4[2] = v4[0]*v[1]-v4[1]*v[0];
			test = v4[0]*v_1[0] + v4[1]*v_1[1] + v4[2]*v_1[2] > v4[0]*v_0[0] + v4[1]*v_0[1] + v4[2]*v_0[2];
		}
		// test MUST be greater than 0, else this is not part of the minimal convex set
		return test > 0 ? TRUE : FALSE;
}

/* void ConvexHull3D(float* p, unsigned int* len_p, float* l, unsigned int* len_l, unsigned int* len_i, int list_type, int debug)
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
 *     p MUST BE PRE-ALLOCATED with size 3*len_p
 *     l MUST BE UNALLOCATED, it is allocated in this function
 *     p, l are NOT FREED IN THIS FUNCTION
 * 
 * PARAMETERS:
 *     float*           p        array of 3D points, size 3*len_p
 *     unsigned int*    len_p    pointer to number of 3D points (actual size of p is 3*len_p)
 *                               it's a pointer so that it can change the original value
 *     unsigned int*    l        array of 3D points, size len_i*len_l
 *     unsigned int*    len_l    pointer to the length of the list (actual size of l is len_i*len_p)
 *                               it's a pointer so that it can change the original value
 *     unsigned int*    len_i    pointer to the length of the list type
 *                               1 if it's points, 2 if it's edges, 3 if it's facets
 *                               it's a pointer so that it can change the original value
 *
 * RETURNS:
 *     void
 * 
 * TODO:
 *     Find names for pointers to sizes. Find better names for sizes. Fix everything.
 * 
 * NOTES:
 *     While not yet implemented, I'm going to write in a section that calculates the normals of each test, and marks points that are in non-triangular facets.
 *     If a new facet is found to have marked points, it's normal is compared to a list of non-triangular facet normals.
 *     The list of points adds n memory, the list of facet normals adds 3*(2*V-4). Still linear.
 */
void ConvexHull3D(float* p, unsigned int size, int facets, int debug) {
	unsigned int* l = NULL;
	unsigned int len_i = 0;
	unsigned int len_l = 0;
	
	// queue
	unsigned int* q = NULL;
	unsigned int q1 = 0;
	unsigned int q2 = 0;
	int* visited = NULL;
	
	float c[3];
	float vii[3];
	float* v = NULL;
	float* theta1 = NULL;
	float* theta2 = NULL;
	unsigned int* s1 = NULL;
	unsigned int* s2 = NULL;
	
	unsigned int i = 0;
	unsigned int ii = 0;
	unsigned int i9 = 0;
	unsigned int i0 = 0;
	unsigned int i1 = 0;
	unsigned int i2 = 0;
	
	float tempf = 0;
	
	int chk_f = FALSE;
	int chk_b = FALSE;
	
	// get rid of pesky repetitions get rid of pesky repetitions
	size = RemoveRepetitions(p3);
	
	// confirm non-linear
	// confirm non-planar
	// vi = pi - p0
	// for i:
	// 	if v0 x vi != 0:
	//		print 'non-linear'
	//	if (v0 x v1) . vi != 0:
	//		print 'non-planar'
	
	// for size < 5, it doesn't seem worth it to check for repeats
	if (size == 1) {
		// point
		len_l = 1;
		len_i = 1;
		l = (float*) malloc(len_i*len_l*sizeof(float));
	} else if (size == 2) {
		// line
		len_l = 1;
		len_i = 2;
		l = (float*) malloc(len_i*len_l*sizeof(float));
		l[0] = 0;
		l[1] = 1;
	} else if (size == 3) {
		// triangle
		if (facets == TRUE) {
			len_l = 1;
			len_i = 3;
			l = (float*) malloc(len_i*len_l*sizeof(float));
			l[0] = 0;
			l[1] = 1;
			l[2] = 2;
		} else {
			len_l = 3;
			len_i = 2;
			l = (float*) malloc(len_i*len_l*sizeof(float));
			l[len_i*0+0] = 0;
			l[len_i*0+1] = 1;
			l[len_i*1+0] = 0;
			l[len_i*1+1] = 2;
			l[len_i*2+0] = 1;
			l[len_i*2+1] = 2;
		}
	} else if (size == 4) {
		// quadralateral or tetrahedron
		// frak it, simple solution works every time
		if (facets == TRUE) {
			len_l = 4;
			len_i = 3;
			l = (float*) malloc(len_i*len_l*sizeof(float));
			l[ 0] = 0;
			l[ 1] = 1;
			l[ 2] = 2;
			l[ 3] = 0;
			l[ 4] = 1;
			l[ 5] = 3;
			l[ 6] = 0;
			l[ 7] = 2;
			l[ 8] = 3;
			l[ 9] = 1;
			l[10] = 2;
			l[11] = 3;
		} else {
			len_l = 3;
			len_i = 2;
			l = (float*) malloc(len_i*len_l*sizeof(float));
			l[0] = 0;
			l[1] = 1;
			l[2] = 0;
			l[3] = 2;
			l[4] = 1;
			l[5] = 2;
		}
		/*
		v1 = [p3[1][0]-p3[0][0], p3[1][1]-p3[0][1], p3[1][2]-p3[0][2]]
		v2 = [p3[2][0]-p3[0][0], p3[2][1]-p3[0][1], p3[2][2]-p3[0][2]]
		v3 = [p3[3][0]-p3[0][0], p3[3][1]-p3[0][1], p3[3][2]-p3[0][2]]
		// check whether all are in the same plane
		// if (v1 x v2) . v3 == 0 then same plane
		if (v1[1]*v2[2]-v1[2]*v2[1])*v3[0] + (v1[2]*v2[0]-v1[0]*v2[2])*v3[1] + (v1[0]*v2[1]-v1[1]*v2[0])*v3[2] == 0:
			// quadralateral
			if facets == True:
				len_l3 = 2
				temp = math.atan2(
					v1[0]*v3[0]+v1[1]*v3[1]+v1[2]*v3[2],
					v2[0]*v3[0]+v2[1]*v3[1]+v2[2]*v3[2]);
				if temp < -math.pi/2:
					l3 = [[1,3,2],[1,0,2]]
				if temp < 0:
					l3 = [[0,3,2],[0,1,2]]
				if temp < math.pi/2:
					l3 = [[1,3,2],[1,0,2]]
				else:
					l3 = [[1,3,0],[1,2,0]]
			} else {
				# frak it, lines are cheap.
				# to do this properly, do a 2D convex hull
				len_l3 = 6
				l3 = [[0,1],[0,2],[0,3],[1,2],[1,3],[2,3]]
		} else {
			if facets == True:
				# tetrahedron
				len_l3 = 4
				l3 = [[0,1,2],[0,1,3],[0,2,3],[1,2,3]]
			} else {
				len_l = 6
				l = [[0,1],[0,2],[0,3],[1,2],[1,3],[2,3]]
			}
		*/
	} else {
		// convex hull
		
		// Determine maximum possible number of edges (E) and facets (F)
		// E <= V + (V-1) + (V-2) + ... + 0 = V*(V-1)/2
		// F = 2 + E - V
		// F <= 2 + V*(V-1)/2 - V
		// F <= 2+V*(V-3)/2
		// This does not account for breaking polygons into triangles (e.g. pyramid), ergo this is too low.
		// But we also know that it can't possibly be be higher than:
		// E <= V*V
		// F = 2 + E - V
		// F <= 2 + V*V - V
		// F <= 2+V*(V-1)
		// Then I did some maths and discovered that if you're building polyhedra from a minimal number of triangles:
		// For V >= 4
		// E = 3*V-6
		// F = 2*V-4
		// Holy frak, that simplifies things.
		if (facets == TRUE) {
			len_i = 3;
			len_l = len_i*(2*size-4);
			l = (float*) malloc(len_i*len_l*sizeof(float));
		} else {
			len_i = 2;
			len_l = len_i*(3*size-6);
			l = (float*) malloc(len_i*len_l*sizeof(float));
		}
		len_l = 0; // becomes incrementor for l
		
		// list of used variables
		visited = (int*) malloc(size*sizeof(int));
		// calculate approximate centroid and get first point
		c[0] = p[0];
		c[1] = p[1];
		c[2] = p[2];
		ii = 0; // index of cx_min, guarrantted to be part of the convex set
		// 
		for (i = 1; i < size; i++) {
			// visited
			visited[i] = UNVISITED;
			// center
			c[0] += p[i*3+0];
			c[1] += p[i*3+1];
			c[2] += p[i*3+2];
			// find min x
			if (p[i*3+0] < p[ii*3+0]) ii = i;
		}
		c[0] /= size;
		c[1] /= size;
		c[2] /= size;
		
		
		// since I know that the queue will contain at most size points, I'll just us an array
		q = (unsigned int*) malloc(size*sizeof(unsigned int)); // queue of points to be visited
		// queue first point
		q[q2++] = ii;
		
		// okay, let's allocate everything here
		v = (float*) malloc(3*size*sizeof(float));
		theta  = (float*) malloc(size*sizeof(float));
		theta2 = (float*) malloc(size*sizeof(float));
		s  = (unsigned int*) malloc(size*sizeof(unsigned int));
		s2 = (unsigned int*) malloc(size*sizeof(unsigned int));
		t  = (int*) malloc(size*sizeof(int));
		
		while (q1 < q2) {
			ii = q[q1++];

			if (debug == TRUE) {
				//print 'POINT LOOP',ii		
			}
			
			// get direction vectors
			// they need not be unit length, this might be wrong
			vii[0] = p3[ii*3+0]-c3[0];
			vii[1] = p3[ii*3+1]-c3[1];
			vii[2] = p3[ii*3+2]-c3[2];
			u0[0] = 0;
			u0[1] = 0;
			u0[2] = 0;
			if (vii[0] == 0) {
				u0[0] = 1;
			} else if (vii[1] == 0) {
				u0[1] = 1;
			} else {
				u0[0] = vii[1];
				u0[1] = -vii[0];
				tempf = sqrt(u0[0]*u0[0] + u0[1]*u0[1]);
				u0[0] /= tempf;
				u0[1] /= tempf;
			}
			// u1 = v x u0
			u1[0] = v[1]*u0[2]-v[2]*u0[1];
			u1[1] = v[2]*u0[0]-v[0]*u0[2];
			u1[2] = v[0]*u0[1]-v[1]*u0[0];
			tempf = sqrt(u1[0]*u1[0] + u1[1]*u1[1] + u1[2]*u1[2]);
			u1[0] /= tempf;
			u1[1] /= tempf;
			u1[2] /= tempf;
			
			// create a list of vectors outward from v
			// calculate their angle around v
			// create a list of test results
			for (i = 0; i < size; i++) {
				v3[i*3+0] = p3[i*3+0] - p3[ii*3+0]
				v3[i*3+1] = p3[i*3+1] - p3[ii*3+1]
				v3[i*3+2] = p3[i*3+2] - p3[ii*3+2]
				theta[i] = atan2(
					u1[0]*v3[i][0] + u1[1]*v3[i][1] + u1[2]*v3[i][2],
					u0[0]*v3[i][0] + u0[1]*v3[i][1] + u0[2]*v3[i][2]);
				t[i] = TRUE;
			}
			// create list of indices of sorted theta
			SortForIndices(&theta, &theta2, &s, &s1);
			// and mark vii as failed
			t[ii] = FALSE;
							
			if (debug == TRUE) {
				//tv = ''
				//ts = ''
				//for k in xrange(size):
				//	tv += '%4d' % k
				//	ts += '%4d' % s[k]
				//print ts + '    | ' + tv
				//print ' ' + '-'*4*size + '---+' + '-'*4*size
			}
				
			// create a set of indices:
			//   i_???? = i9 # mod 10
			// 	i_prev = i0
			// 	i_curr = i1
			// 	i_next = i2 # in and if are already terms
			i0 = -1;
			while (t[s[i0]] == FALSE) i0--;
			i9 = i0-1;
			while (t[s[i9]] == FALSE) i9--;
			i1 = 0;
			while (t[s[i1]] == FALSE) i1++;
			i2 = i1+1
			while (t[s[i2]] == FALSE) i2++;
			// go through the whole list once
			// check forward and backward at each failure
			// (until success in both directions)
			while (i1 < size) {
				// test v_1
				t[s[i1]] = ConvexityTest(&vii[0], &v[s[i0%size]*3], &v[s[i1%size]*3], &v[s[i2%size]*3])

				if (t[s[i1%size]] == FALSE) {
					// concavity detected
					// increment forward
					i1 = i2;
					i2++;
					while (t[s[i2%size]] == FALSE) i2++;
					// check forward and backward until coast is clear
					chk_f = FALSE;
					chk_b = FALSE;
					while (chk_f == FALSE || chk_b == FALSE) {
						// check backward
						chk_b = ConvexityTest(&vii[0], &v[s[i9%size]*3], &v[s[i0%size]*3], &v[s[i1%size]*3]);
						if (chk_b == FALSE) {
							t[s[i0%size]] = FALSE;
							i0 = i9;
							i9--;
							while (t[s[i9%size]] == FALSE) i9--;
						// check forward
						chk_f = ConvexityTest(&vii[0], &v[s[i0%size]], &v[s[i1%size]], &v[s[i2%size]]);
						if (chk_f == FALSE) {
							t[s[i1%size]] = FALSE;
							i1 = i2;
							i2++;
							while (t[s[i2%size]] == FALSE) i2++;
						}
					}
				} else {
					// everything is fine, continue on
					i9 = i0;
					i0 = i1;
					i1 = i2;
					i2++;
					while (t[s[i2%size]] == FALSE) i2++;
				}
				if (debug == TRUE) {
					//tv = ''
					//ts = ''
					//for k in xrange(size):
					//	tv += '%4d' % t[k]
					//	ts += '%4d' % t[s[k]]
					//print ts + '    | ' + tv
				}
		
			// in this version, you can only create lines XOR create triangles
			if (facets == TRUE) {
				// create fan of points that passes
				i1 = 0;
				while (t[s[i1%size]] == FALSE) i1++;
				i2 = i1+1
				while (t[s[i2%size]] == FALSE) i2++;
				while (i1 < size) {
					// if a triangle contains a visited point, do not create
					// and don't queue the other point in the triangle
					// else queue both points
					if (visited[s[i1%size]] != VISITED && visited[s[i2%size]] != VISITED) {
						l[len_l*3+0] = ii;
						l[len_l*3+1] = s[i1%size];
						l[len_l*3+2] = s[i2%size];
						len_l++;
						if (visited[s[i1%size]] != QUEUED) {
							visited[s[i1%size]] = QUEUED;
							q[q2++] = s[i1%size];
						}
						if (visited[s[i2%size]] != QUEUED) {
							visited[s[i2%size]] = QUEUED;
							q[q2++] = s[i2%size];
						}
					}
					i1 = i2;
					i2++;
					while (t[s[i2%size]] == FALSE) i2++;
				}
			} else {
				// create lines of points that passed
				for (i = 0; i < size; i++) {
					// if a triangle contains a visited point, do not create
					// and don't queue the other point in the triangle
					// else queue both points
					if (t[i] == TRUE && visited[i] != VISITED) {
						l[len_l*3+0] = ii;
						l[len_l*3+1] = i;
						len_l += 1;
						if (visited[i] != QUEUED) {
							visited[i] = QUEUED;
							q[q2++] = i;
						}
					}
				}
			}
			if (debug == TRUE) {
				//tv = ''
				//ts = ''
				//for k in xrange(size):
				//	tv += '%4d' % t[k]
				//	ts += '%4d' % t[s[k]]
				//print ts + '    | ' + tv
			}
			visited[ii] = VISITED; // now mark point as off limits
		}
	}
	// resize, in case of non/sub-optimal convex set
	if (facets == TRUE) {
		//del l[len_l:2*2*size-4]
	} else {
		//del l[len_l:2*3*size-6]
	}
	return l;
}
