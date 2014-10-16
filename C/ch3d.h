/* PROGRAM: ch3d.h
 * PURPOSE: Module for 3D convex hull algorithm
 * AUTHOR:  Geoffrey Card
 * DATE:    2014-10-16
 * NOTES:   Time: O(n^2*log(n))    Algorithm requires sort(n) for each n, hence n*n*log(n).
 *          Space: O(n)
 *          Untested.
 * TO DO:   
 */

#ifndef CH3D_H_
#define CH3D_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define PI 3.14159

#define EPSILON 1e-6

#define TRUE  1
#define FALSE 0

#define UNVISITED 0
#define QUEUED    1
#define VISITED   2

#define VERTICES 0
#define EDGES    1
#define FACETS   2

#define CONCAVE 0
#define PLANAR  1
#define CONVEX  2

/* void CrossProduct(float* a, float* b, float* c)
 * 
 * DESCRIPTION:
 *     Cross product: c = a x b
 * 
 * PARAMETERS:
 *     float*           a        3D vector
 *     float*           b        3D vector
 *     float*           c        3D vector, receives result
 *
 * RETURNS:
 *     void
 * 
 * NOTE:
 *     Consider making this inline, due to its small size.
 */
void CrossProduct(float* a, float* b, float* c);

/* float DotProduct(float* a, float* b)
 * 
 * DESCRIPTION:
 *     Dot product: a . b
 * 
 * PARAMETERS:
 *     float*           a        3D vector
 *     float*           b        3D vector
 *
 * RETURNS:
 *     float                     result
 * 
 * NOTE:
 *     Consider making this inline, due to its small size.
 */
float DotProduct(float* a, float* b);

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
void UnitVector(float* u);

/* float TripleProduct(float* a, float* b, float* c)
 * 
 * DESCRIPTION:
 *     Triple product: a . (b x c) == b . (c x a) == c . (a x b)
 * 
 * PARAMETERS:
 *     float*           a        3D vector
 *     float*           b        3D vector
 *     float*           c        3D vector
 *
 * RETURNS:
 *     float                     result
 * 
 * NOTE:
 *     Consider making this inline, due to its small size.
 */
float TripleProduct(float* a, float* b, float* c);

/* void GeneratePoints(float* p, unsigned int size, int convexset)
 * 
 * DESCRIPTION:
 *     Generates a set of size 3D points. If convexset is TRUE, the points are a convex set.
 * 
 * IMPORTANT:
 *     perform srand before this
 *     p MUST BE PRE-ALLOCATED with size 3*size
 *     p IS NOT FREED IN THIS FUNCTION
 * 
 * PARAMETERS:
 *     float*           p        array for storing points
 *     unsigned int     size     number of 3D points (actual size of p is 3*size)
 *
 * RETURNS:
 *     void
 */
void GeneratePoints(float* p, unsigned int size, int convexset);

/* void PrintPoints(float* p, unsigned int lenp)
 * 
 * DESCRIPTION:
 *     Prints an indexed list of all points.
 * 
 * IMPORTANT:
 *     p MUST BE ALLOCATED with size 3*lenp
 *     p IS NOT FREED IN THIS FUNCTION
 * 
 * PARAMETERS:
 *     float*           p        array of 3D points
 *     unsigned int     lenp     number of 3D points (actual size of p is 3*lenp)
 *
 * RETURNS:
 *     void
 */
void PrintPoints(float* p, unsigned int lenp);

/* void PrintResults(float* p, unsigned int lenp, unsigned int* l, unsigned int lenl, unsigned int leni)
 * 
 * DESCRIPTION:
 *     Prints an indexed list of all points.
 *     Prints a list of all vertex/edge/facet indices and points.
 * 
 * IMPORTANT:
 *     p MUST BE ALLOCATED with size 3*lenp
 *     l MUST BE ALLOCATED with size lenl*leni
 *     p, l ARE NOT FREED IN THIS FUNCTION
 * 
 * PARAMETERS:
 *     float*           p        array of 3D points
 *     unsigned int     lenp     number of 3D points (actual size of p is 3*lenp)
 *     float*           l        array of point indices
 *     unsigned int     lenl     number of listed objects (actual size of p is 3*lenl*leni)
 *     unsigned int     leni     number of points per list item (vertex/edge/facet)
 *                               this determines whether the list if of vertices, edges or facets
 *
 * RETURNS:
 *     void
 */
void PrintResults(float* p, unsigned int lenp, unsigned int* l, unsigned int lenl, unsigned int leni);

/* unsigned int RemoveRepetitions (float* p, unsigned int size)
 * 
 * DESCRIPTION:
 *     Removes repeated points (which would cause problems for the 3D convex hull algorithm).
 *     Returns new number of points.
 * 
 * IMPORTANT:
 *     p MUST BE PRE-ALLOCATED with size 3*size
 *     p IS NOT FREED IN THIS FUNCTION
 * 
 * PARAMETERS:
 *     float*           p        array for storing points
 *     unsigned int     size     number of 3D points (actual size of p is size*3)
 *
 * RETURNS:
 *     unsigned int              new value for size
 */
unsigned int RemoveRepetitions (float* p, unsigned int size);

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
void SortForIndices (float* l1, unsigned int* s1, unsigned int* s2, unsigned int size);

/* unsigned int ConvexityTest (float* v, float* v0, float* v1, float* v2)
 * 
 * DESCRIPTION:
 *     Checks whether v1 is above, in, or below the plane of v0 and v2: v1 . (v0 x v2) > 0
 *     If it's in the plane, it checks whether v1 is farther away from v than the line between v0 and v2: 
 * 
 * IMPORTANT:
 * 
 * PARAMETERS:
 *     float*           v        3D vector, vii - c (if you don't understand that, read the convex hull algorithm)
 *     float*           v0       3D vector
 *     float*           v1       3D vector
 *     float*           v2       3D vector
 *
 * RETURNS:
 *     unsigned int              return CONCAVE, PLANAR, or CONVEX
 * 
 * NOTES:
 *     This function is designed for use inside of ConvexHull3D, nowhere else.
 */
unsigned int ConvexityTest (float* v, float* v0, float* v1, float* v2);

/* void ConvexHull3D(float* p, unsigned int* lenp_p, float* l, unsigned int* lenl_p, unsigned int* leni_p, int list_type, int debug)
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
 *     p MUST BE PRE-ALLOCATED with size 3*lenp
 *     l MUST BE UNALLOCATED, it is allocated in this function
 *     p, l are NOT FREED IN THIS FUNCTION
 * 
 * PARAMETERS:
 *     float*           p        array of 3D points, size 3*lenp
 *     unsigned int*    lenp_p   pointer to number of 3D points (actual size of p is 3*lenp)
 *                               it's a pointer so that it can change the original value
 *     unsigned int*    l        array of 3D points, size leni*lenl
 *     unsigned int*    lenl_p   pointer to the length of the list (actual size of l is leni*len_p)
 *                               it's a pointer so that it can change the original value
 *     unsigned int*    leni_p   pointer to the length of the list type
 *                               1 if it's points, 2 if it's edges, 3 if it's facets
 *                               it's a pointer so that it can change the original value
 *
 * RETURNS:
 *     void
 * 
 * NOTES:
 *     While not yet implemented, I'm going to write in a section that calculates the normals of each test, and marks points that are in non-triangular facets.
 *     If a new facet is found to have marked points, it's normal is compared to a list of non-triangular facet normals.
 *     The list of points adds n memory, the list of facet normals adds 3*(2*V-4). Still linear.
 */
void ConvexHull3D(float* p, unsigned int* lenp_p, float* l, unsigned int* lenl_p, unsigned int* leni_p, int list_type, int debug);

#endif /* CH3D_H_ */