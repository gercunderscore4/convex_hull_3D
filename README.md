convex_hull_3D
==============

***IMPORTANT-ISH***
There appears to be a logic error in the MATLAB code. It doesn't actually do anything wrong, it's just not written in a way that makes logical sense.

Basic 3D convex hull algorithm.
Complexity O(n^2 log(n)) because there's a (merge/quick) sort (n log(n)) required per point (n).
Most algorithms I've come across have been very complicated, this one seems pretty simple, I thought it up myself.

How it works (simplified, assumes convex set, wosrt case approaches O(n^3)):

(01)  Take a set of 3D convex points with no repetitions, p3.

(02)  Find the center, c3 (this is my most ill-defined and troublesome step).

(03)  Subtract c3 from p3 to form vectors v3.

(04)  For each vector, v:

(05)    Calculate a 3D basis, u1, u2, u3, using v as the third axis.

(06)    Assign each vector an angle, theta, around u3, using arctan(u2,u1).

(07)    Sort theta, and create a list of relative indices, ind (unsorted -> sorted), and invind (sorted -> unsorted).

(08)    Create a sorted list for test results, t, with all values set to true.

(09)    Set t of v to false.

(10)    While a test has falsified in the past round, for each true vector that tests true, vv:

(11)      Take the vectors ahead, va, and behind, vb.

(12)      t of p gets dot( vv-v, cross(vb-v, va-v) ) >= 0.

(13)    For every vector for which t is true, vv:

(14)      Take the vectors ahead, va, and behind, vb.

(15)      If the unsorted index of vv is smaller than the unsorted indices of va and vb:

(16)        Add to list l3 the triangle of vv with the first true points ahead and behind it, in unsorted form.

(17)  Return the list of all triangles that form the convex hull, l3.

There's no license because I want to let people freely copy, modify, and improve it without any license stuff (and so that they might not put a license on it, so that I can copy it back).
