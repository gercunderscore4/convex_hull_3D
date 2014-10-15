convex_hull_3D
==============

3D convex hull algorithm
Time complexity:  O(n^2*log(n))    Requires sort(n) for every n, in every case.
Space complexity: O(n)
Most algorithms I've come across have been very complicated, this one seems pretty simple, I thought it up myself.

Implemented in MATLAB/Octave and Python (2.7), and soon C/C++.

TODO:
- Write C/C++ implementation.
- Fix floating point errors.
- Correct logic error in MATLAB code (that doesn't actually cause problems, but is still an error).

KNOWN ISSUES:
- The's a problem in dealing with polygonal facets. It'll need something clever to fix it. For now, it creates a valid, but very sub-optimal solution.
- Floating point errors are causing problems. I'll fix it with a minimum value (i.e. abs(val) < epsilon -> val = 0).
- Issues detecting points in a plane (due to the above). I hope to fix this soon.

There's no license because I want to let people freely copy, modify, and improve it without any license. If you do improve it, please let me know.
