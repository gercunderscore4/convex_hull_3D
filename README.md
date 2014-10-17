convex_hull_3D
==============

3D convex hull algorithm

Time complexity:  O(n^2*log(n))

Space complexity: O(n)


Most algorithms I've come across have been very complicated, this one seems pretty simple, I thought it up myself.


Implemented in MATLAB/Octave and Python (2.7), C, and soon C++.

TODO:
- Test new method for detecting polygonal facets.
- Test C implementation.
- Write C++ implementation.

KNOWN ISSUES:
- The's a problem in dealing with non-triangular facets. I have a soluton, I just need to implement it.
- Floating point errors are causing problems. I'll fix it with a minimum value (i.e. abs(val) < epsilon -> val = 0).
- Issues detecting points in a plane (due to the above). I hope to fix this soon.

There's no license because I want to let people freely copy, modify, and improve it without any license. If you do improve it, please let me know.
