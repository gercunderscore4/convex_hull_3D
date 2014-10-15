import time
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
try:
	CH3D
except NameError:
	import CH3D
else:
	reload(CH3D)

sizes = range(1,10)      # How many points in the convex hull (list).
ntest = 1000             # Number of times to test each size.
times = [0.0]*len(sizes) # Records average time for each size.
for j in xrange(len(sizes)):
	timetrial = 0.0
	for i in xrange(ntest):
		p3 = CH3D.GeneratePoints(sizes[j], False) # False -> non-convex set

		t1 = time.clock()
		l3 = CH3D.ConvexHull3D(p3)
		t2 = time.clock()
		
		timetrial += t2-t1 # sum up individual results
	times[j] = timetrial/ntest # record average result
	print '%4d %10.7f' % (sizes[j], times[j]) # print readable results

# plot
fig = plt.figure()
ax = fig.gca()
ax.plot(sizes,times)
plt.xlabel('points (#)')
plt.ylabel('time (s)')
plt.show()
