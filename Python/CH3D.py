# PROGRAM: CH3D.py
# PURPOSE: Module for 3D convex hull algorithm
# AUTHOR:  Geoffrey Card
# DATE:    2014-10-13
# NOTES:   Time: O(n^2*log(n))    Algorithm requires sort(n) for each n, hence n*n*log(n).
#          Space: O(n)
#          Has trouble with non-triangluar facets.
#
# 3D convex hull (written as a low-level implementation)
# I'll be using this to write the C/C++ version.

import math
import Queue
import numpy  # for plotting
import random # for testing

# generate points
# 0. Convex set
# 1. Point
# 2. Line
# 3. Triangle
# 4. Quadralateral
# 5. Tetrahedron
# 6. Pyramid (the optimal base is 2 triangles, not 4)
# 7. Tetrahedron with a co-planar point on one face
#    (The point should either be the starting point, or skipped)
# 8. Cube (all sides should be 2 triangles)
# ?. Random points
def GenerateStandardPoints(sel = 9):
	random.seed()
	if sel == 0:
		# Convex set
		len_p3 = 7
		p3 = [[0 for i in xrange(3)] for i in xrange(len_p3)]
		for i in xrange(len_p3):
			temp1 = 2*math.pi*random.random()
			temp2 = math.asin(2*random.random()-1)
			p3[i][0] = math.sin(temp2)*math.cos(temp1)
			p3[i][1] = math.sin(temp2)*math.sin(temp1)
			p3[i][2] = math.cos(temp2)		
	elif sel == 1:
		# Point
		len_p3 = 1
		p3 = [[random.random() for i in xrange(3)]]
	elif sel == 2:
		# Line
		len_p3 = 2
		p3 = [[random.random() for i in xrange(3)] for i in xrange(len_p3)]
	elif sel == 3:
		# Triangle
		len_p3 = 3
		p3 = [[random.random() for i in xrange(3)] for i in xrange(len_p3)]
	elif sel == 4:
		# Quadralateral
		len_p3 = 4
		p3 = [[random.random() for i in xrange(3)] for i in xrange(len_p3)]
		p3[3] = [p3[0][0]+p3[1][0]-p3[2][0], p3[0][1]+p3[1][1]-p3[2][1], p3[0][2]+p3[1][2]-p3[2][2]]
	elif sel == 5:
		# Tetrahedron
		len_p3 = 4
		p3 = [[random.random() for i in xrange(3)] for i in xrange(len_p3)]
	elif sel == 6:
		# Pyramid, CONVEX HULL FAILED!
		len_p3 = 5
		p3 = [[random.random() for i in xrange(3)] for i in xrange(len_p3)]
		p3[4] = [p3[0][0]+p3[1][0]-p3[2][0], p3[0][1]+p3[1][1]-p3[2][1], p3[0][2]+p3[1][2]-p3[2][2]]
	elif sel == 7:
		# Tetrahedron with co-planar point, CONVEX HULL FAIL
		len_p3 = 5
		p3 = [[random.random() for i in xrange(3)] for i in xrange(len_p3)]
		# p5 = (p3[0]+p3[1])/2
		# p3[4] = (p5+p3[2])/2
		p5 = [(p3[0][0]+p3[1][0])/2, (p3[0][1]+p3[1][1])/2, (p3[0][2]+p3[1][2])/2]
		p3[4] = [(p5[0]+p3[2][0])/2, (p5[1]+p3[2][1])/2, (p5[2]+p3[2][2])/2]
	elif sel == 8:
		# Cube, NOT EVEN CLOSE, TOTALLY FAILED ON THIS ONE
		len_p3 = 8
		p3 = [[0,0,0],[0,0,1],[0,1,0],[0,1,1],[1,0,0],[1,0,1],[1,1,0],[1,1,1]]
	else:
		# Random points
		len_p3 = 12
		p3 = [[random.random() for i in xrange(3)] for i in xrange(len_p3)]
	return p3


def GeneratePoints(count=7,convexset=False):
	random.seed()
	if convexset == True:
		# Convex set
		len_p3 = count
		p3 = [[0 for i in xrange(3)] for i in xrange(len_p3)]
		for i in xrange(len_p3):
			temp1 = 2*math.pi*random.random()
			temp2 = math.asin(2*random.random()-1)
			p3[i][0] = math.sin(temp2)*math.cos(temp1)
			p3[i][1] = math.sin(temp2)*math.sin(temp1)
			p3[i][2] = math.cos(temp2)		
	else:
		# Random points
		len_p3 = count
		p3 = [[random.random() for i in xrange(3)] for i in xrange(len_p3)]
	return p3

################################################################
################################################################
################################################################

# removes repeats from list
def RemoveRepetitions(p3):
	size = len(p3)
	s1 = range(size)
	s2 = [0]*size
	# gc merge(ish) sort
	a = 1
	while a < size*2:
		i = 0
		while i < size:
			lim1 = min(i+1*a,size)
			lim2 = min(i+2*a,size)
			i1 = i
			i2 = lim1
			while i < lim2:
				if i1 == lim1:
					s2[i] = s1[i2]
					i2 += 1
				elif i2 == lim2:
					s2[i] = s1[i1]
					i1 += 1
				elif p3[s1[i1]][0] < p3[s1[i2]][0]: 
					s2[i] = s1[i1]
					i1 += 1
				elif p3[s1[i1]][0] > p3[s1[i2]][0]: 
					s2[i] = s1[i2]
					i2 += 1
				elif p3[s1[i1]][1] < p3[s1[i2]][1]: 
					s2[i] = s1[i1]
					i1 += 1
				elif p3[s1[i1]][1] > p3[s1[i2]][1]: 
					s2[i] = s1[i2]
					i2 += 1
				elif p3[s1[i1]][2] < p3[s1[i2]][2]: 
					s2[i] = s1[i1]
					i1 += 1
				elif p3[s1[i1]][2] > p3[s1[i2]][2]: 
					s2[i] = s1[i2]
					i2 += 1
				else:
					s2[i] = s1[i1]
					i1 += 1
				i += 1
		a *= 2
		temp = s1
		s1 = s2
		s2 = temp
	# repurpose s2 to mark repeats
	s2[s1[size-1]] = 0
	for i in range(size-1):
		if p3[s1[i]][0] == p3[s1[i+1]][0] and p3[s1[i]][1] == p3[s1[i+1]][1] and p3[s1[i]][2] == p3[s1[i+1]][2]:
			s2[s1[i]] = 1
		else:
			s2[s1[i]] = 0
	i = size-1
	while i >= 0:
		if s2[i] == 1:
			del p3[i]
		i -= 1

# convexity test of v_1 around v
def ConvexityTest(v, v_0, v_1, v_2):
		# test = (v_0 x v_2) . v_1
		test = (v_0[1]*v_2[2] - v_0[2]*v_2[1]) * v_1[0] + (v_0[2]*v_2[0] - v_0[0]*v_2[2]) * v_1[1] + (v_0[0]*v_2[1] - v_0[1]*v_2[0]) * v_1[2]
		
		# if it's in the same plane as v_0 and v_2
		# check whether it extends past them
		if test == 0:
			# v4 = (v_2 - v_0) x v
			# if v4 . v_1 > v_4 . v_0 (or v_4 . v_2)
			v4 = [v_2[0]-v_0[0], v_2[1]-v_0[1], v_2[2]-v_0[2]]
			v4 = [v4[1]*v[2]-v4[2]*v[1], v4[2]*v[0]-v4[0]*v[2], v4[0]*v[1]-v4[1]*v[0]]
			test = v4[0]*v_1[0] + v4[1]*v_1[1] + v4[2]*v_1[2] > v4[0]*v_0[0] + v4[1]*v_0[1] + v4[2]*v_0[2]
		
		# test MUST be greater than 0, else this is not part of the minimal convex set
		return test > 0

# sorts and returns (sorted) indices (ascending)
def SortForIndices(l1):
	size = len(l1)
	l2 = [0.0]*size
	s1 = range(size)
	s2 = [0]*size
	a = 1
	while a < size*2:
		i = 0
		while i < size:
			lim1 = min(i+1*a,size)
			lim2 = min(i+2*a,size)
			i1 = i
			i2 = lim1
			while i < lim2:
				if i1 == lim1:
					l2[i] = l1[i2]
					s2[i] = s1[i2]
					i2 += 1
				elif i2 == lim2:
					l2[i] = l1[i1]
					s2[i] = s1[i1]
					i1 += 1
				elif l1[i1] <= l1[i2]:
					l2[i] = l1[i1]
					s2[i] = s1[i1]
					i1 += 1
				else:
					l2[i] = l1[i2]
					s2[i] = s1[i2]
					i2 += 1
				i += 1
		a *= 2
		temp = l1
		l1 = l2
		l2 = temp
		temp = s1
		s1 = s2
		s2 = temp
	return list(s1)

# 3D convex hull
def ConvexHull3D(p3, facets=True, DEBUG=False):
	# get rid of pesky repetitions get rid of pesky repetitions
	RemoveRepetitions(p3)
	
	# confirm non-linear
	# confirm non-planar
	# vi = pi - p0
	# for i:
	# 	if v0 x vi != 0:
	#		print 'non-linear'
	#	if (v0 x v1) . vi != 0:
	#		print 'non-planar'
	
	
	# get size, I don't want to call a function every time I need this
	len_p3 = len(p3)
	
	# for len_p3 < 5, it doesn't seem worth it to check for repeats
	if len_p3 == 1:
		# point
		len_l3 = 1
		l3 = [[0]]
	elif len_p3 == 2:
		# line
		len_l3 = 1
		l3 = [[0,1]]
	elif len_p3 == 3:
		# triangle
		if facets == True:
			len_l3 = 1
			l3 = [[0,1,2]]
		else:
			len_l3 = 3
			l3 = [[0,1],[0,2],[1,2]]
	elif len_p3 == 4:
		# quadralateral or tetrahedron
		v1 = [p3[1][0]-p3[0][0], p3[1][1]-p3[0][1], p3[1][2]-p3[0][2]]
		v2 = [p3[2][0]-p3[0][0], p3[2][1]-p3[0][1], p3[2][2]-p3[0][2]]
		v3 = [p3[3][0]-p3[0][0], p3[3][1]-p3[0][1], p3[3][2]-p3[0][2]]
		# check whether all are in the same plane
		# if (v1 x v2) . v3 == 0 then same plane
		if (v1[1]*v2[2]-v1[2]*v2[1])*v3[0] + (v1[2]*v2[0]-v1[0]*v2[2])*v3[1] + (v1[0]*v2[1]-v1[1]*v2[0])*v3[2] == 0:
			# quadralateral
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
			else:
				# frak it, lines are cheap.
				# to do this properly, do a 2D convex hull
				len_l3 = 6
				l3 = [[0,1],[0,2],[0,3],[1,2],[1,3],[2,3]]
		else:
			if facets == True:
				# tetrahedron
				len_l3 = 4
				l3 = [[0,1,2],[0,1,3],[0,2,3],[1,2,3]]
			else:
				len_l3 = 6
				l3 = [[0,1],[0,2],[0,3],[1,2],[1,3],[2,3]]
	else:
		# convex hull

		# Determine maximum possible number of edges (E) and facets (F)
		# E <= V + (V-1) + (V-2) + ... + 0 = V*(V-1)/2
		# F = 2 + E - V
		# F <= 2 + V*(V-1)/2 - V
		# F <= 2+V*(V-3)/2
		# This does not account for breaking polygons into triangles (e.g. pyramid), ergo this is too low.
		# But we also know that it can't possibly be be higher than:
		# E <= V*V
		# F = 2 + E - V
		# F <= 2 + V*V - V
		# F <= 2+V*(V-1)
		# Then I did some maths and discovered that if you're building polyhedra from a minimal number of triangles:
		# For V >= 4
		# E = 3*V-6
		# F = 2*V-4
		# Holy frak, that simplifies things.
		if facets == True:
			l3 = [[0 for i in xrange(3)] for i in xrange(2*2*len_p3-4)] ###################################################### NOTE, I DOUBLED THIS
			len_l3 = 0
		else:
			l3 = [[0 for i in xrange(2)] for i in xrange(2*3*len_p3-6)]
			len_l3 = 0
	
		# calculate approximate centroid and get first point
		cx_avg = 0.0
		cy_avg = 0.0
		cz_avg = 0.0
		cz_max = p3[0][2]
		ii = 0 # index of cz_max, guarrantted to be part of the convex set
		for i in xrange(len_p3):
			cx_avg += p3[i][0]
			cy_avg += p3[i][1]
			cz_avg += p3[i][2]
			if p3[i][2] > cz_max:
				cz_max = p3[i][2]
				ii = i
		c3 = [cx_avg/len_p3, cy_avg/len_p3, cz_avg/len_p3]
		
		q = Queue.Queue(len_p3); # queue of points to be visited
		queued  = [False for i in xrange(len_p3)] # list of already-queued variables
		visited = [False for i in xrange(len_p3)] # list of used variables
		# queue first point
		q.put(ii)
	
		while q.empty() == False:
			ii = q.get()

			if DEBUG:
				print 'POINT LOOP',ii		
		
			# get direction vectors
			# they need not be unit length, this might be wrong
			v = [p3[ii][0]-c3[0], p3[ii][1]-c3[1], p3[ii][2]-c3[2]]
			u0 = [0, 0, 0]
			u1 = [0, 0, 0]	
			if v[0] == 0:
				u0 = [1, 0, 0]
			elif v[1] == 0:
				u0 = [0, 1, 0]
			else:
				u0 = [v[1], -v[0], 0]
				temp = math.sqrt(u0[0]*u0[0] + u0[1]*u0[1] + u0[2]*u0[2])
				u0[0] /= temp
				u0[1] /= temp
				u0[2] /= temp
			# u1 = v x u0
			u1 = [v[1]*u0[2]-v[2]*u0[1], v[2]*u0[0]-v[0]*u0[2], v[0]*u0[1]-v[1]*u0[0]]
			temp = math.sqrt(u1[0]*u1[0] + u1[1]*u1[1] + u1[2]*u1[2])
			u1[0] /= temp
			u1[1] /= temp
			u1[2] /= temp
			
			# create a list of vectors outward from v
			v3 = [[0.0 for i in xrange(3)] for j in xrange(len_p3)]
			# calculate their angle around v
			theta = [0.0 for j in xrange(len_p3)]
			for i in xrange(len_p3):
				v3[i][0] = p3[i][0] - p3[ii][0]
				v3[i][1] = p3[i][1] - p3[ii][1]
				v3[i][2] = p3[i][2] - p3[ii][2]
				theta[i] = math.atan2(
					u1[0]*v3[i][0] + u1[1]*v3[i][1] + u1[2]*v3[i][2],
					u0[0]*v3[i][0] + u0[1]*v3[i][1] + u0[2]*v3[i][2]);

			# create list of indices of sorted theta
			s = SortForIndices(list(theta))
			# create a list of test results
			t = [True]*len_p3
			# and mark v as failed
			t[ii] = False
							
			if DEBUG:
				tv = ''
				ts = ''
				for k in xrange(len_p3):
					tv += '%4d' % k
					ts += '%4d' % s[k]
				print ts + '    | ' + tv
				print ' ' + '-'*4*len_p3 + '---+' + '-'*4*len_p3
				
			# create a set of indices:
			#   i_???? = i9 # mod 10
			# 	i_prev = i0
			# 	i_curr = i1
			# 	i_next = i2 # in and if are already terms
			i0 = -1
			while t[s[i0]] == False:
				i0 -= 1
			i9 = i0-1
			while t[s[i9]] == False:
				i9 -= 1
			i1 = 0
			while t[s[i1]] == False:
				i1 += 1
			i2 = i1+1
			while t[s[i2]] == False:
				i2 += 1
			# go through the whole list once
			# check forward and backward at each failure
			# (until success in both directions)
			while i1 < len_p3:
				# test v_1
				t[s[i1]] = ConvexityTest(v, v3[s[i0%len_p3]], v3[s[i1%len_p3]], v3[s[i2%len_p3]])

				if t[s[i1%len_p3]] == False:
					# concavity detected
					# increment forward
					i1 = i2;
					i2 += 1
					while t[s[i2%len_p3]] == False:
						i2 += 1
					# check forward and backward until coast is clear
					chk_f = False
					chk_b = False
					while chk_f == False or chk_b == False:
						# check backward
						chk_b = ConvexityTest(v, v3[s[i9%len_p3]], v3[s[i0%len_p3]], v3[s[i1%len_p3]])
						if chk_b == False:
							t[s[i0%len_p3]] = False
							i0 = i9
							i9 -= 1
							while t[s[i9%len_p3]] == False:
								i9 -= 1
						# check forward
						chk_f = ConvexityTest(v, v3[s[i0%len_p3]], v3[s[i1%len_p3]], v3[s[i2%len_p3]])
						if chk_f == False:
							t[s[i1%len_p3]] = False
							i1 = i2
							i2 += 1
							while t[s[i2%len_p3]] == False:
								i2 += 1
				else:
					# everything is fine, continue on
					i9 = i0
					i0 = i1
					i1 = i2
					i2 += 1
					while t[s[i2%len_p3]] == False:
						i2 += 1

				if DEBUG:
					tv = ''
					ts = ''
					for k in xrange(len_p3):
						tv += '%4d' % t[k]
						ts += '%4d' % t[s[k]]
					print ts + '    | ' + tv
		
			# in this version, you can only create lines XOR create triangles
			if facets == True:
				# create fan of points that passes
				i1 = 0
				while t[s[i1%len_p3]] == False:
					i1 += 1
				i2 = i1+1
				while t[s[i2%len_p3]] == False:
					i2 += 1
				while i1 < len_p3:
					# if a triangle contains a visited point, do not create
					# and don't queue the other point in the triangle
					# else queue both points
					if visited[s[i1%len_p3]] == False and visited[s[i2%len_p3]] == False:
						l3[len_l3] = [ii, s[i1%len_p3], s[i2%len_p3]]
						len_l3 += 1
						if queued[s[i1%len_p3]] == False:
							queued[s[i1%len_p3]] = True
							q.put(s[i1%len_p3])
						if queued[s[i2%len_p3]] == False:
							queued[s[i2%len_p3]] = True
							q.put(s[i2%len_p3])
					i1 = i2
					i2 += 1
					while t[s[i2%len_p3]] == False:
						i2 += 1
			else:
				# create lines of points that passed
				for i in xrange(len_p3):
					# if a triangle contains a visited point, do not create
					# and don't queue the other point in the triangle
					# else queue both points
					if t[i] == True and visited[i] == False:
						l3[len_l3] = [ii, i]
						len_l3 += 1
						if queued[i] == False:
							queued[i] = True
							q.put(i)
			if DEBUG:
				tv = ''
				ts = ''
				for k in xrange(len_p3):
					tv += '%4d' % t[k]
					ts += '%4d' % t[s[k]]
				print ts + '    | ' + tv
		
			visited[ii] = True # now mark point as off limits
	
	# resize, in case of non/sub-optimal convex set
	if facets == True:
		del l3[len_l3:2*2*len_p3-4]
	else:
		del l3[len_l3:2*3*len_p3-6]
	
	return l3

################################################################
################################################################
################################################################

import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# plot
def PlotHull(p3, l3, labelled=False, exploded=False):
	len_p3 = len(p3)
	len_l3 = len(l3)
	# convert to array to make plotting easier
	p3 = numpy.array(p3)
	l3 = numpy.array(l3)
	# get center
	c3 = [numpy.mean(p3[:,0]), numpy.mean(p3[:,1]), numpy.mean(p3[:,2])]
	
	# setup plot
	fig = plt.figure()
	ax = fig.gca(projection='3d')
	plt.hold(True)
	
	# point labels
	if labelled == True:
		ax.text(c3[0], c3[1], c3[2], 'c', None)
		for i in xrange(len_p3):
			label = '%d' % (i)
			ax.text(p3[i][0], p3[i][1], p3[i][2], label, None)
			ax.scatter(p3[i][0], p3[i][1], p3[i][2], color=[0.5,0.5,0.5], marker='o')
	
	# draw
	if len(l3[0]) == 1:
		# point
		for i in xrange(len_l3):
			ax.scatter(p3[l3[i][0]][0], p3[l3[i][0]][1], p3[l3[i][0]][2], color=[0.5,0.5,0.5], marker='o')
	elif len(l3[0]) == 2:
		# lines
		if exploded == True:
			i1 = 0;
			i2 = 0;
			while i1 < len_l3:
				i2 = i1+1
				while i2 < len_l3 and l3[i1][0] == l3[i2][0]:
					i2 +=1
				v3 = p3 + 0.1*(p3[l3[i1][0]]-c3)
				for i in range(i1,i2):
					ax.plot([v3[l3[i][0]][0], v3[l3[i][1]][0]], [v3[l3[i][0]][1], v3[l3[i][1]][1]], [v3[l3[i][0]][2], v3[l3[i][1]][2]], label='parametric curve', color=[l3[i1][0]/len_l3, l3[i1][0]/len_l3, l3[i1][0]/len_l3])
				i1 = i2
		else:
			for i in xrange(len_l3):
				ax.plot([p3[l3[i][0]][0], p3[l3[i][1]][0]], [p3[l3[i][0]][1], p3[l3[i][1]][1]], [p3[l3[i][0]][2], p3[l3[i][1]][2]], label='parametric curve', color=[0.5,0.5,0.5])
	elif len(l3[0]) == 3:
		# triangles
		if exploded == True:
			i1 = 0;
			i2 = 0;
			while i1 < len_l3:
				i2 = i1+1
				while i2 < len_l3 and l3[i1][0] == l3[i2][0]:
					i2 +=1
				v3 = p3 + 0.1*(p3[l3[i1][0]]-c3)
				ax.plot_trisurf(v3[:,0], v3[:,1], v3[:,2], triangles=l3[i1:i2], color=[1, l3[i1][0]/len_l3, 0])
				i1 = i2
		else:
			ax.plot_trisurf(p3[:,0], p3[:,1], p3[:,2], triangles=l3, color=[0.5,0.0,0.0])

	# finish and show
	plt.hold(False)
	plt.show()
