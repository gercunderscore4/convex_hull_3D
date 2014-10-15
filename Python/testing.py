try:
	CH3D
except NameError:
	import CH3D
else:
	reload(CH3D)

p3 = CH3D.GenerateStandardPoints(0)                      # vertices
f3 = CH3D.ConvexHull3D(p3, facets=True, DEBUG=False)     # facets
e3 = CH3D.ConvexHull3D(p3, facets=False, DEBUG=False)    # edges
CH3D.PlotHull(p3, f3, labelled=True)
CH3D.PlotHull(p3, e3, labelled=True)
