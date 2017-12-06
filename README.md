Least-Cost Path Networks in QGIS
===================================

**Version: 0.1**

**Supported QGIS version: 2.x**

**Licence: GNU GPLv3**

Method
-------------
Least-Cost Path (LCP) analysis compute optimal ways to move from an origin towards a destination point using a cost surface map.

LCP Networks build on this concept to explore mobility dynamics on a region. They are a popular tool in a diversity of research fields such as archaeology and ecology. This plugin is designed to speed up the process by computing pairwise LCPs between each feature of a point vector layer of <b>origins</b> towards each feature of a point vector layer of <b>destinations</b>

The algorithm currently uses Dijkstra's algorithm to find the optimal paths between every pair of points through a <b>cost surface map</b> (i.e. slope or prominence raster layers). In its current form it performs the most basic analysis using anisotropic costs and the 4 direct neighbors (up, down, left and right)

Output
-------------
The plugin will generate a cost surface map for each of the original points. It will also create a new line vector layer containing all Least-Cost paths.

