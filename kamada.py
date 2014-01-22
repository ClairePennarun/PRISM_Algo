# Powered by Python 2.7

from __future__ import division
from tulip import *
from tulipgui import *
from tulipogl import *

import math
import time

class Prism:

	def applyToGraph(self,graph,epsilon):
		self.graph = graph
		self.createDelaunay()
		nbOverlaps = self.computeOverlapFactors()
		changed = 1
		while nbOverlaps > 0 and changed == 1:
			changed = self.kamada(epsilon)
			self.createDelaunay()
			nbOverlaps = self.computeOverlapFactors()
		changed = 1
		graphCopy = tlp.newGraph()
		tlp.copyToGraph(graphCopy,graph)
		adjustNodesLabels(graphCopy)	
		graphCopy.setName("first phase copy")
		
		self.createDelaunay()
		self.computeBounding()
		self.scanline()
		nbOverlaps = self.computeOverlapFactors()
		while nbOverlaps > 0 and changed == 1:
			changed = self.kamada(epsilon)
			self.createDelaunay()
			self.computeBounding()
			self.scanline()
			nbOverlaps = self.computeOverlapFactors()
		graphCopyFinal = tlp.newGraph()
		tlp.copyToGraph(graphCopyFinal,graph)
		adjustNodesLabels(graphCopyFinal)	
		graphCopyFinal.setName("second phase copy")

	def createDelaunay(self):
		coords = [] # coordonnees des nodes de graph
		nodesOrder = [] # mapping entre les nodes et les indices dans coords
		
		self.viewLayout = self.graph.getLayoutProperty("viewLayout")
		for n in self.graph.getNodes():
			coords.append(self.viewLayout[n])
			nodesOrder.append(n)
		
		(delaunayTuples, triangles) = tlp.delaunayTriangulation(coords)
		self.delaunayEdges = [] # liste d'aretes de DT
		for (s,t) in delaunayTuples:
			self.delaunayEdges.append((nodesOrder[s],nodesOrder[t]))
		
		self.delaunayNeighbors = self.graph.getIntegerVectorProperty("neighbors")
		for n in self.graph.getNodes():
			self.delaunayNeighbors.erase(n)
		for e in self.delaunayEdges:
			self.delaunayNeighbors.pushBackNodeEltValue(e[0],e[1].id)
			self.delaunayNeighbors.pushBackNodeEltValue(e[1],e[0].id)

	def computeOverlapFactors(self):
		self.edgesOverlap = {} # dictionnaire extremites arete -> overlap
		nbOverlaps = 0
		for e in self.delaunayEdges:
			overlap = self.computeOverlap(e[0],e[1])
			self.edgesOverlap[(e[0].id,e[1].id)] = overlap
			self.edgesOverlap[(e[1].id,e[0].id)] = overlap
			if overlap > 1:
				nbOverlaps +=1
		return nbOverlaps
	
	def computeOverlap(self,i,j):
		viewSize = self.graph.getSizeProperty("viewSize")
		sizeI = viewSize[i]
		sizeJ = viewSize[j]
		posI = self.viewLayout[i]
		posJ = self.viewLayout[j]
		firstFactor = 1
		secondFactor = 1
		factor = 1
		a = abs(posI[0] - posJ[0])
		if (a <> 0):
			firstFactor = (sizeI[0]/2 + sizeJ[0]/2)/a
			factor = firstFactor
		b = abs(posI[1] - posJ[1])
		if (b <> 0):
			secondFactor = (sizeI[1]/2 + sizeJ[1]/2)/b
			if (secondFactor < firstFactor):
				factor = secondFactor
		if factor > 1:
			return factor
		return 1

	def findMax(self,deltaM):
		first = self.graph.getOneNode()
		max = deltaM[first]
		maxIndex = first 
		for n in self.graph.getNodes():
			tmp = deltaM[n]
			if tmp > max:
				max = tmp
				maxIndex = n
		return maxIndex
	
	def dist(self,i,j):
		return math.sqrt((self.viewLayout[i][0] - self.viewLayout[j][0])**2 + (self.viewLayout[i][1] - self.viewLayout[j][1])**2)
	
	# resolution de Ax + By = E, Cx + Dy = F
	def solveEquations(self,A,B,C,D,E,F):
		x = (E*D - B*F)/(A*D - B*C)
		y = (1/D)*(F - (C*E*D - B*C*F)/(A*D - B*C))
		return [x,y]
	
	def kamada(self, epsilon):
		hasChanged = 0
		s = 1.5
		eps = epsilon
		layoutPrty = self.viewLayout
		newCoord = layoutPrty
	
		# calcul des dij pour tout couple i,j
		dij = {} # dictionnaire extremites arete -> dij
		for e in self.delaunayEdges:
			i = e[0]
			j = e[1]
			dist = self.dist(i,j)
			if self.edgesOverlap[(i.id,j.id)] > s:
				dij[(i.id,j.id)] = s*dist
				dij[(j.id,i.id)] = s*dist
			else:
				dij[(i.id,j.id)] = self.edgesOverlap[(i.id,j.id)]*dist
				dij[(j.id,i.id)] = self.edgesOverlap[(i.id,j.id)]*dist
		
		# calcul des dExm et dEym de depart
		dExm = self.graph.getDoubleProperty("derivee x")
		dEym = self.graph.getDoubleProperty("derivee y")
		Deltam = self.graph.getDoubleProperty("delta")
		for m in self.graph.getNodes():
			derivx = 0
			derivy = 0
			for index in self.delaunayNeighbors[m]:
				i = tlp.node(index)
				d = dij[(m.id,index)]
				diffx = layoutPrty[m][0] - layoutPrty[i][0]
				factx = 1/(d**2)*(diffx - (d*(diffx)/self.dist(i,m)))
				derivx += factx
				
				diffy = layoutPrty[m][1] - layoutPrty[i][1]
				facty = 1/(d**2)
				facty *= (diffy - (d*(diffy)/self.dist(i,m)))
				derivy += facty
			dExm[m] = derivx
			dEym[m] = derivy
			# recalcul des Deltam
			Deltam[m] = math.sqrt(dExm[m]**2 + dEym[m]**2)
		
		# boucle principale
		while(Deltam[self.findMax(Deltam)] > eps):
			hasChanged = 1
			m = self.findMax(Deltam)
			delta = Deltam[m]
			while (delta > eps):
				A = 0
				B = 0
				D = 0
				for index in self.delaunayNeighbors[m]:
					i = tlp.node(index)
					d = dij[(m.id,i.id)]
					diffx = layoutPrty[m][0] - layoutPrty[i][0]
					diffy = layoutPrty[m][1] - layoutPrty[i][1]
					wim = 1/(dij[(m.id,i.id)]**2)
					den = (diffx**2 + diffy**2)**(3/2)
					A += wim*(1 - (d*(diffy**2))/den)
					B += (wim*d*diffx*diffy)/den
					D += wim*(1 - (d*(diffx**2))/den)
				(dx, dy) = self.solveEquations(A,B,B,D,-dExm[m],-dEym[m])
				oldCoord = layoutPrty[m]
				oldCoord[0] += dx
				oldCoord[1] += dy
				newCoord[m] = oldCoord
				# calcul nouveau deltaM
				derivx = 0
				derivy = 0
				for index in self.delaunayNeighbors[m]:
					i = tlp.node(index)
					diffx = newCoord[m][0] - newCoord[i][0]
					factx = 1/(dij[(m.id,i.id)]**2)
					factx *= (diffx - (dij[(m.id,i.id)]*(diffx)/self.dist(i,m)))
					derivx += factx
				
					diffy = newCoord[m][1] - newCoord[i][1]
					facty = 1/(dij[(m.id,i.id)]**2)
					facty *= (diffy - (dij[(m.id,i.id)]*(diffy)/self.dist(i,m)))
					derivy += facty
				dExm[m] = derivx
				dEym[m] = derivy
				
				# calcul des Deltam
				Deltam[m] = math.sqrt(dExm[m]**2 + dEym[m]**2)
				delta = Deltam[m]
		for n in self.graph.getNodes():
			layoutPrty[n] = newCoord[n]
		return hasChanged
	
	def computeBounding(self):
		self.boundingBox = []
		self.boundingBoxToNode = {}
		self.nodeToBoundingBox = {}
		for n in self.graph.getNodes():
			subg = self.graph.inducedSubGraph([n])
			bb = tlp.computeBoundingBox(subg)
			self.boundingBoxToNode[bb] = n
			self.nodeToBoundingBox[n.id] = bb
			self.boundingBox.append(bb)
			self.graph.delSubGraph(subg)
	
	def scanline(self):
		listOverlapEdges = []
		self.boundingBox = sorted(self.boundingBox, key = lambda b : b[0][0])
		j = 0
		b = self.boundingBox[j]
		while (j < len(self.boundingBox)-2):
			i = self.boundingBox.index(b) +1
			b2 = self.boundingBox[i]
			while (b[1][0] > b2[0][0]) and i < len(self.boundingBox)-1:
				if b.intersect(b2):
					n1 = self.boundingBoxToNode[b]
					n2 = self.boundingBoxToNode[b2]
					listOverlapEdges.append((n1,n2))
				i+=1
				b2 = self.boundingBox[i]
			j+=1
			b = self.boundingBox[j]
		for (n,m) in listOverlapEdges:
			self.delaunayEdges.append((n,m))
			self.delaunayNeighbors.pushBackNodeEltValue(n,m.id)
			self.delaunayNeighbors.pushBackNodeEltValue(m,n.id)

def adjustNodesLabels(graph):
	per = tlpgui.createView("Node Link Diagram view", graph)
	visual = tlpgui.NodeLinkDiagramComponent.getInputData(per)
	labels = visual.getElementLabel()
	lengthLabels = graph.getIntegerProperty("longueur labels")

	for n in graph.getNodes():
		label = labels[n]
		lengthLabels[n] = len(label)
	
	size = visual.getElementFontSize()
	viewSize = graph.getSizeProperty("viewSize")
	for n in graph.getNodes():
		oldSize = viewSize[n]
		newSize = tlp.Vec3f(size[n]*lengthLabels[n],4*size[n],oldSize[2])
		viewSize.setNodeValue(n,newSize)
	visual.setElementSize(viewSize)
	
	rendering = per.getRenderingParameters()
	rendering.setLabelScaled(True)
	rendering.setLabelsDensity(100)
	per.setRenderingParameters(rendering)
	
def computeArea(graph):
	bb = tlp.computeBoundingBox(graph)
	return (bb[1][0] - bb[0][0])*(bb[1][1] - bb[0][1])

def main(graph): 
	viewLabel =  graph.getStringProperty("viewLabel")
	viewLabelBorderWidth =  graph.getDoubleProperty("viewLabelBorderWidth")
	viewLayout =  graph.getLayoutProperty("viewLayout")
	viewShape =  graph.getIntegerProperty("viewShape")
	
	# Ajustement des tailles des noeuds
	for n in graph.getNodes():
		viewShape[n] = tlp.NodeShape.Square
		viewLabelBorderWidth[n] = 0.0
	graphCopy = tlp.newGraph()
	tlp.copyToGraph(graphCopy,graph)
	adjustNodesLabels(graphCopy)
	adjustNodesLabels(graph)
	graph.setName("initial")
	graphCopy.setName("initial copy")
	initialArea = computeArea(graph)
	
	begin = time.time()
	p = Prism()
	p.applyToGraph(graph,0.00005)
	end = time.time()
	
	print "temps : ", end-begin
	finalArea = computeArea(graph)
	print "diff aire : ", finalArea/initialArea
