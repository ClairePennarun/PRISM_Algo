# Powered by Python 2.7

# To cancel the modifications performed by the script
# on the current graph, click on the undo button.

# Some useful keyboards shortcuts : 
#   * Ctrl + D : comment selected lines.
#   * Ctrl + Shift + D  : uncomment selected lines.
#   * Ctrl + I : indent selected lines.
#   * Ctrl + Shift + I  : unindent selected lines.
#   * Ctrl + Return  : run script.
#   * Ctrl + F  : find selected text.
#   * Ctrl + R  : replace selected text.
#   * Ctrl + Space  : show auto-completion dialog.

from tulip import *
from tulipgui import *
from tulipogl import *
import math

# the updateVisualization(centerViews = True) function can be called
# during script execution to update the opened views

# the pauseScript() function can be called to pause the script execution.
# To resume the script execution, you will have to click on the "Run script " button.

# the runGraphScript(scriptFile, graph) function can be called to launch another edited script on a tlp.Graph object.
# The scriptFile parameter defines the script name to call (in the form [a-zA-Z0-9_]+.py)

# the main(graph) function must be defined 
# to run the script on the current graph


def computeOverlapFactor(graph,i,j):
#	print (i,j)
	viewSize = graph.getSizeProperty("viewSize") #tableau des tailles
	viewLayout = graph.getLayoutProperty("viewLayout")
	sizeI = viewSize[i]
	sizeJ = viewSize[j]
	posI = viewLayout[i]
	posJ = viewLayout[j]
	#print posI, posJ
	a = abs(posI[0] - posJ[0])
	firstFactor = 1
	secondFactor = 1
	factor = 1
	#print a
	if (a <> 0):
		firstFactor = (sizeI[0]/2 + sizeJ[0]/2)/a
		factor = firstFactor
	b = abs(posI[1] - posJ[1])
	if (b <> 0):
		secondFactor = (sizeI[1]/2 + sizeJ[1]/2)/b
		if (secondFactor < firstFactor):
			factor = secondFactor
#	print factor
	if factor > 1:
		return factor
	return 1

def adjustNodesLabels(graph):
	per = tlpgui.createView("Node Link Diagram view", graph)
	visual = tlpgui.NodeLinkDiagramComponent.getInputData(per)
	
	labels = visual.getElementLabel()
	lengthLabels = graph.getIntegerProperty("longueur labels")

	for n in graph.getNodes():
		label = labels[n]
		lengthLabels[n] = len(label)
	
#	viewFont =  graph.getStringProperty("viewFont")
#	for n in graph.getNodes():
#		viewFont[n] = "usr/local/share/tulip/bitmaps/fonts/FreeMono/FreeMono.ttf"
#	visual.setElementFont(viewFont)
	
	
	size = visual.getElementFontSize()
	viewSize = graph.getSizeProperty("viewSize")
	for n in graph.getNodes():
		oldSize = viewSize[n]
		newSize = tlp.Vec3f(size[n]*lengthLabels[n],4*size[n],oldSize[2])
		viewSize.setNodeValue(n,newSize)
	visual.setElementSize(viewSize)


def main(graph): 
	externLabel =  graph.getStringProperty("externLabel")
	viewBorderColor =  graph.getColorProperty("viewBorderColor")
	viewBorderWidth =  graph.getDoubleProperty("viewBorderWidth")
	viewColor =  graph.getColorProperty("viewColor")
#	viewFont =  graph.getStringProperty("viewFont")
	viewFontSize =  graph.getIntegerProperty("viewFontSize")
	viewLabel =  graph.getStringProperty("viewLabel")
	viewLabelBorderColor =  graph.getColorProperty("viewLabelBorderColor")
	viewLabelBorderWidth =  graph.getDoubleProperty("viewLabelBorderWidth")
	viewLabelColor =  graph.getColorProperty("viewLabelColor")
	viewLabelPosition =  graph.getIntegerProperty("viewLabelPosition")
	viewLabelRotation =  graph.getDoubleProperty("viewLabelRotation")
	viewLayout =  graph.getLayoutProperty("viewLayout")
	viewMetaGraph =  graph.getGraphProperty("viewMetaGraph")
	viewMetric =  graph.getDoubleProperty("viewMetric")
	viewRotation =  graph.getDoubleProperty("viewRotation")
	viewSelection =  graph.getBooleanProperty("viewSelection")
	viewShape =  graph.getIntegerProperty("viewShape")
	#viewSize =  graph.getSizeProperty("viewSize")
	viewSrcAnchorShape =  graph.getIntegerProperty("viewSrcAnchorShape")
	viewSrcAnchorSize =  graph.getSizeProperty("viewSrcAnchorSize")
	viewTexture =  graph.getStringProperty("viewTexture")
	viewTgtAnchorShape =  graph.getIntegerProperty("viewTgtAnchorShape")
	viewTgtAnchorSize =  graph.getSizeProperty("viewTgtAnchorSize")
	
	adjustNodesLabels(graph)

	coord = []

	for n in graph.getNodes():
		coord.append(viewLayout[n]);
	
	(delaunayEdges, delaunayNodes) = tlp.delaunayTriangulation(coord)
	delaunay = tlp.newGraph()
	tlp.copyToGraph(delaunay,graph)
	delaunay.delEdges(delaunay.getEdges())
	
	for e in delaunayEdges:
		delaunay.addEdge(tlp.node(e[0]), tlp.node(e[1]))

	#tlp.saveGraph(delaunay,"delaunayGraph.tlp")
	
	overlapFactors = [[0 for x in xrange(delaunay.numberOfNodes())] for x in xrange(delaunay.numberOfNodes())] 
	
	for e in delaunay.getEdges():
		src = delaunay.source(e)
		tgt = delaunay.target(e)
		overlapFactors[src.id][tgt.id] = computeOverlapFactor(delaunay, src, tgt)
