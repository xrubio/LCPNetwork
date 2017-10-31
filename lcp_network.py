# -*- coding: utf-8 -*-
"""
/***************************************************************************
 LCPNetwork
                                 A QGIS plugin
 Compute the LCP network from multiple origins to multiple destinations
                              -------------------
        begin                : 2016-11-19
        git sha              : $Format:%H$
        copyright            : (C) 2016 by Xavier Rubio-Campillo
        email                : xavier.rubio@ed.ac.uk
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
"""
from PyQt4.QtCore import QSettings, QTranslator, qVersion, QCoreApplication
from PyQt4.QtGui import QAction, QIcon, QFileDialog, QMessageBox, QProgressBar 

from qgis.core import QgsMapLayer, QgsMapLayerRegistry, QGis, QgsPoint

# Initialize Qt resources from file resources.py
import resources
# Import the code for the dialog
from lcp_network_dialog import LCPNetworkDialog
import os.path
from osgeo import gdal

import numpy as np
from qgis.core import QgsMessageLog
from qgis.core import QgsRasterLayer
from qgis.core import QgsVectorLayer
from qgis.core import QgsContrastEnhancement
from qgis.core import QgsFeature
from qgis.core import QgsGeometry
 
import numpy.ma as ma

import timeit

class LCPNetwork:
    """QGIS Plugin Implementation."""

    def __init__(self, iface):
        """Constructor.

        :param iface: An interface instance that will be passed to this class
            which provides the hook by which you can manipulate the QGIS
            application at run time.
        :type iface: QgsInterface
        """
        # Save reference to the QGIS interface
        self.iface = iface
        # initialize plugin directory
        self.plugin_dir = os.path.dirname(__file__)
        # initialize locale
        locale = QSettings().value('locale/userLocale')[0:2]
        locale_path = os.path.join(
            self.plugin_dir,
            'i18n',
            'LCPNetwork_{}.qm'.format(locale))

        if os.path.exists(locale_path):
            self.translator = QTranslator()
            self.translator.load(locale_path)

            if qVersion() > '4.3.3':
                QCoreApplication.installTranslator(self.translator)


        # Declare instance attributes
        self.actions = []
        self.menu = self.tr(u'&Least-Cost-Paths Network')
        # TODO: We are going to let the user set this up in a future iteration
        self.toolbar = self.iface.addToolBar(u'LCPNetwork')
        self.toolbar.setObjectName(u'LCPNetwork')

        # Create the dialog (after translation) and keep reference
        self.dlg = LCPNetworkDialog()
#        self.dlg.outputFile.clear()
#      self.dlg.browseOutput.clicked.connect(self.selectOutputFile)

    """
    def selectOutputFile(self):
        fileName = QFileDialog.getSaveFileName(self.dlg, "Select output file ","", '*.tif')
        self.dlg.outputFile.setText(fileName)
    """

    # noinspection PyMethodMayBeStatic
    def tr(self, message):
        """Get the translation for a string using Qt translation API.

        We implement this ourselves since we do not inherit QObject.

        :param message: String for translation.
        :type message: str, QString

        :returns: Translated version of message.
        :rtype: QString
        """
        # noinspection PyTypeChecker,PyArgumentList,PyCallByClass
        return QCoreApplication.translate('LCPNetwork', message)


    def add_action(
        self,
        icon_path,
        text,
        callback,
        enabled_flag=True,
        add_to_menu=True,
        add_to_toolbar=True,
        status_tip=None,
        whats_this=None,
        parent=None):
        """Add a toolbar icon to the toolbar.

        :param icon_path: Path to the icon for this action. Can be a resource
            path (e.g. ':/plugins/foo/bar.png') or a normal file system path.
        :type icon_path: str

        :param text: Text that should be shown in menu items for this action.
        :type text: str

        :param callback: Function to be called when the action is triggered.
        :type callback: function

        :param enabled_flag: A flag indicating if the action should be enabled
            by default. Defaults to True.
        :type enabled_flag: bool

        :param add_to_menu: Flag indicating whether the action should also
            be added to the menu. Defaults to True.
        :type add_to_menu: bool

        :param add_to_toolbar: Flag indicating whether the action should also
            be added to the toolbar. Defaults to True.
        :type add_to_toolbar: bool

        :param status_tip: Optional text to show in a popup when mouse pointer
            hovers over the action.
        :type status_tip: str

        :param parent: Parent widget for the new action. Defaults None.
        :type parent: QWidget

        :param whats_this: Optional text to show in the status bar when the
            mouse pointer hovers over the action.

        :returns: The action that was created. Note that the action is also
            added to self.actions list.
        :rtype: QAction
        """

        icon = QIcon(icon_path)
        action = QAction(icon, text, parent)
        action.triggered.connect(callback)
        action.setEnabled(enabled_flag)

        if status_tip is not None:
            action.setStatusTip(status_tip)

        if whats_this is not None:
            action.setWhatsThis(whats_this)

        if add_to_toolbar:
            self.toolbar.addAction(action)

        if add_to_menu:
            self.iface.addPluginToMenu(
                self.menu,
                action)

        self.actions.append(action)

        return action

    def initGui(self):
        """Create the menu entries and toolbar icons inside the QGIS GUI."""

        icon_path = ':/plugins/LCPNetwork/icon.png'
        self.add_action(
            icon_path,
            text=self.tr(u'LCP Network'),
            callback=self.run,
            parent=self.iface.mainWindow())


    def unload(self):
        """Removes the plugin menu item and icon from QGIS GUI."""
        for action in self.actions:
            self.iface.removePluginMenu(
                self.tr(u'&Least-Cost-Paths Network'),
                action)
            self.iface.removeToolBarIcon(action)
        # remove the toolbar
        del self.toolbar


    def loadLayers(self): 
        layers = self.iface.mapCanvas().layers()
        for layer in layers:
            if layer.type() == QgsMapLayer.RasterLayer:
                self.dlg.baseLayer.addItem(layer.name(), layer.id())
            elif layer.geometryType() == QGis.Point:
                self.dlg.origins.addItem(layer.name(), layer.id())
                self.dlg.destinations.addItem(layer.name(), layer.id())
                
    def clearUI(self):
        self.dlg.baseLayer.clear()
        self.dlg.origins.clear()
        self.dlg.destinations.clear()

    def run(self):
        """Run method that performs all the real work"""
        self.clearUI()
        self.loadLayers()
        # show the dialog
        self.dlg.show()
        # Run the dialog event loop
        result = self.dlg.exec_()
        # See if OK was pressed
        if result:
            self.runAlgorithm()

    def loadPoints(self):   
        indexO = self.dlg.origins.currentIndex()
        layerO = self.dlg.origins.itemData(indexO)
        originLayer = QgsMapLayerRegistry.instance().mapLayer(layerO)

        indexD = self.dlg.destinations.currentIndex()
        layerD = self.dlg.destinations.itemData(indexD)
        destinationLayer = QgsMapLayerRegistry.instance().mapLayer(layerD)

        return originLayer,destinationLayer

    def loadBaseRaster(self):
        index = self.dlg.baseLayer.currentIndex()
        layer = self.dlg.baseLayer.itemData(index)
        path= str(QgsMapLayerRegistry.instance().mapLayer(layer).dataProvider().dataSourceUri())
        return gdal.Open(path)

    def getCell( self, point, surface ):
        """ get local coordinates for point in surface """

        transform = surface.GetGeoTransform()
        topLeft = QgsPoint(transform[0], transform[3])

        pointInRaster = QgsPoint(point.x() - topLeft.x(), topLeft.y() - point.y())
        # swap axes
        cell = QgsPoint(int(pointInRaster.y()/-transform[5]), int(pointInRaster.x()/transform[1]))
        return cell

    def getGlobalPos(self, localPos, surface):
        """ get global coordinates for local point in surface """

        transform = surface.GetGeoTransform()
        topLeft = QgsPoint(transform[0], transform[3])

        # swap axes
        pos = QgsPoint(localPos.y()*(-transform[5]), localPos.x()*(-transform[1]))
        globalPoint = QgsPoint(pos.x()+topLeft.x(), pos.y()+topLeft.y())
        
        return globalPoint

    def isInside(self, cell, surface ):
        """ returns true if cellin surface or false if it is not """
        if cell.x() < 0 or cell.x() >= surface.RasterYSize or cell.y() <0 or cell.y() >= surface.RasterXSize :
            return False
        return True


    def getNeighbors(self, point, surface, costMap ):
        """ current: only four direct neighbors """
        neighbors = list()
        nodata = int(surface.GetRasterBand(1).GetNoDataValue())

        candidate = QgsPoint(point.x()-1, point.y())
        if self.isInside(candidate, surface) and costMap[candidate.x(), candidate.y()] != nodata:
            neighbors.append(candidate)

        candidate = QgsPoint(point.x()+1, point.y())
        if self.isInside(candidate, surface) and costMap[candidate.x(), candidate.y()] != nodata:
            neighbors.append(candidate)

        candidate = QgsPoint(point.x(), point.y()-1)
        if self.isInside(candidate, surface) and costMap[candidate.x(), candidate.y()] != nodata:
            neighbors.append(candidate)

        candidate = QgsPoint(point.x(), point.y()+1)
        if self.isInside(candidate, surface) and costMap[candidate.x(), candidate.y()] != nodata:
            neighbors.append(candidate)

#        QgsMessageLog.logMessage("num neighbors of: "+point.toString(0) + " is: " + str(len(neighbors)), "LCPNetwork")
        return neighbors            

    def getMinimumUnvisited(self, visited, distances ):

        # set to nul values of already visited (kind of a mask)
        possibleValues = ma.masked_array(distances, mask=visited)
        candidates = np.where(possibleValues == np.nanmin(possibleValues))  

        if len(candidates[0]) == 0:
            return None

        selected = np.random.randint(len(candidates[0]))
        return QgsPoint(candidates[0][selected], candidates[1][selected])

    def computeCost( self, originGeo, baseRaster):        
        origin = self.getCell(originGeo, baseRaster)
        costValues = baseRaster.GetRasterBand(1).ReadAsArray()

        if not self.isInside(origin, baseRaster):
            return None

        # initialize helper matrices
        width,height = costValues.shape
        visited = np.full([width, height], False, dtype=bool)
        distances = np.full([width,height], np.nan, dtype=np.float32)

        # initialize current
        current = origin

        visited[current.x(), current.y()] = True
        distances[current.x(), current.y()] = 0

        candidates = True
        while candidates: 
            neighbors = self.getNeighbors(current, baseRaster, costValues)
            for neighbor in neighbors:
                
                tentativeDistance = distances[current.x(), current.y()] + costValues[neighbor.x(), neighbor.y()]
                # cost can never be negative
                if tentativeDistance < 0:
                    tentativeDistance = 0
                if np.isnan(distances[neighbor.x(), neighbor.y()]) or distances[neighbor.x(), neighbor.y()] > tentativeDistance:
                    distances[neighbor.x(), neighbor.y()] = tentativeDistance

            visited[current.x(), current.y()] = True

            current = self.getMinimumUnvisited(visited, distances)
            if not current:
                candidates = False

        return distances



    def getPath( self, originGeo, destinationGeo, baseRaster, costMap):
    
        pathLine = []
        origin = self.getCell(originGeo, baseRaster)
        destination = self.getCell(destinationGeo, baseRaster)
            
        if not self.isInside(destination, baseRaster):
            return None

        width,height = costMap.shape
        
        current = destination
    
        while current != origin:
            pathLine.append(current)
            neighbors = self.getNeighbors(current, baseRaster, costMap)
        
            minValue = costMap[current.x(), current.y()]
       
            for neighbor in neighbors:
        
                # if already in path:
                alreadyInPath = False 
                for pathPoint in pathLine:
                    if pathPoint.sqrDist(neighbor)<1.0:
                        alreadyInPath = True
                        break

                if alreadyInPath:
                    continue

                if costMap[neighbor.x(),neighbor.y()] <= minValue:
                    minValue = costMap[neighbor.x(),neighbor.y()]
                    current = neighbor


        pathLine.append(current)

        globalPath = []
        for localPoint in pathLine:
            globalPath.append(self.getGlobalPos(localPoint, baseRaster))

        return globalPath


    def storeCostMap(self, costMap, baseRaster, index):
        outputName = os.path.dirname(__file__)+"/distances_"+str(index)+".tif"

        newRaster = gdal.GetDriverByName('GTiff').Create(outputName, baseRaster.RasterXSize, baseRaster.RasterYSize, 1, gdal.GDT_Float32)
        newRaster.SetProjection(baseRaster.GetProjection())
        newRaster.SetGeoTransform(baseRaster.GetGeoTransform())
        newRaster.GetRasterBand(1).WriteArray(costMap,0,0)
        newRaster.GetRasterBand(1).SetNoDataValue(np.nan)
        newRaster.GetRasterBand(1).FlushCache()
        newRaster = None
        
        newRasterQGIS = QgsRasterLayer(outputName, "distances_"+str(index))
        newRasterQGIS.setContrastEnhancement(QgsContrastEnhancement.StretchToMinimumMaximum)
        QgsMapLayerRegistry.instance().addMapLayer(newRasterQGIS)

    def runAlgorithm(self):
        start = timeit.default_timer()

        QgsMessageLog.logMessage('INIT LCPNetwork plugin - loading points and base raster layers', "LCPNetwork")
        origins,destinations = self.loadPoints()
        baseRaster = self.loadBaseRaster()

        transform = baseRaster.GetGeoTransform()
        nodata = baseRaster.GetRasterBand(1).GetNoDataValue()  
        QgsMessageLog.logMessage("NO DATA: "+str(nodata), "LCPNetwork")

#        baseRaster.GetRasterBand(1).SetNoDataValue(-32768)
#        QgsMessageLog.logMessage('top-left pixel: '+str(transform[0])+'/'+str(transform[3])+' res: '+str(transform[1])+'/'+str(transform[5])+' dims: '+str(baseRaster.RasterXSize)+'/'+str(baseRaster.RasterYSize), "LCPNetwork")
#        QgsMessageLog.logMessage("base surface size: "+str(baseRaster.RasterXSize)+"x"+str(baseRaster.RasterYSize), "LCPNetwork")
        topLeft = QgsPoint(transform[0], transform[3])
    
        pointsListO = []
        for point in origins.getFeatures():
            pointsListO.append(point.geometry().asPoint())
  
        pointsListD = []
        for point in destinations.getFeatures():
            pointsListD.append(point.geometry().asPoint())


        index = 0          

        ## create the shapefile
        lineLayer = self.iface.addVectorLayer("LineString?crs=" + baseRaster.GetProjection(), "least cost path network", "memory")

#        progressBar = QProgressBar()
#        progressBar.setMaximum(len(pointsListO)*len(pointsListD))
#        self.iface.mainWindow().statusBar().addWidget(progressBar)

        for source in pointsListO:

            # compute cost map for the entire area
               
            startCost = timeit.default_timer()
            distances = self.computeCost(source, baseRaster)
            QgsMessageLog.logMessage("seconds to compute cost map: " + str("%.2f"%(timeit.default_timer()-startCost)), "LCPNetwork")

            if distances==None:
                QMessageBox.information(None, "ERROR!", "Cost map could not be computed for source point: "+str(index), "LCPNetwork")
                return 

            self.storeCostMap(distances, baseRaster, index)

            for destination in pointsListD:
                if destination == source:
                    continue
                pathLine = self.getPath(source, destination, baseRaster, distances)
                if pathLine==None:
                    QMessageBox.information(None, "ERROR!", "Route could not be found for destination: "+str(index), "LCPNetwork")
                    return
                
                features = QgsFeature()
                ## set geometry from the list of QgsPoint's to the feature
                features.setGeometry(QgsGeometry.fromPolyline(pathLine))
                lineLayer.dataProvider().addFeatures([features])
            index = index + 1 
#            progressBar.setValue(index+1)
    
        QgsMessageLog.logMessage("seconds for plugin: " + str("%.2f"%(timeit.default_timer()-start)), "LCPNetwork")

