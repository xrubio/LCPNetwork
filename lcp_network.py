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
from PyQt4.QtGui import QAction, QIcon, QFileDialog    

from qgis.core import QgsMapLayer, QgsMapLayerRegistry, QGis, QgsPoint

# Initialize Qt resources from file resources.py
import resources
# Import the code for the dialog
from lcp_network_dialog import LCPNetworkDialog
import os.path
from osgeo import gdal

import numpy


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
        self.dlg.outputFile.clear()
        self.dlg.browseOutput.clicked.connect(self.selectOutputFile)

    def selectOutputFile(self):
        fileName = QFileDialog.getSaveFileName(self.dlg, "Select output file ","", '*.tif')
        self.dlg.outputFile.setText(fileName)

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
                self.dlg.pointsLayer.addItem(layer.name(), layer.id())
                
    def clearUI(self):
        self.dlg.baseLayer.clear()
        self.dlg.pointsLayer.clear()

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


    def runBaseAlgorithm(self, layers):
        outputName = self.dlg.outputFile.text()
        outputFile = open(outputName, 'w')
        
        # select vector layer
        selectedInputIndex = self.dlg.inputLayerList.currentIndex()
        selectedInput = layers[selectedInputIndex]
        fields = selectedInput.pendingFields()
        fieldNames = [field.name() for field in fields]

        for feature in selectedInput.getFeatures():
            line = ';'.join(unicode(feature[x]) for x in fieldNames) + '\n'
            unicodeLine = line.encode('utf-8')
            outputFile.write(unicodeLine)
        outputFile.close()
    


    def loadPoints(self):   
        index = self.dlg.pointsLayer.currentIndex()
        layer = self.dlg.pointsLayer.itemData(index)
        return QgsMapLayerRegistry.instance().mapLayer(layer)

    def loadBaseRaster(self):
        index = self.dlg.baseLayer.currentIndex()
        layer = self.dlg.baseLayer.itemData(index)
        path= str(QgsMapLayerRegistry.instance().mapLayer(layer).dataProvider().dataSourceUri())
        return gdal.Open(path)

    def runAlgorithm(self):
        # TODO 1 - add new raster map to canvas
        # TODO 2 - only selected features

        print('loading points and base raster layers')
        points = self.loadPoints()
        baseRaster = self.loadBaseRaster()

        print('rasterizing '+str(points.featureCount())+' points')
        for feature in points.getFeatures(): 
            point = feature.geometry().asPoint()


        transform = baseRaster.GetGeoTransform()    
        print('top-left pixel: '+str(transform[0])+'/'+str(transform[3])+' res: '+str(transform[1])+'/'+str(transform[5])+' dims: '+str(baseRaster.RasterXSize)+'/'+str(baseRaster.RasterYSize))

        topLeft = QgsPoint(transform[0], transform[3])

        outputValues = numpy.zeros([baseRaster.RasterYSize, baseRaster.RasterXSize])

        for feature in points.getFeatures(): 
            point = feature.geometry().asPoint()
            pointInRaster = QgsPoint(point.x() - topLeft.x(), topLeft.y() - point.y())
            cell = QgsPoint(pointInRaster.x()/transform[1], pointInRaster.y()/-transform[5])
            print('\t'+point.toString()+' -> '+pointInRaster.toString() + '-> ' + cell.toString())

            if cell.x() < 0 or cell.x() >= baseRaster.RasterXSize or cell.y() <0 or cell.y() >= baseRaster.RasterYSize :
                continue

            outputValues[cell.y(),cell.x()] = 1

        outputName = self.dlg.outputFile.currentText()
        newRaster = gdal.GetDriverByName('GTiff').Create(outputName, baseRaster.RasterXSize, baseRaster.RasterYSize, 1, gdal.GDT_Int32)
        newRaster.SetProjection(baseRaster.GetProjection())
        newRaster.SetGeoTransform(baseRaster.GetGeoTransform())
        newRaster.GetRasterBand(1).WriteArray(outputValues,0,0)
                    

