# -*- coding: utf-8 -*-
"""
/***************************************************************************
 LCPNetworkDialog
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

import os

from PyQt5 import uic
from PyQt5.QtWidgets import QDialog

FORM_CLASS, _ = uic.loadUiType(os.path.join(
    os.path.dirname(__file__), 'lcp_network_dialog_base.ui'))


class LCPNetworkDialog(QDialog, FORM_CLASS):
    def __init__(self, parent=None):
        """Constructor."""
        super(LCPNetworkDialog, self).__init__(parent)
        # Set up the user interface from Designer.
        # After setupUI you can access any designer object by doing
        # self.<objectname>, and you can use autoconnect slots - see
        # http://qt-project.org/doc/qt-4.8/designer-using-a-ui-file.html
        # #widgets-and-dialogs-with-auto-connect
        self.setupUi(self)
