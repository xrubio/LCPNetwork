# -*- coding: utf-8 -*-
"""
/***************************************************************************
 LCPNetwork
                                 A QGIS plugin
 Compute the LCP network from multiple origins to multiple destinations
                             -------------------
        begin                : 2016-11-19
        copyright            : (C) 2016 by Xavier Rubio-Campillo
        email                : xavier.rubio@ed.ac.uk
        git sha              : $Format:%H$
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
 This script initializes the plugin, making it known to QGIS.
"""


# noinspection PyPep8Naming
def classFactory(iface):  # pylint: disable=invalid-name
    """Load LCPNetwork class from file LCPNetwork.

    :param iface: A QGIS interface instance.
    :type iface: QgsInterface
    """
    #
    from .lcp_network import LCPNetwork
    return LCPNetwork(iface)
