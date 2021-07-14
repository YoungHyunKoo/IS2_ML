import ee
from ee_plugin import Map
import pandas as pd
import glob, os
import numpy as np

def linear_color_interpolation(c1, c2, t):
    r = int(c1.red() + (c2.red() - c1.red()) * t)
    g = int(c1.green() + (c2.green() - c1.green()) * t)
    b = int(c1.blue() + (c2.blue() - c1.blue()) * t)
    return QColor(r, g, b)

def symbolize_layer(layer):
    single_symbol_renderer = layer.renderer()
    symbol = single_symbol_renderer.symbol()
    symbol.symbolLayer(0).setStrokeWidth(0)
    
    renderer = layer.renderer()
    provider = layer.dataProvider()

    #Specify inputs:
    fld = 'height'
    rampname = 'Magma'
    number_classes = 5

    #COLOR RAMP STYLE:
    myStyle = QgsStyle().defaultStyle()
    ramp = myStyle.colorRamp(rampname)

    vals = []
    for f in layer.getFeatures():
        vals.append(f[fld])
    vals = (list(filter(None,vals))) #REMOVE NONE/NULL from list

    idx = provider.fieldNameIndex(fld)
    max = layer.maximumValue(idx)
    min = layer.minimumValue(idx)
    interval = (max - min)/(number_classes -1 )

    #classes
    sum = min
    classes = [-1000]

    min_val = -0.6
    max_val = 1.0
    interval = 0.2

    for val in np.arange(min_val, max_val, interval):
        classes.append(val)
    classes.append(1000)

    t_list = []
    init = 0
    num_class = len(classes)

    for c in np.linspace(0, 1, num_class-1):
        t_list.append(c)

    myRangeList = []
    for i in range(0, num_class-1):
        symbol = QgsSymbol.defaultSymbol(layer.geometryType())
        symbol.setColor(QColor(ramp.color(t_list[i]).name()))
        myRange = QgsRendererRange(classes[i], classes[i+1], symbol, '%.2f-%.2f'%(classes[i], classes[i+1]))                   
        myRangeList.append(myRange)
        
    myRenderer = QgsGraduatedSymbolRenderer(fld, myRangeList)
    myRenderer.setMode(QgsGraduatedSymbolRenderer.Custom)
    layer.setRenderer(myRenderer)


year = 2020
df = pd.read_csv(glob.glob("F:\\ATL03\\ATL07\\RossSea\\ATL07_S2_{0}*.csv".format(year))[0], index_col = 0)
idx = []

shpdir = "F:\\ATL03\\ATL07\\RossSea\\shpfile\\"

for i in range(0, len(df)):
    shpfiles = glob.glob(shpdir + df['IS2_file'][i].replace("_004_01.h5", "*.shp"))
    
    if len(shpfiles) > 0:
        idx.append(True)
    else:
        idx.append(False)

df = df[idx].reset_index(drop = True)

uname, uidx = np.unique(df['IS2_file'], return_index = True)

for i in uidx:
        
    start_date = df['start_date'][i]
    end_date = df['end_date'][i]

    ross_sea = ee.Geometry.Rectangle([-180, -78, -140, -70])

    image = ee.ImageCollection("COPERNICUS/S2_SR").\
    filterDate(start_date, end_date).filterBounds(ross_sea)

    ids_S2 = df['S2_id'][df['IS2_file'] == df['IS2_file'][i]]
    
    for id in ids_S2:
        image = ee.Image(id)
        visParams = {'bands': ['B5', 'B4', 'B3'], 'min': 0, 'max': 10000}
        Map.addLayer(image, visParams, os.path.basename(id))
        print(id)
        
    shpfiles = glob.glob(shpdir + df['IS2_file'][i].replace("_004_01.h5", "*.shp"))
    for k in range(0, len(shpfiles)):
        layer = iface.addVectorLayer(shpfiles[k], "", "ogr")
        symbolize_layer(layer)
        print(shpfiles[k])

