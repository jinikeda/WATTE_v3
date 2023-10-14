# Wave Attenuation Toolbox (WATTE) version 3.0
# Developed by Center for Computation & Technology and Coastal Ecosystem Design Studio at Louisiana State University (LSU).
# WATTE version 1 was originally developed by M. Foster-Martinez, University of New Orleans: Karim Alizad, University of South Carolina: Scott C. Hagen, LSU (https://digitalcommons.lsu.edu/civil_engineering_data/1/)
# WATTE version 2 was developed using Python3 and QGIS in 2023
# Developer: Jin Ikeda, Shu Gao, Christopher E. Kees, and Peter Bacopoulos
# Last modified 10/14/2023
# Deepest thanks: This software is dedicated to in memory of Dr. Scott C. Hagen. We are truly blessed to have been one of your pupils. We will spread your idea further.
#
# WATTE version 3 is an open source-based toolbox using Python 3 (Verified at version Python 3.10).
# This model estimates and maps wave attenuation along marsh coastlines following an exponential decay.
#
# This software is released under the MIT License, see LICENSE.txt.
#
### Step 1 #############################################################
print("\n#############################################################\n")
print("Step 1: Import modules and inputs parameters")
print("\n#############################################################\n")
########################################################################

### 1.1 Import modules ###
print("Step 1.1: Import modules")

# import os module
import os, sys
import numpy as np
import pandas as pd
import csv
import warnings
import geopandas as gpd
from shapely.geometry import LineString
import glob
import time

# import spatial analysis using Gdal
from osgeo import gdal, ogr, osr
from osgeo.gdalconst import *
gdal.AllRegister()  # Register all of drivers
import rasterio

startTime = time.time()

######################################################################
# Working directory and input parameters
######################################################################

# Print the current working directory
print("\tCurrent working directory: {0}".format(os.getcwd()))

# The target Working directory
Workspace = os.getcwd()
# Workspace = "Z:/CCR_data/ACTIVE/Jin_Projects/WATTE/Shu/Model_development/Test8_Grand_Bay_v3_dev/"

# Change the current working directory
os.chdir(Workspace)

### 1.2 Input parameters and file ###
print("Step 1.2: Input parameters and file")

# Dominant wave direction North True
Wave_direction = float(150)  # North is 0 [degree] and clockwise
Wave_direction_range = 179  # Criterion for consider transects
# Note: Wave_direction = float(input("Dominant wave direction in North True: [degree] between 0-360 type here->"))
assert 0 <= Wave_direction < 360, "\tPlease input the values between 0-360 degree ...{{ (>_<) }}\n"

# Input raster data
Inputspace = os.path.join(Workspace, "Input_data")
Raster_file = os.path.join(Inputspace, "Grand_Bay_2000.tif")

##################################################################################################################
# Input baseline delineation method
##################################################################################################################
# Type 1, Method 1 Land (1): Water (0) isopleth with a moving average method
# Type 2, Method 2 1-way scan based on significant wave direction
# Type 3, Method 3 Baseline delineation based on a manual polyline also provide the path
# Manual_line= "Z:/CCR_data/ACTIVE/Jin_Projects/WATTE/EMBIR/Test_3/Input_data/Polyline.shp"

# Here polyline is estimated to draw from left side to right side, which may affect offset direction (use E2W function)

Baseline_delineation_method = 3
Manual_line = os.path.join(Inputspace, "Polyline.shp")
assert 1 <= Baseline_delineation_method <= 3, "\tPlease input the values between 1-3 ...{{ (>_<) }}\n"
##################################################################################################################

# Classification(s), Water, String
Input_Water_Class = "40"  # Don't use 0 (no data)

# Classification(s), Other, String
Input_Other_Class = "55"

# Classification(s), Marsh, String
Input_Marsh_Class = "16,23,32"
Marsh_Class = [int(x) for x in Input_Marsh_Class.split(',')]

# Decay constant for each classification of water and marshes, String
Input_Decay_Constant = "0.021,0.030,0.090"
Decay_Constant_M = [float(x) for x in Input_Decay_Constant.split(',')]

# Transect length [m]
Input_Transect_length = 1000

# Distance between transects [m]
Input_Distance_Transects = 10

# Filter the wave attenuation values greater than a threshold
Input_Wave_Attenuation_Threshold = 0.02

# Spacing of points along transect [m]
Input_Spacing = 10  # dx = (x_i+1−x_i)

# Interpolation_method = "IDW"          # currently IDW only

# Make process folder
Process_folder = os.path.join(Workspace, 'Process')
os.makedirs(Process_folder, exist_ok=True)

# Make process folder
Outputspace = os.path.join(Workspace, 'Output_data')
os.makedirs(Outputspace, exist_ok=True)

########################### Messages #################################
print("\tDominant wave direction: ", Wave_direction, "[degree]")
print("\tInput raster dataset is " + Raster_file)
print("\tWater classification is ", end=" ")
print(*Input_Water_Class, sep=" ")
Input_Water_Class = int(Input_Water_Class)
print("\tOther classification is ", end=" ")
print(*Input_Other_Class, sep=" ")
Input_Other_Class = int(Input_Other_Class)
print("\tMarsh classification is ", end=" ")
print(*Marsh_Class)
print("\tDecay constant for each classification of marsh is ", end=" ")
print(*Decay_Constant_M)

# print(Decay_Constant_M)
print("\tTransect length is ", Input_Transect_length, "meter")
print("\tDistance between transects is ", Input_Distance_Transects, "meter")

Input_list = [Input_Water_Class] + [Input_Other_Class] + Marsh_Class
print("\tInput_list：", Input_list)

# print ('Wave decay stopping point is', Input_Wave_decay)
print("\tSpacing of points along transect is ", Input_Spacing, "meter")

######################################################################
# Functions
######################################################################

# Function of output x, y coordinates
def getXY(pt):
    return (pt.x, pt.y)

# Function of reordering East to West direction
def E2W(gdf):

    point_start = gdf['East'].iloc[0]
    point_end = gdf['East'].iloc[-1]

    if point_start < point_end:  # The string order is East to West
        pass
    else:
        gdf = gdf.loc[::-1]  # reorder the contents When polyline direction is West to East side

    return gdf

# Function of add East and West coordinate and save as a csv file
def createEW(in_shp):
    gdf = gpd.read_file(in_shp)
    centroidseries = gdf['geometry'].centroid
    East, North = [list(t) for t in zip(*map(getXY, centroidseries))]

    gdf['East'] = East
    gdf['North'] = North

    gdf = E2W(gdf)
    out_csv = in_shp.split(".")[0].split("/")[-1] + '.csv'  # Linux
    Output_point_csv = os.path.join(Process_folder, out_csv)
    gdf.to_csv(Output_point_csv)

    return gdf


# Function of moving average of points
def moving_gdf(gdf,moving_value,transiton_points):
    if moving_value > transiton_points/3:
        print ("check the size of window and transiton")
    # Perform the rolling mean on 'East' and 'North' columns
    gdf[['East', 'North']] = gdf[['East', 'North']].rolling(moving_value,
                                                            min_periods=int(transiton_points/3),
                                                            center=True).mean()
    # Fill NaN values at the beginning and end with original values
    gdf.loc[:int(transiton_points), ['East', 'North']] = gdf.loc[:int(transiton_points), ['East', 'North']].values
    gdf.loc[gdf.shape[0] - int(transiton_points):, ['East', 'North']] = gdf.loc[gdf.shape[0] - int(transiton_points):, ['East','North']].values

    out_csv = "Shoreline_Delineation_Points_MA" + str(int(moving_value)) + ".csv."
    Output_point_csv = os.path.join(Process_folder, out_csv)
    gdf.to_csv(Output_point_csv)

    return gdf


# Function of create a smoothed polyline
def create_smoothed_polyline(gdf,out_shp):
    # Create a LineString geometry from the smoothed points
    smoothed_line = LineString(gdf[['East', 'North']].values)
    smoothed_gdf = gpd.GeoDataFrame(geometry=[smoothed_line])
    smoothed_gdf['ID'] = 1
    smoothed_gdf = smoothed_gdf.set_crs(gdf.crs, allow_override=True)

    smoothed_gdf.to_file(out_shp)

    del smoothed_gdf


# Function to extract point values using a raster file
def extract_raster_value(in_shp, in_raster, out_shp, variable_name, NoData_value):
    NoData_value = NoData_value

    # Open shp file
    pointData = gpd.read_file(in_shp)
    x = pointData['x']
    y = pointData['y']

    with rasterio.open(in_raster) as Raster:
        rows_raster = Raster.height
        cols_raster = Raster.width

        # Vectorized point-to-pixel index conversion
        row, col = Raster.index(x, y)
        # print(row)

        # Convert row and col to NumPy arrays
        row = np.array(row)
        col = np.array(col)

        # Filter out-of-range indices
        valid_indices = (0 <= row) & (row < rows_raster) & (0 <= col) & (col < cols_raster)

        # Create a masked array to handle invalid values
        z = np.ma.masked_array(np.full_like(row, NoData_value), mask=~valid_indices)

        # Read pixel values for valid indices
        z[valid_indices] = Raster.read(1)[row[valid_indices], col[valid_indices]]

    # Reshape the array to match the shape of the original grid
    z = np.array(z)  # Convert the list to a NumPy array
    pointData[variable_name] = z

    # Set the CRS of the output shapefile to match the input raster
    pointData.crs = Raster.crs

    # Save the DataFrame as a shape file
    pointData.to_file(out_shp)

    return pointData


# Function to create points
def create_points(in_gdf, j, distance):  ## add geodataframe
    line = in_gdf['geometry'].iloc[j]
    if j < 1:
        print('\tLine length is:\t', line.length)
    num_points = int(line.length / distance) + 1
    points = [line.interpolate(i * distance) for i in range(num_points)]  # distance is roughly input meter.

    return points


# Function of degree conversion
def cartesian2northtrue(car_deg):  ## input Cartesian system degree

    N_deg = 450 - car_deg  # Convert angles to North True
    # Adjust angles to bring them into [0, 360] range
    N_deg[N_deg >= 360] -= 360
    N_deg[N_deg < 0] += 360

    return N_deg


# Function to create points along a LineString
def points_along_lines(in_shp, out_shp, distance):
    # Initialize lists to store x and y coordinates
    id = []  # id of line
    xx = []  # all x coordinates
    yy = []  # all y coordinates
    dd = []  # A distance from origin in each line
    angles_lines = []

    in_gdf = gpd.read_file(in_shp)
    for j in range(in_gdf['geometry'].shape[0]):

        # Create points along the LineString
        points = create_points(in_gdf, j, distance)
        x = []  # x coordinates of jth line
        y = []  # y coordinates of jth line
        for i, p in enumerate(points):
            # Get the x and y coordinates of the point
            id.append(j), x.append(p.x), y.append(p.y), dd.append(distance * i)
            # print("x coordinate:", p.x, "y coordinate:", p.y, "distance:", distance * i)             # Print the coordinates

        # The angles between consecutive points
        x1 = np.array(x)
        y1 = np.array(y)
        xx.extend(x)
        yy.extend(y)

        # Calculate the angles
        angles = np.degrees(np.arctan2(y1[1:] - y1[:-1], x1[1:] - x1[:-1]))
        angles_N = cartesian2northtrue(angles)
        angles_N = np.append(angles_N[0], angles_N)  # copied first column
        angles_lines.extend(angles_N)

    # Create a DataFrame
    df = pd.DataFrame()
    df['id'] = id
    df['x'] = xx
    df['y'] = yy
    df['distance'] = dd
    df['angle'] = angles_lines

    # Create contour line points
    out_csv = out_shp.split(".")[0].split("/")[-1] + '.csv' # Linux
    #out_csv = out_shp.split(".")[0].split("\\")[-1] + '.csv' # Windows

    print(out_csv)
    Output_point_csv = os.path.join(Process_folder, out_csv)

    # Create a GeoDataFrame from the DataFrame
    out_gdf = gpd.GeoDataFrame(df)
    out_gdf.set_geometry(gpd.points_from_xy(out_gdf['x'], out_gdf['y']), inplace=True, crs=in_gdf.crs)

    # Save the DataFrame as a CSV file
    out_gdf.to_file(out_shp)
    out_gdf.to_csv(Output_point_csv, index=False, float_format='%.12f')

    return Output_point

# Function to offset a polyline
def offset_polyline(in_shp, out_shp, dist_offset, join_style):  # join_style: 1 = round, 2 = mitred
    line = gpd.read_file(in_shp)
    offset = line.geometry.iloc[0].parallel_offset(dist_offset, 'right', join_style=join_style)  # join_style round

    # Create a GeoDataFrame with the offset LineString
    offset_gdf = gpd.GeoDataFrame(geometry=[offset])
    offset_gdf["Length"] = offset_gdf.length

    offset_gdf.crs = line.crs  # gdf_crs
    offset_gdf.to_file(out_shp)

    return out_shp

# Function to create raster file using gdal
def create_raster(Rasterdata, Raster_output, data_array, NoData_value, data_type=gdal.GDT_Float32):

    # data_array = Rasterdata.GetRasterBand(1).ReadAsArray()
    prj = Rasterdata.GetProjection()  # Read projection
    band = Rasterdata.GetRasterBand(1)
    bandnum = Rasterdata.RasterCount
    transform = Rasterdata.GetGeoTransform()

    gtiff_driver = gdal.GetDriverByName('GTiff')  # Use GeoTIFF driver
    out_ds = gtiff_driver.Create(Raster_output,  # Create a output file
                                 band.XSize, band.YSize, bandnum, data_type)
    out_ds.SetProjection(prj)
    out_ds.SetGeoTransform(transform)

    dst_band = out_ds.GetRasterBand(1)
    dst_band.WriteArray(data_array)
    out_ds.GetRasterBand(1).SetNoDataValue(NoData_value)  # Exclude nodata value
    out_ds.GetRasterBand(1).ComputeStatistics(
        0)  # Calculate statistics for Raster pyramids (Pyramids can speed up the display of raster data)
    print('\tMaximum value in this domain is', np.max(data_array))

    del out_ds

    return Raster_output

def read_raster(Rasterdata, GA_stype): #GA_ReadOnly (0) or GA_Update(1)
    Read_raster = gdal.Open(Rasterdata, GA_stype)
    if Rasterdata is None:
        sys.exit('\tCould not open {0}.'.format(Rasterdata))

    return Read_raster

def read_vector(Vectordata, GA_stype): #GA_ReadOnly (0) or GA_Update(1)
    Read_vector = ogr.Open(Vectordata, GA_stype)
    if Vectordata is None:
        sys.exit('\tCould not open {0}.'.format(Vectordata))

    return Read_vector

### Step 2 ###########################################################
print("\n#############################################################\n")
print('Step 2: Outputs')
print("\n#############################################################\n")
######################################################################
Outdir = Process_folder
print("\tOutput Directory is ", Outdir)

### 2.1 Read raster(input) file ###
print("Step 2.1: Read raster file")

######################################################################
# Read raster(input) file
######################################################################

Rasterdata = read_raster(Raster_file, GA_ReadOnly)
print("\tReading raster file (georeference, etc)")

# Coordinate system
prj = Rasterdata.GetProjection()  # Read projection
print("\tProjection:", prj)

# Get raster size and band
rows = Rasterdata.RasterYSize  # number of rows
cols = Rasterdata.RasterXSize  # number of columns
bandnum = Rasterdata.RasterCount  # band number
print("\trows=", rows, "\tcols=", cols)
# print("band=", bandnum)

# Get georeference info
transform = Rasterdata.GetGeoTransform()
xOrigin = transform[0]  # Upperleft x
yOrigin = transform[3]  # Upperleft y
pixelWidth = transform[1]  # cell size x
pixelHeight = transform[5]  # cell size y (value is negative)
print("\txOrigin=", xOrigin, "m", "yOrigin=", yOrigin, "m")
print("\tpixelWidth=", pixelWidth, "m", "pixelHeight=", -pixelHeight, "m")  # pixelHeight is always negative

# Read the raster band
band = Rasterdata.GetRasterBand(1)
# Data type of the values
print('\tdata type is', gdal.GetDataTypeName(band.DataType))  # Each raster file has a different data type
# Get band value info
RV = Rasterdata.GetRasterBand(1).ReadAsArray()  # raster values in the band

RV_1D = RV.reshape(-1)  # 1D array is needed to use set function
Input_list_data = np.unique(RV_1D)
nodata = np.setdiff1d(Input_list_data, np.unique(Input_list))  # find not matched data

#######################################################################################################################
print('\n\tnon-Input value list', nodata, '\n')  # Need to check this
#######################################################################################################################

# Fetch metadata for the band
band.GetMetadata()

# Print only selected metadata:
print("\t[ NO DATA VALUE ] = ", band.GetNoDataValue())  # check nodata
print("\t[ MIN ] = ", band.GetMinimum())
print("\t[ MAX ] = ", band.GetMaximum())

if band.GetMinimum() is None or band.GetMaximum() is None:
    band.ComputeStatistics(0)
    print("\t[ MIN ] = ", band.GetMinimum())
    print("\t[ MAX ] = ", band.GetMaximum())

# 2.2 Make Polygons within the domain
print("Step 2.2: Make Polygons within the domain")

# Convert raster to polygon
Draft_Polygon = os.path.join(Outputspace, 'Draft_Polygon.shp')
drv = ogr.GetDriverByName("ESRI Shapefile")
out_ds = drv.CreateDataSource(Draft_Polygon)
prj2 = osr.SpatialReference()
prj2.ImportFromWkt(Rasterdata.GetProjection())
out_layer = out_ds.CreateLayer(Draft_Polygon, srs=prj2)
New_field = ogr.FieldDefn("Value", ogr.OFTInteger)  # add Value column
out_layer.CreateField(New_field)
gdal.FPolygonize(band, None, out_layer, 0, callback=None)
out_ds.SyncToDisk()

del out_layer
del out_ds

# Classify the polygon
Domain_Polygon = os.path.join(Outputspace, 'Domain_Polygon.shp')  # Polygon within Domain
Output_Land_Polygon = os.path.join(Outputspace, 'Polygon_Land.shp')  # Land Polygon
Output_Water_Polygon = os.path.join(Outputspace, 'Polygon_Water.shp')  # Water Polygon

# Sort the GeoDataFrame in descending order based on the Length column
gdf = gpd.read_file(Draft_Polygon)

gdf_domain = gdf[gdf['Value'] != nodata[0]] ## need to check 09/18/23
gdf_water = gdf[gdf['Value'] == Input_Water_Class]
gdf_land = gdf[gdf['Value'] != Input_Water_Class]

# Create a new GeoDataFrame with only the selected row
gdf_domain = gpd.GeoDataFrame(gdf_domain, geometry='geometry', crs=gdf.crs)
gdf_land = gpd.GeoDataFrame(gdf_land, geometry='geometry', crs=gdf.crs)
gdf_water = gpd.GeoDataFrame(gdf_water, geometry='geometry', crs=gdf.crs)
gdf_domain.to_file(Domain_Polygon)
gdf_land.to_file(Output_Land_Polygon)
gdf_water.to_file(Output_Water_Polygon)

# Convert raster data type to Float32
NoData_value = -99999.
Raster_float = os.path.join(Outputspace, 'Rasterdtype2float.tif')
create_raster(Rasterdata, Raster_float, RV, NoData_value, data_type=gdal.GDT_Float32)

# Open Raster float data and Output land and water classification raster map
Rasterfloat = read_raster(Raster_float, GA_ReadOnly)

# Get band value info
AA = Rasterfloat.GetRasterBand(1).ReadAsArray()  # raster values in the band

WL = AA.copy()
WL[WL == 0.0] = NoData_value
WL[(WL == float(Input_Water_Class))] = 0.0
WL[(WL != 0.0) & (WL != NoData_value)] = 1.0

# Output land and water classification raster map
WL_class = os.path.join(Outputspace, 'WL_class.tif')
create_raster(Rasterfloat, WL_class, WL, NoData_value, data_type=gdal.GDT_Float32)

######################################################################
# Determine the quadrant of wave direction
######################################################################

# Determine the quadrant of dominant wave direction
theta = Wave_direction * np.pi / 180.  # convert degree to radian
# print(theta)
# Set direction regime (qauadrants)
if theta > np.pi / 4. and theta <= 3 * np.pi / 4.:
    q = 41  # Quadrants 4 & 1
elif theta > 3 * np.pi / 4. and theta <= 5 * np.pi / 4.:
    q = 34  # Quadrants 3 & 4
elif theta > 5 * np.pi / 4. and theta <= 7 * np.pi / 4.:
    q = 23  # Quadrants 2 & 3
else:
    q = 12  # Quadrants 1 & 2
print("\tQuadrants is ", q)

### Step 3 ###########################################################
print("\n#############################################################\n")
print("Step 3: Make a baseline")
print("\n#############################################################\n")
######################################################################

##################################################################################################################
### 3.1 Baseline delineation #######

# Method 1 Land (1): Water (0) isopleth with a moving average method
# Method 2 1-way scan based on significant wave direction
# Method 3 Baseline delineation based on a manual polyline
##################################################################################################################
if Baseline_delineation_method == 1:
    ##################################################################################################################
    # Method1 Baseline delineation based on isopleth with a moving average method
    ##################################################################################################################
    print("\tBaseline delineation based on isopleth with a moving average method")

    # Conduct moving average resampling (this method requires large memory. The moving windows shouldn't use a large value on your stand-alone PC)
    Moving_average_size = 31  # User freely change the window size
    ResampleRasterMA = os.path.join(Outputspace, 'ResampleMA.tif')  # Resample file name
    Contour = os.path.join(Outputspace, 'Contour.shp')

    # define number of contours and range
    conNum = 1
    conList = [0.25]  # Extract a contour (water/land isopleth)=[0.25,0.5, 0.75 etc]

    def make_slices(data, win_size):
        # """Return a list of slices given a window size.
        # data - two-dimensional array to get slices from
        # win_size - tuple of (rows, columns) for the moving window
        # """
        row_num = data.shape[0] - win_size[0] + 1
        col_num = data.shape[1] - win_size[1] + 1
        slices = []
        for i in range(win_size[0]):
            for j in range(win_size[1]):
                slices.append(data[i:row_num + i, j:col_num + j])
        return slices


    # Get band value info
    WLclass = gdal.Open(WL_class, GA_ReadOnly)
    indata_pre = WLclass.GetRasterBand(1).ReadAsArray()  # raster values in the band
    indata = np.where((indata_pre == NoData_value), np.nan, indata_pre)
    slices = make_slices(indata, (Moving_average_size, Moving_average_size))
    stacked = np.dstack(slices)
    # outdata = np.full(indata.shape, np.nan)
    outdata = np.zeros(indata.shape, np.float32)
    range1 = int((Moving_average_size - 1) / 2)  # don't use a name "range" because range will reuse in for loops

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="Mean of empty slice")
        outdata[range1:-range1, range1:-range1] = np.nanmean(stacked, 2)

    outdata[0:range1, :] = NoData_value
    outdata[-range1:rows, :] = NoData_value
    outdata[:, 0:range1] = NoData_value
    outdata[:, -range1:cols] = NoData_value
    outdata_mod = np.where(np.isnan(outdata), NoData_value, outdata)  # Change np.nun to numeric value

    create_raster(Rasterdata, ResampleRasterMA, outdata_mod, NoData_value, data_type=gdal.GDT_Float32)

    # Extract contour (Land (1): Water (0) isopleth)

    # Extract contour (Land (1): Water (0) isopleth)
    ResampleRasterdata = gdal.Open(ResampleRasterMA, GA_ReadOnly)
    band =ResampleRasterdata.GetRasterBand(1)

    drv = ogr.GetDriverByName("ESRI Shapefile")  # Set up the shapefile driver
    out_ds = drv.CreateDataSource(Contour)  # Create a data source
    out_layer = out_ds.CreateLayer('Contour', prj2, ogr.wkbLineString)

    # Define fields of id and elev
    fieldDef = ogr.FieldDefn("ID", ogr.OFTInteger)
    out_layer.CreateField(fieldDef)
    fieldDef = ogr.FieldDefn("elev", ogr.OFTReal)
    out_layer.CreateField(fieldDef)

    # Write shapefile
    # ContourGenerate(Band srcBand, double contourInterval, double contourBase, int fixedLevelCount, int useNoData, double noDataValue,
    #                Layer dstLayer, int idField, int elevField
    gdal.ContourGenerate(band, 0, 0, conList, 1, NoData_value, out_layer, 0, 1)

    del out_ds,ResampleRasterdata

    # Add the length in the contour lines (you can also use shapely)
    dst_ds = ogr.Open(Contour, 1)  # 0 means read-only. 1 means writeable.
    layer = dst_ds.GetLayer()
    featureCount = layer.GetFeatureCount()
    New_field = ogr.FieldDefn("Length", ogr.OFTReal)
    layer.CreateField(New_field)

    for i, feat in enumerate(layer):
        Line_length = feat.geometry().Length()
        feat.SetField("Length", Line_length)
        layer.SetFeature(feat)
    # print(feat.geometry().Length())
    del layer
    del dst_ds  # Close the Shapefile

    # Sort the GeoDataFrame in descending order based on the Length column
    gdf = gpd.read_file(Contour)
    gdf = gdf.sort_values('Length', ascending=False)

    if not gdf.empty:
        max_length_gdf = gdf.iloc[[0]]  # Select the first row (the maximum Length value)
        # Rest of your code that uses max_length_gdf
    else:
        print("The DataFrame 'gdf' is empty.")
    # max_length_gdf = gdf.iloc[[0]]  # Select the first row (the maximum Length value)

    # Create a new GeoDataFrame with only the selected row
    max_length_gdf = gpd.GeoDataFrame(max_length_gdf, geometry='geometry', crs=gdf.crs)
    Max_length_contour = os.path.join(Outputspace, 'Max_length_contour.shp')
    max_length_gdf.to_file(Max_length_contour)

    # Create contour line points
    Output_point = os.path.join(Outputspace, 'Contour_points.shp')
    points_along_lines(Max_length_contour, Output_point, Input_Distance_Transects)

    gdf = createEW(Output_point)

    gdf = moving_gdf(gdf, moving_value=31, transiton_points=15) # smooth the line

    Smoothed_polyline = os.path.join(Outputspace, 'Smoothed_line_Contour.shp')  # output file name and location
    create_smoothed_polyline(gdf, Smoothed_polyline) # Create a LineString geometry from the smoothed points

    del gdf  # Close file

    ##################################################################################################################
    # Method2 Baseline delineation based on a 1-way scan based on significant wave direction
    ##################################################################################################################
elif Baseline_delineation_method == 2:
    print("\tBaseline delineation based on a 1-way scan based on significant wave direction")

    ### Find interface (shoreline)points ###
    print("\tFind interface (shoreline) points")

    # Determine scanning direction
    if q == 34:
        RV_ad = np.flip(RV, 0)  # bottom(South) to top(North) Quadrants 3 & 4
    elif q == 12:
        RV_ad = RV  # top(North) to bottom(South) Quadrants 1 & 2
    elif q == 41:
        RV_ad = np.flip(RV, 1)  # right(East) to left(West) Quadrants 4 & 1
    else:
        RV_ad = RV  # left(West) to right(East) Quadrants 2 & 3

    ### Output band value (raw) ###
    print("\tOutput raster values in csv file")

    csv_out = os.path.join(Process_folder, 'raster_values_raw.csv')
    with open(csv_out, "w+", newline='') as csvfile:
        csvwriter = csv.writer(csvfile, delimiter=',')
        csvwriter.writerows(RV)
        csvfile.close()

    # Find interface between water and land)
    print("\tFind interface points between water and land")

    if q == 12 or q == 34:  # Quadrants 1 & 2 or # Quadrants 3 & 4
        fp = np.zeros((cols, 2))  # first point (column, row)
        for j in range(cols):
            #    print(j)
            for i in range(rows):
                if RV_ad[i, j] != 0 and RV_ad[i, j] != nodata and RV_ad[i, j] != Input_Water_Class:
                    break
            fp[j, 0] = j
            fp[j, 1] = i
    # print (fp)
    else:  # Quadrants 2 & 3 or # Quadrants 4 & 1
        fp = np.zeros((rows, 2))  # first point (column, row)
        for i in range(rows):
            for j in range(cols):
                if RV_ad[i, j] != 0 and RV_ad[i, j] != nodata and RV_ad[i, j] != Input_Water_Class:
                    break
            fp[i, 0] = i
            fp[i, 1] = j
    # print (fp)

    csv_out = os.path.join(Process_folder, 'fp.csv')
    with open(csv_out, "w+", newline='') as csvfile:
        csvwriter = csv.writer(csvfile, delimiter=',')
        csvwriter.writerows(fp)  # writing data
        csvfile.close()

    ### Compute coordinates of the points ###
    print("\tCompute and save coordinates of the points (x,y)")

    x = np.zeros(cols)
    y = np.zeros(cols)

    if q == 12:

        for j in range(cols):
            x[j] = xOrigin + pixelWidth / 2 + fp[j, 0] * pixelWidth
            y[j] = yOrigin - (fp[j, 1]) * (-pixelHeight)  # pixel height is negative

    elif q == 34:
        for j in range(cols):
            x[j] = xOrigin + pixelWidth / 2 + fp[j, 0] * pixelWidth
            y[j] = yOrigin - (rows - fp[j, 1]) * (-pixelHeight)  # pixel height is negative

    elif q == 23:
        for i in range(rows):
            x[i] = xOrigin + fp[i, 1] * pixelWidth
            y[i] = yOrigin - (-pixelHeight) / 2 - (fp[i, 0]) * (-pixelHeight)  # pixel height is negative

    else:
        for i in range(rows):
            x[i] = xOrigin + (cols - fp[i, 1]) * pixelWidth
            y[i] = yOrigin - (-pixelHeight) / 2 - (fp[i, 0]) * (-pixelHeight)  # pixel height is negative

    csv_out = os.path.join(Process_folder, 'coordinate_W.csv')
    Input_point_W = os.path.join(Outputspace, 'Input_points_W.shp')


    gdf = gpd.GeoDataFrame({'East': x, 'North': y, 'FID': range(len(x))})
    gdf.set_geometry(gpd.points_from_xy(gdf['East'], gdf['North']), inplace=True, crs=gdf_domain.crs)
    gdf.to_file(Input_point_W)

    del gdf, gdf_domain

    # Extract Points within water region
    Output_point = os.path.join(Outputspace, 'Points.shp')

    points = gpd.read_file(Input_point_W)
    polygon = gpd.read_file(Output_Water_Polygon)
    intersect = gpd.sjoin(points, polygon, how='inner', predicate='intersects')

    intersect.crs = points.crs
    intersect.to_file(Output_point)

    gdf = createEW(Output_point)
    gdf = moving_gdf(gdf, moving_value=31, transiton_points=15) # smooth the line

    Smoothed_polyline = os.path.join(Outputspace, 'Smoothed_line_W.shp')  # output file name and location
    create_smoothed_polyline(gdf, Smoothed_polyline) # Create a LineString geometry from the smoothed points

    del gdf, intersect  # Close file

    ##################################################################################################################
    # Method3 Baseline delineation based on a manual input
    ##################################################################################################################
else:
    print("\tBaseline delineation based on a manual polyline")
    # Create manual line points
    Output_point = os.path.join(Outputspace, 'Manual_points.shp')

    points_along_lines(Manual_line, Output_point, Input_Distance_Transects)
    gdf = createEW(Output_point)
    gdf = moving_gdf(gdf, moving_value=31, transiton_points=15) # smooth the line

    Smoothed_polyline = os.path.join(Outputspace, 'Smoothed_line_Manual.shp')  # output file name and location
    create_smoothed_polyline(gdf,Smoothed_polyline) # Create a LineString geometry from the smoothed points

    del gdf # Close file

###########################################################################################
# Offset to offshore direction and Make a baseline
###########################################################################################

Offset_line_pre = os.path.join(Outputspace, 'Offset_Line_pre.shp')

Dist = (30 * ((abs(pixelWidth) + abs(pixelHeight)) / 2))  # 30* average cell size
if q == 23 or q == 34:
    Dist_offset = Dist # right side will be positive
else:
    Dist_offset = -Dist

Offset_line = os.path.join(Outputspace, 'Offset_Line.shp')
offset_polyline(Smoothed_polyline, Offset_line, Dist_offset,1)
print('\tDone offset')

# Sort the GeoDataFrame in descending order based on the Length column
gdf = gpd.read_file(Offset_line)
gdf = gdf.sort_values('Length', ascending=False)
max_length_gdf = gdf.iloc[[0]]  # Select the first row (the maximum Length value)

# Create a new GeoDataFrame with only the selected row
max_length_gdf = gpd.GeoDataFrame(max_length_gdf, geometry='geometry', crs=gdf.crs)
Baseline_pre = os.path.join(Outputspace, 'Baseline_pre.shp')
max_length_gdf.to_file(Baseline_pre)
del gdf, max_length_gdf  # geopandas is still open or blocked to delete

# Create baseline_pre points
Baseline_pre_points=os.path.join(Outputspace,'Baseline_pre_points.shp')

Output_point = points_along_lines(Baseline_pre, Baseline_pre_points, Input_Distance_Transects)
gdf = createEW(Baseline_pre_points)
gdf = moving_gdf(gdf, moving_value=31, transiton_points=15)  # smooth the line

Baseline = os.path.join(Outputspace,'Baseline.shp')   # output file name and location
create_smoothed_polyline(gdf, Baseline)  # Create a LineString geometry from the smoothed points

del gdf  # Close file

print('\tDone baseline')

# #Create baseline points
Baseline_points = os.path.join(Outputspace, 'Baseline_points.shp')
Output_point = points_along_lines(Baseline, Baseline_points, Input_Distance_Transects)

### Step 4 ###########################################################
print("\n#############################################################\n")
print("Step 4: Make transect")
print("\n#############################################################\n")
######################################################################
### 4.1 Compute coordinates of the mid/end points for two transects ###
print("Step 4.1: Compute coordinates of the mid/end points for two transects")

dis = Input_Transect_length + 1  # Avoid python computation errors, add aditional distance e.g. +1 m
NP = int(np.ceil(dis / Input_Spacing))  # number of points on the transects, don't use np becasue we alreaqdy use np as numpy shortcut

dst_ds = read_vector(Baseline_points, 0)  # 0 means read-only. 1 means writeable.
layer = dst_ds.GetLayer()
featureCount = layer.GetFeatureCount()

print("\tNumber of features in %s: %d" % (os.path.basename(Baseline_points), featureCount))

fp = np.zeros((featureCount, 10))  # feature input (column, row)
for i, feat in enumerate(layer):
    # feat = layer.GetFeature(i)
    point = feat.GetGeometryRef()
    x, y = point.GetX(), point.GetY()
    azimuth = feat["angle"]  # Get the angle from the shapefile feat[5] is feat["angle"]
    az1 = (azimuth + 90) % 360  # azimuth1 positive direction
    az2 = (azimuth - 90) % 360  # azimuth2 negative direction

    fp[i, 0] = x
    fp[i, 1] = y
    fp[i, 2] = az1
    fp[i, 3] = az2

mx = np.zeros((featureCount - 1, 1))  # mid point x  (column, row)
my = np.zeros((featureCount - 1, 1))  # mid point y (column, row)
# print((fp[i, 0] + fp[i+1, 0])/2)

for i in range(0, featureCount - 1):
    mx[i, 0] = (fp[i, 0] + fp[i + 1, 0]) / 2
    my[i, 0] = (fp[i, 1] + fp[i + 1, 1]) / 2
    fp[i, 0] = mx[i, 0]
    fp[i, 1] = my[i, 0]

    # Convert to Cartesian coordinate
    if -(fp[i, 2] - 90) < 0:
        az1_Cartesian = 360 - (fp[i, 2] - 90)
    else:
        az1_Cartesian = -(fp[i, 2] - 90)

    fp[i, 4] = az1_Cartesian
    # print(az1_Cartesian)

    dx2 = dis * np.cos(az1_Cartesian * np.pi / 180.)
    x2 = mx[i, 0] + dx2
    dy2 = dis * np.sin(az1_Cartesian * np.pi / 180.)
    y2 = my[i, 0] + dy2

    if -(fp[i, 3] - 90) < 0:
        az2_Cartesian = 360 - (fp[i, 3] - 90)
    else:
        az2_Cartesian = -(fp[i, 3] - 90)

    fp[i, 5] = az2_Cartesian

    dx3 = dis * np.cos(az2_Cartesian * np.pi / 180.)
    x3 = mx[i, 0] + dx3
    dy3 = dis * np.sin(az2_Cartesian * np.pi / 180.)
    y3 = my[i, 0] + dy3

    fp[i, 6] = x2
    fp[i, 7] = y2
    fp[i, 8] = x3
    fp[i, 9] = y3

fp.resize(featureCount - 1, 10)  # delete last row
# print(fp.shape)

# Output coordinates of the points as a csv file
csv_out = os.path.join(Process_folder, 'Azimuth.csv')
with open(csv_out, "w+", newline='') as csvfile:
    csvwriter = csv.writer(csvfile, delimiter=',')
    fields = ["mx", "my", "az1", "az2", "az1_c", "az2_c", "x2", "y2", "x3", "y3"]  # writing field names
    csvwriter.writerow(fields)
    csvwriter.writerows(fp)  # writing data
    csvfile.close()

### 4.2 Make two transects
print("Step 4.2: Make two transects")

df = pd.read_csv(csv_out)  # Read Azimuth.csv using pandas

paths_a = ogr.Geometry(ogr.wkbMultiLineString)  # azimuth1 on the right side facing forward
paths_b = ogr.Geometry(ogr.wkbMultiLineString)  # azimuth2 on the left side facing forward

###### paths_a and paths_b are sensitively affected by Wave_direction_range ##########################################################################
###### Users can switch the logical statemenet depending on the situations ###########################
for i in range(df.shape[0]):
    if abs(Wave_direction - df.iloc[i, 3]) <= Wave_direction_range:  # abs(Wave_direction-df.iloc[i,3]) != None:
        path_a = ogr.Geometry(ogr.wkbLineString)
        path_a.AddPoint(df.iloc[i, 0], df.iloc[i, 1])
        path_a.AddPoint(df.iloc[i, 6], df.iloc[i, 7])
        paths_a.AddGeometry(path_a)

    if abs(Wave_direction - df.iloc[i, 2]) <= Wave_direction_range:  # abs(Wave_direction-df.iloc[i,2]) != None:
        path_b = ogr.Geometry(ogr.wkbLineString)
        path_b.AddPoint(df.iloc[i, 0], df.iloc[i, 1])
        path_b.AddPoint(df.iloc[i, 8], df.iloc[i, 9])
        paths_b.AddGeometry(path_b)
#######################################################################################################################################################

# azimuth1 on the right side facing forward
Transect_a = os.path.join(Outputspace, 'Transect_a.shp')
drv = ogr.GetDriverByName("ESRI Shapefile")  # Set up the shapefile driver
out_ds = drv.CreateDataSource(Transect_a)  # Create the data source

layer = out_ds.CreateLayer("Transect_a", prj2, ogr.wkbLineString)  # Create a layer
idField = ogr.FieldDefn("id", ogr.OFTInteger)  # Add an ID field
layer.CreateField(idField)

for i, Ft in enumerate(paths_a):
    featureDefn = layer.GetLayerDefn()  # Create the feature and set values
    feature = ogr.Feature(featureDefn)
    feature.SetGeometry(Ft)
    feature.SetField("id", i)
    layer.CreateFeature(feature)
    del feature
na_transect = layer.GetFeatureCount()  # Number of transect1
del out_ds  # Save and close the data source

# azimuth2 on the left side facing forward
Transect_b = os.path.join(Outputspace, 'Transect_b.shp')
drv = ogr.GetDriverByName("ESRI Shapefile")  # Set up the shapefile driver
out_ds = drv.CreateDataSource(Transect_b)  # Create the data source

layer = out_ds.CreateLayer("Transect_b", prj2, ogr.wkbLineString)  # Create a layer
idField = ogr.FieldDefn("id", ogr.OFTInteger)  # Add an ID field
layer.CreateField(idField)

for i, xb in enumerate(paths_b):
    # Create the feature and set values
    featureDefn = layer.GetLayerDefn()
    feature = ogr.Feature(featureDefn)
    feature.SetGeometry(xb)
    feature.SetField("id", i)
    layer.CreateFeature(feature)
    del feature
nb_transect = layer.GetFeatureCount()  # Number of transect2
del out_ds  # Save and close the data source

print("\tNumber of transect1 =", na_transect, ",Number of transect2 =", nb_transect)
print("\tNeed to check the number of transect1 or 2 is zero")
print("\tSuccessfully made transects (* ^ ω ^)")

### 4.3 Convert polyline to points and read raster values ###
print("Step 4.3: Convert polyline to points and copy raster values")

###Select correct transect ###
if na_transect >= nb_transect:
    transect = Transect_a
    n_transect = na_transect
else:
    transect = Transect_b
    n_transect = nb_transect

# Create transect points
Transect_points = os.path.join(Outputspace, 'Transect_points.shp')
Output_point = points_along_lines(transect, Transect_points, Input_Spacing)

# Copy the raster values
RV_transect_points = os.path.join(Outputspace, 'RV_transect_points.shp')
gdf = extract_raster_value(Transect_points, Raster_file, RV_transect_points, 'Raster_val', NoData_value)

### Step 5 ###########################################################
print("\n#############################################################\n")
print("Step 5: Calculate wave attenuation")
print("\n#############################################################\n")
######################################################################
### 5.1 Add decay constants ###
print("Step 5.1: Add decay constants")

# Read Raster values and deploy its associated decay constants in RV_transect_points.shp
# Create a dictionary first
marsh_decay_constants = {
    class_val: decay_val
    for class_val, decay_val in zip(Marsh_Class, Decay_Constant_M)
}
# print(marsh_decay_constants)

# Determine  'DecayC' based on class values
gdf['DecayC'] = 9999 # Set initial values to 9999
gdf.loc[gdf['Raster_val'].isin(Marsh_Class), 'DecayC'] = gdf['Raster_val'].map(marsh_decay_constants)
gdf.loc[gdf['Raster_val'] == Input_Water_Class, 'DecayC'] = 0.0

### 5.2 Calculate wave attenuation ###
print("Step 5.2: Calculate wave attenuation")

DC = gdf['DecayC'].values.tolist()
DC = np.array(DC)

# convert to size to calculate the wave attenuation value on each transect
k_const = np.reshape(DC, (n_transect, NP))  # n_transect: number of transect, NP: number of points on the transect

csv_out = os.path.join(Process_folder, 'Decay_Const.csv')
np.savetxt(csv_out, k_const, delimiter=',')

#WT = np.zeros((k_const.shape[0], k_const.shape[1]), 'd')  # output double
WT = np.full_like(k_const, np.nan, dtype=float)
##### probably np.nan should be better. Future modifications ######

for i in range(k_const.shape[0]):
    WT_pre = 1.
    for j in range(k_const.shape[1]):
        k = k_const[i, j]
        if k == 9999:  # 9999 means land region and no water propagate downstream
            break
        else:
            WT_Temp = WT_pre * np.exp(-k * int(float(Input_Spacing)))
            WT[i, j] = WT_Temp
            WT_pre = WT_Temp
            if WT_pre < Input_Wave_Attenuation_Threshold: # no more calculation
                break

csv_out = os.path.join(Process_folder, 'WT.csv')
with open(csv_out, "w+", newline='') as csvfile:
    csvwriter = csv.writer(csvfile, delimiter=',')
    csvwriter.writerows(WT)
    csvfile.close()

# Reshape the size (2D array to 1D)
WT_reshape = WT.flatten()

try:
    gdf['WT'] = WT_reshape
except:
  print("The size does not matched Check the calculation")

gdf.to_file(RV_transect_points) # Overwrite the shapefile

### 5.3 delete unnecessary wave attenuation points (small wt value affects the results.) ###
print("Step 5.3: Delete unnecessary wave attenuation points")

### #1 Filter by Input_Wave_Attenuation_Threshold using geopandas (can skip this part)
gdf_new = gdf[gdf['WT'] > Input_Wave_Attenuation_Threshold] # Create a new geodataframe with only records where WT > no more calculation 0.02

# Save the new shapefile
fn = os.path.join(Outputspace, 'RV_transect_points_threshold.shp')
gdf_new.to_file(fn)

del gdf_new  # geopandas is still open or blocked to delete
del gdf  # Close the Shapefile

### #2 Create a Fishnet ###
spacing = round(float(Input_Spacing) * 3 / 4, 2)

# Still considering to change raster extent or not # future modification Jin 09/25/2023
dst_ds = ogr.Open(fn, 0)
layer = dst_ds.GetLayer()
extent = layer.GetExtent()
print("\tExtent of points id", extent)
del dst_ds

x_width = extent[1] - extent[0]
y_height = extent[3] - extent[2]
transform2 = (extent[0] - spacing / 2, spacing, 0.0, extent[3] + spacing / 2, 0.0, -spacing)

row_num = int(np.ceil(y_height / spacing) + 1)
col_num = int(np.ceil(x_width / spacing) + 1)

AA = np.arange(1, row_num * col_num + 1, dtype=int).reshape((row_num, col_num)) # ID for each cell

Fishnet_file = os.path.join(Outputspace,'Fishnet.tif')           # output filename
gtiff_driver = gdal.GetDriverByName('GTiff')                     # Use GeoTIFF driver
out_ds = gtiff_driver.Create(Fishnet_file,                       # Create a output file
col_num, row_num, bandnum, GDT_UInt32)

out_ds.SetProjection(prj)
out_ds.SetGeoTransform(transform2)                               # determine the position

out_band = out_ds.GetRasterBand(1)
out_band.WriteArray(AA)
out_band = out_ds.GetRasterBand(1).ComputeStatistics(0)          # calculate statistics for Raster pyramids (Pyramids can speed up the display of raster data)
print('maximum id in this fishnet is', np.max(AA))

del out_ds

RV_transect_points_fishnet = os.path.join(Outputspace,'RV_transect_points_fishnet.shp')

# overwrite  RV_transect_points
gdf = extract_raster_value(fn, Fishnet_file, RV_transect_points_fishnet, 'FishnetID_', -99999.)

# Filtering the WT points (can be replaced by geopandas)
gdf2 = gdf.sort_values(by=['FishnetID_','WT'],ascending=False).groupby('FishnetID_').head(1) # only select maximum value
gdf3 = gdf2.dropna(subset=['WT']).sort_values(by=['id','distance'], ascending=True)  # drop NaN value and resorted again.

RV_transect_points_filter = os.path.join(Outputspace, 'RV_transect_points_filtered.shp')  # output file name and location

# Save the DataFrame as a CSV file
gdf3.to_file(RV_transect_points_filter)

### Step 6 ###########################################################
print("\n#############################################################\n")
print("Step 6: Interpolate points to create a raster")
print("\n#############################################################\n")
######################################################################
WT_raster_gdal = os.path.join(Outputspace, 'WT_raster.tif')
print("\tEXtent is: ",extent)

alg_setting = "invdist:power=2.0:smoothing=0.0:radius1=20*Input_Spacing:radius2=20*Input_Spacing:angle=0.0:max_points=12:min_points=2:nodata=0.0"
idw = gdal.Grid(WT_raster_gdal, RV_transect_points_filter, format="GTiff",
                outputBounds=[extent[0], extent[3], extent[0] + spacing * col_num, extent[3] - spacing * row_num],
                outputSRS=prj, width=col_num, height=row_num, outputType=gdal.GDT_Float32, algorithm=alg_setting,
                zfield='WT')
# caution for the assigned output bounds: [ulx, uly, lrx, lry]
del idw

# Extract IDW data only land regions.
Output_raster = os.path.join(Outputspace,'WT_raster_extracted.tif')

Input_raster = gdal.Open(WT_raster_gdal, 0)
dst_ds = read_vector(Output_Land_Polygon, 0)  # 0 means read-only. 1 means writeable.
layer = dst_ds.GetLayer()

out_ds = gdal.Warp(Output_raster, WT_raster_gdal, cutlineDSName=Output_Land_Polygon, cropToCutline=True,
                   dstNodata=np.nan)
out_ds.SetProjection(prj)
band = out_ds.GetRasterBand(1)
band.GetMetadata()
band.SetNoDataValue(NoData_value)  # exclude nodata value
band.ComputeStatistics(0)  # calculate statistics

del out_ds  # Close file
del Input_raster
del Rasterdata  # Close the raster image

# list of files to be deleted (geopandas products cannot delete now... not sure the reasons.)
# file_list = [file.split(".")[0] + '*' for file in [Baseline_points, Baseline_pre, Offset_line_pre, fn3,RV_transect_points_fishnet, RV_transect_points_fishnet, Smoothed_polyline]]
file_list = [file.split(".")[0] + '*' for file in
             [Baseline_points, Baseline_pre_points, Offset_line_pre, RV_transect_points_fishnet,
              RV_transect_points_fishnet, Smoothed_polyline]]

# Delete each file
print("\nDelete files\n")

for pattern in file_list:
    for file in glob.glob(pattern):
        try:
            os.remove(file)
            print(f"{file} has been deleted.")
        except FileNotFoundError:
            print(f"{file} not found.")
        except PermissionError:
            print(f"{file} cannot be deleted due to permission error.")

print('\n#################################################\n')
print("Job Finished ʕ •ᴥ•ʔ")
print("--- %s seconds ---" % (time.time() - startTime))
print('\n#################################################\n')
