# Wave Attenuation Toolbox (WATTE) version 3.1
# Developed by Center for Computation & Technology and Coastal Ecosystem Design Studio at Louisiana State University (LSU).
# WATTE version 1 was originally developed by M. Foster-Martinez, University of New Orleans: Karim Alizad, University of South Carolina: Scott C. Hagen, LSU (https://digitalcommons.lsu.edu/civil_engineering_data/1/)
# WATTE version 2 was developed using Python3 and QGIS in 2023
# Developer: Jin Ikeda, Shu Gao, Christopher E. Kees, and Peter Bacopoulos
# Last modified 04/16/2024
# Deepest thanks: This software is dedicated to in memory of Dr. Scott C. Hagen. We are truly blessed to have been one of your pupils. We will spread your idea further.
#
# WATTE version 3 is an open source-based toolbox using Python 3 (Verified at version Python 3.10).
# This model estimates and maps wave attenuation along marsh coastlines following an exponential decay.
# version 3.1: This version includes set decay constants from a combination of Inun_level and biomass_density classifications
#
# This software is released under the MIT License, see LICENSE.txt.
#
### Step 1 #############################################################################################################
print("\n#############################################################\n")
print("Step 1: Import modules and inputs parameters")
print("\n#############################################################\n")
########################################################################################################################

### 1.1 Import modules ###
print("Step 1.1: Import modules")

from WATTE_functions import *

# import spatial analysis using Gdal
from osgeo import gdal, ogr, osr
from osgeo.gdalconst import *
gdal.AllRegister()  # Register all of drivers
import rasterio

startTime = time.time()

########################################################################################################################
# Working directory and input parameters
########################################################################################################################

# Print the current working directory
print("\tCurrent working directory: {0}".format(os.getcwd()))

# The target Working directory
Workspace = os.getcwd()

# Change the current working directory
os.chdir(Workspace)

### 1.2 Input parameters and file ###
print("Step 1.2: Input parameters and file")

# Dominant wave direction North True
Wave_direction = float(300)  # North is 0 [degree] and clockwise
Wave_direction_range = 179  # Criterion for consider transects
# Note: Wave_direction = float(input("Dominant wave direction in North True: [degree] between 0-360 type here->"))
assert 0 <= Wave_direction < 360, "\tPlease input the values between 0-360 degree ...{{ (>_<) }}\n"

# Input raster data
Inputspace = os.path.join(Workspace, "Input_data")
Raster_file = os.path.join(Inputspace, "Productivity_TCB.tif")

# Make process folder
Process_folder = os.path.join(Workspace, 'Process')
os.makedirs(Process_folder, exist_ok=True)

# Make process folder
Outputspace = os.path.join(Workspace, 'Output_data')
os.makedirs(Outputspace, exist_ok=True)

########################################################################################################################
# Input baseline delineation method
########################################################################################################################
# Type 1, Method 1 Land (1): Water (0) isopleth with a moving average method
# Type 2, Method 2 1-way scan based on significant wave direction
# Type 3, Method 3 Baseline delineation based on a manual polyline also provide the path
# Manual_line= "Z:/CCR_data/ACTIVE/Jin_Projects/WATTE/EMBIR/Test_3/Input_data/Polyline.shp"

# Here polyline is estimated to draw from left side to right side, which may affect offset direction (use E2W function)

Baseline_delineation_method = 2
Manual_line = os.path.join(Inputspace, "Polyline.shp")
assert 1 <= Baseline_delineation_method <= 3, "\tPlease input the values between 1-3 ...{{ (>_<) }}\n"

########################################################################################################################
# Calculate inundation depth when Inundation_flag is True
########################################################################################################################
# import argparse, sys
# def main(argv):
#     parser = argparse.ArgumentParser(description='Reading inundation depth') # Create the parser
#     parser.add_argument("--inundationFile", type=str,default=None, help="Path to inundation file <*.tif>")
#     args = parser.parse_args(argv) # Parse the arguments
#
#     print('Inundation file:', args.inundationFile)
#
#     return args.inundationFile
#
# if __name__ == "__main__":
#     main(sys.argv[1:])
#     print('Find inundation file')

# Swith to True if you want to use inundation file
Inundation_flag = True
Inundation_file = os.path.join(Inputspace, "Inundation_depth_TCB.tif")

if Inundation_flag:
    check_flag_file(Inundation_flag, Inundation_file)
    print("\tInundation file is set to be used")
    marsh_height = 0.25  # Assumption of vegetation (marsh) height [m]
    ndv = -99999  # No data value (ndv) using ADCIRC conversion

    Inun_depth = read_raster_info(Inundation_file,1)
    print ('Inundation depth \t',Inun_depth.min(),Inun_depth.max())
    reference_raster = read_raster(Inundation_file,GA_ReadOnly)

    min_depth = 1e-6  # Minimum depth for wave attenuation calculation [m]
    Inun_depth = np.where(Inun_depth < min_depth, 0, Inun_depth)
    Inun_level = np.where((0 < Inun_depth) & (marsh_height != 0), Inun_depth / marsh_height, ndv)

    # Output inundation level
    Raster_output = os.path.join(Outputspace, 'Inundation_level.tif')
    create_raster(reference_raster, Raster_output, Inun_level, ndv, data_type=gdal.GDT_Float32)
    print ('Inundation depth \t',Inun_depth.max(),'and level\t', Inun_level.max())

########################################################################################################################
# Read marsh classifications
########################################################################################################################

# Classification(s), Water, String
Input_Water_Class = "40"  # Don't use 0 (no data)

# Classification(s), Other, String
Input_Other_Class = "55"

# Classification(s), Marsh, String
Input_Marsh_Class = "16,23,32"
Marsh_Class = [int(x) for x in Input_Marsh_Class.split(',')]

if Inundation_flag == False:
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

no_decay_value = 9999  # No decay value

########################### Messages ###################################################################################
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

if Inundation_flag == False:
    print("\tDecay constant for each classification of marsh is ", end=" ")
    print(*Decay_Constant_M)

# print(Decay_Constant_M)
print("\tTransect length is ", Input_Transect_length, "meter")
print("\tDistance between transects is ", Input_Distance_Transects, "meter")

Input_list = [Input_Water_Class] + [Input_Other_Class] + Marsh_Class
print("\tInput_list：", Input_list)

# print ('Wave decay stopping point is', Input_Wave_decay)
print("\tSpacing of points along transect is ", Input_Spacing, "meter")

### Step 2 #############################################################################################################
print("\n#############################################################\n")
print('Step 2: Outputs')
print("\n#############################################################\n")
########################################################################################################################
Outdir = Process_folder
print("\tOutput Directory is ", Outdir)

### 2.1 Read raster(input) file ###
print("Step 2.1: Read raster file")

########################################################################################################################
# Read raster(input) file
########################################################################################################################

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

########################################################################################################################
# Set decay constant values to a tiff file when Inundation_flag is True

# Set decay constants for marsh from a combination of Inun_level and biomass_density classifications
# The details, please see set_decay_values2tiff function and guidance_decay_constants_table
########################################################################################################################

if Inundation_flag:
    decay_tiff_values = set_decay_values2tiff(Inun_level,RV,Input_Water_Class,Input_Other_Class, guidance_decay_constants_table) #,Marsh_Class

    decay_tiff = os.path.join(Outputspace, 'decay_tiff_values.tif')
    create_raster(Rasterdata, decay_tiff, decay_tiff_values, no_decay_value, data_type=gdal.GDT_Float32)

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
# Set direction regime (quadrants)
if theta > np.pi / 4. and theta <= 3 * np.pi / 4.:
    q = 41  # Quadrants 4 & 1
elif theta > 3 * np.pi / 4. and theta <= 5 * np.pi / 4.:
    q = 34  # Quadrants 3 & 4
elif theta > 5 * np.pi / 4. and theta <= 7 * np.pi / 4.:
    q = 23  # Quadrants 2 & 3
else:
    q = 12  # Quadrants 1 & 2
print("\tQuadrants is ", q)

### Step 3 #############################################################################################################
print("\n#############################################################\n")
print("Step 3: Make a baseline")
print("\n#############################################################\n")

### 3.1 Baseline delineation #######

# Method 1 Land (1): Water (0) isopleth with a moving average method
# Method 2 1-way scan based on significant wave direction
# Method 3 Baseline delineation based on a manual polyline
########################################################################################################################
Moving_average_size = 3  # User freely change the window size
MA_range = int((Moving_average_size - 1) / 2)  # don't use a name "range" because range will reuse in for loops

if Baseline_delineation_method == 1:
    ####################################################################################################################
    # Method1 Baseline delineation based on isopleth with a moving average method
    ####################################################################################################################
    print("\tBaseline delineation based on isopleth with a moving average method")

    # Conduct moving average resampling (this method requires large memory. The moving windows shouldn't use a large value on your stand-alone PC)

    ResampleRasterMA = os.path.join(Outputspace, 'ResampleMA.tif')  # Resample file name
    Contour = os.path.join(Outputspace, 'Contour.shp')

    # define number of contours and range
    conNum = 1
    conList = [0.20]  # Extract a contour (water/land isopleth)=[0.25,0.5, 0.75 etc]

    # Get band value info
    WLclass = gdal.Open(WL_class, GA_ReadOnly)
    indata_pre = WLclass.GetRasterBand(1).ReadAsArray()  # raster values in the band
    indata = np.where((indata_pre == NoData_value), np.nan, indata_pre)
    slices = make_slices(indata, (Moving_average_size, Moving_average_size))
    stacked = np.dstack(slices)
    # outdata = np.full(indata.shape, np.nan)
    outdata = np.zeros(indata.shape, np.float32)


    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="Mean of empty slice")
        outdata[MA_range:-MA_range, MA_range:-MA_range] = np.nanmean(stacked, 2)

    outdata[0:MA_range, :] = NoData_value
    outdata[-MA_range:rows, :] = NoData_value
    outdata[:, 0:MA_range] = NoData_value
    outdata[:, -MA_range:cols] = NoData_value
    outdata_mod = np.where(np.isnan(outdata), NoData_value, outdata)  # Change np.nun to numeric value

    create_raster(Rasterdata, ResampleRasterMA, outdata_mod, NoData_value, data_type=gdal.GDT_Float32)

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
    points_along_lines(Max_length_contour, Output_point, Input_Distance_Transects, Process_folder)

    gdf = sort2polyline(Output_point, Process_folder,q) # modify on 04/15/24
    gdf = moving_gdf(gdf, Moving_average_size, MA_range, Process_folder) # smooth the line moving_gdf(gdf,moving_value,transition points,target_folder)

    Smoothed_polyline = os.path.join(Outputspace, 'Smoothed_line_Contour.shp')  # output file name and location
    create_smoothed_polyline(gdf, Smoothed_polyline) # Create a LineString geometry from the smoothed points

    del gdf  # Close file

    ####################################################################################################################
    # Method2 Baseline delineation based on a 1-way scan based on significant wave direction
    ####################################################################################################################
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

    if q == 12 or q == 34:
        x = np.zeros(cols)
        y = np.zeros(cols)
    else:
        x = np.zeros(rows)  # modify the bug 04/15/24
        y = np.zeros(rows)

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

    gdf = sort2polyline(Output_point, Process_folder,q)
    gdf = moving_gdf(gdf, Moving_average_size, MA_range, Process_folder) # smooth the line moving_gdf(gdf,moving_value,transition points,target_folder)

    Smoothed_polyline = os.path.join(Outputspace, 'Smoothed_line_W.shp')  # output file name and location
    create_smoothed_polyline(gdf, Smoothed_polyline) # Create a LineString geometry from the smoothed points

    del gdf, intersect  # Close file

    ####################################################################################################################
    # Method3 Baseline delineation based on a manual input
    ####################################################################################################################
else:
    print("\tBaseline delineation based on a manual polyline")
    # Create manual line points
    Output_point = os.path.join(Outputspace, 'Manual_points.shp')

    points_along_lines(Manual_line, Output_point, Input_Distance_Transects, Process_folder)
    gdf = sort2polyline(Output_point, Process_folder,q)
    gdf = moving_gdf(gdf, Moving_average_size, MA_range, Process_folder) # smooth the line moving_gdf(gdf,moving_value,transition points,target_folder)

    Smoothed_polyline = os.path.join(Outputspace, 'Smoothed_line_Manual.shp')  # output file name and location
    create_smoothed_polyline(gdf,Smoothed_polyline) # Create a LineString geometry from the smoothed points

    del gdf # Close file

########################################################################################################################
# Offset to offshore direction and Make a baseline
########################################################################################################################

Offset_line_pre = os.path.join(Outputspace, 'Offset_Line_pre.shp')
offset_cell = 3  # Default 30 cells

Dist = (offset_cell * ((abs(pixelWidth) + abs(pixelHeight)) / 2))  # 30* average cell size
if q == 23 or q == 34:
    Dist_offset = Dist # right side will be positive. For positive distance the offset will be at the left side of the input line. For a negative distance it will be at the right side. In general, this function tries to preserve the direction of the input.
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

Output_point = points_along_lines(Baseline_pre, Baseline_pre_points, Input_Distance_Transects, Process_folder)
gdf = sort2polyline(Baseline_pre_points, Process_folder,q)
gdf = moving_gdf(gdf, Moving_average_size, MA_range, Process_folder)  # smooth the line moving_gdf(gdf,moving_value,transition points,target_folder)

Baseline = os.path.join(Outputspace,'Baseline.shp')   # output file name and location
create_smoothed_polyline(gdf, Baseline)  # Create a LineString geometry from the smoothed points

del gdf  # Close file

print('\tDone baseline')

# #Create baseline points
Baseline_points = os.path.join(Outputspace, 'Baseline_points.shp')
Output_point = points_along_lines(Baseline, Baseline_points, Input_Distance_Transects, Process_folder)

### Step 4 #############################################################################################################
print("\n#############################################################\n")
print("Step 4: Make transect")
print("\n#############################################################\n")
########################################################################################################################
### 4.1 Compute coordinates of the mid/end points for two transects ###
print("Step 4.1: Compute coordinates of the mid/end points for two transects")

dis = Input_Transect_length + 1  # Avoid python computation errors, add additional distance e.g. +1 m
NP = int(np.ceil(dis / Input_Spacing))  # number of points on the transects, don't use np because we already use np as numpy shortcut

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
    if 90 <= Wave_direction and Wave_direction <= 270:
        az1 = (azimuth - 90) % 360  # azimuth1 positive (downstream) direction
        az2 = (azimuth + 90) % 360  # azimuth2 negative (upstream) direction
    else:
        az1 = (azimuth + 90) % 360 # azimuth1 positive (downstream) direction
        az2 = (azimuth - 90) % 360 # azimuth2 negative (upstream) direction

    fp[i, 0] = x
    fp[i, 1] = y
    fp[i, 2] = az1
    fp[i, 3] = az2

for i in range(0, featureCount - 1):
    fp[i, 0] = (fp[i, 0] + fp[i + 1, 0]) / 2 # Overwrite to mid point x  (column, row)
    fp[i, 1] = (fp[i, 1] + fp[i + 1, 1]) / 2 # Overwrite to mid point y (column, row)

    # Convert to Cartesian coordinate
    az1_Cartesian = (-(fp[i, 2] - 90)) % 360
    fp[i, 4] = az1_Cartesian
    # print(az1_Cartesian)

    dx2 = dis * np.cos(az1_Cartesian * np.pi / 180.)
    x2 = fp[i, 0] + dx2
    dy2 = dis * np.sin(az1_Cartesian * np.pi / 180.)
    y2 = fp[i, 1] + dy2

    az2_Cartesian = (-(fp[i, 3] - 90)) % 360
    fp[i, 5] = az2_Cartesian

    dx3 = dis * np.cos(az2_Cartesian * np.pi / 180.)
    x3 = fp[i, 0] + dx3
    dy3 = dis * np.sin(az2_Cartesian * np.pi / 180.)
    y3 = fp[i, 1] + dy3

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

###### paths_a and paths_b are sensitively affected by Wave_direction_range ############################################
###### Users can switch the logical statement depending on the situations ##############################################
for i in range(df.shape[0]):

    # Calculate the upper and lower bounds of the range
    W1 = (Wave_direction + Wave_direction_range) % 360
    W2 = (Wave_direction - Wave_direction_range) % 360

    if W2 <= fp[i, 3] <= W1 or (W1 < W2 and (fp[i, 3] >= W2 or fp[i, 3] <= W1)):  # downstream is opposite (add +180 degree or use fp[i,3] to switch)
        path_a = ogr.Geometry(ogr.wkbLineString)
        path_a.AddPoint(df.iloc[i, 0], df.iloc[i, 1])
        path_a.AddPoint(df.iloc[i, 6], df.iloc[i, 7])
        paths_a.AddGeometry(path_a)

    if W2 <= fp[i, 2] <= W1 or (W1 < W2 and (fp[i, 2] >= W2 or fp[i, 2] <= W1)):  # downstream is opposite (add +180 degree)
        path_b = ogr.Geometry(ogr.wkbLineString)
        path_b.AddPoint(df.iloc[i, 0], df.iloc[i, 1])
        path_b.AddPoint(df.iloc[i, 8], df.iloc[i, 9])
        paths_b.AddGeometry(path_b)
########################################################################################################################

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
Output_point = points_along_lines(transect, Transect_points, Input_Spacing, Process_folder)

# Copy the raster values
RV_transect_points = os.path.join(Outputspace, 'RV_transect_points.shp')
gdf = extract_raster_value(Transect_points, Raster_file, RV_transect_points, 'Raster_val', NoData_value)

### Step 5 #############################################################################################################
print("\n#############################################################\n")
print("Step 5: Calculate wave attenuation")
print("\n#############################################################\n")
########################################################################################################################
### 5.1 Add decay constants ###
print("Step 5.1: Add decay constants")

# Read Raster values and deploy its associated decay constants in RV_transect_points.shp
if Inundation_flag:
    # Extract the decay constant from the inundation raster
    indices = gdf['Raster_val'].isin(Marsh_Class)
    extract_point_values(decay_tiff, RV_transect_points,indices,no_decay_value)
    # Open the shapefile
    gdf = gpd.read_file(RV_transect_points)
    gdf.loc[gdf['Raster_val'] == Input_Water_Class, 'DecayC'] = 0.0 # Set the decay constant for water to 0
else:
    # Create a dictionary first
    marsh_decay_constants = {
        class_val: decay_val
        for class_val, decay_val in zip(Marsh_Class, Decay_Constant_M)
    }
    # print(marsh_decay_constants)

    # Determine  'DecayC' based on class values
    gdf['DecayC'] = no_decay_value # Set initial values to 9999
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
        if k == no_decay_value:  # 9999 means land region and no water propagate downstream
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
print("\tExtent is: ",extent)

alg_setting = "invdist:power=2.0:smoothing=0.0:radius1=20*Input_Spacing:radius2=20*Input_Spacing:angle=0.0:max_points=12:min_points=2:nodata=0.0"
idw = gdal.Grid(WT_raster_gdal, RV_transect_points_filter, format="GTiff",
                outputBounds=[extent[0], extent[3], extent[0] + spacing * col_num, extent[3] - spacing * row_num],
                outputSRS=prj, width=col_num, height=row_num, outputType=gdal.GDT_Float32, algorithm=alg_setting,
                zfield='WT')
# caution for the assigned output bounds: [ulx, uly, lrx, lry]
print ("Spatial interpolation is done.")
del idw  # Close file

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
