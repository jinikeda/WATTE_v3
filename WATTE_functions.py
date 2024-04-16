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

def N2S(gdf):
    point_start = gdf['North'].iloc[0]
    point_end = gdf['North'].iloc[-1]

    if point_start > point_end:  # The string order is North to South
        pass
    else:
        gdf = gdf.loc[::-1]  # reorder the contents When polyline direction is South to North side

    return gdf

# Function of add East and North coordinates and reorder the direction of the polyline
def sort2polyline(in_shp, target_folder,q): # add quadrant information to modify the direction of the polyline (04/15/2024)
    gdf = gpd.read_file(in_shp)
    centroidseries = gdf['geometry'].centroid
    East, North = [list(t) for t in zip(*map(getXY, centroidseries))]

    gdf['East'] = East
    gdf['North'] = North

    if q == 12 or q == 34:
        gdf = E2W(gdf)
    else:
        gdf = N2S(gdf)

    out_csv = in_shp.split(".")[0].split("/")[-1] + '.csv'  # Linux
    Output_point_csv = os.path.join(target_folder, out_csv)
    gdf.to_csv(Output_point_csv)

    return gdf


# Function of moving average of points
def moving_gdf(gdf,moving_value,transiton_points,target_folder):
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
    Output_point_csv = os.path.join(target_folder, out_csv)
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

def extract_point_values(raster_path, points_path, indices, no_decay_value):
    # Load the shapefile of points
    points = gpd.read_file(points_path)
    # Load the DEM raster
    dem = rasterio.open(raster_path)

    # extract xy from point geometry
    raster_values = []
    array = dem.read(1)

    for index, point in zip(indices, points.geometry):
        if index: # If the point is within the marsh area, extract the z value
            x = point.xy[0][0]
            y = point.xy[1][0]
            row, col = dem.index(x, y)

            # Append the z value to the list of z values
            raster_values.append(array[row, col])
        else:
            # Append no_decay_value to the list of z values
            raster_values.append(no_decay_value)

    points['DecayC'] = raster_values
    points.to_file(points_path, driver='ESRI Shapefile')
    del points

    return raster_values


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
    N_deg = (450 - car_deg) % 360  # Convert angles to North True

    return N_deg


# Function to create points along a LineString
def points_along_lines(in_shp, out_shp, distance, target_folder):
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
    Output_point_csv = os.path.join(target_folder, out_csv)

    # Create a GeoDataFrame from the DataFrame
    out_gdf = gpd.GeoDataFrame(df)
    out_gdf.set_geometry(gpd.points_from_xy(out_gdf['x'], out_gdf['y']), inplace=True, crs=in_gdf.crs)

    # Save the DataFrame as a CSV file
    out_gdf.to_file(out_shp)
    out_gdf.to_csv(Output_point_csv, index=False, float_format='%.12f')

    return out_shp

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

    if np.any(np.isnan(data_array)):
        # Replace None with np.nan in data_array
        data_array = np.where(data_array == None, np.nan, data_array)
        print('\tMaximum value in this domain without nan is', np.nanmax(data_array))
    else:
        print('\tMaximum value in this domain is', np.max(data_array))

    del out_ds

    return Raster_output

def read_raster(Rasterdata, GA_stype): #GA_ReadOnly (0) or GA_Update(1)
    Read_raster = gdal.Open(Rasterdata, GA_stype)
    if Rasterdata is None:
        sys.exit('\tCould not open {0}.'.format(Rasterdata))

    return Read_raster


def read_raster_info(raster_file, band_num, read_info = False): # Future modify the function read transform and projection information March 11, 2024 (Jin)
    raster = gdal.Open(raster_file)

    if raster is None:
        sys.exit('\tCould not open {0}.'.format(raster_file))

    band_values = raster.GetRasterBand(band_num).ReadAsArray()

    if read_info:
        transforms = raster.GetGeoTransform()
        X = raster.RasterXSize        # Read column size
        Y = raster.RasterYSize        # Read row size (negative value)
        prj = raster.GetProjection()  # Read projection
        return transforms, X, Y, prj, band_values

        # transforms[0] # Origin x coordinate
        # transforms[1] # Pixel width
        # transforms[2] # x pixel rotation (0° if image is north up)
        # transforms[3] # Origin y coordinate
        # transforms[4] # y pixel rotation (0° if image is north up)
        # transforms[5] # Pixel height (negative)

    else:
        return band_values


def read_vector(Vectordata, GA_stype): #GA_ReadOnly (0) or GA_Update(1)
    Read_vector = ogr.Open(Vectordata, GA_stype)
    if Vectordata is None:
        sys.exit('\tCould not open {0}.'.format(Vectordata))

    return Read_vector

def degree2radian(degree):
    radian = degree * np.pi / 180
    return radian

def radian2degree(radian):
    degree = radian * 180 / np.pi
    return degree

def make_slices(data, win_size):
    # """Return a list of slices given a window size.
    # data - two-dimensional array to get slices from
    # win_size - tuple of (rows, columns) for the moving window
    # """
    row_num = data.shape[0] - win_size[0] + 1 # Number of rows for the slice
    col_num = data.shape[1] - win_size[1] + 1 # Number of columns for the slice
    slices = []
    for i in range(win_size[0]):
        for j in range(win_size[1]):
            slices.append(data[i:row_num + i, j:col_num + j])
    return slices


# Define the dictionary of marsh decay constants first one is inundation depth and the second one is marsh density
guidance_decay_constants_table = {
    ("low_low"): 0.035,
    ("medium_low"): 0.021,
    ("high_low"): 0.001,
    ("low_medium"): 0.055,
    ("medium_medium"): 0.030,
    ("high_medium"): 0.006,
    ("low_high"): 0.107,
    ("medium_high"): 0.090,
    ("high_high"): 0.015,
    ("water"): 0.0,
    ("land"): 9999
}

def set_decay_values2tiff(Inun_level, bio_density, Input_Water_Class,Input_Other_Class, guidance_decay_constants_table):

    # Still Considering use Marsh_Class instead of 16, 23, 32. Here, 16 is low, 23 is medium, 32 is high

    print("\n#############################################################")
    print("Set dacay constant values to tiff file using inundation level and marsh density")
    print("\n#############################################################")

    # Define the masks
    marsh_mask_low_low = (Inun_level < 1) & (bio_density == 16)
    marsh_mask_medium_low = ((1 <= Inun_level) & (Inun_level < 2)) & (bio_density == 16)
    marsh_mask_high_low = (2 <= Inun_level) & (bio_density == 16)

    marsh_mask_low_medium = (Inun_level < 1) & (bio_density == 23)
    marsh_mask_medium_medium = ((1 <= Inun_level) & (Inun_level < 2)) & (bio_density == 23)
    marsh_mask_high_medium = (2 <= Inun_level) & (bio_density == 23)

    marsh_mask_low_high = (Inun_level < 1) & (bio_density == 32)
    marsh_mask_medium_high = ((1 <= Inun_level) & (Inun_level < 2)) & (bio_density == 32)
    marsh_mask_high_high = (2 <= Inun_level) & (bio_density == 32)

    mask_water = (bio_density == Input_Water_Class)
    mask_other = (bio_density == Input_Other_Class)

    # Retrieve the corresponding values from the dictionary
    key_list =["water", "low_low", "medium_low", "high_low", "low_medium", "medium_medium", "high_medium", "low_high", "medium_high", "high_high", "water"]
    value_list = [guidance_decay_constants_table.get(key) for key in key_list]
    value_other = guidance_decay_constants_table.get("other")

    # Create a list of masks
    mask_list = [mask_water, marsh_mask_low_low, marsh_mask_medium_low, marsh_mask_high_low,
                 marsh_mask_low_medium, marsh_mask_medium_medium, marsh_mask_high_medium,
                 marsh_mask_low_high, marsh_mask_medium_high, marsh_mask_high_high,
                 ]

    # Initialize decay_tiff_values
    decay_tiff_values = np.full(bio_density.shape, value_other, dtype=np.float32)

    # Iterate over the masks and values
    for mask, key, value in zip(mask_list, key_list, value_list):
        # Apply the mask to the decay_tiff_values array and set the corresponding elements to the value
        print(f"Value for {key}:", value)
        decay_tiff_values[mask] = value

    return decay_tiff_values

# def create_raster(file, reference_raster, z_array, dtype, no_data_value,stats_flag=False):
#     # Create the output raster dataset
#     gtiff_driver = gdal.GetDriverByName('GTiff')
#     out_ds = gtiff_driver.Create(file, reference_raster.RasterXSize, reference_raster.RasterYSize, reference_raster.RasterCount, dtype) # dtype is e.g. gdal.GDT_Int32 and gdal.GDT_Float32
#     out_ds.SetProjection(reference_raster.GetProjection())
#     out_ds.SetGeoTransform(reference_raster.GetGeoTransform())
#     dst_band = out_ds.GetRasterBand(1)
#     dst_band.WriteArray(z_array)
#     dst_band.SetNoDataValue(no_data_value)  # Exclude nodata value
#     stats = dst_band.ComputeStatistics(0)
#     min_val, max_val, mean_val, std_dev_val = stats
#     if stats_flag:
#         print(
#         f'Made a raster file. Statistics:\n\tMinimum: {min_val}, Maximum: {max_val}, Mean: {mean_val}, Standard Deviation: {std_dev_val}')
#     else:
#         print('Made a raster file')
#     out_ds = None
#     return