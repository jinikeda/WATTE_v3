# import os module
import os, sys
import numpy as np
import pandas as pd
import csv
import warnings
import geopandas as gpd
from shapely.geometry import LineString, MultiLineString
from shapely.geometry import Point
import glob
import time
import xarray as xr
from scipy.spatial import cKDTree as KDTree

# import spatial analysis using Gdal
from osgeo import gdal, ogr, osr
from osgeo.gdalconst import *
gdal.AllRegister()  # Register all of drivers
gdal.UseExceptions()  # Enable exceptions
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

# Function to get the longest length of a LineString
def get_max_length_row(gdf):
    gdf = gdf.sort_values('Length', ascending=False)

    if not gdf.empty:
        max_length_gdf = gdf.iloc[[0]]  # Select the first row (the maximum Length value)
    else:
        print("The DataFrame 'gdf' is empty.")
        return None

    # Create a new GeoDataFrame with only the selected row
    max_length_gdf = gpd.GeoDataFrame(max_length_gdf, geometry='geometry', crs=gdf.crs)

    return max_length_gdf


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

def MultiLineString2LineString(shapefile_path, number_cells,raster_size): # number_cells: 10, cell_size:
    # Open the shapefile
    gdf = gpd.read_file(shapefile_path)
    crs = gdf.crs


    # Iterate over each geometry in the GeoSeries
    for geometry in gdf.geometry:
        # Check if the geometry is a MultiLineString
        if isinstance(geometry, MultiLineString):
            # Iterate over each LineString in the MultiLineString

            x_coords = []    # Initialize empty lists for the x and y coordinates
            y_coords = []
            for line in geometry.geoms:
                # Iterate over the coordinates of the line
                for x, y in line.coords:
                    # Add the coordinates to the respective lists
                    x_coords.append(x)
                    y_coords.append(y)

            # Calculate the Euclidean distance between the points
            distances = [np.sqrt((x_coords[i + 1] - x_coords[i]) ** 2 + (y_coords[i + 1] - y_coords[i]) ** 2) for i
                         in range(len(x_coords) - 1)]
            x_coords = np.array(x_coords)
            y_coords = np.array(y_coords)

            # Sort the distances in descending order and select the top 10
            top_10_distances = np.sort(distances)[-10:]
            print("Top 10 MultiLineString distances:", top_10_distances)  # to split the line
            np.savetxt('top_10_distances.txt', top_10_distances, fmt='%f')

            # Find the indices exceeding the threshold
            threshold = number_cells * raster_size

            # Convert distances to a numpy array
            distances = np.array(distances)
            # Find the indices exceeding the threshold
            indices = np.where(distances > threshold)[0]
            print("Indices exceeding the threshold:", indices)

            # Initialize an empty list for the lines
            lines = []

            # Iterate over the indices
            for num, i in enumerate(indices):
                if num == 0:
                    # Create a LineString from the x and y coordinates
                    line = LineString(zip(x_coords[:i + 1], y_coords[:i + 1]))
                    lines.append(line)
                else:
                    line = LineString(
                        zip(x_coords[indices[num - 1] + 1:i + 1], y_coords[indices[num - 1] + 1:i + 1]))
                    lines.append(line)

            # Create a GeoDataFrame from the LineStrings
            gdf = gpd.GeoDataFrame({"ID": range(len(lines))}, geometry=lines, crs=crs)
            gdf['Length'] = gdf.length

        else:
            # # The geometry is not a MultiLineString, so we can proceed as before
            # for x, y in geometry.coords:
            #     # Add the coordinates to the respective lists
            #     x_coords.append(x)
            #     y_coords.append(y)
            pass

    return gdf


# Function to create a buffer around a polyline
def get_extent(shapefile_path):
    gdf = gpd.read_file(shapefile_path)
    minx, miny, maxx, maxy = gdf.geometry.total_bounds
    return minx, miny, maxx, maxy


def create_buffer(shapefile_path, output_file_path, distance):
    # Read the shapefile
    gdf = gpd.read_file(shapefile_path)

    # Create a buffer around the geometry
    gdf['geometry'] = gdf['geometry'].buffer(distance) # Ensure the geometry is a LineString

    # Save the GeoDataFrame with the buffered geometry to a new shapefile
    gdf.to_file(output_file_path, driver='ESRI Shapefile')


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


def get_epsg_code(prj):
    projection = osr.SpatialReference(wkt=prj) # Read projection here prj = raster.GetProjection()
    EPSG_code = projection.GetAttrValue('AUTHORITY', 1)  # Get the EPSG code
    return EPSG_code


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

    print("\n#############################################################\n")
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

def check_flag_file(flag_file, file_path):
    if flag_file:
        try:
            assert os.path.exists(file_path), "\tPlease input the valid file path ...{{ (>_<) }}\n"
        except:
            print('Flag file is not found')
    else:
        print('Flag file is not set to be used')

def process_local_wave_file(local_wave_file, output_shp_file,buffer_file, ESPG): # 'EPSG:26916'
    with xr.open_dataset(local_wave_file) as ds:
        print ('reading the local wave file')

    print (ds.info())

    grid_x_values = ds['grid_x'].values # longitude
    grid_y_values = ds['grid_y'].values # latitude
    grid_Hs_values = ds['Hs'].values # significant wave height
    grid_Tp_values = ds['Tp'].values.astype('float64')/(1e+9) # peak period
    grid_Dir_values = ds['Dir'].values # wave direction

    grid_Hs_min, grid_Hs_max, grid_Hs_mean, grid_Hs_median, grid_Hs_std = calculate_array_stats(grid_Hs_values)
    grid_Tp_min, grid_Tp_max, grid_Tp_mean, grid_Tp_median,grid_Tp_std = calculate_array_stats(grid_Tp_values)
    grid_Dir_min, grid_Dir_max, grid_Dir_mean, grid_Dir_median, grid_Dir_std = calculate_array_stats(grid_Dir_values)

    print(f"Significant Wave Height: Min={grid_Hs_min}, Max={grid_Hs_max}, Mean={grid_Hs_mean}, Median={grid_Hs_median}, Std={grid_Hs_std}")
    print(f"Peak Period: Min={grid_Tp_min}, Max={grid_Tp_max}, Mean={grid_Tp_mean},Median={grid_Tp_median}, Std={grid_Tp_std}")
    print(f"Wave Direction: Min={grid_Dir_min}, Max={grid_Dir_max}, Mean={grid_Dir_mean},Median={grid_Dir_median}, Std={grid_Dir_std}")

    df = pd.DataFrame({
        'x': grid_x_values.flatten(),
        'y': grid_y_values.flatten(),
        'Hs': grid_Hs_values.flatten(),
        'Tp': grid_Tp_values.flatten(),
        'Dir': grid_Dir_values.flatten()
    })


    # Create a GeoDataFrame with the points
    geometry = [Point(xy) for xy in zip(df['x'], df['y'])]
    gdf = gpd.GeoDataFrame(df, geometry=geometry)
    gdf = gdf.set_crs(ESPG)

    # Perform a spatial join between the points and the polygon
    gdf_buffer = gpd.read_file(buffer_file)
    filter_gdf = gpd.sjoin(gdf, gdf_buffer, predicate='within')
    filter_gdf.drop(['index_right', 'ID'], axis=1, inplace=True, errors='ignore') # The dropped columns may differ depending on the input shapefile
    filter_gdf = filter_gdf.reset_index(drop=True)

    print(f"Output file: {output_shp_file}")
    filter_gdf.to_file(output_shp_file)

    del ds

    return filter_gdf

# internal function of process_local_wave_file
def calculate_array_stats(array):
    min_value = np.min(array)
    max_value = np.max(array)
    mean_value = np.mean(array)
    median_value = np.median(array)
    std_value = np.std(array)

    return min_value, max_value, mean_value, median_value, std_value

def add_wave_climate2transect_points(Transect_points_file,wave_shp_file,out_shp_file, nearest_wave_file, no_decay_value,   scale_factor, leafsize=10):
    gdf = gpd.read_file(Transect_points_file)
    gdf['Hs'] = float(no_decay_value)
    gdf['Dir'] = float(no_decay_value)
    gdf['Tp'] = float(no_decay_value)
    base_index = gdf[gdf['distance'].values == 0].index
    gdf_xy_point = gdf.loc[base_index, ['x', 'y', 'geometry']]

    Wave_base_gdf = get_nearest_variable_index(gdf_xy_point, wave_shp_file, nearest_wave_file, 'd3d', scale_factor, leafsize=10, xy_list=['x', 'y'])

    gdf.loc[base_index, 'Hs'] = Wave_base_gdf['Hs'].values
    gdf.loc[base_index, 'Dir'] = (Wave_base_gdf['Dir'].values + 180) % 360 # use same direction as the wave file (convert)
    gdf.loc[base_index, 'Tp'] = Wave_base_gdf['Tp'].values


    gdf.to_file(out_shp_file)


# Function to pick the index of the nearest variable point
def get_nearest_variable_index(target_gdf, reference_shp, new_shp, str_variable, scale_factor, leafsize=10, xy_list=['x', 'y']):

    # Read the shapefiles
    reference_gdf = gpd.read_file(reference_shp)

    # Use a projected system
    xy_reference = reference_gdf[xy_list].to_numpy()
    xy_reference = xy_reference / scale_factor

    xy_target = target_gdf[xy_list].to_numpy()
    xy_target = xy_target / scale_factor

    dist_tree = KDTree(xy_reference, leafsize=leafsize)
    dist, ix = dist_tree.query(xy_target, 1)  # find the nearest neighbor for each point


    new_gdf = target_gdf.copy()
    new_gdf[str_variable + '_dis'] = dist
    new_gdf[str_variable + '_index'] = ix
    new_gdf['Hs'] = reference_gdf['Hs'].iloc[ix].values
    new_gdf['Dir'] = reference_gdf['Dir'].iloc[ix].values
    new_gdf['Tp'] = reference_gdf['Tp'].iloc[ix].values

    new_gdf.to_file(new_shp)

    return new_gdf

# function for filtering (reliable) wave data
def filter_wave_data(gdf, output_shp_file):
    # Read the shapefile into a GeoDataFrame
    #gdf = gpd.read_file(wave_shp_file)
    print('Before filtering: \t gdf.shape', gdf.shape)

    # Get the indices of the rows where 'Hs' is not 0
    indices = gdf[gdf['Hs'] != 0.0].index

    # Select specific columns
    columns = ['x','y','Hs', 'Dir', 'Tp','geometry'] # This part may need to be modified depending on the input shapefile
    gdf_filter = gdf.loc[indices, columns]

    print('After filtering: \t gdf.shape', gdf_filter.shape)
    gdf_filter.to_file(output_shp_file)

    return gdf_filter


def calculate_wave_refraction(gdf, n_transect, NP, check_diff_file):
    # Convert to list and then to numpy array
    Hs_base = np.array(gdf['Hs'].tolist())
    Hs_dir_base = np.array(gdf['Dir'].tolist())
    Transect_ang = np.array(gdf['angle'].tolist())

    # Reshape to calculate the wave attenuation value on each transect
    Hs_transect = np.reshape(Hs_base, (n_transect, NP))
    Hs_dir_transect = np.reshape(Hs_dir_base, (n_transect, NP))
    Transect_ang = np.reshape(Transect_ang, (n_transect, NP))

    # Calculate the difference between the first column of Hs_dir_transect and Transect_ang for each row
    Hs_transect_base = Hs_transect[:, 0]

    # Calculate the absolute difference
    raw_diff = np.abs(Hs_dir_transect[:, 0] - Transect_ang[:, 0])

    # Use modulo operator to handle the wrap-around at 360 degrees
    difference = np.where(raw_diff > 180, 360 - raw_diff, raw_diff)
    cos_diff = np.cos(difference * np.pi / 180)

    # If the angle difference is greater than 90 degrees, set the wave height to zero
    Hs_modify = np.where(cos_diff < 0, 0, Hs_transect_base * cos_diff)

    # Save the difference to a CSV file
    np.savetxt(check_diff_file, difference, delimiter=',')

    return Hs_modify



## use the class method. remove function as its slower
class Invdisttree:

    """ inverse-distance-weighted interpolation using KDTree:
invdisttree = Invdisttree( X, z )  -- data points, values
interpol = invdisttree( q, nnear=3, eps=0, p=1, weights=None, stat=0 )
    interpolates z from the 3 points nearest each query point q;
    For example, interpol[ a query point q ]
    finds the 3 data points nearest q, at distances d1 d2 d3
    and returns the IDW average of the values z1 z2 z3
        (z1/d1 + z2/d2 + z3/d3)
        / (1/d1 + 1/d2 + 1/d3)
        = .55 z1 + .27 z2 + .18 z3  for distances 1 2 3
    q may be one point, or a batch of points.
    eps: approximate nearest, dist <= (1 + eps) * true nearest
    p: use 1 / distance**p
    weights: optional multipliers for 1 / distance**p, of the same shape as q
    stat: accumulate wsum, wn for average weights
How many nearest neighbors should one take ?
a) start with 8 11 14 .. 28 in 2d 3d 4d .. 10d; see Wendel's formula
b) make 3 runs with nnear= e.g. 6 8 10, and look at the results --
    |interpol 6 - interpol 8| etc., or |f - interpol*| if you have f(q).
    I find that runtimes don't increase much at all with nnear -- ymmv.
p=1, p=2 ?
    p=2 weights nearer points more, farther points less.
    In 2d, the circles around query points have areas ~ distance**2,
    so p=2 is inverse-area weighting. For example,
        (z1/area1 + z2/area2 + z3/area3)
        / (1/area1 + 1/area2 + 1/area3)
        = .74 z1 + .18 z2 + .08 z3  for distances 1 2 3
    Similarly, in 3d, p=3 is inverse-volume weighting.
Scaling:
    if different X coordinates measure different things, Euclidean distance
    can be way off.  For example, if X0 is in the range 0 to 1
    but X1 0 to 1000, the X1 distances will swamp X0;
    rescale the data, i.e. make X0.std() ~= X1.std() .
A nice property of IDW is that it's scale-free around query points:
if I have values z1 z2 z3 from 3 points at distances d1 d2 d3,
the IDW average
    (z1/d1 + z2/d2 + z3/d3)
    / (1/d1 + 1/d2 + 1/d3)
is the same for distances 1 2 3, or 10 20 30 -- only the ratios matter.
In contrast, the commonly-used Gaussian kernel exp( - (distance/h)**2 )
is exceedingly sensitive to distance and to h.
    """

#leafsize: The number of points at which the algorithm switches over to brute-force. Default: 10.

    def __init__(self, X, z, leafsize=10, stat=0): # scale=scaling factor for distances
        assert len(X) == len(z), "len(X) %d != len(z) %d" % (len(X), len(z))
        self.tree = KDTree(X, leafsize=leafsize)  # build the tree
        self.z = z
        self.stat = stat
        self.wn = 0
        self.wsum = None;

    def __call__(self, q, nnear=6, eps=0, p=1, weights=None):
            # nnear nearest neighbours of each query point --
        q = np.asarray(q)
        qdim = q.ndim
        if qdim == 1:
            q = np.array([q])
        if self.wsum is None:
            self.wsum = np.zeros(nnear)

        self.distances, self.ix = self.tree.query( q, k=nnear, eps=eps ) #scipy.spatial.KDTree.query
        interpol = np.zeros((len(self.distances),) + np.shape(self.z[0]))
        weight_factors = []  # list to store weight factors
        jinterpol = 0
        for dist, ix in zip(self.distances, self.ix):
            if nnear == 1:
                weight_factors.append(1.0)
                wz = self.z[ix]
            elif dist[0] < 1e-10:
                wz = self.z[ix[0]]
            else:  # weight z s by 1/dist --
                w = 1 / dist ** p
                if weights is not None:
                    w *= weights[ix]  # >= 0
                w /= np.sum(w)
                weight_factors.append(w)  # store the weight factors for each query point
                wz = np.dot(w, self.z[ix])
                if self.stat:
                    self.wn += 1
                    self.wsum += w
            interpol[jinterpol] = wz
            jinterpol += 1
        #return interpol, dist if qdim > 1 else interpol[0]
        return interpol, weight_factors if qdim > 1 else interpol[0]

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