# WATTE (Wave ATTEnuation Model) ver3.2: Last updates on April 18th, 2024

### § 1. WATTE toolbox
WATTE is a one-dimensional model based on an exponential decay formula for marsh coastlines. The functional basis of WATTE is wave energy dissipation by vegetation: 

$$W_T(x) = e^{-kx}$$

where $W_T$ is the fraction of wave height transmission along the distance $x$ into a marsh coastline, and $k$ is a decay constant that depends on inundation level and biomass density


<p style="text-align: center;"><strong>Table 1.</strong> Decay constants as a function of inundation level and biomass density</p>

|                 |               |      **Inundation level**     |                              |                             |         
|-----------------|---------------|---------------------------|------------------------------|-----------------------------|
|                 |               |      Low ($\alpha^{a}$ < 1)         |      Medium (1 ≤ $\alpha$ < 2)     |      High (2 ≤ $\alpha$ < 4 $^{b}$)     |
| **Biomass density** |     Low       |     0.035                 |     0.021                    |     0.001                   |
|                 |     Medium    |     0.055                 |     0.030                    |     0.006                   |
|                 |     High      |     0.107                 |     0.090                    |     0.015                   |

<p style="text-align: center;">$^{a}$ Ratio of water (inundation) depth to vegetation height, $^{b}$ Upper limit due to measurements</p> 

<p style="text-align: right;"><a href="https://www.sciencedirect.com/science/article/abs/pii/S1364815219304530">Foster-Martinez et al. (2020)</a></p>

<br />


WATTE version 3 is an **open-source toolbox that uses Python 3 only**, enabling it to run on **high-performance computing**.
It estimates and maps wave attenuation along marsh coastlines following an exponential decay. \

ver3.0: Release Fall 2023
ver3.1: Add the command of set decay constants from a combination of Inun_level and biomass_density classifications, Spring 2024
ver3.2: Add the command of local wave evaluation from a Delft3d file, Spring 2024

**Developer: Jin Ikeda, Shu Gao, Christopher E. Kees, and Peter Bacopoulos**


### § 2. Software Requirements
See: WATTE_environment.yml 
#### Require Python packages: ####
* Python (tested, 3.10.12) 
* GDAL (tested, 3.7.2) 
* geopandas (tested, 0.14.0)
* numpy (tested, 1.26.0)
* pandas (tested, 2.1.0)
* rasterio (tested, 1.3.8)
* scipy (tested, 1.11.2

#### Python code ####
&ensp;**WATTE_v3.py** 

* Run the code:
  ##### Virtual environment setting
    * Step 1: Download a Python code and input data. Example, *git clone git@github.com:jinikeda/WATTE_v3.git*    
    * Step 2: *conda env create -f WATTE_environment.yml*
    * Step 3  *conda activate WATTE*

  ##### Running Python code
    * Step 1: Change the path of the working directory in the Python code
    * Step 2: *python WATTE_v3.py*

### § 3. Flowchart
![Flowchart](https://github.com/jinikeda/WATTE_v3/blob/main/Photo/Figure1.png "")


### § 4. Input data

#### Must ####
* Classification map in a coordinate system (land, water, and marshes in raster format)

#### User's optional ####
* Manual shoreline definition in vector format when the user selects a manual baseline delineation
* Inundation depth to calculate inundation level

### § 5. Input parameters in WATTE ver3

#### Must ####
* Raster values (Water, Marshes, and Land or others)
* Baseline delineation method (select between 1-3). The user has three options to provide flexibility for the given domain of interest: 
  * Method 1: Isopleth Water/Land with a moving average window (auto baseline delineation) 
  * Method 2: 1-way scanning based on a dominant wave direction (auto baseline delineation) 
  * Method 3: Manual definition
* Dominant wave direction (0-360 [degrees] clockwise from true North) and evaluate wave angle range (e.g., 90 degrees for each side)

#### User's optional ####
* Wave angle (default = 90 [degrees] from a dominant wave direction, which evaluates whole transects)
* Transect length [m] (default = 1000 [m])
* Distance between transects [m] (default = 10 [m])
* Spacing of points along transect [m] (default = 10 [m])
* Filter the wave attenuation values greater than a threshold (default = 0.02)
* Offset distance for a baseline delineation (default = 30 * raster cell size [m])
* Moving average size for isopleth Water/Land (Isopleth with moving average method requires large memory. The moving windows shouldn't use a large value on your stand-alone PC. default = 31)
* Extract contour for isopleth Water/Land (default = 0.25))

<p style="text-align: center;"><strong>Table 2.</strong> ADCIRC and Hydro-MEM outputs can be used for parameterizing WATTE</p>

|      HydroMEM output              |      WATTE parameterization                 |                                           |
|-----------------------------------|---------------------------------------------|-------------------------------------------|
|                                   |      Generic utility                        |      Specific utility                     |
|     Wave climate (ADCIRC)                  |     Wave direction and angle range          | Baseline definition (one-way scan) and offset    |
|     Elevation and water levels (ADCIRC)    |     Inundation depth (low, medium, high)    |     Selection of k-value                  |
|     Biomass distribution (Hydro-MEM)          |     Biomass density (low, medium, high)     |     Selection of k-value                  |

Here, the user manually determines decay constants from  **Table 1** when the user doesn't have ADCIRC and Hydro-MEM results. But, a classification map (Land, water, and marshes) is a must. \
&ensp; [Information about Hydro-MEM](https://www.arcgis.com/apps/MapJournal/index.html?appid=85242c8a228945f3b943f3ec7f01e035)


### § 6. Major outputs

* 2D mapping of wave attenuation in raster format (WT_raster_extracted.tif)
* Baseline (Baseline.shp)
* Transect points (Transect_points.shp)
* Wave attenuation ratios over the transects (RV_transect_points_filtered.shp)

### § 7. Output examples
#### Baseline delineation ####
![Baseline](https://github.com/jinikeda/WATTE_v3/blob/main/Photo/Figure2.png "")

<p style="text-align: left;"><strong>Figure.1</strong> Baseline delineation: three methods in Grand Bay, MS and AL</p>

Grand Bay has a complex marine-dominated embayment system. There is a general consistency between the three baseline definitions; however, there are also local differences. The based on moving-average isopleth overlooks along the western edge of the protruding landform in the center of the image. The method based on a one-way scan doesn’t capture the local embayment in the center of the image due to the northward scanning terminating before reaching the local embayment. The one-way scan captures the two small islands in the eastern portion of the image but not the back-barrier shoreline because of the terminating northward scan. The two automated methods cannot fully capture such complex geometry. The user can adjust the input parameters or use a manual definition shoreline (baseline) at that time. 

#### Calculate attenuation over transects ####
![Transects](https://github.com/jinikeda/WATTE_v3/blob/main/Photo/Figure3.png "")

<p style="text-align: left;"><strong>Figure.2</strong> Example demonstration of attenuation calculation over transect points</p>

After obtaining the baseline, transects normal to baseline extend inshore every 10 m along the baseline. In this case, the baseline is offset to 25 m offshore, as an example. The underlying raster contains five bands: low, medium, and high biomass density, water; and upland. Transect points are assigned k-values (Table 1) based on the underlying raster. The coloration of the transect points corresponds to the calculated wave height transmission, which ranges from one (i.e., full wave height) at the baseline to zero (or near zero) (i.e., fully dissipated wave) over the marsh landscape or at the upland.

#### Map attenuation ####
![Map](https://github.com/jinikeda/WATTE_v3/blob/main/Photo/Figure4.png "")

<p style="text-align: left;"><strong>Figure.3</strong> $W_T$ at the isopleth W/L = 0.25 with a moving average window = 31 cells. Here, the available wave angle is set at 179.99 so that evaluate whole wave directions in this case</p>



This software is released under the MIT License; see LICENSE.txt.

WATTE version 1 was developed by M. Foster-Martinez, Karim Alizad, and Scott C. Hagen (2020). WATTE version 1 needs an **ArcGIS** license (refer to https://coastalscience.noaa.gov/news/gis-toolbox-for-estimating-wave-attenuation-by-coastal-marshes-developed-by-noaa-funded-study/) \
<br />&ensp; [Download WATTE ver1](https://digitalcommons.lsu.edu/civil_engineering_data/1/)

WATTE version 2 is an **open source-based** toolbox using Python 3 and QGIS (Verified at version **Python 3.10** and **QGIS 3.22**, https://github.com/jinikeda/WATTE_v2).


