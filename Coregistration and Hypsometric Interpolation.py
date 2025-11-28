import xdem
from rasterio.enums import Resampling
from pathlib import Path
import geoutils as gu
import numpy as np
import matplotlib.pyplot as plt

# functions
def clean_array(r):
    a = r.data
    return a.filled(np.nan) if np.ma.isMaskedArray(a) else a

def nmad(x):
    x = x[np.isfinite(x)]
    if x.size == 0:
        return np.nan
    med = np.median(x)
    return 1.4826 * np.median(np.abs(x - med))

# Import DEMs
dem_ref = xdem.DEM("BEV/24 - ALS_DTM_5m_begrenzt - v2.tif")
#__________ change the year to import a different DEM
#dem_tba = xdem.DEM("BEV/Output - v2/2006 - DEM_004 - 5m - loch.tif")
dem_tba = xdem.DEM("Daten_Grossglockner/GROSSGLOCKNER PLEIADES STEREO PRIMARY/2022 - Pleiades-DEM - 5m_begrenzt - v2.tif")


# Align dem_tba to dem_ref
dem_tba = dem_tba.reproject(ref = dem_ref, resampling = Resampling.bilinear)


# Import all glacier polygons into list for hypsometric interpolation
#__________change glacier inventory to fit year of DEM
GI_folder = Path("Daten_Grossglockner/Glocknergruppe/GI_3")
glacier_files = sorted(GI_folder.glob("*.gpkg"))
glacier_polygons = [gu.Vector(f) for f in glacier_files]

# Import all outlines at once
#__________ change glacier inventory
glacier_outlines = gu.Vector("Daten_Grossglockner/Glocknergruppe/GI_3.shp")
lake_outlines = gu.Vector("BEV/Sonstiges/Seen Maske.shp")
cloud_outlines = gu.Vector("Daten_Grossglockner/GROSSGLOCKNER PLEIADES STEREO PRIMARY/Wolken Maske Pleiades.shp")

# create cleared arrays for statistical calculations
arr_dem_ref = clean_array(dem_ref)
arr_dem_tba = clean_array(dem_tba)

# Clip glacier outlines to extent of dem_ref
glacier_outlines = glacier_outlines.crop(dem_ref, clip = True)
lake_outlines = lake_outlines.crop(dem_ref, clip = True)
cloud_outlines = cloud_outlines.crop(dem_ref, clip = True)

# Create a mask from glacier outlines
glacier_mask = glacier_outlines.create_mask(dem_ref).data.astype(bool)
lake_mask = lake_outlines.create_mask(dem_ref).data.astype(bool)
cloud_mask = cloud_outlines.create_mask(dem_ref).data.astype(bool)

# extract the slope and create a mask
slope = xdem.terrain.slope(dem_ref)
slope_mask = slope > 40

# create a mask of all usable datapoints 
dh_before = arr_dem_ref - arr_dem_tba   # create DEM difference 
hard_p_outliers = np.abs(dh_before) > 100   # find hard outliers above +100 m elevation change
hard_n_outliers = np.abs(dh_before) < 100   # find hard outliers above -100 m elevation change
both_valid = np.isfinite(arr_dem_ref) & np.isfinite(arr_dem_tba)   # find values that are not NaN and not infinite in both DEMs
outlier_mask = hard_p_outliers & hard_n_outliers & both_valid   # create a mask of outliers

# create a mask of invalid areas consisting of glacier, lakes, outliers, high degree slopes, and clouds
#__________ add cloud mask, if working with pleiades data with clouds
invalid_mask = glacier_mask | lake_mask | outlier_mask | slope_mask #| cloud_mask

# create a mask of valid datapoints without the invalid mask
inlier_mask = np.asarray(both_valid & (~invalid_mask).data, dtype = bool)

# Print statistical status before the coregistration
dhb_in = dh_before[inlier_mask]
print("Before coregistration:   mean = {:.6f} m | median = {:.6f} m | std = {:.6f} m | NMAD = {:.6f} m | n = {:,}"
      .format(np.nanmean(dhb_in), np.nanmedian(dhb_in), np.nanstd(dhb_in), nmad(dhb_in), dhb_in.size))

# Coregistration: 3D shift + 2nd-order polynomial
coreg = xdem.coreg.NuthKaab() + xdem.coreg.Deramp(poly_order = 2)
coreg.fit(dem_ref, dem_tba, inlier_mask = inlier_mask)

# apply the changes
dem_aligned = coreg.apply(dem_tba)

# create DEM difference
arr_dem_a = clean_array(dem_aligned)
dh_after = arr_dem_ref - arr_dem_a

# Print statistical status after the coregistration
dha_in = dh_after[inlier_mask]
print("After coregistration:  mean = {:.3f} m | median = {:.3f} m | std = {:.3f} m | NMAD = {:.3f} m | n = {:,}"
      .format(np.nanmean(dha_in), np.nanmedian(dha_in), np.nanstd(dha_in), nmad(dha_in), dha_in.size))

#__________ activate cloud mask, if working with Pleiades data
#dem_aligned.data[cloud_mask] = np.nan
dem_aligned.data[lake_mask] = np.nan

#__________ change the year before saving DEM without interpolation
dem_aligned.save("BEV/Output - v2/coreg_hyps/2006_DEM_5m - coreg - v1.tif")

# create DEM difference with xdem
#__________ change dates for respective data
ddem = xdem.dDEM(
    raster = dem_ref - dem_aligned,
    start_time=np.datetime64("2006-09-21"),   # 1969-10-12   # 1991-09-04   # 2006-09-21   # 2022-08-04
    end_time=np.datetime64("2022-09-22")
)

combined = np.zeros_like(ddem.data.squeeze(), dtype = "float32")

# perform regional hypsometric interpolation for each glacier polygon, as local hypsometric interpolation is not working
for outline in glacier_polygons:
    print(outline["Gletschern"])
    interpolated = ddem.interpolate(
        method = "regional_hypsometric",
        reference_elevation = dem_ref,
        mask = outline
    )
    combined = np.where(~np.isnan(interpolated), interpolated, combined)
    print("END")

# create DEM from previous interpolation
ddem_interpolated = xdem.DEM.from_array(
    combined,
    transform = ddem.transform,
    crs = ddem.crs
)

# only for Pleiades!!!

# add study area that was cropped out of original Pleiades data
#include_mask = gu.Vector("BEV/Sonstiges/Inklusionsmaske.shp")
#include_mask = include_mask.create_mask(dem_ref).data.astype(bool)

# set all data below and above 500 m elevation change to NaN
#ddem_interpolated.data[(ddem_interpolated.data > 500) | (ddem_interpolated.data < -500)] = np.nan
# set datapoints in additional area to 0 (when reconstructing DEM, zero values get values from reference DEM)
#ddem_interpolated.data[include_mask] = 0

# Plot before and after
f, ax = plt.subplots(1, 2, figsize = (18, 18))
ax[0].set_title("before")
ddem.plot(cmap='RdBu', vmin=-20, vmax=20, ax=ax[0])
glacier_outlines.plot(ddem, fc='none', ec='k', lw=0.5)
ax[1].set_title(" after regional hypsometric gap-filling for each glacier")
ddem_interpolated.plot(cmap='RdBu', vmin=-20, vmax=20, ax=ax[1], cbar_title="Elevation differences (m)")
glacier_outlines.plot(ddem_interpolated, fc='none', ec='k', lw=0.5)
_ = ax[1].set_yticklabels([])

# Beaware of overwriting already saved DEMs!!!
# DEM reconstruction
dem_final = dem_ref - ddem_interpolated

#__________ change the year before saving DEM with interpolation
dem_final.save("BEV/Output - v2/coreg_hyps/2006_DEM_5m - coreg_hyps - v1.tif")
