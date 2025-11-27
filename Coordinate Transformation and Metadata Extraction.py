# Code to transform metadata files for aerial images of Bundesamt für Eich- und Vermessungswesen
# Transform coordinate system
# Merge files into one csv-file

import os
import csv
from pyproj import Transformer

# filepath - where files are at
filepath = "./2006 - Luftbilder"
files = [f for f in os.listdir(filepath) if f.endswith(".txt")]

# get foldername
foldername = os.path.basename(filepath.strip("/\\"))

# filepath for final csv-file
final_csv = os.path.join(filepath, f"02 - metadaten_EPSG-3035 - {foldername}.csv")

# setting of columns
columns = ["image_number", "Y_Lon", "X_Lat", "Altitude"]

# Coordinate transformation setting
transformer = Transformer.from_crs("EPSG:31255", "EPSG:3035", always_xy = True)

# list for all data
data = []

# loop to extract all data and merge it again
for file_name in files:
    path = os.path.join(filepath, file_name)
    with open(path, "r", encoding="utf-8") as f:
        row = f.readline().strip()
        split = row.split()

        if len(split) == 9:
            
            # original cooridinates
            y = float(split[6])
            x = float(split[7])
            
            # transformation
            lon, lat = transformer.transform(y, x)
            
            # extract wanted data
            entry = {
                "image_number": split[0] + ".tif",
                "Y_Lon": lon,
                "X_Lat": lat,
                "Altitude": split[8],
            }
            data.append(entry)
        else:
            print(f"❗ Format error in file: {file_name} → {row}")
            
            
# write csv-file
with open(final_csv, "w", newline="", encoding="utf-8") as csvfile:
    writer = csv.DictWriter(csvfile, fieldnames=columns)
    writer.writeheader()
    writer.writerows(data)
    
print(f"✅ CSV saved at: {final_csv}")
print(data[0])
