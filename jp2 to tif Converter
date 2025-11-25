import os
import processing
from qgis.core import QgsRasterLayer, QgsProject, QgsRasterFileWriter
from qgis.analysis import QgsNativeAlgorithms
QgsApplication.processingRegistry().addProvider(QgsNativeAlgorithms())

input_folder = "C:/Users/Martin/Documents/UNI - Master USW_Umweltmonitoring/Masterarbeit/Daten/BEV/2006 - Luftbilder"
output_folder = "C:/Users/Martin/Documents/UNI - Master USW_Umweltmonitoring/Masterarbeit/Daten/BEV/2006 - Luftbilder - Tiff"

for filename in os.listdir(input_folder):
    if filename.endswith(".jp2"):
        jp2_path = os.path.join(input_folder, filename)
        tif_name = os.path.splitext(filename)[0] + ".tif"
        tif_path = os.path.join(output_folder, tif_name)

        layer = QgsRasterLayer(jp2_path, filename)
        if not layer.isValid():
            print(f"Fehler beim Laden: {filename}")
            continue

        #Export with Processing Tool (GDAL Translate)
        params = {
        'INPUT': jp2_path,
        'TARGET_CRS': layer.crs().authid(),
        'OUTPUT': tif_path
        }

        processing.run("gdal:translate", params)
        print(f"Exported to: {tif_path}")
