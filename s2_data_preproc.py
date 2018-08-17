"""IDEAM datacube preprocessing script. Includes:
jp2 conversion to geotiff using old GDAL
Duplicate product cleaner
By John Roberts, jfr10@le.ac.uk"""

from collections import deque
from xml.etree import ElementTree
import subprocess
import gdal, osr
import re
import argparse
import os
import glob
import logging

import pdb

gdal.UseExceptions()

def init_log():
    global log
    log = logging.getLogger(__name__)
    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(logging.DEBUG)
    log.addHandler(logging.StreamHandler())
    log.setLevel(logging.DEBUG)

    
def apply_sen2cor(image_path, L2A_path, delete_unprocessed_image=False):
    """Applies sen2cor to the SAFE file at image_path. Returns the path to the new product. From github.com/clcr/pyeo"""
    # Here be OS magic. Since sen2cor runs in its own process, Python has to spin around and wait
    # for it; since it's doing that, it may as well be logging the output from sen2cor. This
    # approatch can be multithreaded in future to process multiple image (1 per core) but that
    # will take some work to make sure they all finish before the program moves on.
    sen2cor_proc = subprocess.Popen([L2A_path, image_path],
                                    stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                    universal_newlines=True)
    try:
        while True:
            nextline = sen2cor_proc.stdout.readline()
            if nextline == '' and sen2cor_proc.poll() is not None:
                break
            if "CRITICAL" in nextline:
                log.error(nextline)
                break
            log.info(nextline)
    except subprocess.CalledProcessError:
        log.error("Sen2Cor failed")
        raise
    log.info("sen2cor processing finished for {}".format(image_path))
    if delete_unprocessed_image:
        os.rmdir(image_path)
    return image_path.replace("MSIL1C", "MSIL2A")


def atmospheric_correction(image_directory, out_directory, L2A_path, delete_unprocessed_image=False):
    """Applies Sen2cor cloud correction to level 1C images. From github.com/clcr/pyeo"""
    images = [image for image in os.listdir(image_directory)
              if image.startswith('MSIL1C', 4)]
    log.info(images)
    # Opportunity for multithreading here
    for image in images:
        image_path = os.path.join(image_directory, image)
        image_timestamp = get_sen_2_image_timestamp(image)
        if glob.glob(os.path.join(out_directory, r"*_{}_*".format(image_timestamp))):
            log.warning("{} exists. Skipping.".format(image))
            continue
        try:
            l2_path = apply_sen2cor(image_path, L2A_path, delete_unprocessed_image=delete_unprocessed_image)
        except subprocess.CalledProcessError:
            log.error("Atmospheric correction failed for {}. Moving on to next image.".format(image))
            pass
        else:
            l2_name = os.path.basename(l2_path)
            os.rename(l2_path, os.path.join(out_directory, l2_name))    
    
    
def remove_extra_s2_images(safe_file):
    """Removes images stored at a higher resolution than required"""
    log.info("Removing redundant imagery from {}".format(safe_file))
    # Two regexes that match every band except those listed after the ?!
    re_20_m=r"L2A_\S+_B(?!05|06|07|8A|11|12)\S+(.jp2|.tif)"
    re_60_m=r"L2A_\S+_B(?!01|09|10)\S+(.jp2|.tif)"
    
    base_glob = os.path.join(safe_file, "GRANULE", '*', "IMG_DATA")
    glob_20_m = os.path.join(base_glob, "R20m")
    glob_60_m = os.path.join(base_glob, "R60m")
    path_20_m = glob.glob(glob_20_m)[0]
    path_60_m = glob.glob(glob_60_m)[0]
    list_20_m = os.listdir(path_20_m)
    list_60_m = os.listdir(path_60_m)

    
    regex_output_20 = list(re.match(re_20_m, path) for path in list_20_m)
    to_delete_20 = list(os.path.join(path_20_m,match.group(0))
                        for match in regex_output_20 if match)
    
    regex_output_60 = list(re.match(re_60_m, path) for path in list_60_m)
    to_delete_60 = list(os.path.join(path_60_m,match.group(0))
                        for match in regex_output_60 if match)
    
    for tbd in to_delete_20+to_delete_60:
        log.info("Deleting "+tbd)
        os.remove(tbd)

        
def extract_metadata(jp2_path):
    """Extracts the XML metadata in a jp2 geo-image as a bytestring. No GDAL."""
    with open(jp2_path, "rb") as jp2_bytes:
        # This iterates through the jp2 file byte-by-byte, looking 
        # for the start tag FeatureCollection. Once found, it reads
        # bytes into output string until the matching closing tag
        # is found. Finally, it strips out extra characters
        log.info("Extracting metadata from {}".format(jp2_path))
        testqueue = deque(maxlen=24)
        testqueue.append(jp2_bytes.read(1))
        out = b""
        # Look for first xml element
        while b"<gml:FeatureCollection x" not in b"".join(testqueue):
            testqueue.append(jp2_bytes.read(1))
        out += b"".join(testqueue)
        # Look for the bytes that explicitly are NOT xml
        while b'\x00\x00\x00\x00' not in b"".join(testqueue):
            testqueue.append(jp2_bytes.read(1))
            out += testqueue[-1]
        # Removing newline characters and stripping extra bytes at end
        # out = out.replace("\n", "")
        out = out[:out.rfind(b">")+1]
    return out


def get_geotransform_from_metadata(metadata_string):
    """Returns the gdal geotransform, raster size and CRS string from metadata."""
    metadata = ElementTree.fromstring(metadata_string)
    namespaces = {"ns0":r"http://www.opengis.net/gml"}

    origin_point = metadata.find(r".//ns0:origin/ns0:Point", namespaces)
    projection = origin_point.attrib.get("srsName")

    origin_pos = origin_point.find(r".//ns0:pos", namespaces)
    origin_pos_tuple = origin_pos.text.split(" ")
    origin_pos_tuple = tuple(int(coord) for coord in origin_pos_tuple)

    size_node = metadata.find(r".//ns0:GridEnvelope/ns0:high", namespaces)
    size_tuple = size_node.text.split(" ")
    size_tuple = tuple(int(coord) for coord in size_tuple)

    res_nodes = metadata.findall(r".//ns0:RectifiedGrid/ns0:offsetVector", namespaces)
    x_res = int(res_nodes[0].text.split(" ")[0])
    y_res = int(res_nodes[1].text.split(" ")[1])

    # build the geotransform
    gt = (origin_pos_tuple[0], x_res, 0, origin_pos_tuple[1], 0, y_res)
    log.debug("Metadata extracted: {}".format(gt))
    return gt, size_tuple, projection


def add_projection_and_extent_to_tiff(tif_path, metadata_tuple, gdal_home=r"/usr/share/gdal/1.11"):
    """Restores geospatial data to a tiff. Takes output from
    get_geotransform_from_metadata. If not set, sets the gdal
    data path to gdal home"""
    log.info("Adding {} to {}".format(metadata_tuple, tif_path))
    if not os.getenv("GDAL_DATA"):
        os.putenv("GDAL_DATA", gdal_home)
    gt, size, projection = metadata_tuple
    EPSG = int(projection.split("::")[-1])
    log.debug("EPSG number: {}".format(EPSG))
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(EPSG)
    image = gdal.Open(tif_path, gdal.GA_Update)
    image.SetGeoTransform(gt)
    image.SetProjection(srs.ExportToWkt())
    log.info("Projection set to: {}".format(srs.ExportToWkt()))
    image = None
    

def convert_jp2_to_tiff(jp2_path, tif_path):
    """Converts a jp2 to a tif (NOT Geotiff) using Jpeg2000 tools"""
    log.info("Converting {} to {}".format(os.path.basename(jp2_path), os.path.basename(tif_path)))
    subprocess.call(["opj_decompress", "-i", jp2_path, "-o", tif_path])
    return tif_path


def replace_jp2_with_geotiff(jp2_path, remove_old=False):
    """Will replace a jp2 with a geotiff of the same name.
    Can remove old jp2"""
    log.info("Preparing to convert {}".format(jp2_path))
    tif_path = jp2_path[:-4]+".tif"
    metadata_string = extract_metadata(jp2_path)
    metadata = get_geotransform_from_metadata(metadata_string)
    convert_jp2_to_tiff(jp2_path, tif_path)
    add_projection_and_extent_to_tiff(tif_path, metadata)
    if remove_old:
        os.rm(jp2_path)

        
def convert_dir_to_geotiff(jp2_dir, remove_old=False):
    """Converts an entire directory of jp2 to geotiff"""
    log.info("Converting {}".format(jp2_dir))
    images = glob.glob(os.path.join(jp2_dir, "*.jp2"))
    for image in images:
        replace_jp2_with_geotiff(image, remove_old)

        
def convert_safe_to_geotiff(safe_dir, remove_old=False):
    """Converts an entire .SAFE file to geotiff"""
    log.info("Converting {}".format(safe_dir))
    image_glob = os.path.join(safe_dir, "GRANULE", '*', "IMG_DATA", "R*")
    image_dirs = glob.glob(image_glob)
    for image_dir in image_dirs:
        convert_dir_to_geotiff(image_dir, remove_old)



        

if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Remove duplicated data and convert jp2 to geotiff")
    parser.add_argument("safe_path", help="Path to the .SAFE file to preprocess")
    args = parser.parse_args()
    init_log()
    safe_path = args.safe_path
    remove_extra_s2_images(safe_path)
    convert_safe_to_geotiff(safe_path)
    
    