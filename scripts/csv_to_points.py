#!/usr/bin/env python

import pandas
from collections import defaultdict
import tempfile
import os
import decimal

import omero.clients
import omero.cli
import omero
from omero.rtypes import rint, rdouble, rstring, unwrap
from omero_metadata.populate import ParsingContext
from omero.util.metadata_utils import NSBULKANNOTATIONSRAW

"""
This script parses Dots_3Dcoordinates_miFISH_chr2_rep1.csv file, to create
a Point for each row, on one of the 25 Dataset1 images
"""

# project_name = "idr0123-mota-mifish/experimentA"
DATASET_NAMES = [
    "Dataset1_miFISH chr2 Replicate 1",
    "Dataset2_miFISH chr2 Replicate 2",
    "Dataset3_miFISH with one dual-color probe AF594-AF488 and two single-color probes AT542 and AT647N",
    "Dataset4_miFISH with one dual-color probe AT647N-AT542 and two single-color probes AF594 and AF488",
    "Dataset5_iFISH chr1 spotting Replicate 1",
    "Dataset6_iFISH chr2 spotting Replicate 1",
    "Dataset7_iFISH chr2 spotting Replicate 2",
    "Dataset8_iFISH chr10 spotting Replicate 1",
    # "Dataset9_iFISH chr10 spotting Replicate 2"
]

CSV_NAMES = [
    "experimentA/coordinates/Dataset1_Dots_3Dcoordinates_miFISH_chr2_rep1.csv",
    "experimentA/coordinates/Dataset2_Dots_3Dcoordinates_miFISH_chr2_rep2.csv",
    "experimentA/coordinates/Dataset3_Dots_3Dcoordinates_miFISH_dual_colour_probe_AF594-AF488.csv",
    "experimentA/coordinates/Dataset4_Dots_3Dcoordinates_miFISH_dual_colour_probe_AT647N-AT546.csv",
    "experimentA/coordinates/Dataset5_Dots_3Dcoordinates_iFISH_chr1.csv",
    "experimentA/coordinates/Dataset6_Dots_3Dcoordinates_iFISH_chr2_rep1.csv",
    "experimentA/coordinates/Dataset7_Dots_3Dcoordinates_iFISH_chr2_rep2.csv",
    "experimentA/coordinates/Dataset8_Dots_3Dcoordinates_iFISH_chr10_rep1.csv",
    "experimentA/coordinates/Dataset9_Dots_3Dcoordinates_iFISH_chr10_rep2.csv"
]

BED_NAME = "experimentA/BED_files/Dataset%s.bed"

probes = {
    "a488": {"label": "AF488", "color": (255, 0, 0)},
    "tmr": {"label": "AT542", "color": (0, 255, 0)},
    "a594": {"label": "AF594", "color": (255, 255, 0)},
    "Cy5": {"label": "AT647N", "color": (0, 255, 255)},
    "a700": {"label": "AF700", "color": (255, 0, 255)},
    "ir800": {"label": "AF790", "color": (50, 50, 255)}
}

def create_roi(updateService, image, shapes, name=None):
    roi = omero.model.RoiI()
    if name:
        roi.name = rstring(name)
    roi.setImage(image._obj)
    for shape in shapes:
        roi.addShape(shape)
    return updateService.saveAndReturnObject(roi)


def delete_rois(conn, im):
    result = conn.getRoiService().findByImage(im.id, None)
    to_delete = []
    for roi in result.rois:
        to_delete.append(roi.getId().getValue())
    if to_delete:
        print("Deleting existing {} rois".format(len(to_delete)))
        conn.deleteObjects("Roi", to_delete, deleteChildren=True, wait=True)


def rgba_to_int(red, green, blue, alpha=255):
    return int.from_bytes([red, green, blue, alpha], byteorder="big", signed=True)


def get_omero_col_type(dtype):
    """Returns s for string, d for double, l for long/int"""
    if dtype == "int":
        return "l"
    elif dtype == "float":
        return "d"
    return "s"


def populate_metadata(image, file_path, file_name):
    """Parses csv to create OMERO.table"""
    client = image._conn.c
    ctx = ParsingContext(
        client, image._obj, file=file_path, allow_nan=True
    )
    ctx.parse()


def process_image(conn, image, df, bed_file):

    updateService = conn.getUpdateService()

    print(image.name)
    cell_id = image.name.split("[")[1].replace("]", "")
    cell_id = int(cell_id)
    print("Cell ID", cell_id)

    # collect rows by Nuclei ID
    rows_by_roi = defaultdict(list)

    # Process just rows for this Image (Field of View)
    bed_rows = bed_file[(bed_file["FoV"] == cell_id)]
    for index, row in bed_rows.iterrows():
        x = row['x']
        y = row['y']
        z = row['z']
        cell_id = row["FoV"]
        matching_rows = df[(df["File"] == cell_id) & (df["x"] > x-0.01) & (df["x"] < x+0.01) & (df["y"] > y-0.01) & (df["y"] < y+0.01) & (df["z"] > z-0.01) & (df["z"] < z+0.01)]
        if len(matching_rows) == 0:
            # Bed row has no match in Dot file - Don't know 'Label' - use 3 for now.... 
            roi_key = (100 * row['Nucleus_ID']) + 3
            row["Label"] = "_"   # used for ROI name
            rows_by_roi[roi_key].append(row)
            print("NO MATCH!")
            continue

        # Need to duplicate the bed rows
        for r_index, dot_row in matching_rows.iterrows():
            for col in df.columns:
                if col not in ('x', 'y', 'z', 'File', 'Nuclei'):
                    value = dot_row[col]
                    row[col] = value
            roi_key = (100 * row['Nucleus_ID']) + row['Label']       # e.g. 120 (Label is 0, 1 or 2)
            rows_by_roi[roi_key].append(row)

    roi_names = list(rows_by_roi.keys())
    if len(roi_names) == 0:
        print("No Rows found for Image")
        return

    roi_names.sort()
    new_rows = []
    for roi_key in roi_names:
        rows = rows_by_roi[roi_key]
        points = []
        roi_name = None
        for row in rows:
            # create an ROI with single Point for each row
            point = omero.model.PointI()
            # Try switching x and y...
            point.y = rdouble(row['x'])
            point.x = rdouble(row['y'])
            # We don't want Python3 behaviour of rounding .5 to even number - always round up
            point.theZ = rint(int(decimal.Decimal(row['z']).quantize(decimal.Decimal('1'), rounding=decimal.ROUND_HALF_UP)))

            if row['Channel'] in probes:
                probe = probes[row['Channel']]
                point.textValue = rstring(probe["label"])
                point.strokeColor = rint(rgba_to_int(*probe["color"]))
            points.append(point)

            # these are the same for all Points in an ROI
            nuclei_id = row["Nucleus_ID"]
            cluster = row["Label"]  # 0, 1 or 2
            roi_name = f'Nuclei_{nuclei_id}_{cluster}'
        roi = create_roi(updateService, image, points, roi_name)

        # Need to get newly saved shape IDs
        shapes = list(roi.copyShapes())
        print("Nuclei: %s, added %s shapes" % (roi_key, len(shapes)))
        for row, shape in zip(rows, shapes):
            # checks that the order of shapes is same as order of rows
            assert shape.theZ.val == decimal.Decimal(row['z']).quantize(decimal.Decimal('1'), rounding=decimal.ROUND_HALF_UP)
            row["roi"] = roi.id.val
            row["shape"] = shape.id.val
            new_rows.append(row)

    return new_rows


def main(conn):

    for count, dataset_name in enumerate(DATASET_NAMES):
        csv_name = CSV_NAMES[count]
        dataset = conn.getObject("Dataset", attributes={"name": dataset_name})
        if dataset is None:
            continue
        bed_file_name = BED_NAME % (count + 1)
        df = pandas.read_csv(csv_name, delimiter=",")
        bed_file = pandas.read_csv(bed_file_name, delimiter="\t")

        col_types = [get_omero_col_type(t) for t in bed_file.dtypes]
        col_names = list(bed_file.columns)
        # combine column names from bed file and dot file
        for col_name, type in zip(df.columns, df.dtypes):
            # ignore x,y,z duplicates, File (Fov) and Nuclei (Nucleus_ID)
            if col_name not in ('x', 'y', 'z', 'File', 'Nuclei'):
                col_names.append(col_name)
                col_types.append(get_omero_col_type(type))

        # Create output table with extra columns
        df2 = pandas.DataFrame(columns=(["roi", "shape"] + col_names))

        for image in dataset.listChildren():
            rows = process_image(conn, image, df, bed_file)
            print("Created ROWS", len(rows))
            for row in rows:
                df2 = df2.append(row)

        # data has been added to df2, so we can write it...
        csv_name = dataset_name + ".csv"
        csv_path = os.path.join(tempfile.gettempdir(), csv_name)
        print("writing CSV to ", csv_path)
        # Add # header roi, shape, other-col-types...
        with open(csv_path, "w") as csv_out:
            csv_out.write("# header roi,l," + ",".join(col_types) + "\n")

        df2.to_csv(csv_path, mode="a", index=False)

        # Create OMERO.table from csv
        populate_metadata(dataset, csv_path, csv_name)

# Usage:
# cd idr0123-mota-mifish
# python scripts/csv_to_points.py

if __name__ == "__main__":
    with omero.cli.cli_login() as c:
        conn = omero.gateway.BlitzGateway(client_obj=c.get_client())
        main(conn)
        # conn.close()
