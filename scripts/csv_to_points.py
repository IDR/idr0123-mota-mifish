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
dataset_name = "Dataset1_miFISH chr2 Replicate 1"
CSV_NAME = "experimentA/coordinates/Dots_3Dcoordinates_miFISH_chr2_rep1.csv"


def create_roi(updateService, image, shapes):
    roi = omero.model.RoiI()
    # Give ROI a name if first shape has one
    name = unwrap(shapes[0].textValue)
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


def process_image(conn, image, df):


    updateService = conn.getUpdateService()

    col_types = [get_omero_col_type(t) for t in df.dtypes]
    col_names = list(df.columns)

    # Create output table with extra columns
    df2 = pandas.DataFrame(columns=(["roi", "shape"] + col_names))

    # rows_by_chr = defaultdict(list)
    # max_chr = 0

    print(image.name)
    cell_id = image.name.split("[")[1].replace("]", "")
    cell_id = int(cell_id)
    print("Cell ID", cell_id)

    for index, row in df.loc[df["File"] == cell_id].iterrows():
        # create an ROI with single Point for each row
        point = omero.model.PointI()
        # Try switching x and y...
        point.y = rdouble(row['x'])
        point.x = rdouble(row['y'])
        # We don't want Python3 behaviour of rounding .5 to even number - always round up
        point.theZ = rint(int(decimal.Decimal(row['z']).quantize(decimal.Decimal('1'), rounding=decimal.ROUND_HALF_UP)))

        # TODO: better textValue
        point.textValue = rstring(str(row['Channel']))
        roi = create_roi(updateService, image, [point])

        # Need to get newly saved shape IDs
        shapes = list(roi.copyShapes())
        print("Row: %s, added shape ID: %s" % (index, shapes[0].id.val))
        row["roi"] = roi.id.val
        row["shape"] = shapes[0].id.val
        df2 = df2.append(row)

    csv_name = dataset_name + ".csv"
    csv_path = os.path.join(tempfile.gettempdir(), csv_name)
    # Add # header roi, shape, other-col-types...
    with open(csv_path, "w") as csv_out:
        csv_out.write("# header roi,l," + ",".join(col_types) + "\n")

    df2.to_csv(csv_path, mode="a", index=False)

    # Create OMERO.table from csv
    populate_metadata(image, csv_path, csv_name)


def main(conn):

    dataset = conn.getObject("Dataset", attributes={"name": dataset_name})
    df = pandas.read_csv(CSV_NAME, delimiter=",")

    for image in dataset.listChildren():
        process_image(conn, image, df)

# Usage:
# cd idr0123-mota-mifish
# python scripts/csv_to_points.py

if __name__ == "__main__":
    with omero.cli.cli_login() as c:
        conn = omero.gateway.BlitzGateway(client_obj=c.get_client())
        main(conn)
        conn.close()
