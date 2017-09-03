#!/usr/bin/env python

import argparse
import simi
from os.path import splitext
from mamut_xml_templates import *

'''Convert Simi data into MaMuT XML format.

Attention! This script is an example to be used as a template. It will not work
as is on your lineage file. For instance, it assumes that the lineage has two
founder cells named specifically "AB" and "CD" (see around line 120). Each
founder cell will make a MaMuT track. If your lineage begins with four cells,
you need to edit the code below to create four tracks.

A required parameter is the z_calibration. You will know this value once you
convert your image data into HDF5 for the BigDataViewer. The value is found in
the XML file that pairs with the HDF5 file (".h5"). It is the second to last
value in the <affine> tag: 8.0 in this case:

    <ViewRegistration timepoint="2220" setup="0">
      <ViewTransform type="affine">
        <affine>1.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 8.0 0.0</affine>
      </ViewTransform>
    </ViewRegistration>

In most cases the z_calibration value will be 1.0.

Here are a few examples on how to run with different parameters. The simplest:

    python simi2mamut.py --sbc lineage.sbc --sbd lineage.sbd --z_calibration 1.0 --out lineage-simi.xml

If your lineage is too large and you just want to test if it works, you can limit the number of frames:

    python simi2mamut.py --sbc lineage.sbc --sbd lineage.sbd --z_calibration 1.0 --frame_limit 10 --out lineage-simi.xml

In Simi, you can skip time points when marking a cell. When you playback the
movie, the 3D renderer will interpolate the positions between these two points.
MaMuT does not do that. To get a smooth movie without gaps, every timepoint
must have a marked spot. You can use the interpolate parameter to create the
intermediate spots in your lineage:

    python simi2mamut.py --sbc lineage.sbc --sbd lineage.sbd --z_calibration 1.0 --interpolate --out lineage-simi.xml

You can also control how much interpolation you want by using the fraction
parameter where 1.0 means all intermediate spots. Here I just want 75% of the
intermediate spots to be created:

    python simi2mamut.py --sbc lineage.sbc --sbd lineage.sbd --z_calibration 1.0 --interpolate --fraction 0.75 --out lineage-simi.xml

After creating your MaMuT file, make sure that the <ImageData> field is
actually pointing to the correct image data XML file (the one where you got the
z_calibration parameter). You need to edit the filename and folder manually,
but also check that the dimensions match that of your image data (i.e. width,
height, nslices and nframes):

    <ImageData filename="lineage.xml" folder="/media/nelas/Actinotroch/4d" width="1376" height="1040" nslices="40" nframes="908" pixelwidth="1.0" pixelheight="1.0" voxeldepth="1.0" timeinterval="1.0" />
'''

def main():
    '''Convert Simi to MaMuT.'''

    # Parse arguments.
    parser = argparse.ArgumentParser(description='Convert Simi BioCell cell lineage data to MaMuT XML.')
    parser.add_argument('-c', '--sbc', help='.sbc file', required=True)
    parser.add_argument('-d', '--sbd', help='.sbd file', required=True)
    parser.add_argument('-o', '--out', help='output file', required=True)
    parser.add_argument('-z', '--z_calibration', help='calibration for z (check HDF5 XML)', required=True)
    parser.add_argument('-i', '--interpolate', help='interpolate spots', action='store_true')
    parser.add_argument('-f', '--fraction', help='interpolate fraction (all=1.0, half=0.5)', default=1.0)
    parser.add_argument('-l', '--frame_limit', help='limit number of frames', default=None)
    args = parser.parse_args()

    # Parse a Simi BioCell .sbd file.
    s = simi.SimiProject(args.sbc, args.sbd)

    # Declare initial variables.
    output = open(args.out, 'w')
    spot_id = 1
    if args.frame_limit:
        last_frame = int(args.frame_limit)
    else:
        last_frame = s.sbd.last_frame

    # Lists aggregating spots and edges.
    cell_edges = []

    # Build template for spots per frame.
    spots_per_frame = []
    for f in range(0, last_frame + 1):
        spots_per_frame.append([])

    # Get calibration factor.
    calibration = s.sbd.get_calibration_factor()

    # Iterate through cells.
    for key, cell in s.sbd.valid_cells.items():
        # Define a list of edges.
        cell.spot_edges = []

        # Get interpolated spots for every frame (best for MaMuT).
        if args.interpolate:
            all_spots = cell.interpolate_spots(float(args.fraction))
        # Or not.
        else:
            all_spots = cell.spots

        if cell.last_frame > last_frame:
            continue

        # Iterate through cell interpolated spots.
        for spot in all_spots:
            # Define new id variable.
            spot.id = spot_id
            # Define new cell variable.
            spot.cell = key
            # Fix X value to MaMuT (based on CALIBRATION field of .sbc).
            spot.x = spot.x * calibration
            # Fix Y value to MaMuT.
            spot.y = spot.y * calibration
            # Fix Z value to MaMuT (multiply by 10).
            spot.z = spot.z * float(args.z_calibration)
            # Append spot to the list of his frame.
            spots_per_frame[spot.frame].append(spot)
            # Is this the first spot?
            spot_index = all_spots.index(spot)
            # If not, create an edge (first spot is skipped).
            if spot_index != 0:
                # Create an edge using the previous spot as source and current spot as target.
                cell.spot_edges.append(edge_template.format(source_id=spot_id-1, target_id=spot_id))
            # Increment unique spot id.
            spot_id += 1
        # Define cell's source_id == the id of the last spot.
        cell.source_id = all_spots[-1].id
        # Define cell's target_id == the id of the first spot.
        cell.target_id = all_spots[0].id
        # Append cells to generate cell edges.
        cell_edges.append(cell)

    # Begin XML file.
    output.write(begin_template)

    # Begin AllSpots.
    output.write(allspots_template.format(nspots=spot_id))

    # Loop through lists of spots.
    for frame, spots in enumerate(spots_per_frame):
        if spots:
            output.write(inframe_template.format(frame=frame))
            for mamut_spot in spots:
                output.write(spot_template.format(id=mamut_spot.id, name=mamut_spot.cell, frame=mamut_spot.frame, x=mamut_spot.x, y=mamut_spot.y, z=mamut_spot.z))
            output.write(inframe_end_template)
        else:
            output.write(inframe_empty_template.format(frame=frame))

    # End AllSpots.
    output.write(allspots_end_template)

    # Begin AllTracks.
    output.write(alltracks_template)

    # Begin Track.
    output.write(track_template.format(id=1, duration=last_frame, stop=last_frame, nspots=spot_id))

    # Get list of CD descendants.
    cd = s.sbd.cells['CD'].get_descendants()

    # Loop through cells printing cell and spot edges.
    for cell in cell_edges:
        if cell.generic_name in cd.keys():
            if cell.parent:
                try:
                    output.write(edge_template.format(source_id=cell.parent.source_id, target_id=cell.target_id))
                except:
                    pass
            for edge in cell.spot_edges:
                output.write(edge)

    # End Track.
    output.write(track_end_template)

    # Begin Track.
    output.write(track_template.format(id=2, duration=last_frame, stop=last_frame, nspots=spot_id))

    # Get list of AB descendants.
    ab = s.sbd.cells['AB'].get_descendants()

    # Loop through cells printing cell and spot edges.
    for cell in cell_edges:
        if cell.generic_name in ab.keys():
            if cell.parent:
                try:
                    output.write(edge_template.format(source_id=cell.parent.source_id, target_id=cell.target_id))
                except:
                    pass
            for edge in cell.spot_edges:
                output.write(edge)

    # End Track.
    output.write(track_end_template)

    # End AllTracks.
    output.write(alltracks_end_template)

    # Filtered tracks.
    output.write(filteredtracks_template)

    # Get some variables from the .sbc file.
    n_slices = s.sbc.settings['DISC']['LEVELCOUNT']
    file_name = splitext(args.sbc)[0]

    # End XML file.
    output.write(end_template.format(filename=file_name, nslices=n_slices, nframes=last_frame))


if __name__ == '__main__':
    main()
