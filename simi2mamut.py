#!/usr/bin/env python

import argparse
import simi
from os.path import splitext
from mamut_xml_templates import *

# Convert Simi data into MaMuT XML format.

# Example:
# ./simi2mamut.py --sbc wt1.sbc --sbd wt1.sbd --out wt1-inter.xml --interpolate

# TODO:
#   - Set voxel size from wt2.xml

def main():
    '''Convert Simi to MaMuT.'''

    # Parse arguments.
    parser = argparse.ArgumentParser(description='Convert Simi BioCell cell lineage data to MaMuT XML.')
    parser.add_argument('-c', '--sbc', help='.sbc file', required=True)
    parser.add_argument('-d', '--sbd', help='.sbd file', required=True)
    parser.add_argument('-o', '--out', help='output file', required=True)
    parser.add_argument('-i', '--interpolate', help='interpolate spots', action='store_true')
    args = parser.parse_args()

    # Parse a Simi BioCell .sbd file.
    s = simi.SimiProject(args.sbc, args.sbd)

    # Declare initial variables.
    output = open(args.out, 'w')
    spot_id = 1
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
            all_spots = cell.interpolate_spots()
        # Or not.
        else:
            all_spots = cell.spots

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
            spot.z = spot.z * 10.0
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
