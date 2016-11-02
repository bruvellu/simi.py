#!/usr/bin/env python

'''Read and access cell lineage data from Simi BioCell.

This is a small library to parse SBD files from Simi*BioCell. It collects
relevant data such as cell name, time point and coordinates for all the cell
lineage data points. The goal is to provide an easy interface to access Simi
data and allow better interoperability with more recent software for cell
lineage data analyses.

Usage:

    # Load library.
    import simi

    # Parse a Simi BioCell .sbd file.
    s = simi.Sbd('lineage.sbd')

    # Access all cells.
    s.cells
    [<cell1>, <cell2>, <cell3>]

    # Access cell data.
    cell = s.cells[0]
    cell.generic_name
    '1b1'

    cell.frame
    198

    cell.spots
    [<spot1>, <spot2>, <spot3>]

    # Print out cell data.
    cell.print_data()

    Name: 1b1 (1b1)
    Frame: 198
    Spots: 7
            frame	x	y	z
            198	390	323	16
            213	407	332	19
            228	403	332	16
            243	394	326	14
            252	394	332	15
            258	394	331	13
            273	393	333	14

'''

class Sbd:
    '''Database file for Simi BioCell.'''
    def __init__(self, sbd_file):
        self.sbd_file = self.open_sbd(sbd_file)
        self.cells = {}
        self.last_frame = 0

        # Parse file.
        self.parse_sbd()

    def open_sbd(self, filepath):
        '''Get user input, save abspath and try to read file.'''
        try:
            sbd_file = open(filepath)
            print('Loaded "{0}"'.format(sbd_file.name))
            return sbd_file
        except:
            print('Could not load the file "{0}"'.format(filepath))

    def parse_sbd(self):
        '''Parse SBD information.'''
        print('Parsing "{0}"'.format(self.sbd_file.name))
        # Temporary buffer for each line of a cell record.
        tmp_cell = ''
        # Read SBD file line by line.
        for line in self.sbd_file.readlines():
            if line.startswith('---'):
                if tmp_cell:
                    # Create cell instance with raw data.
                    cell = Cell(tmp_cell)
                    # If cell is valid, add it to list.
                    if cell.valid:
                        self.add_cell(cell)
                    # Clean temporary cell.
                    tmp_cell = ''
            else:
                # Add next line to cell data.
                tmp_cell = tmp_cell + line

    def add_cell(self, cell):
        '''Add new cell to list of cells.'''
        self.cells[cell.generic_name] = cell
        self.update_last_frame(cell)

    def update_last_frame(self, cell):
        '''Updates the last_frame value for the simi instance.'''
        if cell.last_frame > self.last_frame:
            self.last_frame = cell.last_frame

    def write_matrix(self, outfile):
        '''Output flat matrix files with all cells.'''
        matrix = open(outfile, 'w')
        # Write header.
        matrix.write('{0},{1},{2},{3},{4},{5}\n'.format('generic_name', 'frame', 'mitosis', 'color', 'n_spots', 'custom_name'))
        for cell_id, cell in self.cells.items():
            matrix.write('{0},{1},{2},{3},{4},{5}\n'.format(cell.generic_name, cell.frame, cell.mitosis, cell.color, cell.n_spots, cell.custom_name))
        matrix.close()


class Cell:
    '''Store all cell-related information'''
    def __init__(self, raw_data):
        self.raw_data = raw_data
        self.valid = False

        # Cell attributes.
        self.generic_name = ''
        self.mitosis = 0
        self.frame = 0
        self.last_frame = 0
        self.color = 0
        self.custom_name = ''
        self.n_spots = 0
        self.spots = []

        # TODO: Each cell has a parent and a sister. Extract such information
        # from the SBD file and add to cell attributes.

        # Parse data, any error returns False (invalid).
        self.valid = self.parse_data()

    def parse_data(self):
        '''Extract attributes from raw data.'''
        split_lines = self.raw_data.split('\n')

        # Generic name.
        try:
            line_one = split_lines[0]
            self.generic_name = line_one.split()[4]
        except:
            return False

        # Mitosis ID.
        try:
            line_two = split_lines[1]
            self.mitosis = int(line_two.split()[0])
        except:
            return False

        # Frame, color and custom name.
        try:
            line_three = split_lines[2]
            self.frame = int(line_three.split()[0])
            self.color = int(line_three.split()[5])
            self.custom_name = line_three.split()[6]
        except:
            return False

        # Number of spots.
        try:
            line_four = split_lines[3]
            self.n_spots = int(line_four.split()[0])
        except:
            return False

        # If there are no spots, invalid cell.
        if self.n_spots == 0:
            return False
        else:
            # Position of last spot in the list.
            last_spot_index = self.n_spots + 4
            # Add spots one by one.
            for spot_line in split_lines[4:last_spot_index]:
                new_spot = Spot(spot_line)
                self.spots.append(new_spot)
                self.update_last_frame(new_spot)
            # Change status to valid.
            return True

    def update_last_frame(self, spot):
        '''Updates the last_frame value for the cell.'''
        if spot.frame > self.last_frame:
            self.last_frame = spot.frame

    def print_data(self):
        '''Print out cell data.'''
        print('\nName: {0} ({1})\nFrame: {2}\nSpots: {3}'.format(
            self.custom_name, self.generic_name, self.frame, self.n_spots))
        print('\tframe\tx\ty\tz')
        for spot in self.spots:
            print('\t{0}\t{1}\t{2}\t{3}'.format(spot.frame, spot.x, spot.y, spot.z))


class Spot:
    '''A spot is a manually tracked point with x, y, z, t dimensions.'''
    def __init__(self, raw_data):
        self.raw_data = raw_data
        self.valid = False

        # Spot attributes.
        self.frame = 0
        self.x = 0
        self.y = 0
        self.z = 0

        # Parse and validate data.
        self.valid = self.parse_data()

    def parse_data(self):
        '''Parse spot coordinates.'''
        coordinates = self.raw_data.split()

        # Get attributes.
        self.frame = int(coordinates[0])
        self.x = int(coordinates[1])
        self.y = int(coordinates[2])
        self.z = int(coordinates[3])

        if self.frame:
            return True

