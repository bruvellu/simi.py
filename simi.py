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

# TODO:
#
#   - Read and parse .sbc information.
#   - SimiProject class to integrate data from .sbc and .sbd files.
#   - Identify and write parsers for all the .sbd fields.
#   - Reconstruct parent children hierarchy.

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
            return sbd_file
        except:
            print('Could not load the file "{0}"'.format(filepath))

    def parse_sbd(self):
        '''Parse SBD information.'''
        # Temporary buffer for each line of a cell record.
        tmp_cell = ''
        # Read SBD file line by line.
        for line in self.sbd_file.readlines()[7:]:  # skip headers
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
        matrix.write('{0},{1},{2},{3},{4},{5}\n'.format('generic_name', 'birth_frame', 'gen_birth_time', 'color', 'n_spots', 'name'))
        for cell_id, cell in self.cells.items():
            matrix.write('{0},{1},{2},{3},{4},{5}\n'.format(cell.generic_name, cell.birth_frame, cell.generation_birth_time, cell.color, cell.n_spots, cell.name))
        matrix.close()


class Cell:
    '''Store all cell-related information'''
    def __init__(self, raw_data):
        self.raw_data = raw_data
        self.valid = False

        # Cell attributes line 1. Setting the default values to None because
        # the format uses 0 and -1 switches.
        self.cells_left = None
        self.cells_right = None
        self.active_cells_left = None
        self.active_cells_right = None
        self.generic_name = u''

        # Cell attributes line 2. Not sure what "generation" means, but it's
        # not really used anyway.
        self.generation_birth_time = None  # seconds (ms?)
        self.generation_level = None
        self.generation_wildtype = None
        self.generation_color = None
        self.generation_name = u''

        # Cell attributes line 3.
        self.birth_frame = None  # frame
        self.birth_level = None  # level (z) at mitosis
        self.wildtype = None
        self.size = None  # integer from BioCell
        self.shape = None  # integer from BioCell
        self.color = None  # decimal value #TODO convert to HEX.
        self.name = u''  # custom name defined by user

        # Cell attributes line 4.
        self.n_spots = None
        self.comment = u''

        # Cell attributes line 5 and beyond: please see class Spot. Spots are
        # stored in a list.
        self.spots = []

        # Additional attributes.
        self.last_frame = None
        self.parent = None
        self.daughter = None

        # TODO: Each cell has a parent and a sister. Extract such information
        # from the SBD file and add to cell attributes.

        # Parse data, any error returns False (invalid).
        self.valid = self.parse_data()

        if self.cells_right > 1:
            print(self.generic_name, len(self.spots), self.valid)
            print(self.cells_left, self.cells_right, self.active_cells_left, self.active_cells_right)
            print


    def parse_data(self):
        '''Extract attributes from raw data.'''
        split_lines = self.raw_data.split('\n')

        # Line one values.
        try:
            line_one = split_lines[0]
            line_one_split = line_one.split()
            # Loop to capture values.
            for index, value in enumerate(line_one_split):
                if index == 0:
                    self.cells_left = int(value)
                elif index == 1:
                    self.cells_right = int(value)
                elif index == 2:
                    self.active_cells_left = int(value)
                elif index == 3:
                    self.active_cells_right = int(value)
                elif index == 4:
                    self.generic_name = value
        except:
            print('Error parsing line 1!')
            print(self.raw_data)
            return False

        # Line two values.
        try:
            line_two = split_lines[1]
            line_two_split = line_two.split()
            # Loop to capture values.
            for index, value in enumerate(line_two_split):
                if index == 0:
                    self.generation_birth_time = int(value)
                elif index == 1:
                    self.generation_level = int(value)
                elif index == 2:
                    self.generation_wildtype = int(value)
                elif index == 3:
                    self.generation_color = int(value)
                elif index == 4:
                    self.generation_name = value
        except:
            print('Error parsing line 2!')
            print(self.raw_data)
            return False

        # Line three values.
        try:
            line_three = split_lines[2]
            line_three_split = line_three.split()
            # Loop to capture values.
            for index, value in enumerate(line_three_split):
                if index == 0:
                    self.birth_frame = int(value)
                elif index == 1:
                    self.birth_level = int(value)
                elif index == 2:
                    self.wildtype = int(value)
                elif index == 3:
                    self.size = int(value)
                elif index == 4:
                    self.shape = int(value)
                elif index == 5:
                    self.color = int(value)
                elif index == 6:
                    self.name = value
        except:
            print('Error parsing line 3!')
            print(self.raw_data)
            return False

        # Line four values.
        try:
            line_four = split_lines[3]
            line_four_split = line_four.split()
            # Define variables.
            self.n_spots = int(line_four_split[0])
            # Remove first item to join comment string.
            line_four_split.pop(0)
            self.comment = ' '.join(line_four_split)
        except:
            print('Error parsing line 4!')
            print(self.raw_data)
            return False

        # If there are no spots, invalid cell.
        if self.n_spots == 0:
            return False
        else:
            # Get the list index of the last spot in split_lines.
            last_spot_index = 4 + self.n_spots
            # Use the index to add the exact number of spots.
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
            self.name, self.generic_name, self.birth_frame, self.n_spots))
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

