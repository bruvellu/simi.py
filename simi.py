#!/usr/bin/env python

'''Read and access cell lineage data from Simi BioCell.

Small library to parse Simi*BioCell project files (.sbc and .sbd). It collects
all the relevant cell lineage data such as cell name, frame and coordinates,
and provides a programmatic interface to access the data. In this manner the
data can be analyzed by custom functions, or exported to different formats for
data analyses in more up-to-date cell lineage software.

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

'''

from collections import OrderedDict

# TODO:
#
#   - Read and parse .sbc information.
#   - SimiProject class to integrate data from .sbc and .sbd files.
#   - Identify and write parsers for all the .sbd fields.


class SimiProject:
    '''Simi BioCell project file.'''
    def __init__(self, sbc_file, sbd_file):
        self.sbc = Sbc(sbc_file)
        self.sbd = Sbd(sbd_file)


class Sbc:
    '''Settings file for Simi BioCell.'''
    def __init__(self, sbc_file):
        self.sbc_file = self.open_sbc(sbc_file)

        # Main dictionary with all settings.
        self.settings = {}

        # Parse data.
        self.parse_sbc()

    def open_sbc(self, filepath):
        '''Save abspath and try to read file.'''
        try:
            sbc_file = open(filepath)
            return sbc_file
        except:
            print('Could not load the file "{0}"'.format(filepath))

    def parse_sbc(self):
        '''Parse SBD information.'''
        # Stores the current setting.
        current_setting = u''
        # Read SBC file line by line.
        for line in self.sbc_file.readlines():
            if line.startswith(';') or line.startswith('\r'):
                continue
            else:
                if line.startswith('['):
                    current_setting = line.strip('[]\r\n')
                    self.settings[current_setting] = {}
                else:
                    splitted = line.split('=')
                    key = splitted[0]
                    value = splitted[1].strip('\r\n')
                    self.settings[current_setting][key] = value



class Sbd:
    '''Database file for Simi BioCell.'''
    def __init__(self, sbd_file):
        self.sbd_file = self.open_sbd(sbd_file)

        # Main dictionary with all cells.
        self.cells = OrderedDict()

        # Main dictionary with valid cells.
        self.valid_cells = OrderedDict()

        # Dictionary for invalid cells (=without spots).
        self.invalid_cells = OrderedDict()

        # Last frame of the recording.
        self.last_frame = 0  # TODO: get from sbc?

        # Parse file.
        self.parse_sbd()

        # Generate descendants.
        for name, cell in self.cells.iteritems():
            cell.get_descendants()

    def __str__(self):
        return self.sbd_file.name

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
        # Temporary switches for parenthood.
        has_parent = False
        parent_cell = None
        # Store for generation birth time of right sibling cell.
        sister_cells = {}

        # Read SBD file line by line.
        for line in self.sbd_file.readlines()[7:]:  # skip headers
            if line.startswith('---'):
                if tmp_cell:
                    # Create cell instance with raw data.
                    cell = Cell(tmp_cell)
                    cell.sbd = self

                    # Define parent cell.
                    if has_parent:
                        cell.parent = parent_cell
                        cell.parent.daughters.append(cell)
                    else:
                        # If cell has a parent, but it's upstream in the lineage.
                        if cell.generation_birth_time in sister_cells.keys():
                            # Get common parent of sibling cells.
                            common_parent = sister_cells[cell.generation_birth_time].parent
                            cell.parent = common_parent
                            cell.parent.daughters.append(cell)

                    # If cell has a left daughter, turn switch on and define parent.
                    if cell.cells_left == 1:
                        has_parent = True
                        parent_cell = cell
                    elif cell.cells_left == 0:
                        has_parent = False
                        parent_cell = None

                    # If cell has a right daughter, add to the sibling dictionary.
                    if cell.cells_right == 1:
                        sister_cells[cell.generation_birth_time] = cell

                    # Add cell to main dictionary.
                    self.cells[cell.generic_name] = cell

                    # If cell is valid, add it to list.
                    if cell.valid:
                        self.valid_cells[cell.generic_name] = cell
                        self.update_last_frame(cell)
                    else:
                        self.invalid_cells[cell.generic_name] = cell

                    # Clean temporary cell.
                    tmp_cell = ''
            else:
                # Add next line to cell data.
                tmp_cell = tmp_cell + line

    def update_last_frame(self, cell):
        '''Updates the last_frame value for the simi instance.'''
        if cell.last_frame > self.last_frame:
            self.last_frame = cell.last_frame

    def get_cells_without_parent(self):
        '''Returns a list of valid cells that have no parent.'''
        no_parent = []
        for cell in self.list_of_cells:
            if not cell.parent:
                no_parent.append(cell)
        return no_parent

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
        self.sbd = None
        self.sbc = None  # TODO: not defined yet.
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
        self.wildtype = None  # fate of the cell, check .sbc for list
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
        self.daughters = []

        # Traverse values.
        self.parents = {}
        self.descendants = {}

        # Parse data, any error returns False (invalid).
        self.valid = self.parse_data()

    def __str__(self):
        return 'CELL={name}'.format(name=self.generic_name)

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
                new_spot.cell = self
                self.spots.append(new_spot)
                self.update_last_frame(new_spot)
            # Change status to valid.
            return True

    def update_last_frame(self, spot):
        '''Updates the last_frame value for the cell.'''
        if spot.frame > self.last_frame:
            self.last_frame = spot.frame

    def get_descendants(self):
        '''Returns the descendants.

        This function calculates the descendants for all the descendants
        recursively.  It can be used as an initial method to populate the
        descendants or to get the current descendants of a cell.

        A dictionary is returned.
        '''
        # Iterate through daughter cells.
        for child in self.daughters:
            # Add child to dictionary.
            self.descendants[child.generic_name] = child
            # Descend in the tree.
            self.descendants.update(child.get_descendants())
        return self.descendants

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
        self.cell = None
        self.valid = False

        # Spot attributes.
        self.frame = 0
        self.x = 0
        self.y = 0
        self.z = 0

        # Parse and validate data.
        self.valid = self.parse_data()

    def __str__(self):
        return 'CELL={cell} FRAME={frame} X={x} Y={y} Z={z}'.format(
                cell=self.cell.generic_name, frame=self.frame, x=self.x, y=self.y, z=self.z)

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

