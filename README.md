# simi.py

Parse and access cell lineage data from Simi BioCell using Python.

## Usage example

First go into Python prompt:

```shell
user@computer:~/path/to/simi.py$ python

Python 2.7.13 (default, Jan 19 2017, 14:48:08) 
[GCC 6.3.0 20170118] on linux2
Type "help", "copyright", "credits" or "license" for more information.
>>> 
```

Now import the `simi.py` library to load and parse `.sbc` and `.sbd` files.

```python
# Load library.
import simi

# Parse .sbc and .sbd files of a SimiProject.
s = simi.SimiProject(sbc_file="lineage.sbc", sbd_file="lineage.sbd")
```

Access cell lineage data.

```python
# The .sbc file contains all settings of the project.
print s.sbc.settings
{'3DCENTER': {'AUTO': '1',
  'CHECK': '1',
  'CX': '458',
  'CY': '262',
  'DX': '0',
  'DY': '0',
  'LBX': '0',
  ...

# The .sbd file contains all the tracked cells.
print s.sbd.cells
OrderedDict([('ABCD', <simi.Cell instance at 0x7fb27c5803f8>),
             ('CD', <simi.Cell instance at 0x7fb27c5806c8>),
             ('D', <simi.Cell instance at 0x7fb27c5805a8>),
             ('1D', <simi.Cell instance at 0x7fb27c580fc8>),
             ('2D', <simi.Cell instance at 0x7fb27c580368>),
             ('3D', <simi.Cell instance at 0x7fb27c580e18>),
             ('4D', <simi.Cell instance at 0x7fb27c5c6d88>),
             ...
```

Get information from every individual cell and tracked spots.

```python
# Access each cell individually by name.
c = s.sbd.cells['3D']
print c.parent
CELL=2D

# Get all spots with xyz coordinates:
for spot in c.spots:
    print spot.frame, spot.x, spot.y, spot.z

329 241 34 348
320 256 39 369
324 256 40 381
326 251 42 409
326 248 40 426
322 248 41 430
...
```

## Convert Simi to MaMuT

[simi2mamut.py](simi2mamut.py) is the script I use to convert my Simi lineage
files to MaMuT. It will _not work_ as is on your own lineage files, but it
should work after adaptating the code. Please check out the
[script documentation](simi2mamut.py#L8) and let me know if you run into
[issues](https://github.com/nelas/simi.py/issues).

## Disclaimer

The code has only been tested in my own lineage files. I had to manually fix
some issues with Simi format before parsing. Most crucial is to make sure each
cell has a unique name. The library also cannot handle very well alternative
lineages, if present. Always double check the data before analysing!

Please report any [issues](https://github.com/nelas/simi.py/issues).

## Citation

The code herein is part of the following publication:

Vellutini, B. C., Martín-Durán, J. M. & Hejnol, A. Cleavage modification did
not alter blastomere fates during bryozoan evolution. BMC Biol. 15, 33 (2017).
doi: [10.1186/s12915-017-0371-9](http://dx.doi.org/10.1186/s12915-017-0371-9)

