Python library to parse and access cell lineage data from Simi BioCell

Example usage:

```
# Load library.
import simi

# Parse .sbc and .sbd files of a SimiProject.
s = simi.SimiProject(sbc_file="lineage.sbc", sbd_file="lineage.sbd")

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

# Access each cell individually.
c = s.sbd.cells['3D']
print c.parent
CELL=2D

# Get all spots with xyz coordinates:
for spot in c.spots:
    print spot.frame spot.x, spot.y, spot.z

329 241 34 348
320 256 39 369
324 256 40 381
326 251 42 409
326 248 40 426
322 248 41 430
...
```

