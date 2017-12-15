General description:
-------------------
The script create_contact_maps_from_pdb_files.py uses the distances program to retrieve contacts maps
from pdb files, and saves them in CSV and numpy format.
All of the pdb files in the specified folder will be processed.

Usage:
-----
On linux/mac:
$ python3 create_contact_maps_from_pdb_files.py

Use the following command to view the arguments that the script can parse:
$ python3 create_contact_maps_from_pdb_files.py -h

Default arguments are shown in angle brackets <>.

Defaults:
--------
pdb files folder:     pdb_files
csv output folder:    CSV_files
numput output folder: numpy_files
