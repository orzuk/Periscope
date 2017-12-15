General description:
-------------------
The script `create_contact_maps_from_pdb_files.py` uses the distances program to retrieve contacts maps
from pdb files, and saves them in CSV and numpy format.
All of the pdb files in the specified folder will be processed.

**Note:** If you wish to use the default folders, make sure you first create the missing folders (CSV_files, numpy_files) in the same directory as the `create_contact_maps_from_pdb_files.py` script. The default folders are listed at the end.

Usage:
-----
On Linux/Mac:
```
$ python3 create_contact_maps_from_pdb_files.py
```

Use the following command to view the arguments that the script can parse:
```
$ python3 create_contact_maps_from_pdb_files.py -h
```

Default arguments are shown in angle brackets <>.

An example of running the script and passing all possible arguments to it:
```
$ python3 create_contact_maps_from_pdb_files.py -pdb_dir /data/my_pdb_dir \
  -csv_out_dir /data/my_csv_output_dir -np_out_dir /data/my_numpy_output_dir
```



Default folders:
--------
| folder | default folder name |
| --- | --- |
| pdb files folder | PDB_files |
| csv output folder |    CSV_files |
| numpy output folder |  numpy_files |
