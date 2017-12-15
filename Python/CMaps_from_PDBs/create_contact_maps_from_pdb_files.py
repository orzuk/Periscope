import numpy as np
import glob
import os
import argparse
import subprocess


# ############# #
#   FUNCTIONS   #
# ############# #
def parse_arguments():
    parser = argparse.ArgumentParser()

    parser.add_argument('-pdb_dir', action="store", dest="pdb_dir", default="PDB_files",
                        help="The directory which contains the PDB files. <PDB_files>")

    parser.add_argument('-csv_out_dir', action="store", dest="csv_out_dir", default="CSV_files",
                        help="The directory which will contain the output in CSV format. <CSV_files>")

    parser.add_argument('-np_out_dir', action="store", dest="np_out_dir", default="numpy_files",
                        help="The directory which will contain the contact maps in numpy format. <numpy_files>")

    return parser.parse_args()


def validate_args(args):
    if not os.path.isdir(args.pdb_dir):
        msg = "Invalid value for -pdb_dir: cannot find the folder '" + args.pdb_dir + "'."
        raise Exception(msg)

    if not os.path.isdir(args.csv_out_dir):
        msg = "Invalid value for -csv_out_dir: cannot find the folder '" + args.csv_out_dir + "'."
        raise Exception(msg)

    if not os.path.isdir(args.np_out_dir):
        msg = "Invalid value for -np_out_dir: cannot find the folder '" + args.np_out_dir + "'."
        raise Exception(msg)

    pdb_files = glob.glob(os.path.join(args.pdb_dir, "*.pdb"))

    if len(pdb_files) == 0:
        msg = "No PDB files found in folder '" + args.pdb_dir + "'."
        raise Exception(msg)


def run(args):
    # Get the names of the PDB files.
    pdb_files = glob.glob(os.path.join(args.pdb_dir, "*.pdb"))

    for idx, f in enumerate(pdb_files):
        print("processing file " + str(idx + 1) + " of " + str(len(pdb_files)))

        cmd = os.path.join(".", "distances") + " " + f

        # Run the distances program.
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)

        # Get the distances program output.
        stdout, err = p.communicate()

        stdout = stdout.decode()
        lines = stdout.splitlines()
        # Create a numpy 2d martix from the distances program output.
        raw_data = np.array([ x.split(' ') for x in lines ], dtype="float64")

        # print("matrix shape: " + str(raw_data.shape))

        raw_data_output_fname = os.path.join(args.csv_out_dir, os.path.basename(f).split('.')[0] + "_raw_data.csv")

        # Save the raw output to CSV format.
        np.savetxt(raw_data_output_fname, raw_data, delimiter=",", fmt="%.6f")

        nrows = len(np.unique(raw_data[:, 0]))
        ncols = nrows

        contact_mat = np.zeros(shape = (nrows, ncols))

        # Create the contact map.
        for row in raw_data:
            i_idx = int(row[0]) - 1
            j_idx = int(row[1]) - 1
            dist = row[4]
            contact_mat[i_idx, j_idx] = dist
            contact_mat[j_idx, i_idx] = dist

        contact_mat_csv_output_fname = os.path.join(args.csv_out_dir, os.path.basename(f).split('.')[0] + ".csv")
        # Save the contact map in CSV format.
        np.savetxt(contact_mat_csv_output_fname, contact_mat, delimiter=",", fmt="%.6f")


        contact_mat_np_output_fname = os.path.join(args.np_out_dir, os.path.basename(f).split('.')[0] + ".npy")
        # Save the contact map in numpy format.
        np.save(contact_mat_np_output_fname, contact_mat)



# ######## #
#   MAIN   #
# ######## #
if __name__ == "__main__":
    args = parse_arguments()
    validate_args(args)
    run(args)
