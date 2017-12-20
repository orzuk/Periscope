import numpy as np
import glob
import os
import argparse
import traceback
import subprocess
import sys


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
        try:
            print("processing file " + f + " - " + str(idx + 1) + " of " + str(len(pdb_files)))

            cmd = os.path.join(".", "distances") + " " + f

            # Run the distances program.
            p = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)

            # Get the distances program output.
            stdout, err = p.communicate()

            stdout = stdout.decode()
            lines = stdout.splitlines()

            i_indices_vec = np.empty(len(lines), dtype = np.int32)
            j_indices_vec = np.empty(len(lines), dtype = np.int32)

            for idx, line in enumerate(lines):
                i_indices_vec[idx] = line.split(' ')[0]

            for idx, line in enumerate(lines):
                j_indices_vec[idx] = line.split(' ')[1]

            unique_i_indices = np.unique(i_indices_vec)
            unique_j_indices = np.unique(j_indices_vec)

            # Row indices should be 0 to n, where n = len(unique_i_indices) - 1
            if not (unique_i_indices == np.arange(len(unique_i_indices))).all():
                raise Exception("Non consecutive row indices for contact matrix.")

            # Column indices should be 1 to n, where n = len(unique_j_indices)
            if not ( unique_j_indices == (np.arange(len(unique_j_indices)) + 1) ).all():
                raise Exception("Non consecutive column indices for contact matrix.")

            # We use +1 since the raw data indices are:
            # i=0...n-1 and j=1...n, so we have n+1 indices in total
            dim_size = len(np.unique(unique_i_indices)) + 1
            contact_mat = np.zeros(shape = (dim_size, dim_size), dtype = np.float16)

            # Create the contact map.
            for line in lines:
                vals = line.split(' ')
                i_idx = int(vals[0])
                j_idx = int(vals[1])
                dist = float(vals[4])
                contact_mat[i_idx, j_idx] = dist
                contact_mat[j_idx, i_idx] = dist

            # Sanity check to ensure the matrix is symmetric
            if not (contact_mat.transpose() == contact_mat).all():
                raise Exception("Contact matrix for " + f + " is not symmetric.")

            contact_mat_csv_output_fname = os.path.join(args.csv_out_dir, os.path.basename(f).split('.')[0] + ".csv")
            # Save the contact map in CSV format.
            np.savetxt(contact_mat_csv_output_fname, contact_mat, delimiter=",", fmt="%.6f")


            contact_mat_np_output_fname = os.path.join(args.np_out_dir, os.path.basename(f).split('.')[0] + ".npy")
            # Save the contact map in numpy format.
            np.save(contact_mat_np_output_fname, contact_mat)

        except KeyboardInterrupt:
                raise
        except:
            print("ERROR: " + str(sys.exc_info()[1]))



# ######## #
#   MAIN   #
# ######## #
if __name__ == "__main__":
    args = parse_arguments()
    validate_args(args)
    run(args)
