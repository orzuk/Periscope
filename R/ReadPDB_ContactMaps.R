#This alphabet contains all letters from the Single-Letter Amino Acid Code (see ?AMINO_ACID_CODE) plus "*" (the stop letter), "-" (the gap letter), "+" (the hard masking letter), and "." (the not a letter or not available letter). It is stored in the AA_ALPHABET predefined constant (character vector).
setwd("C:\\Users\\oorzu\\Google Drive\\CoEvolution\\Code") # setwd("D:\\Yair\\Proteins\\")  # set path

library(Rpdb);

source("ContactMap.R")


data_dir <- "..\\Data\\CPD_ARATH_O04147_1-181_e-4_m30_evfold_results.tar\\"; 
#pdb_file <- "CPD_ARATH_O04147_1-181_e-4_m30_190_20_hMIN.pdb" # predicted structure 
pdb_file <- "CPD_ARATH_O04147_1-181_e-4_m30_pdb1fsi.ent" # true structure (?)

pdb_file_full <- paste(data_dir, pdb_file, sep="");
dist_thresh <- 5; # threshold for declaring contact (in Angstram)

contact_map <- PDB_to_contact_map(pdb_file_full, dist_thresh)

image(1-contact_map$C); # Display contact map for atoms
Cij <- which(contact_map$C != 0, arr.ind=T); # get pairs which are in conctact 
contact_file <- paste(data_dir, 'Output\\\\', pdb_file, '.contact_atom', sep='')
dev.print(pdf, paste(contact_file, '.pdf', sep='')); # save to file
write.table(Cij, file=contact_file, row.names=FALSE, col.names=FALSE);

image(1-contact_map$C_AA); # Display contact map for amino acids 
Cij <- which(contact_map$C_AA != 0, arr.ind=T); # get pairs which are in conctact 
contact_file <- paste(data_dir, 'Output\\\\', pdb_file, '.contact_AA', sep='')
dev.print(pdf, paste(contact_file, '.pdf', sep='')); # save to file
write.table(Cij, file=contact_file, row.names=FALSE, col.names=FALSE);


#P <- read.pdb(system.file(pdb_file,package="Rpdb")); # read file


P <- read.pdb(pdb_file_full) # read file
n <- length(P$atoms$x1); # 999; # get protein length in atoms 
D <- sqrt( kronecker(matrix(1,1,n),P$atoms$x1^2) +t(kronecker(matrix(1,1,n),P$atoms$x1^2))  -2*P$atoms$x1 %*% t(P$atoms$x1) + 
             kronecker(matrix(1,1,n),P$atoms$x2^2) +t(kronecker(matrix(1,1,n),P$atoms$x2^2))  -2*P$atoms$x2 %*% t(P$atoms$x2) + 
             kronecker(matrix(1,1,n),P$atoms$x3^2) +t(kronecker(matrix(1,1,n),P$atoms$x3^2))  -2*P$atoms$x3 %*% t(P$atoms$x3) ); 
C <- matrix(0, n, n); C[which(D<dist_thresh)] <- 1; 

n_AA <- tail(P$atoms$resid, n=1); # get  protein length in amino-acids

C_AA <- matrix(0, n_AA, n_AA); # Create contact map for amino acids 
for (i in 1:n_AA)
{
    i_AA <- which(P$atoms$resid == i); # sometimes residues are skipped (e.g. 112 in example protein) - why? 
    for (j in 1:n_AA)
    {
      j_AA <- which(P$atoms$resid == j);
      C_AA[i,j] <- max(C[i_AA, j_AA]); # set AA contact if at least one atom is in contact  
    }    
}
  
image(1-C_AA); # Display contact map for amino acids 

