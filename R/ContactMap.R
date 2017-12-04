library(Rpdb);

# Convert 3D structure to residue contact map
PDB_to_contact_map <- function(pdb_file, dist_thresh)
{

  P <- read.pdb(pdb_file); # read file
  n <- length(P$atoms$x1); # 999; # get protein length 
  D <- sqrt( kronecker(matrix(1,1,n),P$atoms$x1^2) +t(kronecker(matrix(1,1,n),P$atoms$x1^2))  -2*P$atoms$x1 %*% t(P$atoms$x1) + 
               kronecker(matrix(1,1,n),P$atoms$x2^2) +t(kronecker(matrix(1,1,n),P$atoms$x2^2))  -2*P$atoms$x2 %*% t(P$atoms$x2) + 
               kronecker(matrix(1,1,n),P$atoms$x3^2) +t(kronecker(matrix(1,1,n),P$atoms$x3^2))  -2*P$atoms$x3 %*% t(P$atoms$x3) ); 
  C <- matrix(0, n, n); C[which(D<dist_thresh)] <- 1; # contact map for atoms
  
  
  n_AA <- tail(P$atoms$resid, n=1); # get  protein length in amino-acids
  
  C_AA <- matrix(0, n_AA, n_AA); D_AA <- matrix(0, n_AA, n_AA); # Create contact map for amino acids 
  for (i in 1:n_AA)
  {
    i_AA <- which(P$atoms$resid == i); # sometimes residues are skipped (e.g. 112 in example protein) - why? 
    for (j in 1:n_AA)
    {
      j_AA <- which(P$atoms$resid == j);
      C_AA[i,j] <- max(C[i_AA, j_AA]); # set AA contact if at least one atom is in contact  
      D_AA[i,j] <- min(D[i_AA, j_AA]); # set AA distance
    }    
  }

  return(list(C=C, D=D, C_AA=C_AA, D_AA=D_AA)); 
}


