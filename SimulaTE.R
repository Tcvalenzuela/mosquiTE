args = commandArgs(trailingOnly = TRUE)

# test if there are at least seven arguments: if not, return an error
if (length(args) < 7) {
  stop("Six arguments must be supplied
       - Genome to analyze as Genome.fasta
       - Index of genome as Genome.fai
       - Repeat library as TE.lib
       - Percentage of insertion on the genome as number
       - Consensus sequences as a multi-fasta file
       - Maximum level of divergence (k)
       - A prefix for the output files", call. = FALSE)
} else {
  fasta = args[1]
  fai = args[2]
  lib = args[3]
  percentage_te = as.numeric(args[4])
  consensus_file = args[5]
  k = as.numeric(args[6])
  prefix = args[7]
}

'''
In this section, the script checks if the required number of arguments is provided. 
It expects a total of seven arguments to be passed to the script: the genome file, 
genome index file, repeat library file, percentage of insertion, consensus sequences file,
maximum level of divergence (k), and a prefix for the output files.
If the number of arguments is incorrect, an error message is displayed.
'''
library(dplyr)
library(data.table)
library(seqinr)

# Read index file
index = fread(fai)

# Read TE monomers
cmd = paste0("grep '>' ", lib, " | cut -c2- ") # Keep the name of the TE lib
names = read.table(text = system(cmd, intern = TRUE), stringsAsFactors = FALSE, comment = "")
cmd = paste0("grep -v '>' ", lib)
seqs = read.table(text = system(cmd, intern = TRUE), stringsAsFactors = FALSE)

# Data frame with TE names and TE monomers
tes = data.frame(Repeat = names$V1, Sequence = seqs$V1)
rm(names, seqs)

'''
This section loads the required libraries (dplyr, data.table, seqinr)
and reads the genome index file (fai) using the fread function from the data.table package.
It also reads the TE monomers by extracting the names and sequences from the repeat library file (lib) 
using shell commands and storing them in a data frame called tes.
'''

# Set total percentage of TEs to add
genomeSize = sum(index$V2)
bps_te = floor((genomeSize / 100) * percentage_te)

# Counter that needs to be updated each cycle
counter = 1

# Number of TE insertions and base pairs occupied
n_insertions = 0
bps_occupied = 0

# Name for the new set of TE insertions .fasta
te_filename = paste0(prefix, "_te_insertions.fasta")

# BED file with all the insertions to add
master_bed = data.frame(Chromosome = index$V1, Start = 0, End = index$V2)

'''
In this section, the total number of base pairs to insert (bps_te) is calculated 
based on the genome size and the specified percentage of insertion (percentage_te).
The counters (counter, n_insertions, bps_occupied) are initialized to keep track of the TE insertions. 
The te_filename variable is set to store the name of the file that will contain the TE insertions. 
The master_bed data frame is created to hold the information about the TE insertions in BED format.
'''


# Read consensus sequences
consensus_seqs <- read.fasta(consensus_file)

# Function to generate a diverged sequence based on the consensus sequence and divergence level
diverge_sequence <- function(consensus, divergence_level) {
  divergence_sites <- sample(1:length(consensus), divergence_level, replace = FALSE)
  diverged_seq <- consensus
  for (site in divergence_sites) {
    base <- consensus[site]
    new_base <- sample(setdiff(c("A", "C", "G", "T"), base), 1)
    diverged_seq[site] <- new_base
  }
  return(paste(diverged_seq, collapse = ""))
}

'''
In this section, the consensus sequences are read from the multi-fasta file
specified by consensus_file using the read.fasta function from the seqinr package. 
The consensus_seqs object contains the consensus sequences, with each sequence associated with its repeat name.
The diverge_sequence function is defined to generate diverged sequences based on a given 
consensus sequence and a divergence level. It randomly selects positions in the sequence and replaces
the bases at those positions with a different randomly chosen base. The function returns the diverged sequence as a string.
'''


while (bps_occupied <= bps_te) {
  # Get a random chromosome
  chr <- sample(index$V1, 1)
  
  # Check if the chromosome has not been used for TE insertion
  if (!grepl(paste0(chr, "_"), master_bed$Chromosome)) {
    # Get a random location on chr
    temp <- master_bed[master_bed$Chromosome == chr, ][sample(1:nrow(master_bed[master_bed$Chromosome == chr, ]), 1), ]
    new_locus <- sample(temp[, 2]:temp[, 3], 1)
    
    # Create TE sequence to insert in the new locus
    te <- sample(1:nrow(tes), 1)
    n_monomers <- sample(1:100, 1)
    consensus_seq <- consensus_seqs[[tes$Repeat[te]]]
    divergence <- sample(0:k, 1)
    diverged_seq <- diverge_sequence(consensus_seq, divergence)
    new_insertion_name <- paste(chr, tes$Repeat[te], n_monomers, divergence, counter, sep = "_")
    
    # Divide the chromosome where the new_locus is
    # Locate the right interval to cut
    row <- which(master_bed$Chromosome == chr & master_bed$Start <= new_locus & master_bed$End >= new_locus)
    new_5end <- data.frame(Chromosome = chr, Start = master_bed$Start[row], End = new_locus)
    new_3end <- data.frame(Chromosome = chr, Start = new_locus, End = master_bed$End[row])
    new_locus_line <- data.frame(Chromosome = new_insertion_name, Start = 0, End = nchar(diverged_seq))
    
    # Add the rows to the data.frame
    master_bed <- add_row(master_bed[-row, ], .after = row - 1, new_5end)
    master_bed <- add_row(master_bed, .after = row, new_locus_line)
    master_bed <- add_row(master_bed, .after = row + 1, new_3end)
    
    counter <- counter + 1
    n_insertions <- n_insertions + 1
    bps_occupied <- bps_occupied + nchar(diverged_seq)
    
    sink(te_filename, append = TRUE)
    cat(paste0(">", new_insertion_name, "\n", diverged_seq, "\n"))
    sink()
  }
}

print(paste("Number of TE insertions:", n_insertions))
print(paste("Base pairs occupied:", bps_occupied))

# Create the new genome one chromosome at a time
# ...
