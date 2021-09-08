\name{gen.simuHaplo}
\alias{gen.simuHaplo}
\title{Gene dropping simulations - haplotypes}
\description{Gene dropping from founders to specified probands. If sequence data is provided the function can produce simulated sequences, otherwise will produce proband haplotypes with labelled w.r.t. the founder the segment came from}
\usage{gen.simuHaplo = function (gen, pro, ancestors, simulNo = 1, RecombRate=c(0,0), 
                                    Reconstruction =0, BP=0, Hapfile=NULL, Mapfile=NULL, seed= 0)}

\arguments{
\item{gen}{An object of class GLgen obtained with gen.genealogy, gen.lineages or gen.branching. Required.}
\item{pro}{Vector of proband id numbers to be included. Required.}
\item{ancestors}{Vector of ancestors id numbers to be included}
\item{simulNo}{Number of simulations to run. Results from extra simulations will be included in the same output text file}
\item{RecombRate}{The expected number of recombination events in the simulated region per meiosis. E.g.: if simulating a region of 1Centimorgan, recombination rate should be 0.01
    Vector of 2 elements, the first is recombination rate for male offspring, the second element specifies rate for female offspring}
\item{Reconstruction}{Indicator specifying whether the simulation results should be converted into sequence data.
    Default is Reconstruction = 0, which is no sequence data Reconstruction
    If Reconstruction = 1 then user must provide BP, Hapfile, and Mapfile}
\item{BP}{The length of the simulated region in basepairs. Only required if Reconstruction = 1}
\item{Hapfile}{The name of a file in the working directory  containing sequence data for the founder haplotypes. Each line specifies a haplotype, and should be the same format as the example}
\item{Mapfile}{The name of a file in the working directory which specifies the locations of the SNPs of Hapfile. This file should contain only one number per each line, each line specifies the position (in BP) of the corresponding SNP in the hapfile }
\item{seed}{this will seed the random number gennerator (Mersenne twister) for the simulation. Default is seed = 0 which will choose a random seed. Use the seed for replication purposes}
}

\value{A list containing the following:

 	the seed number of the simulation,

 	the number of meiosis in the simulation,

 	the number of recombination events in the simulation.
}
\seealso{
\code{\link{gen.genealogy}}
\code{\link{gen.simuSample}}
\code{\link{gen.simuSet}}
\code{\link{gen.simuSampleFreq}}
}
\examples{

#load the genealogy table
genealogy_table<-read.csv("~/sample_gen.csv")
View(genealogy_table)

#convert the genealogy table into a GENLIB object
sample_gen<-gen.genealogy(genealogy_table, autoComplete = TRUE) 

#visualize the genealogy
gen.graph(sample_gen)

probands <- gen.pro(sample_gen)
founders <- gen.founder(sample_gen)

# Simulate the haplotypes of all probands, with segments labelled corresponding to the founder it is inherited from
# we let "Reconstruction = 0"  as the default so no sequence data is required, and the generated haplotypes will be in terms of ancestral haplotypes with no sequence info

gen.simuHaplo(sample_gen,pro = probands, ancestors= founders, simulNo=1, RecombRate=c(0.01,0.02), seed=0)

#The output files will be written to the working directory. If you do not wish to give access to this directory, 
# then change the working directory with the command setwd() before running the function

# Read the first 10 lines of the output file "Proband_Haplotypes.txt"
# the first line of the output file is "x;y" where x is the number of simulations and y is the number of probands
# each of the following lines has 3 elements contained in curly braces {1}{2}{3}
# the first curly brace specifies the simulation number and proband number 
# the next two braces specify the two haplotypes for that proband in the following format: {position; founder hap; position; founder hap; position...}
# Where the positions are between 0 and 1 and specify the relative length along the segment being simulated, and founder hap is the ID of the founder haplotype the segment is inherited from

results<-readLines("~/Proband_Haplotypes.txt", n =10)
writeLines(results)

#read the first 10 lines of the output file "all_nodes_haplotypes.txt"
# this output file contains information on all the nodes, not just the probands
# the output format is the same for specifying the haplotypes, except there is additional information contained within the first set of curly braces to allow tracing back recombinations up the genealogy

results <- readLines ("~/all_nodes_haplotypes.txt", n=10)
writeLines(results)

#We can also re-run the function with the Reconstruction = 1 option to obtain sequence data
# But we must provide the function with a Hapfile and Mapfile.

#Hapfile contains the sequence data for each founder haplotype. Each line consists of only 2 elements, the founder haplotype ID #, and the sequence, seperated by a single whitespace
# the sequence data can use any characters, but must be the same length between all founders, and each position in the sequence corresponds to a BP position in the haplotype so any missing data should be encoded using a missing data character and not skipped
# example:
writeLines(readLines("~/sequences.txt", n = 10))

#Mapfile contains the base pair positions of all the characters of the sequences in the Hapfile. Each line describes the position of the corresponding character in the sequence
# there should be no empty lines or whitespaces. Example:
writeLines(readLines("~/SNP_positions.txt", n =10))

gen.simuHaplo(sample_gen,pro = probands, ancestors= founders, simulNo=2, RecombRate=c(0.01,0.02), Reconstruction = 1, BP= 1000000, Hapfile= "sequences.txt", Mapfile = "SNP_positions.txt", seed=777)

#Since we used the reconstruction option, we will now get an additional output file with sequence data for the proband haplotypes:
results<-readLines("~/reconstructed_sequences.txt", n = 10)
writeLines(results)

}
\keyword{datagen}
