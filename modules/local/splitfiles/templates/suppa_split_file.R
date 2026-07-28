#!/usr/bin/env Rscript
# Author: Zifo Bioinformatics
# Email: bioinformatics@zifornd.com
# License: MIT

# Splits any input file (e.g. tpm) by column using samplesheet information
# Utility script for Suppa processing
# Also takes note of how many samples for each condition for down stream Clustering module

#' Parse out options from a string without recourse to optparse
#'
#' @param x Long-form argument list like --opt1 val1 --opt2 val2
#'
#' @return named list of options and values similar to optparse

parse_args <- function(x){
    args_list <- unlist(strsplit(x, ' ?--')[[1]])[-1]
    args_vals <- lapply(args_list, function(x) scan(text=x, what='character', quiet = TRUE))

    # Ensure the option vectors are length 2 (key/ value) to catch empty ones
    args_vals <- lapply(args_vals, function(z){ length(z) <- 2; z})

    parsed_args <- structure(lapply(args_vals, function(x) x[2]), names = lapply(args_vals, function(x) x[1]))
    parsed_args[! is.na(parsed_args)]
}

#########################################################
####### Define function for splitting input files #######
#########################################################

# Function for taking all sample names associated with a given condition
split_files <- function(condition, samplesheet, input_data, output_file_suffix, prefix, calculate_ranges) {

    # Get indices of rows which cover given condition for ranges
    indices <- which(samplesheet\$condition == condition)

    # Get sample names for given condition
    sample_names <- samplesheet[samplesheet\$condition == condition, ]\$sample

    # Check header of input_data contains all samples from processed samplesheet
    if (!all(samplesheet\$sample %in% colnames(input_data))) {
        err_msg <- "suppa_split_file.R Input_file must contain samplesheet samples."
        stop(err_msg, call.=FALSE)
    }

    if (prefix == ""){
        output_file <- paste0(condition, output_file_suffix)
    } else {
        output_file <- paste0(prefix, "_" , condition, output_file_suffix)
    }

    # Subset input files and save out as new file
    write.table(
        input_data[, sample_names, drop=F],
        file = output_file,
        quote = FALSE,
        sep = "\t"
    )

    # Get Cluster ranges which match the tpm and psi files above (1-3 4-6)
    # Column numbers have to be continuous, with no overlapping or missing columns between them. Ex:1-3,4-6
    if (calculate_ranges) {
        ranges_file <- paste0(prefix, "_ranges.txt")
        if (indices[1] == 1) {
            range <- paste0(as.character(indices[1]), "-" , as.character(indices[length(indices)]))
            cat(range, file = ranges_file, sep = "")
        } else {
            range <- paste0(as.character(indices[1]), "-", as.character(indices[length(indices)]))
            cat(c(",", range), file = ranges_file, sep = "", append = TRUE)
        }
    }
}


######################################
########### Collect inputs ###########
######################################

input_file         <- "${tpm_psi}"
samplesheet        <- "${samplesheet}"
output_file_suffix <- "${output_type}"
calculate_ranges   <- as.logical("${calc_ranges}" == "true") # TRUE if we want to return cluster ranges
output_prefix = ifelse('$task.ext.prefix' == 'null', '$meta.id', '$task.ext.prefix')
args_opt <- parse_args('$task.ext.args')

######################################
######### Read in input file #########
######################################

input_data <- read.csv(input_file, sep = "\t", header = TRUE, check.names = FALSE)

######################################
####### Process samplesheet ##########
######################################

# Read in samplesheet
samplesheet <- read.csv(samplesheet, header = TRUE, check.names = FALSE)

# check header of sample sheet
if (!"sample" %in% colnames(samplesheet) | !"condition" %in% colnames(samplesheet)) {
    err_msg <- "suppa_split_file.R samplesheet must contain 'sample' and 'condition' column headers."
    stop(err_msg, call.=FALSE)
}

# Take only sample and condition columns
samplesheet <- samplesheet[, c("sample", "condition")]
stopifnot(colnames(samplesheet) == c("sample", "condition"))

# filter for unique rows based on sample name
samplesheet <- samplesheet[!duplicated(samplesheet[,"sample"]), ]

# Take unique conditions
conditions <- unique(samplesheet[, "condition"])

####################################################################
##### Iterate through conditions - split files and get ranges ######
####################################################################

for (cond in conditions) {
    split_files(
        cond,
        samplesheet,
        input_data,
        output_file_suffix,
        output_prefix,
        calculate_ranges
    )
}

####################################
########## Session info ############
####################################

sink(paste0(output_prefix, "_R_sessionInfo.log"))
print(sessionInfo())
sink()


#########################################
######## Write out versions.yml #########
#########################################

r.version <- strsplit(R.version\$version.string, " ")[[1]][3]
writeLines(
    c(
        "${task.process}:",
        paste0('    r-base: "', r.version, '"')
    ),
    con = "versions.yml"
)
