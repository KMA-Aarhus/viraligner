__author__ = "Rikke Myrhøj Jensen, Tine Sneibjerg Ebsen"
__version__ = "0.1"

import os
import glob
import pandas as pd
import numpy as np
from pathlib import Path
from os.path import isfile, isdir, join, exists, expanduser
from datetime import datetime
import time
import atexit


configfile: "config.yaml"

print(f"         Viral genome aligner v{__version__}  -  Aarhus Universityhospital  -  Department of Clinical Microbiology")
print()
print("                    ██╗   ██╗██╗██████╗  █████╗ ██╗     ██╗ ██████╗ ███╗   ██╗███████╗██████╗")
print("                    ██║   ██║██║██╔══██╗██╔══██╗██║     ██║██╔════╝ ████╗  ██║██╔════╝██╔══██╗")
print("                    ██║   ██║██║██████╔╝███████║██║     ██║██║  ███╗██╔██╗ ██║█████╗  ██████╔╝")
print("                    ╚██╗ ██╔╝██║██╔══██╗██╔══██║██║     ██║██║   ██║██║╚██╗██║██╔══╝  ██╔══██╗")
print("                     ╚████╔╝ ██║██║  ██║██║  ██║███████╗██║╚██████╔╝██║ ╚████║███████╗██║  ██║")
print("                      ╚═══╝  ╚═╝╚═╝  ╚═╝╚═╝  ╚═╝╚══════╝╚═╝ ╚═════╝ ╚═╝  ╚═══╝╚══════╝╚═╝  ╚═╝")
print()
print("                                 Press ctrl+c at any time to stop this pipeline.")
print()

#TO RUN THIS WORKFLOW ACTIVATE CONDA ENVIRONMENT WITH SNAKEMAKE INSTALLED AND ENTER: "snakemake --profile profile"
#PRIMER SEQUENCES AND REFERENCE CAN BE CHANGED IN THE CONFIG
#THIS WORKFLOW IS CURRENTLY NOT USING THE PRIMER SEQUENCES

####################
# INPUT VALIDATION #
####################

#MISSING

#########
# SETUP #
#########

tab = "\t"
nl = "\n"

if config["samplesheet"] == "NA":
    raise Exception("No samplesheet file was given. Please specify a samplesheet by appending --config rundir=\"path/to/samplesheet/\" to the command line call.")
if config["rundir"] == "NA":
    raise Exception("No rundir path was given. Please specify a rundir by appending --config rundir=\"path/to/rundir/\" to the command line call.")

#DIRECTORIES
samplesheet = config["samplesheet"]
rundir = config["rundir"]
scheme_directory = config["scheme_directory"]
scheme_version = config["scheme_version"]
scheme_name = config["scheme_name"]

# For CoverMon
reference = config["reference"]
regions = config["regions"]
threshold = config["threshold"]
maxDepth = config["maxDepth"]

if rundir[-1] == "/":
    print("Removing trailing slash from rundir")
    rundir = rundir[0:-1]

print("Checking that the rundir exists ...                    ", end = "", flush = True)
if not os.path.isdir(rundir):
    raise Exception(f"The rundir does not exist.")
print("✓")

###################
# Validate rundir #
###################

print(f"Looking for MinKNOW-characteristic output:")

for i in range(200):
    print("  Looking ... ", end = "", flush = True)
    fastq_pass_bases = glob.glob(rundir + "/**/fastq_pass", recursive = True) # Find any occurrence of the wanted path
    if len(fastq_pass_bases) == 0:
        print("nothing found yet, waiting 10 secs ...")
        time.sleep(10) # Wait 10 seconds.
    elif(i == 10):
        print() # clean newline
        raise Exception("nothing found after 10 tries. Aborting.")
    else: 
        print(f"Found                                    ✓")
        break

if not len(fastq_pass_bases) == 1:
    raise Exception(f"There seems to be more than one fastq_pass sub-directory beneath the given rundir. These paths were found:{nl} {str(nl + ' ').join(fastq_pass_bases)}{nl}Please specify a more specific rundir.")


fastq_pass_base = fastq_pass_bases[0]
del fastq_pass_bases
print(f"Found the following fastq_pass base which will be given to CoverMon: {nl}  {fastq_pass_base}{nl}")

base_dir = os.path.dirname(fastq_pass_base) # This only works because there is NOT a trailing slash on the fastq_pass_base
print(f"Base directory: {base_dir}")

out_base = "output" # out_base is the directory where the pipeline will write its output to.
print(f"Output directory: {out_base}")


###################################
# Validate and format samplesheet #
###################################

print(f"samplesheet: {samplesheet} \n")

samplesheet_extension = samplesheet.split(".")[-1]
print(f"Reading .{samplesheet_extension}-type sample sheet \"{samplesheet}\"")

# Read spreadsheet
if samplesheet_extension == "xlsx":
    # Uses openpyxl
    df = pd.read_excel(samplesheet, dtype = str)

elif samplesheet_extension == "xls":
    df = pd.read_excel(samplesheet)

elif samplesheet_extension == "tsv":
	df = pd.read_csv(samplesheet, sep="\t")

else:
    raise Exception(f"The spreadsheet file extension {samplesheet_extension} is not yet implemented.")

#######################################################

# Clean up the spreadsheet
print("Cleaning sample sheet ...                              ", end = "", flush = True)
df.columns = map(str.lower, df.columns) # Lowercase
df.columns = map(str.strip, df.columns) # Remove edge-spaces
df.columns = map(lambda x: str(x).replace(" ", "_"), df.columns) # Replace spaces with underscore
df["barcode"] = df["barcode"].apply(np.vectorize(lambda x: str(x).strip().replace(" ", ""))) # Because we are later going to join using this column, it is necessary to strip it for spaces.
df = df.dropna(subset = ["sample_id"])# remove rows not containing a sample ID
print("✓")

# Check that the spreadsheet complies
print("Checking that the necessary columns exist ...          ", end = "", flush = True)
for i in ["barcode", "sample_id"]:
    if not i in df.columns:
        raise Exception(f"The sample sheet is missing a necessary column. The sample sheet must contain the column {i}, but it only contains {df.columns.tolist()}")
print("✓")

acceptable_barcodes = [f"NB{i:02d}" for i in range(1,97)]

print("Checking that the barcodes are correctly formatted ... ", end = "", flush = True)
for i in df["barcode"]:
    if not i in acceptable_barcodes: 
        raise Exception(f"The given barcode \'{i}\' is not an acceptable barcode. Here is a list of acceptable barcodes for inspiration:{nl} {' '.join(acceptable_barcodes)}")
print("✓")

print("Checking that the barcodes are unique ...              ", end = "", flush = True)
if not len(df["barcode"]) == len(set(df["barcode"])):
    bc_counts = pd.DataFrame(df['barcode'].value_counts())
    bc_counts.columns = ["count"]
    bc_counts = bc_counts[bc_counts["count"] > 1]
    raise Exception(f"{nl}One or more barcodes are duplicated. Each barcode may only be used once:{nl}{bc_counts}")
print("✓")

print("Checking that the sample id's are unique ...           ", end = "", flush = True)
if not len(df["sample_id"]) == len(set(df["sample_id"])):
    sid_counts = pd.DataFrame(df['sample_id'].value_counts())
    sid_counts.columns = ["count"]
    sid_counts = sid_counts[sid_counts["count"] > 1]
    raise Exception(f"{nl}One or more sample_id's are duplicated. Each sample_id may only be used once:{nl}{sid_counts}")
print("✓")

#######################################################

# Assign controls
print("Marking sample-types following these definitions:")
print("  positive_control: the sample_id must start with \"SEQPOS\" (case insensitive).")
print("  negative_control: the sample_id must end with \"NEG\" (case insensitive).")

df = df.assign(type = ['positive_control' if a.lower().startswith("seqpos") else ('negative_control' if a.lower().endswith("neg") else 'sample') for a in df['sample_id']])

print()
print("These are the samples from the samplesheet you have given:")
print(df.to_string())
print("//")
print()

######################
# Set up directories #
######################

base_dir = os.path.dirname(fastq_pass_base) # This only works because there is NOT a trailing slash on the fastq_pass_base
print(f"This is the batch base directory:{nl}  {base_dir}") # base_dir is the place where fastq_pass, fast5_pass and the sequencing summary resides.

very_long_batch_id = base_dir.split("/")[-1]
print(f"This is the very long batch id:", very_long_batch_id)

date_parse, time_parse, minion_parse, flowcell_parse, arbhash_parse = very_long_batch_id.split("_")

print("date:    ", date_parse)
print("time:    ", time_parse)
print("minion:  ", minion_parse)
print("flowcell:", flowcell_parse)
print("arbhash: ", arbhash_parse)

batch_id = ".".join(very_long_batch_id.split("_")[0:2]) # The first two words (date, time), joined by a dot.
print(f"This is the parsed batch_id:", batch_id)

out_base = os.path.join(base_dir, "viralign_output") # out_base is the directory where the pipeline will write its output to.

##################
# Start CoverMon #
##################


if config["run_monitoring"] and exists(expanduser("~/viraligner/CoverMon.flag")):
    print(f"Starting covermon in the background ... ")
    # Now we have all resources to start monitoring in the background
    os.system("rm ~/viraligner/CoverMon.flag") 
    print(f"Starting monitoring in the background ... ")
    cov_mon_sh = open("start_mon.sh", "w")
    command = f"#!/bin/bash{nl}source ~/miniconda3/etc/profile.d/conda.sh{nl}cd CoverMon{nl}conda activate covermon {nl}python seq_mon.py '{samplesheet}' {rundir} {reference} {threshold} {maxDepth} {regions}"
    cov_mon_sh.write(command)
    cov_mon_sh.close()
    # Start the sequence monitoring in a new terminal
    os.system("gnome-terminal --tab -- bash start_mon.sh")

##########################
# Wait for run to finish #
##########################

minutes_wait = 10
print("Checking that the sequencing_summary_*.txt-file has been written to disk ...")
while True:
    sequencing_summary_file = glob.glob(base_dir + "/sequencing_summary_*.txt")
    if len(sequencing_summary_file) == 0:
        print(f"  Still sequencing/basecalling; waiting {minutes_wait} minutes ...")
        time.sleep(60*minutes_wait)
    else:
        break

sequencing_summary_file = sequencing_summary_file[0]
print("  The sequencing summary has been found                ✓")
print(f"  This is the sequencing_summary_*.txt-file (full): \"{sequencing_summary_file}\"")

#######################################################

# Backup samplesheet
sample_sheet_given_file = f"{fastq_pass_base}/../sample_sheet_given.tsv"
print(f"Backing up the original sample sheet ...               ", end = "", flush = True)
df.to_csv(sample_sheet_given_file, sep = "\t")
print("✓")

print()

#######################################################


###Create df to match barcode directories to sample id
disk_barcodes_list  = sorted(glob.glob(fastq_pass_base + "/barcode*")) # Find all fastq_pass/barcode* directories
fast5_barcodes_list = sorted(glob.glob(rundir + "/**/fast5_pass/barcode*", recursive = True)) # Find all fast5_pass/barcode* directories

disk_barcodes_df = pd.DataFrame({'barcode_path': disk_barcodes_list, 'fast5_path': fast5_barcodes_list})
disk_barcodes_df = disk_barcodes_df.assign(barcode = ["NB" + i[-2:] for i in disk_barcodes_df["barcode_path"]]) #add barcode

# the workflow_table is the table that contains the records where the barcode could be found on the disk.
workflow_table = disk_barcodes_df.merge(df, how='left', on='barcode') # left join (merge) the present barcodes onto the df table.
workflow_table = workflow_table.dropna(subset = ["sample_id"])

#workflow_table = disk_barcodes_df.merge(df, how='right', on='barcode') # Merge barcodes directories with the samplesheet



#############################
# GENOME ALIGNMENT WORKFLOW #
#############################

rule all:
	input: 
		expand(["{out_base}/{batch_id}_{sample_id}/consensus/{batch_id}_{sample_id}.consensus.fasta", \
				"{out_base}/{batch_id}_{sample_id}/read_filtering/{batch_id}_{sample_id}.fastq", \
				"{out_base}/collected/{batch_id}_collected_input_long.tsv", \
				"{out_base}/flags/{batch_id}_clean_ready.flag.ok", \
				"{out_base}/{batch_id}/{batch_id}_metadata_init.tsv", \
				"{out_base}/{batch_id}/{batch_id}_batch_report.html"],
				out_base = out_base, sample_id = df["sample_id"], batch_id = batch_id, scheme_directory=scheme_directory,scheme_version=scheme_version,scheme_name=scheme_name) 


rule read_filtering:
    input:
        barcode_dir = directory(lambda wildcards: workflow_table[workflow_table["sample_id"] == wildcards.sample_id]["barcode_path"].values[0])
    output: "{out_base}/{batch_id}_{sample_id}/read_filtering/{batch_id}_{sample_id}.fastq" #_barcode00.fastq
    conda: "envs/viraligner.yaml"
    shell: """


    artic guppyplex --skip-quality-check --min-length 100 --directory {input.barcode_dir} --output {output}


    """


rule minion:
    input:
        fastq ="{out_base}/{batch_id}_{sample_id}/read_filtering/{batch_id}_{sample_id}.fastq",
    output:
        consensus = "{out_base}/{batch_id}_{sample_id}/consensus/{batch_id}_{sample_id}.consensus.fasta",
        depths = ["{out_base}/{batch_id}_{sample_id}/consensus/{batch_id}_{sample_id}.coverage_mask.txt.1.depths",
                  "{out_base}/{batch_id}_{sample_id}/consensus/{batch_id}_{sample_id}.coverage_mask.txt.2.depths"]
    params:
        base_dir = base_dir,
        output_dir = "{out_base}/{batch_id}_{sample_id}/consensus/",
        sequencing_summary_file = sequencing_summary_file
    threads: 4
    conda: "envs/minion.yaml"
    shell: """

    # Check that artic minion can run, so we are sure that we are not creating a blank output due to dependency errors.
    artic minion -h > /dev/null 2> /dev/null && echo "artic minion in itself runs fine."


    # Original nanopolish
    artic minion \
        --normalise 200 \
        --threads 4 \
        --scheme-directory {scheme_directory} \
        --scheme-version {scheme_version} \
        --read-file {input.fastq} \
        --fast5-directory {params.base_dir}/fast5_pass \
        --sequencing-summary {params.sequencing_summary_file} \
        {scheme_name} {wildcards.batch_id}_{wildcards.sample_id} \
        || echo ">{wildcards.batch_id}_{wildcards.sample_id}_notenoughdata" > {output.consensus} \
            && cat scripts/197233N.txt >> {output.consensus} 

    # If there is not enough data, the job should exit gracefully and create a blank output file with 197233 N's


    # I have considered that it should only or-exit gracefully when the sample type is "positive_control" or "negative_control". But I think it is very possible that normal samples can also fail, which should not halt the complete batch in terms of the pipeline.
    mv {wildcards.batch_id}_{wildcards.sample_id}.* {params.output_dir}
    """

rule depth:
    input:
        depths = ["{out_base}/{batch_id}_{sample_id}/consensus/{batch_id}_{sample_id}.coverage_mask.txt.1.depths",
                  "{out_base}/{batch_id}_{sample_id}/consensus/{batch_id}_{sample_id}.coverage_mask.txt.2.depths"]
    output: "{out_base}/{batch_id}_{sample_id}/consensus/{batch_id}_{sample_id}_depths.tsv"
    shell: """

        # hattespase

        cat {input.depths} \
        | awk -v sample={wildcards.batch_id}_{wildcards.sample_id} '{{ print sample "\\t" $0 }}' \
        > {output}

    """

rule cat_depth:
	input:
		expand("{out_base}/{batch_id}_{sample_id}/consensus/{batch_id}_{sample_id}_depths.tsv", batch_id = batch_id, out_base = out_base, sample_id = df["sample_id"])
	output: 
		"{out_base}/collected/depths.tsv"
	shell: """

	cat {input} > {output}

	"""


rule merge_variant_data:
    input:
        "{out_base}/collected/depths.tsv"
    output:
        "{out_base}/collected/{batch_id}_collected_input_long.tsv"
    run:
        wft_out = workflow_table
        wft_out = wft_out[['barcode', 'sample_id', 'barcode_path', 'type']]

        wft_out = wft_out.assign(batch_id = batch_id,
            very_long_batch_id = very_long_batch_id)
        wft_out['upload_consensus_file'] = wft_out.apply(lambda row: row.batch_id + "_" + row.sample_id + ".consensus.fasta", axis=1)
        
        wft_out = wft_out.rename(columns = {"batch_id": "#batch_id"})
        
        wft_out.to_csv(str(output[0]), index = False, sep = "\t")
        
rule final_merge:
    input: 
        consensuses = expand("{out_base}/{batch_id}_{sample_id}/consensus/{batch_id}_{sample_id}.consensus.fasta", out_base = out_base, batch_id = batch_id, sample_id = workflow_table["sample_id"]),
        depths = expand("{out_base}/{batch_id}_{sample_id}/consensus/{batch_id}_{sample_id}_depths.tsv", out_base = out_base, batch_id = batch_id, sample_id = workflow_table["sample_id"]),
        collected_input = "{out_base}/collected/{batch_id}_collected_input_long.tsv"
    output:
        file = "{out_base}/collected/{batch_id}_all.tsv",
        dir = directory("{out_base}/clean_upload_{batch_id}"),
        clean_ready = "{out_base}/flags/{batch_id}_clean_ready.flag.ok"
    shell: """

        # Collect all long metadata files together in the collected-directory.
        cat {input.collected_input} | grep -vE "^#" >> {output.file}


        # Make a copy of clean outputs
        mkdir -p {output.dir}
        rm -rf {output.dir}/* # Delete old content if any.
        

        # Copy consensus-files and metadata to the upload directory
        cp {input.consensuses} {output.dir}
        cat {input.depths} > {output.dir}/depths.tsv
        cp {output.file} {output.dir}
        cp {base_dir}/barcode_alignment*.tsv {output.dir}/{batch_id}_barcode_alignment.tsv
        cp {base_dir}/barcode_alignment*.tsv {out_base}/collected/{batch_id}_barcode_alignment.tsv
        cp {base_dir}/final_summary_*.txt {output.dir}/{batch_id}_final_summary.txt
        echo "machine_hostname=$(hostname)" >> {output.dir}/{batch_id}_final_summary.txt
        echo "final_merge_date=$(date --iso-8601=s)" >> {output.dir}/{batch_id}_final_summary.txt

        head -n 1 {base_dir}/throughput_*.csv > {output.dir}/{batch_id}_final_throughput.txt
        tail -n 1 {base_dir}/throughput_*.csv >> {output.dir}/{batch_id}_final_throughput.txt


        # Touch a flag to show that the clean files are OK.
        touch {output.clean_ready}


        # The output.dir now contains the essential files to backup and keep for eternity.


        """
        
# Read the long metadata from the clean upload.
rule metadata_init:
    input: 
        "{out_base}/flags/{batch_id}_clean_ready.flag.ok"
    output: 
        file = "{out_base}/{batch_id}/{batch_id}_metadata_init.tsv"
    conda: "envs/tidyverse.yaml"

    shell: """
        Rscript scripts/metadata_init.r {wildcards.batch_id} {out_base}/clean_upload_{batch_id} {output.file}

    """

rule batch_report:
    input: 
        metadata_init="{out_base}/{batch_id}/{batch_id}_metadata_init.tsv"
    output: "{out_base}/{batch_id}/{batch_id}_batch_report.html"
    params:
        markdown_template_rmd = "scripts/batch_report.rmd",
        markdown_template_html = "scripts/batch_report.html"
    conda: "envs/tidyverse.yaml"
    shell: """

        # The R-markdown template runs with the same working dir as where it lies.
        # It is also hard to give command line arguments to the template.
        # Therefore we need to move the template and the data to a predictable location:
        cp {input.metadata_init} scripts/metadata_init.tsv
        cp {out_base}/clean_upload_{batch_id}/depths.tsv scripts/depths.tsv

        # Generate the batch report
        Rscript -e 'library(rmarkdown); rmarkdown::render("scripts/batch_report.rmd", "html_document")'

        # Clean up temporary files
        mv {params.markdown_template_html} {output}
        rm scripts/metadata_init.tsv
        rm scripts/depths.tsv

        """



rule plot_depth:
	input: 
		depth = "{out_base}/depths.tsv",
		samplesheet = config["samplesheet"]
	output: "{out_base}/depths.png"
	conda: "envs/tidyverse.yaml"
	script: "scripts/depths.R"


