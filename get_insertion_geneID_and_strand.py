import os
import subprocess

import pandas as pd

# ============================================================================================================
# get_insertion_geneID_and_strand.py
#
# Sebastian Mueller at Oregon State University Department of Chemical, Biological, Environmental Engineering
# email: muellese@oregonstate.edu
# 1/22/2024
#
# The goal of this program is to take mutation sequence locations and locate what genes overlap with the
# mutation sequence location. This will provide information about which genes have a mutation insertion and
# which mutation it is.

# Specifically, This python file will take a Ds Gfp insertion gff file and a genome anotation gff file.
# Overlapping regions between Ds Gfp insertions file and the genome annotation file are found and program
# returns the cooresponding gene ID and the orientation/strand of the cooresponding overlap gene.
#
# Lastly, the program will insert the gene ID and strand information into a copy of template CSV and provide
# a new CSV named "new_insertion_features.csv". Template CSV can be blank, otherwise it should have same key
# used for merge strategy at end of file.
#
# Required files include: insertion mutation gff named "insertion_location.gff3", refrence genome annotation
# named "genome_annotation.gff3", a template csvnamed "template.csv", a copy of Dr. Molly Megraw's binary
# search perl program named "gff_genomics.pl". Place these files in same directory of bash script and python
# file

# ===========================================================================================================

# ============================ Check for needed files and run binary search =================================
print("checking for needed files:")

cwd = os.getcwd()
valid_files = {
    "genome_annotation.gff3",
    "insertion_location.gff3",
    "template.csv",
    "gff_genomics.pl",
}

# check to make sure files are present
test_set = (
    set()
)  # empty set to populate with matching file names to check all files present
for filename in os.listdir(cwd):
    if filename in valid_files:
        test_set.add(filename)

if valid_files == test_set:
    print("All files present and named correctly!")

# initalize pathname variables for all files
annotation_path = os.path.join(cwd, "genome_annotation.gff3")
insertion_path = os.path.join(cwd, "insertion_location.gff3")
template_csv = os.path.join(cwd, "template.csv")
perl_path = os.path.join(cwd, "gff_genomics.pl")

# produce the binary search string to pass to binary search program
binary_search_command = (
    f"perl {perl_path} -G {annotation_path} {insertion_path} > bisearch_overlaps.tab"
)

# run binary search
subprocess.run(binary_search_command, shell=True)

# ================ produce binary search df and clean ========================================================

# produce pd df with the info produced by the pl bisearch program
bisearch_output_df = pd.read_csv(
    "bisearch_overlaps.tab",
    sep="\t",
    comment="#",
    names=["Ds_GFP_allele", "v5_gene", "Gene_ID_strand"],
)

# Molly's pl program produced an output file with CSV headers, want to drop this to keep DF clean
bisearch_output_df.drop(0, inplace=True)

# clean 'Ds_GFP_allele' column to only contain insertion ID string
bisearch_output_df["Ds_GFP_allele"] = bisearch_output_df["Ds_GFP_allele"].apply(
    lambda x: x.split("=")[1].split(";")[0] if "=" in x and ";" in x else None
)

# clean 'v5_column' to contain gene ID and populate 'Gene_ID_strand' column
# key word and delimiter variables to use in string splits for finiding just the gene ID and strand
keyword = "logic_name=cshl_gene"
delimiter = "|"


def extract_ID_and_strand(x: str) -> list:
    """Function takes a string does two string splits to pull gene ID and Strand"""

    # this will pull only the gene id that is in a deliniated section with keyword indicating it is gene feature
    gene_ID = (
        [pair.split("=")[1] for pair in x.split(delimiter) if keyword in pair][0].split(
            ";"
        )[0]
        if any(keyword in pair for pair in x.split(delimiter))
        else None
    )

    # this will pull cooresponding strand info from the gene id pull above
    strand = (
        [pair.split("Strand:")[1] for pair in x.split(delimiter) if keyword in pair][0]
        if any(keyword in pair for pair in x.split(delimiter))
        else None
    )
    return [gene_ID, strand]


# apply the function to all pandas columns
bisearch_output_df[["v5_gene", "Gene_ID_strand"]] = bisearch_output_df["v5_gene"].apply(
    lambda x: pd.Series(extract_ID_and_strand(x))
)

# ============================ produce final CSV with new data ==============================================

# produce new df with the template and wranggled overlap dfs (can manually chage key for merge)
current_insertion_features_df = pd.read_csv(template_csv)
new_insertion_features_df = pd.merge(
    current_insertion_features_df, bisearch_output_df, how="outer", on="Ds_GFP_allele"
)

# produce the new overlap features csv
new_insertion_features_df.to_csv("new_insertion_features.csv", index=False)
