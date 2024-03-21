import argparse
import logging
import csv
import subprocess

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("append_amino_acids")
logger.setLevel(logging.INFO)

'''
requires three parts:
    - gene
    - transcript
    - exon and CDS info
'''
def get_info(in_file_path: str):
    """
    Generator function that reads in the bed file.

    :param in_file_path: String containing file path to bed file.
    """
    with open(in_file_path, newline="") as infile:
        reader = csv.reader(infile, delimiter="\t")
        #bed row format: chr start stop geneID score strand thick-start thick-end itemRgb block(exon)-count blockSizes blockStarts
        for row in reader:
            #rows_to_write = []
            #try:
            #    exon_length = row[10].split(",")
            #    exon_start = [int(x)+1 for x in row[10].split(",")] #convert from 0 based to 1 based coordinates
            #    exons = zip(exon_start, exon_length)
            #except IndexError:
            #    logger.error(f"No exon information found. Exon information for {row[3]} will not be written.")
            #gtf row format: chr location feature start stop score strand frame attribute
            #gene_row = [row[0], row[]]
            #transcript_row = 
            #rows_to_write.append([])
            if row[2]=="transcript":
                gene_row = row.copy()
                gene_row[2] = "gene"
                yield [ gene_row, row ]
            elif row[2]=="CDS" or row[2]=="start_codon" or row[2]=="stop_codon":
                continue
            else:
                yield [row]

def main(args):
    with open(args.out, "w", newline="") as outfile:
        for x in get_info(args.bed):
            writer = csv.writer(outfile, delimiter="\t", lineterminator="\n", quoting = csv.QUOTE_NONE, quotechar="")
            writer.writerows(x)
    logger.info(f"File written to {args.out}.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Convert bed file to gtf."
    )
    parser.add_argument(
        "-b",
        "--bed",
        required=True,
        help="Bed file",
    )
    parser.add_argument(
        "-o",
        "--out",
        required=True,
        help="Name for output gtf file. Include the .gtf extension when inputing.",
    )
    args = parser.parse_args()
    main(args)