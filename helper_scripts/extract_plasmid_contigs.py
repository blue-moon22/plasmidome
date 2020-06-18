from Bio import SeqIO
import os, sys

def main(argv):
    input_fasta = argv[0]
    scaffold_res = argv[1]
    output_fasta = argv[2]

    headers = []
    with open(scaffold_res, "r") as f:
        for line in f:
            if line.split(",")[1] == "Plasmid":
                headers.append(line.split(",")[0])

    records = []
    for record in SeqIO.parse(input_fasta, "fasta"):
        if record.id in headers:
            records.append(record)

    SeqIO.write(records, output_fasta, "fasta")


if __name__ == "__main__":
    main(sys.argv[-3:])
