import sys, getopt

def write_reads(sam_file, arg_blast, out_file):

    arg_list = []
    node_list = []
    with open(arg_blast, "r") as bl:
        for line in bl:
            arg_list.append(line)
            node_list.append(line.split("\t")[0])

    contig = ""
    with open(out_file, "w") as out:
        with open(sam_file, "r") as sam:
            for line in sam:
                curr_contig = line.split("\t")[2]
                if curr_contig in node_list and not curr_contig == contig:
                    for i, node in enumerate(node_list):
                        if curr_contig == node:
                            start_pos = int(arg_list[i].split("\t")[6])
                            end_pos = int(arg_list[i].split("\t")[7]) - 100
                            pos = int(line.split("\t")[3])
                            if pos > start_pos and pos < end_pos:
                                out.write(arg_list[i])
                                contig = node


def main(argv):
    try:
        opts, args = getopt.getopt(argv, "h")
    except getopt.GetoptError:
        print("get_reads_mapped_to_ARGs.py <sam_file> <blast_arg_out> <sam_out_file>")
        sys.exit(2)
    for opt, arg in opts:
        if opt == "-h":
            print("get_reads_mapped_to_ARGs.py <sam_file> <blast_arg_out> <sam_out_file>")
            sys.exit()
    sam_file = sys.argv[1]
    arg_blast = sys.argv[2]
    out_file = sys.argv[3]

    write_reads(sam_file, arg_blast, out_file)


if __name__ == "__main__":
    main(sys.argv[1:])
