import gzip
min_len = -1
max_len = 0
with gzip.open(snakemake.input.input_fa, "rt") as infile:
    for line in infile:
        if line.startswith(">"):
            continue
        curr_len = int(len(line)) - 1
        if min_len == -1:
            min_len = curr_len
        elif curr_len < min_len:
            min_len = curr_len
        if curr_len > max_len:
            max_len = curr_len
with open(snakemake.output.raw_len, "w") as out, open(snakemake.input.raw_cnt, "r") as cnt:
    out.write("%s" % cnt.read() )
    out.write("reads.raw.length.min\t%i\n" % min_len)
    out.write("reads.raw.length.max\t%i\n" % max_len)
