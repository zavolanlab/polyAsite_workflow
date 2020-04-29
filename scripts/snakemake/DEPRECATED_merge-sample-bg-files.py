import gzip

sites = []
reads = {}
header_lines = []

with gzip.open(snakemake.input.raw_table, "rt") as infile:
    for line in infile:
        if line.startswith("#"):
            header_lines.append(line)
            continue
        F = line.rstrip().split("\t")
        site_id = ":".join(F[0:3])
        sites.append(site_id)
        reads[site_id] = [F[0], F[1], F[2], F[-2], F[-1]]

for in_f in snakemake.input.filtered:
    with open(in_f, "r") as ifile:
        for line in ifile:
            line_list = line.rstrip().split("\t")
            curr_id = ":".join(line_list[0:3])
            reads[curr_id].insert(-2, line_list[3])

with gzip.open(snakemake.output.table_adjusted, "wt") as out_file:
    for h in header_lines:
        out_file.write("%s" % h)
    for s in sites:
        out_file.write("%s\n" % "\t".join( reads[s] ) )
