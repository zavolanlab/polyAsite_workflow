import gzip

map={}
sites = []
reads = {}
header_lines = []

# Get info from original 3pSites file
with gzip.open(snakemake.input.raw_table, "rt") as infile:
    for line in infile:
        # Header lines
        if line.startswith("#"):
            l= line.rstrip().split(";")
            col=int(l[0].lstrip("#"))
            sample=l[1].split("/")[-1].split(".")[0]
            # map column to sample
            map[sample]=col
            header_lines.append(line)
            continue
        F = line.rstrip().split("\t")
        site_id = ":".join(F[0:3])
        sites.append(site_id)
        # For each site, store list with appropriate length to accomodate all samples
        reads[site_id] = [F[0], F[1], F[2]] + [None]*len(map) + [F[-2], F[-1]]

for in_f in snakemake.input.filtered:
    # Find sample id
    sample=in_f.rstrip().split("/")[-1].split(".")[0]
    with open(in_f, "r") as ifile:
        for line in ifile:
            line_list = line.rstrip().split("\t")
            curr_id = ":".join(line_list[0:3])
            reads[curr_id][map[sample]]=line_list[3]

with gzip.open(snakemake.output.table_adjusted, "wt") as out_file:
    for h in header_lines:
        out_file.write("%s" % h)
    for s in sites:
        out_file.write("%s\n" % "\t".join( reads[s] ) )
