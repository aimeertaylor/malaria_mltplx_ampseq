#!/usr/bin/env python

def read_gatk_table(file):
    table = {}
    with open(file) as f:
        f.readline() # header
        for line in f:
            line = line.strip().split("\t")
            chrom = line[0]
            if chrom not in table:
                table[chrom] = {}
            pos = int(line[1])
            table[chrom][pos] = line[2:]
    return table

def get_amp_frqs(table, filter_missing=True):
    amp_gts = {}
    for chrom in table:
        if chrom not in amplicons:
            continue
        for pos in table[chrom]:
            if pos not in amplicons[chrom]:
                continue
            amp = amplicons[chrom][pos]
            if amp in amp_gts:
                for i, gt in enumerate(table[chrom][pos]):
                    amp_gts[amp][i] += gt
            else:
                amp_gts[amp] = table[chrom][pos]
    
    if filter_missing:
        for amp in list(amp_gts.keys()):
            amp_gts[amp] = [gt for gt in amp_gts[amp] if '.' not in gt]
    
    frqs = {}
    for amp in amp_gts:
        frqs[amp] = {}
        for gt in set(amp_gts[amp]):
            frqs[amp][gt] = float(amp_gts[amp].count(gt)) / len(amp_gts[amp])
    
    return frqs

amplicons = {}
with open("../amplicons.txt") as f:
    f.readline()
    for line in f:
        (chrom, start, stop, name) = line.strip().split("\t")
        if chrom not in amplicons:
            amplicons[chrom] = {}
        start = int(start)
        stop = int(stop) + 1
        if start < stop:
            for i in range(start, stop):
                amplicons[chrom][i] = name
        else:
            for i in range(stop, start):
                amplicons[chrom][i] = name


senegal = read_gatk_table("Senegal.table")
senegal_amplicon_frqs = get_amp_frqs(senegal)

columbia = read_gatk_table("Columbia.table")
columbia_amplicon_frqs = get_amp_frqs(columbia)

french_guiana = read_gatk_table("FrenchGuiana.table")
french_guiana_amplicon_frqs = get_amp_frqs(french_guiana)


all_amps = set(list(senegal_amplicon_frqs.keys())).union(list(columbia_amplicon_frqs.keys())).union(list(french_guiana_amplicon_frqs.keys()))

combined_amp_frqs = {}
countries = ['Senegal', 'Columbia', 'FrenchGuiana']
for amp in all_amps:
    combined_amp_frqs[amp] = {}
    for i, amp_frqs in enumerate([senegal_amplicon_frqs, columbia_amplicon_frqs, french_guiana_amplicon_frqs]):
        if amp in amp_frqs:
            for gt in amp_frqs[amp]:
                if gt not in combined_amp_frqs[amp]:
                    combined_amp_frqs[amp][gt] = {}
                combined_amp_frqs[amp][gt][countries[i]] = amp_frqs[amp][gt]

max_gts = max([len(combined_amp_frqs[amp]) for amp in combined_amp_frqs])

with open("combined_amp_frqs.txt", "wt") as w:
    w.write("Amplicon\tGenotype\tSenegal\tColumbia\tFrenchGuiana\n")
    for amp in sorted(combined_amp_frqs):
        for gt in sorted(combined_amp_frqs[amp]):
            s_frq = combined_amp_frqs[amp][gt].get('Senegal', 0)
            c_frq = combined_amp_frqs[amp][gt].get('Columbia', 0)
            f_frq = combined_amp_frqs[amp][gt].get('FrenchGuiana', 0)
            w.write("{}\t{}\t{}\t{}\t{}\n".format(amp, gt, s_frq, c_frq, f_frq))



with open("combined_amp_frqs_sparse_matrix.txt", "wt") as w:
    w.write("Amplicon\tCountry\t{}\n".format("\t".join(["Genotype.{}".format(i+1) for i in range(max_gts)])))
    for country in sorted(countries):
        for amp in sorted(combined_amp_frqs):
            w.write("{}\t{}\t".format(amp, country))
            w.write("\t".join(['{:g}'.format(combined_amp_frqs[amp][gt].get(country, 0)) for gt in sorted(combined_amp_frqs[amp])]))
            fill_zeros = max_gts - len(combined_amp_frqs[amp])
            if fill_zeros > 0:
                w.write("\t{}".format("\t".join(["0"]*fill_zeros)))
            w.write("\n")
                


