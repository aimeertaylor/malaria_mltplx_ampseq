#!/usr/bin/env python

# fill in variables here:
# path to file of genomic coordinates for panel:
coords = "gtseq_coords.txt"
# dictionary of country names and their GATK table file paths:
country_paths = {"colombia": "colombia.table",
                 "frenchguiana": "frenchguiana.table"}


def read_gatk_table(file):
    """ Reads a single file in GATK table format into memory.

    Args:
      file: file path for GATK table for a single country/population.

    Returns:
      table: a nested dictionary, with a key for each chromosome. each
      of these points to a nested dictionary with keys of coordinates pointing
      to lists of variants per sample.
    """
    table = {}
    with open(file) as f:
        f.readline()  # header
        for line in f:
            line = line.strip().split("\t")
            chrom = line[0]
            if chrom not in table:
                table[chrom] = {}
            pos = int(line[1])
            table[chrom][pos] = [x[0] for x in line[2:]]
    return table


def read_all_tables(country_paths_dict):
    """ Wrapper function to read in any number of GATK tables.

    Args:
      country_paths_dict: dictionary with keys of country/population names and
      values of the paths to the GATK table per country.

    Returns:
        dictionary where a key is a country/population name and a value is a
        table object, from read_gatk_table, for that country.
    """
    country_gatk_table_dict = {}
    for key, value in country_paths_dict.items():
        country_gatk_table_dict[key] = read_gatk_table(value)
    return country_gatk_table_dict


def read_amplicon_coordinates(coord_file):
    """ Reads in genomic coordinates from a BED file.

    Args:
      coord_file: path to tab-separated file with columns chromosome, start
      coordinate, stop coordinate, and name of amplicon.

     Returns:
       dictionary where a key is a chromosome name and a value is a nested
       dictionary with a key for each position covered by an amplicon and a
       value of the amplicon name.
    """
    amplicons = {}
    with open(coord_file) as f:
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
    return amplicons


def get_amp_freqs(table, amplicon_var, filter_missing=True):
    """ Calculates haplotype frequencies for a list of window coordinates.

    Args:
      table: a processed GATK table, output from read_gatk_table.
      amplicon_var:

    Returns:
      table: a nested dictionary, with a key for each chromosome. each
      of these points to a nested dictionary with keys of coordinates pointing
      to lists of variants per sample.
    """
    amp_gts = {}
    for chrom in table:
        if chrom not in amplicon_var:
            continue
        for pos in table[chrom]:
            if pos not in amplicon_var[chrom]:
                continue
            amp = amplicon_var[chrom][pos]
            if amp in amp_gts:
                for i, gt in enumerate(table[chrom][pos]):
                    amp_gts[amp][i] += gt
            else:
                amp_gts[amp] = table[chrom][pos]
    if filter_missing:
        for amp in list(amp_gts.keys()):
            amp_gts[amp] = [gt for gt in amp_gts[amp] if '.' not in gt]
    freqs = {}
    for amp in amp_gts:
        freqs[amp] = {}
        for gt in set(amp_gts[amp]):
            freqs[amp][gt] = float(amp_gts[amp].count(gt)) / len(amp_gts[amp])
    return freqs


def get_freqs_per_country_one_panel(coord_file, country_gatk_table_dict,
                                    output_prefix, output_matrix=True):
    """ Wrapper function to write files with amplicon frequencies in all
    countries listed for a single amplicon panel. Will output a tsv of
    haplotype frequencies; can also output in sparse matrix format.

    Args:
      coord_file: tsv of chromosome, start coordinate, stop coordinate, and
        amplicon name, called by read_amplicon_coordinates.
      country_gatk_table_dict: output from read_all_tables, dictionary with
        country names linked to dictionaries of variant data.
      output_prefix: string with text to prefix output filenames.
      output_matrix: output in sparse matrix format.

     Returns:
       writes haplotype frequencies in tsv format.
    """
    amplicons = read_amplicon_coordinates(coord_file)
    countries = sorted(country_gatk_table_dict.keys())

    # get amp freqs per country
    all_amps = {}
    amp_freq_dict = {}
    for key, value in country_gatk_table_dict.items():
        amp_freq_dict[key] = get_amp_freqs(value, amplicons)
        if len(all_amps) == 0:
            all_amps = set(list(amp_freq_dict[key].keys()))
        else:
            all_amps = all_amps.union(list(amp_freq_dict[key].keys()))

    # combine all country amp freqs into single data structure
    comb_amp_freqs = {}
    for amp in all_amps:
        comb_amp_freqs[amp] = {}
        for i, amp_freqs in enumerate(amp_freq_dict.values()):
            if amp in amp_freqs:
                for j, gt in enumerate(sorted(amp_freqs[amp],
                                              key=lambda x: amp_freqs[amp][x],
                                              reverse=True)):
                    if j not in comb_amp_freqs[amp]:
                        comb_amp_freqs[amp][j] = {}
                    comb_amp_freqs[amp][j][countries[i]] = amp_freqs[amp][gt]

    if max([len(comb_amp_freqs[amp]) for amp in comb_amp_freqs]) == 0:
        max_gts = 1
    else:
        max_gts = max([len(comb_amp_freqs[amp]) for amp in comb_amp_freqs])

    # write tsv output
    header_dict = {}
    header_string = "Amplicon\tGenotype"
    for country in countries:
        header_dict.update({country + "_freq": {}})
        header_string = header_string + "\t" + country
    header_string = header_string + "\n"
    header_dict_key_list = list(header_dict.keys())

    with open(''.join([output_prefix, "_comb_amp_freqs.txt"]), "wt") as w:
        w.write(header_string)
        for amp in sorted(comb_amp_freqs):
            to_write = str(amp)
            for gt in sorted(comb_amp_freqs[amp]):
                to_write = to_write + "\t" + str(gt)
                for key in header_dict_key_list:
                    country_name = key[:-5]
                    header_dict.update(
                        {key: comb_amp_freqs[amp][gt].get(country_name, 0)})
                    to_write = to_write + "\t" + str(header_dict[key])
                to_write = to_write + "\n"
                w.write(to_write)

    # write output: sparse matrix
    if output_matrix:
        with open(''.join([output_prefix,
                          "_comb_amp_freqs_sparse_matrix.txt"]), "wt") as w:
            w.write("Amplicon\tCountry\t{}\n".format(
                    "\t".join(["Genotype.{}".format(
                        i+1) for i in range(max_gts)])))
            for amp in sorted(comb_amp_freqs):
                for country in countries:
                    w.write("{}\t{}\t".format(amp, country))
                    w.write("\t".join(['{:g}'.format(
                        comb_amp_freqs[amp][gt].get(country, 0)) for gt in
                            sorted(comb_amp_freqs[amp])]))
                    fill_zeros = max_gts - len(comb_amp_freqs[amp])

                    if fill_zeros > 0:
                        w.write("\t{}".format("\t".join(["0"]*fill_zeros)))
                    w.write("\n")


country_dict = read_all_tables(country_paths)
get_freqs_per_country_one_panel(coords, country_dict, "test_out",
                                output_matrix=True)
