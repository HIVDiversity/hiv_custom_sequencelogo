import argparse
import sys
import pathlib
from itertools import groupby
import collections
import subprocess
# external libraries
import pandas as pd


__author__ = 'Colin Anthony'


def py3_fasta_iter(fasta_name):
    """
    modified from Brent Pedersen: https://www.biostars.org/p/710/#1412
    given a fasta file. yield tuples of header, sequence
    """
    fh = open(str(fasta_name), 'r')
    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
    for header in faiter:
        # drop the ">"
        header_str = header.__next__()[1:].strip()
        # join all sequence lines to one.
        seq = "".join(s.strip() for s in faiter.__next__())
        yield (header_str, seq)


def fasta_to_dct(file_name):
    """
    :param file_name: The fasta formatted file to read from.
    :return: a dictionary of the contents of the file name given. Dictionary in the format:
    {sequence_id: sequence_string, id_2: sequence_2, etc.}
    """
    dct = collections.defaultdict(str)
    my_gen = py3_fasta_iter(file_name)
    for i, (k, v) in enumerate(my_gen):
        # resolve for duplicate sequence headers
        new_key = k.replace(" ", "_") + str(i).zfill(4)
        dct[new_key] = v.upper()

    return dct


def gethxb2(dict):
    """
    :param dict: a dictionary of your aligned input sequences. Must contain HXB2, with HXB2 in the header
    :return: the HXB2 sequence as a string
    """
    found = False
    hxb2_seq = None
    hxb2_key = None
    for k in dict.keys():
        if "HXB2" in k.upper():
            found = True
            hxb2_key = k
            hxb2_seq = dict[k]
            print(f"Found hxb2 ref. seq. Its full name is: {hxb2_key}")
            break
    if not found:
        print("We could not find a sequence with 'HXB2' in its name. "
              "Please make sure you have an HXB2 ref. seq. included")
    return str(hxb2_key), str(hxb2_seq)


def posnumcalc(hxb2seq, start):
    pos_num = []
    n = start
    s = 0.01
    m = len(hxb2seq) - 1
    for i, resi in enumerate(hxb2seq):
        if i == 0 and resi == '-':
            print("Can't start numbering. HXB2 starts with a gap. Use a longer HXB2 sequence for the numbering")
        if i != m:
            if resi != '-' and hxb2seq[i+1] != '-':
                pos_num.append(n)
                n += 1
            elif resi != '-' and hxb2seq[i+1] == '-':
                g = n
                pos_num.append(g)
            elif resi == '-' and hxb2seq[i+1] == '-':
                g = n + s
                pos_num.append(g)
                s += 0.01
            elif resi == '-' and hxb2seq[i+1] != '-':
                g = n + s
                pos_num.append(g)
                s = 0.01
                n += 1
        else:
            if resi != '-':
                pos_num.append(n)
            elif resi == '-':
                g = n + s
                pos_num.append(g)
    return pos_num


def main(path_to_weblogo, path, sites_of_interest, fasta_file, start, color_scheme, title, x_label):
    """
    make a sequence logo from an aligned fasta file (with HXB2) and a csv file with a list of sites
    :param path_to_weblogo: (str) The path to where weblogo was installed, eg "~/anaconda3/bin/weblogo"
    :param path: (str) path to where the output should be written
    :param sites_of_interest: (str) csv file. column heading 'sites' with the positions to use in this column
    :param fasta_file: (str) the path and name of the fasta file with your virus sequences (include hxb2)
    :param start: (int) hxb2 start position of alignment
    :param color_scheme: (str) the name of the color scheme to use for the logo
    :param title: (str) the title for the sequence logo
    :param x_label: (str) the x-axis label for the sequence logo
    :returns: a fasta file with the sites used to make the logo, a png file of the sequence logo
    """
    weblogo_path = pathlib.Path(path_to_weblogo).absolute()
    path = pathlib.Path(path).absolute()
    fasta_file = pathlib.Path(fasta_file).absolute()
    sites_of_interest = pathlib.Path(sites_of_interest).absolute()
    name = fasta_file.stem + f"_for_logo.fasta"
    fasta_outfile = pathlib.Path(outpath, name)
    logo_outfile = pathlib.Path(outpath, name.replace(".fasta", ""))
    if path_to_weblogo == "":
        print("no path specified for weblogo, will assume that it is in your $PATH\n")
        weblogo_path = ""

    elif not weblogo_path.is_dir():
        print("path to weblogo executalbe not found\n")
        sys.exit("exiting")

    if not path.is_dir():
        print("path to project folder is not a directory\n")
        sys.exit("exiting")

    if not fasta_file.is_file():
        print("fasta file does not exist")
        sys.exit("exiting")

    if not sites_of_interest.is_file():
        print("sites of interest csv file does not exist")
        sys.exit("exiting")

    # read in CSV with sites to extract in one column
    sites_of_interest_df = pd.read_csv(sites_of_interest, sep=None, engine='python')

    # get sites to extract into a list
    try:
        sites_list = sorted((sites_of_interest_df["sites"].dropna().unique().tolist()))
    except KeyError as e:
        print("csv file did not contain the heading 'sites'\n")
        print(e)
        sys.exit("exiting")

    # read in fasta file with hxb2
    virus_d = fasta_to_dct(fasta_file)
    try:
        hxb2_name, hxb2_seq = gethxb2(virus_d)
    except KeyError:
        print("HXB2 not found in alignment. Make sure it is in the alignment with 'HXB2' in the name")
        sys.exit("exiting")
    else:
        del virus_d[hxb2_name]

    # get hxb2 numbering
    pos_num = posnumcalc(hxb2_seq, start)

    # for each sequence, collect the required sites
    collected_sites_d = collections.defaultdict(list)
    for seq_name, seq in virus_d.items():
        for site in sites_list:
            idx = pos_num.index(site)
            collected_sites_d[seq_name].append(seq[idx])

    # write out the collected sites as a new fasta file
    with open(fasta_outfile, 'w') as fh:
        for collected_name, collected_sites_list in collected_sites_d.items():
            sites_to_use = "".join(collected_sites_list)
            fh.write(f">{collected_name}\n{sites_to_use}\n")

    sites = ','.join([str(int(x)) for x in sites_list])
    path_to_weblogo_install = pathlib.Path(weblogo_path, "weblogo")

    web_logo_cmd = f"{path_to_weblogo_install} --fin {fasta_outfile} --datatype fasta " \
        f"--fout {logo_outfile}.png --format png --sequence-type protein --units probability --size large " \
        f"--title '{title}' --xlabel '{x_label}' --annotate '{sites}' --ylabel 'Frequency' " \
        f"--color-scheme {color_scheme} --fontsize 12  --title-fontsize 14 --text-font ArialMT " \
        f"--logo-font ArialMT --title-font ArialMT --resolution 600"
    try:
        subprocess.call(web_logo_cmd, shell=True)
    except subprocess.CalledProcessError as e:
        print("error with the run: ", e)
        sys.exit("exiting")

    print("done")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Create a sequence logo from an aligned fasta file and a cvs file '
                                                 'with a list of sites. HXB2 must be in the sequence alignment',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-p', '--path_to_weblogo', default="", type=str, required=False,
                        help='The path to where weblogo was installed, eg "~/anaconda3/bin/weblogo". '
                             'Default assumes weblogo is in your $PATH')
    parser.add_argument('-o', '--outpath', type=str, required=True,
                        help='The path to the project folder, where the folders will be created')
    parser.add_argument('-in', '--sites_of_interest', type=str, default=None, required=False,
                        help='The path and file name of the csv file with neut data from lanl')
    parser.add_argument('-f', '--fasta_file', required=False, type=str,
                        help='The path and name of the aligned and translated fasta file of virus sequences')
    parser.add_argument('-s', '--start', default=1, type=int,
                        help='the HXB2 gene start position of your alignment', required=False)
    parser.add_argument('-c', '--color_scheme', default="chemistry", choices=["chemistry", "charge", "monochrome",
                                                                              "classic", "hydrophobicity",
                                                                              "base pairing"], type=str,
                        help='The color scheme for the sequence logo', required=False)
    parser.add_argument('-t', '--title', default="", type=str,
                        help='the title for the sequence logo. eg: "Global HIV ENV C3"', required=False)
    parser.add_argument('-x', '--x_label', default="Position (HXB2 numbering)", type=str,
                        help='the x-axis label for the sequence logo. eg: "N332-V3 Supersite"', required=False)

    args = parser.parse_args()
    path_to_weblogo = args.path_to_weblogo
    outpath = args.outpath
    sites_of_interest = args.sites_of_interest
    fasta_file = args.fasta_file
    start = args.start
    color_scheme = args.color_scheme
    title = args.title
    x_label = args.x_label

    main(path_to_weblogo, outpath, sites_of_interest, fasta_file, start, color_scheme, title, x_label)
