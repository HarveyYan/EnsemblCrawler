import threading
import csv
import os
import urllib
import urllib.request as request
from bs4 import BeautifulSoup
import warnings
import numpy as np

OUTPUT_CDNA = './genes/cDNA/'
OUTPUT_UTR3 = './genes/utr3/'
if not os.path.exists(OUTPUT_CDNA):
    os.makedirs(OUTPUT_CDNA)
if not os.path.exists(OUTPUT_UTR3):
    os.makedirs(OUTPUT_UTR3)
supp_file_path = './Supplemental_File_3.tsv'
base_link = 'https://www.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g='
download_link_cdna = 'https://www.ensembl.org/Homo_sapiens/Export/Output/Transcript?db=core;flank3_display=0;flank5_display=0;g=g_id;output=fasta;r_value;strand=feature;t=t_value;param=cdna;genomic=off;_format=Text'
download_link_3utr = 'https://www.ensembl.org/Homo_sapiens/Export/Output/Transcript?db=core;flank3_display=0;flank5_display=0;g=g_id;output=fasta;r_value;strand=feature;t=t_value;param=utr3;genomic=off;_format=Text'
output = open("./cDNA", "w")
cdna_missing = open("./cDNA_missing", "w")
output = open("./utr3", "w")
utr3_missing = open("./utr3_missing", "w")
file = open('./example.tsv')
# 0: ensembl id
# 1: unique identifier
# 4: longest isoform
# 5 ~ 12: distributions

def crawler(gene):
    # start download from website
    id = gene[0]
    longest_isoform = gene[1]
    res = request.urlopen(base_link + id)
    r_value = res.geturl().split(';')[-1]
    t_value = None
    site = BeautifulSoup(res.read())
    if site.find('table', {'id': 'transcripts_table'}) is None:
        # gene not in ensembl database
        print(id + " not found")
        cdna_missing.write(id)
        cdna_missing.write('\n')
        cdna_missing.flush()
        utr3_missing.write(id)
        utr3_missing.write('\n')
        utr3_missing.flush()
        return
    transcripts = site.find('table', {'id': 'transcripts_table'}).find_all('tr')[1:]
    for t in transcripts:
        if int(t.find_all('td')[2].getText()) == longest_isoform:
            t_value = t.find('a').getText().split('.')[0]
            break
    # supplemental file isoform length not found
    if t_value == None:
        longest = 0
        for t in transcripts:
            if int(t.find_all('td')[2].getText()) > longest:
                longest = int(t.find_all('td')[2].getText())
                t_value = t.find('a').getText().split('.')[0]

    # now assemble r_value and t_value

    # for cDNA
    link = download_link_cdna.replace('g_id', id)
    link = link.replace('r_value', r_value)
    link = link.replace('t_value', t_value)
    # now go for the download!
    count = 0
    while True:
        try:
            seq = request.urlopen(link).read()
            output = open(OUTPUT_CDNA+id,'w')
            print(id, "cDNA downloaded")
            seqs = seq.decode('utf-8').split("\r\n")
            output.write('>gene_id:{0} trans_id:{1}\n'.format(id, seqs[0][1:]))
            output.writelines(seqs[1:])
            output.flush()
            output.close()
        except (urllib.error.HTTPError, urllib.error.URLError) as e:
            print(id," error")
            count += 1
            if count > 5:
                cdna_missing.write(id)
                cdna_missing.write('\n')
                cdna_missing.flush()
                break
            else:
                continue

    # for utr_3
    link = download_link_3utr.replace('g_id', id)
    link = link.replace('r_value', r_value)
    link = link.replace('t_value', t_value)
    # now go for the download!
    count = 0
    while True:
        try:
            utr3 = request.urlopen(link).read()
            output = open(OUTPUT_UTR3 + id, 'w')
            print(utr3, "cDNA downloaded")
            utr = utr3.decode('utf-8').split("\r\n")
            output.write('>gene_id:{0} trans_id:{1}\n'.format(id, utr[0][1:]))
            output.writelines(utr[1:])
            output.flush()
        except (urllib.error.HTTPError, urllib.error.URLError) as e:
            print(id, " error")
            count += 1
            if count > 5:
                utr3_missing.write(id)
                utr3_missing.write('\n')
                utr3_missing.flush()
                break
            else:
                continue

def download_all():
    genes = []
    reader = csv.reader(file, delimiter='\t')
    next(reader)

    for line in reader:
        id = line[0]
        longest_isoform = int(line[4])
        genes.append([id, longest_isoform])
    genes = sorted(genes)

    threads = [threading.Thread(target=crawler, args=(gene,)) for gene in genes]
    for thread in threads:
        thread.start()
    for thread in threads:
        thread.join()


if __name__ == "__main__":
    download_all()