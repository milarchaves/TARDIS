import sys
import os
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

for out_dir in os.listdir(os.getcwd()):
    if os.path.isdir(out_dir):
        os.getcwd(out_dir)
        sequence_data = open(out_dir + '_CHOIR_MonomerSequence.fasta').read()
        result_handle = NCBIWWW.qblast("blastp", "nr", sequence_data, entrez_query='txid9606[ORGN]')
        with open('results.xml', 'a') as save_file:
            blast_results = result_handle.read()
            save_file.write(blast_results)

        E_VALUE_THRESH = 1e-20
        for record in NCBIXML.parse(open("results.xml")):
            if record.alignments:
                print("\n")
                print("query: %s" % record.query[:100])
                for align in record.alignments:
                    for hsp in align.hsps:
                        if hsp.expect < E_VALUE_THRESH:
                            print("match: %s " % align.title[:100])
