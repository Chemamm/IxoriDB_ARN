import random
import os
from IxoriDB.settings import BLAST_DIR, UPLOAD_DIR
from Bio import SeqIO


class Blast:
    def __init__(self,file):
        with open(file) as blast:
            self.lines = list(blast)
        self.hits = self.getHits()

    class Hit:
        def __init__(self, line):
            self.row = line.split("\t")
            self.ID = self.row[0]
            self.hitID = self.row[1]
            self.identity = self.row[2]
            self.length = self.row[3]
            self.mismatch = self.row[4]
            self.gapopen = self.row[5]
            self.qstart = self.row[6]
            self.qend = self.row[7]
            self.sstart = self.row[8]
            self.send = self.row[9]
            self.evalue = self.row[10]
            self.bitscore = self.row[11]

    def getHits(self):
        hits = []
        for line in self.lines:
            hit = Blast.Hit(line)
            hits.append(hit)
        return hits


def blastCommandLine(seq, type_of_seq, evalue, word_size, max_target_seq):
    work = random.randint(0, 10000000)
    if seq:
        with open(UPLOAD_DIR + "tmp/%i.fa" %work, "wt") as fa_tmp:
            fa_tmp.write(seq)

    fasta_check = is_fasta(UPLOAD_DIR + "tmp/%i.fa" %work)
    if fasta_check:
        if type_of_seq == "Peptide":
            cmd = "%s/blastp -db blast_db/Trinity_MGSG.pep -query %stmp/%i.fa -max_target_seqs %i -word_size " \
                  "%i -evalue %f -out %stmp/%i.blast -outfmt 6" %(BLAST_DIR, UPLOAD_DIR, work, int(max_target_seq), int(word_size), float(evalue), UPLOAD_DIR, work)
            os.system(cmd)

        elif type_of_seq == "Nucleotide":
            os.system("%s/blastx -db Trinity_MGSG.pep -query %stmp/%i.fa -max_target_seqs %i -word_size "
                      "%i -evalue %f -out %stmp/%i.blast -outfmt 6"%(BLAST_DIR, UPLOAD_DIR, work, int(max_target_seq), int(word_size), float(evalue), UPLOAD_DIR, work))

    return work, fasta_check


def is_fasta(filename):
    with open(filename, "r") as handle:
        fasta = SeqIO.parse(handle, "fasta")
        return any(fasta)  # False when `fasta` is empty, i.e. wasn't a FASTA file


def intersection(lst1, lst2):
    return list(set(lst1) & set(lst2))


