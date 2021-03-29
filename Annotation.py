#!/usr/bin/env/python3
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pylab import *
import os


class CompleteExcel:
    def __init__(self, file):
        with open(file) as ex:
            self.lines = list(ex)
        self.IDs, self.transcriptIDs = self.getIDs()


    class Line:
        def __init__(self, line):
            self.row = line.split("\t")
            self.cdsID = self.row[0]
            self.transcriptID = self.cdsID.split(".")[0]
            self.transcript_seq = self.row[1]
            self.cds_seq = self.row[2]
            self.seq = self.row[3]
            self.family = self.row[4]
            self.hitID = self.row[5]
            self.identity = self.row[6]
            self.coverage = self.row[7]
            self.bitscore = self.row[8]
            self.description = self.row[9]
            self.db = self.row[10]
            self.GOterm = self.row[11]
            self.GOterm_des = self.row[12]
            self.KW = self.row[13]
            self.KW_des = self.row[14]
            self.ipro_id = self.row[15]
            self.ipro_des = self.row[16]
            self.ipro_goterm = self.row[17]
            self.ipro_goterm_des = self.row[18]
            self.merged = self.row[19]
            self.merged_des = self.row[20]
            self.sigp = self.row[21]
            self.tm_domain_all = self.row[22]
            self.tm_domain = self.row[23]
            self.classification = self.row[24]
            self.eclass = self.row[25]
            self.SG_FC = self.row[26]
            self.MG_FC = self.row[27]
            self.total_FPKM = self.row[28]
            self.DE_str = self.row[29]
            self.DE_FC = self.row[30]
            self.DE_PPDE = self.row[31]
            self.evalues = self.row[32::]

    def getIDs(self):
        IDs = []
        transcriptIDs = []
        for n in range(1,len(self.lines)):
            line = CompleteExcel.Line(self.lines[n])
            IDs.append(line.cdsID)
            transcriptIDs.append(line.transcriptID)
        return IDs, transcriptIDs

    def splitSecreted(self, output_label):
        with open("%s_secreted.fa" %output_label, "wt") as secreted_out:
            with open("%s_non-secreted.fa" %output_label, "wt") as nosecreted_out:
                for n in range(1,len(self.lines)):
                    line = CompleteExcel.Line(self.lines[n])
                    if line.classification == "Putatively Secreted, Annotated" or line.classification == "Putatively Secreted, not Annotated":
                        secreted_out.write(">%s\n%s\n" %(line.cdsID,line.seq))
                    else:
                        nosecreted_out.write(">%s\n%s\n" %(line.cdsID,line.seq))


    def getstats(self, output):
        GO = 0
        notGO = 0
        GO_SA = 0
        GO_SNA = 0
        GO_NSA = 0
        GO_NSNA = 0
        notGO_SA = 0
        notGO_SNA = 0
        notGO_NSA = 0
        notGO_NSNA = 0
        GO_interpro = 0
        notGO_interpro = 0
        GO_SA_SG = 0
        GO_SA_MG = 0
        GO_SA_N = 0
        GO_SNA_SG = 0
        GO_SNA_MG = 0
        GO_SNA_N = 0
        GO_NSA_SG = 0
        GO_NSA_MG = 0
        GO_NSA_N = 0
        GO_NSNA_SG = 0
        GO_NSNA_MG = 0
        GO_NSNA_N = 0
        notGO_SA_SG = 0
        notGO_SA_MG = 0
        notGO_SA_N = 0
        notGO_SNA_SG = 0
        notGO_SNA_MG = 0
        notGO_SNA_N = 0
        notGO_NSA_SG = 0
        notGO_NSA_MG = 0
        notGO_NSA_N = 0
        notGO_NSNA_SG = 0
        notGO_NSNA_MG = 0
        notGO_NSNA_N = 0
        for n in range(1, len(self.lines)):
            line = CompleteExcel.Line(self.lines[n])
            if line.GOterm != "None":
                GO += 1
                if line.classification == "Putatively Secreted, Annotated":
                    GO_SA += 1
                    if line.eclass == "SG specific":
                        GO_SA_SG += 1
                    elif line.eclass == "MG specific":
                        GO_SA_MG += 1
                    elif line.eclass == "Neutral":
                        GO_SA_N += 1
                elif line.classification == "Putatively Secreted, not Annotated":
                    GO_SNA += 1
                    if line.eclass == "SG specific":
                        GO_SNA_SG += 1
                    elif line.eclass == "MG specific":
                        GO_SNA_MG += 1
                    elif line.eclass == "Neutral":
                        GO_SNA_N += 1
                elif line.classification == "Putatively Non-Secreted, Annotated":
                    GO_NSA += 1
                    if line.eclass == "SG specific":
                        GO_NSA_SG += 1
                    elif line.eclass == "MG specific":
                        GO_NSA_MG += 1
                    elif line.eclass == "Neutral":
                        GO_NSA_N += 1
                elif line.classification == "Putatively Non-Secreted, not Annotated":
                    GO_NSNA += 1
                    if line.eclass == "SG specific":
                        GO_NSNA_SG += 1
                    elif line.eclass == "MG specific":
                        GO_NSNA_MG += 1
                    elif line.eclass == "Neutral":
                        GO_NSNA_N += 1
                if line.ipro_id != "None":
                    GO_interpro += 1

            else:
                notGO += 1
                if line.classification == "Putatively Secreted, Annotated":
                    notGO_SA += 1
                    if line.eclass == "SG specific":
                        notGO_SA_SG += 1
                    elif line.eclass == "MG specific":
                        notGO_SA_MG += 1
                    elif line.eclass == "Neutral":
                        notGO_SA_N += 1
                elif line.classification == "Putatively Secreted, not Annotated":
                    notGO_SNA += 1
                    if line.eclass == "SG specific":
                        notGO_SNA_SG += 1
                    elif line.eclass == "MG specific":
                        notGO_SNA_MG += 1
                    elif line.eclass == "Neutral":
                        notGO_SNA_N += 1
                elif line.classification == "Putatively Non-Secreted, Annotated":
                    notGO_NSA += 1
                    if line.eclass == "SG specific":
                        notGO_NSA_SG += 1
                    elif line.eclass == "MG specific":
                        notGO_NSA_MG += 1
                    elif line.eclass == "Neutral":
                        notGO_NSA_N += 1
                elif line.classification == "Putatively Non-Secreted, not Annotated":
                    notGO_NSNA += 1
                    if line.eclass == "SG specific":
                        notGO_NSNA_SG += 1
                    elif line.eclass == "MG specific":
                        notGO_NSNA_MG += 1
                    elif line.eclass == "Neutral":
                        notGO_NSNA_N += 1
                if line.ipro_id != "None":
                    notGO_interpro += 1

        with open("%s_uniprotGO.stats" % output, "wt") as out:
            out.write("Total of CDS with GO annotation: %s. (%s of them have interpro domains)\nFrom them:\n\t%s Putatively Secreted, Annotated\n\t\t%s SG specific\n\t\t%s MG specific\n\t\t%s Neutral\n\t%s Putatively Secreted, not Annotated\n\t\t%s SG specific\n\t\t%s MG specific\n\t\t%s Neutral\n\t%s Putatively Non-Secreted, Annotated\n\t\t%s SG specific\n\t\t%s MG specific\n\t\t%s Neutral\n\t%s Putatively Non-Secreted, not Annotated\n\t\t%s SG specific\n\t\t%s MG specific\n\t\t%s Neutral\n\n" % (GO, GO_interpro, GO_SA, GO_SA_SG, GO_SA_MG, GO_SA_N, GO_SNA, GO_SNA_SG, GO_SNA_MG, GO_SNA_N, GO_NSA, GO_NSA_SG, GO_NSA_MG, GO_NSA_N, GO_NSNA, GO_NSNA_SG, GO_NSNA_MG, GO_NSNA_N))
            out.write("Total of CDS without GO annotation: %s. (%s of them have interpro domains)\nFrom them:\n\t%s Putatively Secreted, Annotated\n\t\t%s SG specific\n\t\t%s MG specific\n\t\t%s Neutral\n\t%s Putatively Secreted, not Annotated\n\t\t%s SG specific\n\t\t%s MG specific\n\t\t%s Neutral\n\t%s Putatively Non-Secreted, Annotated\n\t\t%s SG specific\n\t\t%s MG specific\n\t\t%s Neutral\n\t%s Putatively Non-Secreted, not Annotated\n\t\t%s SG specific\n\t\t%s MG specific\n\t\t%s Neutral\n\n" % (notGO, notGO_interpro, notGO_SA, notGO_SA_SG, notGO_SA_MG, notGO_SA_N, notGO_SNA, notGO_SNA_SG, notGO_SNA_MG, notGO_SNA_N, notGO_NSA, notGO_NSA_SG, notGO_NSA_MG, notGO_NSA_N, notGO_NSNA, notGO_NSNA_SG, notGO_NSNA_MG, notGO_NSNA_N))

        #Interpro
        ipro = 0
        notipro = 0
        ipro_SA = 0
        ipro_SNA = 0
        ipro_NSA = 0
        ipro_NSNA = 0
        notipro_SA = 0
        notipro_SNA = 0
        notipro_NSA = 0
        notipro_NSNA = 0
        ipro_GO = 0
        notipro_GO = 0
        ipro_SA_SG = 0
        ipro_SA_MG = 0
        ipro_SA_N = 0
        ipro_SNA_SG = 0
        ipro_SNA_MG = 0
        ipro_SNA_N = 0
        ipro_NSA_SG = 0
        ipro_NSA_MG = 0
        ipro_NSA_N = 0
        ipro_NSNA_SG = 0
        ipro_NSNA_MG = 0
        ipro_NSNA_N = 0
        notipro_SA_SG = 0
        notipro_SA_MG = 0
        notipro_SA_N = 0
        notipro_SNA_SG = 0
        notipro_SNA_MG = 0
        notipro_SNA_N = 0
        notipro_NSA_SG = 0
        notipro_NSA_MG = 0
        notipro_NSA_N = 0
        notipro_NSNA_SG = 0
        notipro_NSNA_MG = 0
        notipro_NSNA_N = 0
        for n in range(1, len(self.lines)):
            line = CompleteExcel.Line(self.lines[n])
            if line.ipro_id != "None":
                ipro += 1
                if line.classification == "Putatively Secreted, Annotated":
                    ipro_SA += 1
                    if line.eclass == "SG specific":
                        ipro_SA_SG += 1
                    elif line.eclass == "MG specific":
                        ipro_SA_MG += 1
                    elif line.eclass == "Neutral":
                        ipro_SA_N += 1
                elif line.classification == "Putatively Secreted, not Annotated":
                    ipro_SNA += 1
                    if line.eclass == "SG specific":
                        ipro_SNA_SG += 1
                    elif line.eclass == "MG specific":
                        ipro_SNA_MG += 1
                    elif line.eclass == "Neutral":
                        ipro_SNA_N += 1
                elif line.classification == "Putatively Non-Secreted, Annotated":
                    ipro_NSA += 1
                    if line.eclass == "SG specific":
                        ipro_NSA_SG += 1
                    elif line.eclass == "MG specific":
                        ipro_NSA_MG += 1
                    elif line.eclass == "Neutral":
                        ipro_NSA_N += 1
                elif line.classification == "Putatively Non-Secreted, not Annotated":
                    ipro_NSNA += 1
                    if line.eclass == "SG specific":
                        ipro_NSNA_SG += 1
                    elif line.eclass == "MG specific":
                        ipro_NSNA_MG += 1
                    elif line.eclass == "Neutral":
                        ipro_NSNA_N += 1
                if line.GOterm != "None":
                    ipro_GO += 1

            else:
                notipro += 1
                if line.classification == "Putatively Secreted, Annotated":
                    notipro_SA += 1
                    if line.eclass == "SG specific":
                        notipro_SA_SG += 1
                    elif line.eclass == "MG specific":
                        notipro_SA_MG += 1
                    elif line.eclass == "Neutral":
                        notipro_SA_N += 1
                elif line.classification == "Putatively Secreted, not Annotated":
                    notipro_SNA += 1
                    if line.eclass == "SG specific":
                        notipro_SNA_SG += 1
                    elif line.eclass == "MG specific":
                        notipro_SNA_MG += 1
                    elif line.eclass == "Neutral":
                        notipro_SNA_N += 1
                elif line.classification == "Putatively Non-Secreted, Annotated":
                    notipro_NSA += 1
                    if line.eclass == "SG specific":
                        notipro_NSA_SG += 1
                    elif line.eclass == "MG specific":
                        notipro_NSA_MG += 1
                    elif line.eclass == "Neutral":
                        notipro_NSA_N += 1
                elif line.classification == "Putatively Non-Secreted, not Annotated":
                    notipro_NSNA += 1
                    if line.eclass == "SG specific":
                        notipro_NSNA_SG += 1
                    elif line.eclass == "MG specific":
                        notipro_NSNA_MG += 1
                    elif line.eclass == "Neutral":
                        notipro_NSNA_N += 1
                if line.GOterm != "None":
                    notipro_GO += 1

        with open("%s_interpro.stats" % output, "wt") as out:
            out.write("Total of CDS with ipro annotation: %s. (%s of them have GO terms)\nFrom them:\n\t%s Putatively Secreted, Annotated\n\t\t%s SG specific\n\t\t%s MG specific\n\t\t%s Neutral\n\t%s Putatively Secreted, not Annotated\n\t\t%s SG specific\n\t\t%s MG specific\n\t\t%s Neutral\n\t%s Putatively Non-Secreted, Annotated\n\t\t%s SG specific\n\t\t%s MG specific\n\t\t%s Neutral\n\t%s Putatively Non-Secreted, not Annotated\n\t\t%s SG specific\n\t\t%s MG specific\n\t\t%s Neutral\n\n" % (ipro, ipro_GO, ipro_SA, ipro_SA_SG, ipro_SA_MG, ipro_SA_N, ipro_SNA, ipro_SNA_SG, ipro_SNA_MG, ipro_SNA_N, ipro_NSA, ipro_NSA_SG, ipro_NSA_MG, ipro_NSA_N, ipro_NSNA, ipro_NSNA_SG, ipro_NSNA_MG, ipro_NSNA_N))
            out.write("Total of CDS without ipro annotation: %s. (%s of them have GO terms)\nFrom them:\n\t%s Putatively Secreted, Annotated\n\t\t%s SG specific\n\t\t%s MG specific\n\t\t%s Neutral\n\t%s Putatively Secreted, not Annotated\n\t\t%s SG specific\n\t\t%s MG specific\n\t\t%s Neutral\n\t%s Putatively Non-Secreted, Annotated\n\t\t%s SG specific\n\t\t%s MG specific\n\t\t%s Neutral\n\t%s Putatively Non-Secreted, not Annotated\n\t\t%s SG specific\n\t\t%s MG specific\n\t\t%s Neutral\n\n" % (notipro, notipro_GO, notipro_SA, notipro_SA_SG, notipro_SA_MG, notipro_SA_N, notipro_SNA, notipro_SNA_SG, notipro_SNA_MG, notipro_SNA_N, notipro_NSA, notipro_NSA_SG, notipro_NSA_MG, notipro_NSA_N, notipro_NSNA, notipro_NSNA_SG, notipro_NSNA_MG, notipro_NSNA_N))
    #Merged GO
        merged = 0
        notmerged = 0
        merged_SA = 0
        merged_SNA = 0
        merged_NSA = 0
        merged_NSNA = 0
        notmerged_SA = 0
        notmerged_SNA = 0
        notmerged_NSA = 0
        notmerged_NSNA = 0
        merged_interpro = 0
        notmerged_interpro = 0
        merged_SA_SG = 0
        merged_SA_MG = 0
        merged_SA_N = 0
        merged_SNA_SG = 0
        merged_SNA_MG = 0
        merged_SNA_N = 0
        merged_NSA_SG = 0
        merged_NSA_MG = 0
        merged_NSA_N = 0
        merged_NSNA_SG = 0
        merged_NSNA_MG = 0
        merged_NSNA_N = 0
        notmerged_SA_SG = 0
        notmerged_SA_MG = 0
        notmerged_SA_N = 0
        notmerged_SNA_SG = 0
        notmerged_SNA_MG = 0
        notmerged_SNA_N = 0
        notmerged_NSA_SG = 0
        notmerged_NSA_MG = 0
        notmerged_NSA_N = 0
        notmerged_NSNA_SG = 0
        notmerged_NSNA_MG = 0
        notmerged_NSNA_N = 0
        for n in range(1, len(self.lines)):
            line = CompleteExcel.Line(self.lines[n])
            if line.merged != "None":
                merged += 1
                if line.classification == "Putatively Secreted, Annotated":
                    merged_SA += 1
                    if line.eclass == "SG specific":
                        merged_SA_SG += 1
                    elif line.eclass == "MG specific":
                        merged_SA_MG += 1
                    elif line.eclass == "Neutral":
                        merged_SA_N += 1
                elif line.classification == "Putatively Secreted, not Annotated":
                    merged_SNA += 1
                    if line.eclass == "SG specific":
                        merged_SNA_SG += 1
                    elif line.eclass == "MG specific":
                        merged_SNA_MG += 1
                    elif line.eclass == "Neutral":
                        merged_SNA_N += 1
                elif line.classification == "Putatively Non-Secreted, Annotated":
                    merged_NSA += 1
                    if line.eclass == "SG specific":
                        merged_NSA_SG += 1
                    elif line.eclass == "MG specific":
                        merged_NSA_MG += 1
                    elif line.eclass == "Neutral":
                        merged_NSA_N += 1
                elif line.classification == "Putatively Non-Secreted, not Annotated":
                    merged_NSNA += 1
                    if line.eclass == "SG specific":
                        merged_NSNA_SG += 1
                    elif line.eclass == "MG specific":
                        merged_NSNA_MG += 1
                    elif line.eclass == "Neutral":
                        merged_NSNA_N += 1
                if line.ipro_id != "None":
                    merged_interpro += 1

            else:
                notmerged += 1
                if line.classification == "Putatively Secreted, Annotated":
                    notmerged_SA += 1
                    if line.eclass == "SG specific":
                        notmerged_SA_SG += 1
                    elif line.eclass == "MG specific":
                        notmerged_SA_MG += 1
                    elif line.eclass == "Neutral":
                        notmerged_SA_N += 1
                elif line.classification == "Putatively Secreted, not Annotated":
                    notmerged_SNA += 1
                    if line.eclass == "SG specific":
                        notmerged_SNA_SG += 1
                    elif line.eclass == "MG specific":
                        notmerged_SNA_MG += 1
                    elif line.eclass == "Neutral":
                        notmerged_SNA_N += 1
                elif line.classification == "Putatively Non-Secreted, Annotated":
                    notmerged_NSA += 1
                    if line.eclass == "SG specific":
                        notmerged_NSA_SG += 1
                    elif line.eclass == "MG specific":
                        notmerged_NSA_MG += 1
                    elif line.eclass == "Neutral":
                        notmerged_NSA_N += 1
                elif line.classification == "Putatively Non-Secreted, not Annotated":
                    notmerged_NSNA += 1
                    if line.eclass == "SG specific":
                        notmerged_NSNA_SG += 1
                    elif line.eclass == "MG specific":
                        notmerged_NSNA_MG += 1
                    elif line.eclass == "Neutral":
                        notmerged_NSNA_N += 1
                if line.ipro_id != "None":
                    notmerged_interpro += 1

        with open("%s_merged.stats" % output, "wt") as out:
            out.write("Total of CDS with merged annotation: %s. (%s of them have interpro domains)\nFrom them:\n\t%s Putatively Secreted, Annotated\n\t\t%s SG specific\n\t\t%s MG specific\n\t\t%s Neutral\n\t%s Putatively Secreted, not Annotated\n\t\t%s SG specific\n\t\t%s MG specific\n\t\t%s Neutral\n\t%s Putatively Non-Secreted, Annotated\n\t\t%s SG specific\n\t\t%s MG specific\n\t\t%s Neutral\n\t%s Putatively Non-Secreted, not Annotated\n\t\t%s SG specific\n\t\t%s MG specific\n\t\t%s Neutral\n\n" % (merged, merged_interpro, merged_SA, merged_SA_SG, merged_SA_MG, merged_SA_N, merged_SNA, merged_SNA_SG, merged_SNA_MG, merged_SNA_N, merged_NSA, merged_NSA_SG, merged_NSA_MG, merged_NSA_N, merged_NSNA, merged_NSNA_SG, merged_NSNA_MG, merged_NSNA_N))
            out.write("Total of CDS without merged annotation: %s. (%s of them have interpro domains)\nFrom them:\n\t%s Putatively Secreted, Annotated\n\t\t%s SG specific\n\t\t%s MG specific\n\t\t%s Neutral\n\t%s Putatively Secreted, not Annotated\n\t\t%s SG specific\n\t\t%s MG specific\n\t\t%s Neutral\n\t%s Putatively Non-Secreted, Annotated\n\t\t%s SG specific\n\t\t%s MG specific\n\t\t%s Neutral\n\t%s Putatively Non-Secreted, not Annotated\n\t\t%s SG specific\n\t\t%s MG specific\n\t\t%s Neutral\n\n" % (notmerged, notmerged_interpro, notmerged_SA, notmerged_SA_SG, notmerged_SA_MG, notmerged_SA_N, notmerged_SNA, notmerged_SNA_SG, notmerged_SNA_MG, notmerged_SNA_N, notmerged_NSA, notmerged_NSA_SG, notmerged_NSA_MG, notmerged_NSA_N, notmerged_NSNA, notmerged_NSNA_SG, notmerged_NSNA_MG, notmerged_NSNA_N))

    def length_distribution(self, output):
        GOseq_dict = {}
        notGOseq_dict = {}
        iproseq_dict = {}
        notiproseq_dict = {}
        mergedseq_dict = {}
        notmergedseq_dict = {}
        for n in range(1, len(self.lines)):
            line = CompleteExcel.Line(self.lines[n])
            if line.GOterm != "None":
                GOseq_dict[line.cdsID] = line.seq
            else:
                notGOseq_dict[line.cdsID] = line.seq
            if line.ipro_id != "None":
                iproseq_dict[line.cdsID] = line.seq
            else:
                notiproseq_dict[line.cdsID] = line.seq
            if line.merged != "None":
                mergedseq_dict[line.cdsID] = line.seq
            else:
                notmergedseq_dict[line.cdsID] = line.seq

        GOlength_dict = getlength(GOseq_dict)
        notGOlength_dict = getlength(notGOseq_dict)
        iprolength_dict = getlength(iproseq_dict)
        notiprolength_dict = getlength(notiproseq_dict)
        mergedlength_dict = getlength(mergedseq_dict)
        notmergedlength_dict = getlength(notmergedseq_dict)

        sns.set_theme(style="darkgrid")
        subplot(1,2,1)
        ax = sns.histplot(list(GOlength_dict.values()))
        ax.set(xlabel='With Uniprot/B2GO GO terms')
        plt.xlim(0, 3000)
        plt.ylim(0, 900)
        subplot(1,2,2)
        ax = sns.histplot(list(notGOlength_dict.values()))
        ax.set(xlabel='Without Uniprot/B2GO GO terms')
        plt.xlim(0, 3000)
        plt.ylim(0, 900)
        fig = ax.get_figure()
        fig.savefig("%s_uniprotGO_length_distribution.pdf" %output)
        plt.close()

        sns.set_theme(style="darkgrid")
        subplot(1,2,1)
        ax = sns.histplot(list(iprolength_dict.values()))
        ax.set(xlabel='With Interpro Annotation')
        plt.xlim(0, 3000)
        plt.ylim(0, 900)
        subplot(1,2,2)
        ax = sns.histplot(list(notiprolength_dict.values()))
        ax.set(xlabel='Without Interpro Annotation')
        plt.xlim(0, 3000)
        plt.ylim(0, 900)
        fig = ax.get_figure()
        fig.savefig("%s_interpro_length_distribution.pdf" %output)
        plt.close()

        sns.set_theme(style="darkgrid")
        subplot(1,2,1)
        ax = sns.histplot(list(mergedlength_dict.values()))
        ax.set(xlabel='With GO terms')
        plt.xlim(0, 3000)
        plt.ylim(0, 900)
        subplot(1,2,2)
        ax = sns.histplot(list(notmergedlength_dict.values()))
        ax.set(xlabel='Without GO terms')
        plt.xlim(0, 3000)
        plt.ylim(0, 900)
        fig = ax.get_figure()
        fig.savefig("%s_mergedGO_length_distribution.pdf" %output)
        plt.close()

    def expression_distribution(self, output):
        self.evalues_secreted, self.evalues_nonsecreted = get_expression_dicts(self.lines)
        sns.set_theme(style="darkgrid")
        subplot(1,2,1)
        ax = sns.histplot(data=self.evalues_secreted, log_scale=True)
        ax.set(xlabel='Secreted FPKM')
        subplot(1,2,2)
        ax = sns.histplot(data=self.evalues_nonsecreted, log_scale=True)
        ax.set(xlabel='Non Secreted FPKM')
        fig = ax.get_figure()
        fig.savefig(output)
        plt.close()

    def getBackgrounds(self,output_label):
        GO_dict = {}
        KW_dict = {}
        IP_dict = {}
        family_dict = {}
        for n in range(1,len(self.lines)):
            line = CompleteExcel.Line(self.lines[n])
            #GO
            if line.merged != "None":
                GO_dict[line.transcriptID] = line.merged.replace("; ",",")
            #KW
            if line.KW != "None":
                KW_dict[line.transcriptID] = line.KW.replace("; ",",")
            #IP
            if line.ipro_id != "None":
                IP_dict[line.transcriptID] = line.ipro_id.replace("; ",",")
            #family
            if line.family != "None":
                family_dict[line.transcriptID] = line.family.replace("; ",",")

        with open("%s_GO_all.txt" %output_label,"wt") as out:
            for key in GO_dict:
                out.write("%s\t%s\n" %(key, GO_dict[key]))
        with open("%s_KW_all.txt" %output_label,"wt") as out:
            for key in KW_dict:
                out.write("%s\t%s\n" %(key, KW_dict[key]))
        with open("%s_IP_all.txt" %output_label,"wt") as out:
            for key in IP_dict:
                out.write("%s\t%s\n" %(key, IP_dict[key]))
        with open("%s_family_all.txt" %output_label,"wt") as out:
            for key in family_dict:
                out.write("%s\t%s\n" %(key, family_dict[key]))

    def getannotationdicts(self):
        go_dict = {}
        kw_dict = {}
        ip_dict = {}
        for n in range(1, len(self.lines)):
            line = CompleteExcel.Line(self.lines[n])
            go = line.GOterm.split("; ")
            go_des = line.GOterm_des.split("; ")
            ip_go = line.ipro_goterm.split("; ")
            ip_go_des = line.ipro_goterm_des.split("; ")
            ip = line.ipro_id.split("; ")
            ip_des = line.ipro_des.split(";")
            kw = line.KW.split("; ")
            kw_des = line.KW_des.split(";")
            for n in range(0,len(go)):
                if len(go) == len(go_des):
                    go_dict[go[n]] = go_des[n]
                else:
                    print(go)
                    print(go_des)
            for n in range(0, len(ip_go)):
                go_dict[ip_go[n]] = ip_go_des[n]
            for n in range(0, len(kw)):
                kw_dict[kw[n]] = kw_des[n]
            for n in range(0,len(ip)):
                ip_dict[ip[n]] = ip_des[n]
        return go_dict, kw_dict, ip_dict



class Ematrix:
    def __init__(self, file):
        with open(file) as ex:
            self.lines = list(ex)
        self.df = pd.read_csv(file, sep="\t", index_col=0)
        self.header, self.dict = self.matrix2dict()
        # Renaming
        rename_dict = dict()
        rename_dict["pre_2_SG_S235"] = "SG_unfed_1,5mg_pre2"
        rename_dict["pre_4_SG_S267"] = "SG_unfed_1,7mg_pre4"
        rename_dict["pre_6_SG_S236"] = "SG_unfed_1,4mg_pre6"
        rename_dict["pre_8_SG_S274"] = "SG_unfed_1,6mg_pre8"
        rename_dict["pre_10_SG_S237"] = "SG_unfed_1,4mg_pre10"
        rename_dict["1_12h_SG_S268"] = "SG_1_12h_1,5mg_1"
        rename_dict["105_12h_SG_S1"] = "SG_1_12h_2mg_105"
        rename_dict["5_12h_SG_S224"] = "SG_1_12h_1,8mg_5"
        rename_dict["101_12h_SG_S225"] = "SG_1_12h_1,7mg_101"
        rename_dict["103_12h_SG_S275"] = "SG_1_12h_2mg_103"
        rename_dict["7_24h_SG_S270"] = "SG_1_24h_1,5mg_7"
        rename_dict["9_24h_SG_S226"] = "SG_1_24h_1,7mg_9"
        rename_dict["109_24h_SG_S312"] = "SG_1_24h_2,3mg_109"
        rename_dict["106_24h_SG_S227"] = "SG_1_24h_1,6mg_106"
        rename_dict["8_24h_SG_S269"] = "SG_1_24h_2,3mg_8"
        rename_dict["13_48h_SG_S255"] = "SG_1_48h_3,2mg_13"
        rename_dict["15_48h_SG_S256"] = "SG_1_48h_3mg_15"
        rename_dict["111_48h_SG_S228"] = "SG_1_48h_2,9mg_111"
        rename_dict["113_48h_SG_S229"] = "SG_1_48h_3,1mg_113"
        rename_dict["114_48h_SG_S271"] = "SG_1_48h_2,9mg_114"
        rename_dict["18_72h_SG_S257"] = "SG_1_72h_5,2mg_18"
        rename_dict["19_72h_SG_S258"] = "SG_1_72h_7,4mg_19"
        rename_dict["20_72h_SG_S259"] = "SG_1_72h_5,5mg_20"
        rename_dict["116_72h_SG_S307"] = "SG_1_72h_6,2mg_116"
        rename_dict["118_72h_SG_S230"] = "SG_1_72h_5,9mg_118"
        rename_dict["21_96h_SG_S239"] = "SG_1_96h_10,3mg_21"
        rename_dict["22_96h_SG_S240"] = "SG_1_96h_10,3mg_22"
        rename_dict["24_96h_SG_S260"] = "SG_1_96h_10,2mg_24"
        rename_dict["121_96h_SG_S231"] = "SG_1_96h_9,1mg_121"
        rename_dict["122_96h_SG_S261"] = "SG_1_96h_9,5mg_122"
        rename_dict["52_12h_SG_S280"] = "SG_2_12h_2,2mg_52"
        rename_dict["53_12h_SG_S232"] = "SG_2_12h_1,9mg_53"
        rename_dict["152_12h_SG_S233"] = "SG_2_12h_2mg_152"
        rename_dict["153_12h_SG_S234"] = "SG_2_12h_1,7mg_153"
        rename_dict["155_12h_SG_S272"] = "SG_2_12h_1,8mg_155"
        rename_dict["56_1d_SG_S262"] = "SG_2_24h_2,2mg_56"
        rename_dict["59_1d_SG_S263"] = "SG_2_24h_2,4mg_59"
        rename_dict["156_1d_SG_S241"] = "SG_2_24h_2,2mg_156"
        rename_dict["157_1d_SG_S273"] = "SG_2_24h_2,3mg_157"
        rename_dict["159_1d_SG_S242"] = "SG_2_24h_2,5mg_159"
        rename_dict["62_2d_SG_S243"] = "SG_2_48h_4,6mg_62"
        rename_dict["63_2d_SG_S244"] = "SG_2_48h_4,5mg_63"
        rename_dict["65_2d_SG_S245"] = "SG_2_48h_4,2mg_65"
        rename_dict["162_2d_SG_S246"] = "SG_2_48h_5mg_162"
        rename_dict["163_2d_SG_S264"] = "SG_2_48h_4,3mg_163"
        rename_dict["67_3d_SG_S247"] = "SG_2_72h_8,7mg_67"
        rename_dict["68_3d_SG_S265"] = "SG_2_72h_8,2mg_68"
        rename_dict["69_3d_SG_S248"] = "SG_2_72h_9,3mg_69"
        rename_dict["166_3d_SG_S249"] = "SG_2_72h_6,7mg_166"
        rename_dict["167_3d_SG_S250"] = "SG_2_72h_7,2mg_167"
        rename_dict["71_4d_SG_S266"] = "SG_2_96h_11,1mg_71"
        rename_dict["74_4d_SG_S251"] = "SG_2_96h_12,9mg_74"
        rename_dict["75_4d_SG_S252"] = "SG_2_96h_15,1mg_75"
        rename_dict["171_4d_SG_S253"] = "SG_2_96h_10,2mg_171"
        rename_dict["172_4d_SG_S254"] = "SG_2_96h_10,1mg_172"
        rename_dict["pre4_MG_S24"] = "MG_unfed_1,7mg_pre4"
        rename_dict["pre10_MG_S2"] = "MG_unfed_1,4mg_pre10"
        rename_dict["pre_8_MG_S279"] = "MG_unfed_1,6mg_pre8"
        rename_dict["4_12h_MG_S310"] = "MG_1_12h_1,8mg_4"
        rename_dict["5_12h_MG_S311"] = "MG_1_12h_1,8mg_5"
        rename_dict["104_12h_MG_S238"] = "MG_1_12h_1,8mg_104"
        rename_dict["9_24h_MG_S297"] = "MG_1_24h_1,7mg_9"
        rename_dict["10_24h_MG_S298"] = "MG_1_24h_1,6mg_10"
        rename_dict["107_24h_MG_S299"] = "MG_1_24h_1,9mg_107"
        rename_dict["13_48h_MG_S281"] = "MG_1_48h_3,2mg_13"
        rename_dict["15_48h_MG_S282"] = "MG_1_48h_3mg_15"
        rename_dict["114_48h_MG_S283"] = "MG_1_48h_2,9mg_114"
        rename_dict["20_72h_MG_S300"] = "MG_1_72h_5,5mg_20"
        rename_dict["116_72h_MG_S284"] = "MG_1_72h_6,2mg_116"
        rename_dict["118_72h_MG_S285"] = "MG_1_72h_5,9mg_118"
        rename_dict["22_96h_MG_S286"] = "MG_1_96h_10,3mg_22"
        rename_dict["24_96h_MG_S287"] = "MG_1_96h_10,2mg_24"
        rename_dict["122_96h_MG_S288"] = "MG_1_96h_9,5mg_122"
        rename_dict["155_12h_MG_S308"] = "MG_2_12h_1,8mg_155"
        rename_dict["54_12_MG_S309"] = "MG_2_12h_2,2mg_54"
        rename_dict["152_12h_MG_S276"] = "MG_2_12h_2mg_152"
        rename_dict["56_1d_MG_S277"] = "MG_2_24h_2,2mg_56"
        rename_dict["59_1d_MG_S301"] = "MG_2_24h_2,4mg_59"
        rename_dict["157_1d_MG_S278"] = "MG_2_24h_2,3mg_157"
        rename_dict["62_2d_MG_S289"] = "MG_2_48h_4,6mg_62"
        rename_dict["63_2d_MG_S290"] = "MG_2_48h_4,5mg_63"
        rename_dict["163_2d_MG_S302"] = "MG_2_48h_4,3mg_163"
        rename_dict["67_3d_MG_S291"] = "MG_2_72h_8,7mg_67"
        rename_dict["68_3d_MG_S292"] = "MG_2_72h_8,2mg_68"
        rename_dict["167_3d_MG_S293"] = "MG_2_72h_7,2mg_167"
        rename_dict["71_4d_MG_S294"] = "MG_2_96h_11,1mg_71"
        rename_dict["74_4d_MG_S295"] = "MG_2_96h_12,9mg_74"
        rename_dict["171_4d_MG_S296"] = "MG_2_96h_10,2mg_171"
        self.rename_dict = rename_dict
        self.order = []
        for key in self.rename_dict:
            self.order.append("%s_FPKM" % rename_dict[key])

    def matrix2dict(self):
        matrix_dict = {}
        header = ""
        for line in self.lines:
            if line.startswith("gene") == False:
                ID = line.split("\t")[0]
                evalues = line.replace("\n", "").split("\t")[1::]
                matrix_dict[ID] = evalues
            else:
                header = header + line
        header = header.split("\n")[0].split("\t")[1::]
        return header, matrix_dict

    def orderandrename(self, output):
        print(self.df.head())
        order = []
        rename_df = self.df
        for key in self.rename_dict:
            order.append("%s_FPKM" % rename_dict[key])
            rename_df.columns = rename_df.columns.str.replace(key, rename_dict[key])

        print(rename_df.head())
        print(order)
        order_df = rename_df[order]
        order_df.to_csv(output, sep="\t")

class DE:
    def __init__(self, file, excel_transcriptIDs=False):
        with open(file) as de:
            self.lines = list(de)
        self.IDsUPfc = {}
        self.IDsDOWNfc = {}
        self.IDsUPpvalue = {}
        self.IDsDOWNpvalue = {}
        for n in range(1, len(self.lines)):
            line = DE.Line(self.lines[n])
            if excel_transcriptIDs == False:
                if line.fc >= 0:
                    self.IDsUPfc[line.ID] = "%.4f" %line.fc
                    self.IDsUPpvalue[line.ID] = "%.4f" %line.pvalue
                else:
                    self.IDsDOWNfc[line.ID] = "%.4f" %line.fc
                    self.IDsDOWNpvalue[line.ID] = "%.4f" %line.pvalue
            else:
                if line.ID in excel_transcriptIDs:
                    if line.fc >= 0:
                        self.IDsUPfc[line.ID] = "%.4f" %line.fc
                        self.IDsUPpvalue[line.ID] = "%.4f" %line.pvalue
                    else:
                        self.IDsDOWNfc[line.ID] = "%.4f" %line.fc
                        self.IDsDOWNpvalue[line.ID] = "%.4f" %line.pvalue

        self.IDsUP = list(self.IDsUPfc.keys())
        self.IDsDOWN = list(self.IDsDOWNfc.keys())

    class Line:
        def __init__(self, line):
            self.row = line.split("\t")
            self.ID = self.row[0].replace('"','')
            self.pvalue = float(self.row[2])
            self.fc = np.log2(float(self.row[3]))



class Excel:
    def __init__(self, file):
        with open(file) as ex:
            self.lines = list(ex)
        self.IDs = self.getIDlist()
        self.dict = self.excel2dict()

    class Line:
        def __init__(self, line):
            self.row = line.replace("\n", "").split("\t")
            self.cdsID = self.row[0]
            self.seq = self.row[1]
            self.hitID = gethitID(self.row)
            self.identity = self.row[3]
            self.coverage = self.row[4]
            self.bitscore = self.row[5]
            self.description = self.row[6]
            self.db = self.row[7]
            self.ex_count = self.row[8:96]
            self.fpkm = self.row[96:184]

    def getnames(self):
        names = []
        for n in range(1, len(self.lines)):
            line = Excel.Line(self.lines[n])
            names.append(line.hitID)
        names = deduplicate_list(names)
        return names

    def getIDlist(self):
        IDlist = []
        for n in range(1, len(self.lines)):
            line = Excel.Line(self.lines[n])
            IDlist.append(line.cdsID)
        return IDlist

    def excel2dict(self):
        ex_dict = {}
        for n in range(1, len(self.lines)):
            line = Excel.Line(self.lines[n])
            ex_dict[line.cdsID] = self.lines[n]
        return ex_dict


class Uniprot:
    def __init__(self, file, obo_dict):
        with open(file) as uni:
            self.lines = list(uni)
        self.goterm_dict, self.goterm_des_dict, self.keyword_dict, self.keyword_des_dict, self.family_dict = self.uniprot2dict(obo_dict)

    class Line:
        def __init__(self, line):
            self.row = line.replace("\n", "").split("\t")
            self.hitID = self.row[-1]
            self.entry = self.row[0]
            self.entry_name = self.row[1]
            self.protein_name = self.row[2]
            self.organism = self.row[3]
            self.GObp = self.row[4]
            self.GOcc = self.row[5]
            self.GO = self.row[6]
            self.GOmf = self.row[7]
            self.GOterm = self.row[8]
            self.keyword = self.row[9]
            self.keywordID = self.row[10]
            self.family = self.row[11]


    def uniprot2dict(self,obo_dict):
        goterm_dict = {}
        goterm_des_dict = {}
        kw_dict = {}
        kw_des_dict = {}
        family_dict = {}
        for n in range(1,len(self.lines)):
            uniprot_line = Uniprot.Line(self.lines[n])
            if uniprot_line.GOterm != "":
                goterm_dict[uniprot_line.hitID] = uniprot_line.GOterm
                goterm_des_dict[uniprot_line.hitID] = "; ".join([obo_dict[x] for x in uniprot_line.GOterm.split("; ")])
            else:
                goterm_dict[uniprot_line.hitID] = "None"
                goterm_des_dict[uniprot_line.hitID] = "None"
            if uniprot_line.keywordID != "":
                kw_dict[uniprot_line.hitID] = uniprot_line.keywordID
                kw_des_dict[uniprot_line.hitID] = uniprot_line.keyword
            else:
                kw_dict[uniprot_line.hitID] = "None"
                kw_des_dict[uniprot_line.hitID] = "None"
            if uniprot_line.family != "":
                family_dict[uniprot_line.hitID] = uniprot_line.family

        return goterm_dict, goterm_des_dict, kw_dict, kw_des_dict, family_dict

    def getAnnotated(self):
        hitID_list = []
        for n in range(1, len(self.lines)):
            line = Uniprot.Line(self.lines[n])
            if line.GOterm != "":
                hitID_list.append(line.hitID)
        return hitID_list


class Blast2GO:
    def __init__(self, file, obo_dict):
        with open(file) as bg:
            self.lines = list(bg)
        self.dict, self.des_dict = self.blast2GO2dict(obo_dict)

    def blast2GO2dict(self, obo_dict):
        goterm_dict = {}
        goterm_des_dict = {}
        for n in range(1, len(self.lines)):
            ID = self.lines[n].split("\t")[2]
            goterms = self.lines[n].split("\t")[9].replace("C:", "").replace("P:", "").replace("F:", "")
            if goterms != "" and goterms != "no GO terms":
                goterm_dict[ID] = goterms
                goterm_des_dict[ID] = "; ".join([obo_dict[x] for x in goterms.split("; ")])
            else:
                goterm_dict[ID] = "None"
                goterm_des_dict[ID] = "None"
        return goterm_dict, goterm_des_dict


class Fasta:
    def __init__(self, file):
        with open(file) as ft:
            self.lines = list(ft)
        self.dict = self.fasta2dict()

    def split(self, n, output):
        file_number = 1
        count = 0
        out = open("%s_%s.fa" % (output, file_number), "wt")
        for key in self.dict:
            out.write(">%s\n%s\n" % (key, self.dict[key]))
            count += 1
            if count >= int(n):
                out.close()
                file_number += 1
                count = 0
                out = open("%s_%s.fa" % (output, file_number), "wt")

    def fasta2dict(self):
        fasta_dict = dict()
        seqid = ""
        for i in range(0, len(self.lines)):
            if self.lines[i].startswith('>'):
                if seqid != "":
                    fasta_dict[seqid] = seq
                seqid = self.lines[i].split(' ')[0].replace(">", "").replace("\n", "")
                seq = ""
            else:
                seq += self.lines[i].replace("\n", "")
        if seqid != "":
            fasta_dict[seqid] = seq
        return fasta_dict

    def filterlongestisoform(self,output):
        form_dict = {}
        for key in self.dict:
            form_ID = "_".join(key.split("_")[0:-1])
            if form_ID not in form_dict:
                form_dict[form_ID] = key
            else:
                if len(self.dict[key]) > len(self.dict[form_dict[form_ID]]):
                    form_dict[form_ID] = key
        with open(output,"wt") as out:
            for key in form_dict:
                out.write(">%s\n%s\n" %(form_dict[key], self.dict[form_dict[key]]))




class Interpro:
    def __init__(self, file, obo_dict):
        with open(file) as ip:
            self.lines = list(ip)
        self.ipro_dict, self.ipro_des_dict, self.goterm_dict, self.goterm_des_dict = self.interpro2dict(obo_dict)

    class Line:
        def __init__(self, line):
            fields = line.replace("\n", "").split("\t")
            self.ID = fields[0]
            self.Analysis = fields[3]
            self.sign_acc = fields[4]
            self.sign_des = fields[5]
            self.ipro_annotation = fields[11]
            self.ipro_description = fields[12]
            if len(fields) > 13:
                if fields[13] != "-":
                    self.goterm = fields[13].replace("|","; ")
                else:
                    self.goterm = "None"
            else:
                self.goterm = "None"

    def interpro2dict(self, obo_dict):
        ipro_dict = {}
        ipro_des_dict = {}
        goterm_dict = {}
        goterm_des_dict = {}
        for line in self.lines:
            ipro_line = Interpro.Line(line)
            if ipro_line.ID not in ipro_dict and ipro_line.ipro_annotation != "-":
                ipro_dict[ipro_line.ID] = [ipro_line.ipro_annotation]
                ipro_des_dict[ipro_line.ID] = [ipro_line.ipro_description]
            elif ipro_line.ID in ipro_dict and ipro_line.ipro_annotation != "-":
                if ipro_line.ipro_annotation not in ipro_dict[ipro_line.ID]:
                    ipro_dict[ipro_line.ID].append(ipro_line.ipro_annotation)
                    ipro_des_dict[ipro_line.ID].append(ipro_line.ipro_description)
            if ipro_line.ID not in goterm_dict and ipro_line.goterm != "None":
                goterm_dict[ipro_line.ID] = [ipro_line.goterm]
                goterm_des_dict[ipro_line.ID] = [obo_dict[x] for x in ipro_line.goterm.split("; ")]
            elif ipro_line.ID in goterm_dict and ipro_line.goterm != "None":
                if ipro_line.goterm not in goterm_dict[ipro_line.ID]:
                    goterm_dict[ipro_line.ID].append(ipro_line.goterm)
                    goterm_des_dict[ipro_line.ID] += [obo_dict[x] for x in ipro_line.goterm.split("; ")]

        for key in ipro_dict:
            ipro_dict[key] = "; ".join(deduplicate_list("; ".join(ipro_dict[key]).split("; ")))
            ipro_des_dict[key] = "; ".join(deduplicate_list("; ".join(ipro_des_dict[key]).split("; ")))
        for key in goterm_dict:
            goterm_dict[key] = "; ".join(deduplicate_list("; ".join(goterm_dict[key]).split("; ")))
            goterm_des_dict[key] = "; ".join(deduplicate_list("; ".join(goterm_des_dict[key]).split("; ")))

        return ipro_dict, ipro_des_dict, goterm_dict, goterm_des_dict

class SignalP:
    def __init__(self, file):
        with open(file) as sp:
            self.lines = list(sp)
        self.header = self.lines[0:2]
        self.dict = self.signalP2dict()

    def signalP2dict(self):
        sp_dict = {}
        for n in range(2, len(self.lines)):
            ID = self.lines[n].split("\t")[0]
            signal = self.lines[n].split("\t")[1]
            if signal == "OTHER":
                signal = "None"
            sp_dict[ID] = signal
        return sp_dict

class SignalP4:
    def __init__(self, file, tmhmm_mature):
        with open(file) as sp:
            self.lines = list(sp)
        self.dict = self.signalP42dict()

    def signalP42dict(self):
        sp_dict = {}
        count = 0
        for n in range(2,len(self.lines)):
            ID = self.lines[n].split("  ")[0]
            sp = self.lines[n].split("  ")[8][-1]
            if sp == "Y":
                count += 1
                sp = "Yes"
            else:
                sp = "None"
            sp_dict[ID] = sp
        print("%s total of signal peptides" %count)
        return sp_dict

class Phobius:
    def __init__(self, file):
        with open(file) as ph:
            self.lines = list(ph)
        self.dict = self.getSecreted()

    class Line:
        def __init__(self,line):
            self.row = "\t".join(line.split()).split("\t")
            self.ID = self.row[0]
            self.tmhmm = self.row[1]
            self.sp = self.row[2]
            self.prediction = self.row[3]

    def getSecreted(self):
        secreted_dict = {}
        count = 0
        for n in range(1,len(self.lines)):
            line = Phobius.Line(self.lines[n])
            if line.tmhmm == "0" and line.sp == "Y":
                count += 1
                secreted_dict[line.ID] = "Yes"
            else:
                secreted_dict[line.ID] = "None"
        perc = count/(len(self.lines) - 1) * 100
        #print("%i Phobius secreted proteins; %f of the total of %i" %(count,perc,len(self.lines)-1))
        return secreted_dict


class TMHMM:
    def __init__(self, file):
        with open(file) as tm:
            self.lines = list(tm)
        self.TMHlines = self.getTMHlines()
        self.dict = self.tmh2dict()

    def getTMHlines(self):
        TMHlines = []
        for line in self.lines:
            if "Number of predicted TMHs" in line:
                TMHlines.append(line)
        return TMHlines

    def tmh2dict(self):
        tmh_dict = {}
        for line in self.TMHlines:
            ID = line.split(" Number of predicted TMHs:  ")[0].replace("# ", "")
            tmh_n = line.split(" Number of predicted TMHs:  ")[-1].replace("\n", "")
            if tmh_n == "0":
                tmh_n = "None"
            tmh_dict[ID] = tmh_n
        return tmh_dict

class GOobo:
    def __init__(self, file):
        with open(file) as obo:
            self.lines = list(obo)
        self.dict = {}
        for n in range(0,len(self.lines)):
            if "[Term]" in self.lines[n]:
                ID = self.lines[n+1].split(": ")[1].replace("\n","")
                name = self.lines[n+2].split(": ")[1].replace("\n","")
                self.dict[ID] = name
            elif "alt_id:" in self.lines[n]:
                ID = self.lines[n].split(": ")[1].replace("\n","")
                self.dict[ID] = name

class EnrichmentOut:
    def __init__(self, file):
        with open(file) as en:
            self.lines = list(en)

    class Line:
        def __init__(self,line):
            self.row = line.split("\t")
            self.ID = self.row[0]
            self.des = self.row[1]
            self.inbackground = self.row[2]
            self.ingenelist = self.row[3]
            self.notingenelist = self.row[4]
            self.pvalue = self.row[5]
            self.fdr = self.row[6]
            self.re = self.row[7]
            self.type = self.row[8]
            self.genes = self.row[9].replace("\n","")

    def parsedescription(self, des_dict, output):
        with open(output,"wt") as out:
            out.write(self.lines[0])
            for n in range(1,len(self.lines)):
                line = EnrichmentOut.Line(self.lines[n])
                des = des_dict[line.ID]
                out.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(line.ID,des,line.inbackground,line.ingenelist,line.notingenelist,line.pvalue,line.fdr,line.re,line.type,line.genes))



def intersectfasta(fasta, excel, output, reverse=False):
    ft = Fasta(fasta)
    ex = Excel(excel)
    with open(output, "wt") as out:
        for key in ft.dict:
            if not reverse:
                if key in ex.IDs:
                    out.write(">%s\n%s\n" % (key, ft.dict[key]))
            elif reverse:
                if key not in ex.IDs:
                    out.write(">%s\n%s\n" % (key, ft.dict[key]))

def intersectexcel(excel, fasta, output, reverse=False):
    ft = Fasta(fasta)
    ex = CompleteExcel(excel)
    with open(output, "wt") as out:
        out.write(ex.lines[0])
        for n in range(1,len(ex.lines)):
            line = CompleteExcel.Line(ex.lines[n])
            if not reverse:
                if line.cdsID in ft.dict:
                    out.write(ex.lines[n])
            elif reverse:
                if line.cdsID not in ft.dict:
                    out.write(ex.lines[n])

def gethitID(row):
    db = row[7]
    if db == "Uniref90":
        hit = row[2].split("_")[-1]
    elif db == "TickSialoFam":
        hit = row[2]
    elif db == "ArachnidaDB" or "SwissProt":
        hit = row[2].split("|")[-1]
    return hit

def get_expression_dicts(lines):
    evalues_secreted = []
    evalues_nonsecreted = []
    for n in range(1,len(lines)):
        line = CompleteExcel.Line(lines[n])
        if line.classification == "Putatively Secreted":
            evalues_secreted.append(float(line.total_FPKM))
        else:
            evalues_nonsecreted.append(float(line.total_FPKM))
    return evalues_secreted, evalues_nonsecreted



def writenames(names, output):
    with open(output, "wt") as out:
        out.write("\n".join(names))


def deduplicate_list(mylist):
    mylist = list(dict.fromkeys(mylist))
    return mylist


def getlength(dictionary):
    length_dict = dict()
    for key in dictionary:
        length_dict[key] = len(dictionary[key])
    return length_dict


def getUnknown(excel, uniprot, output):
    ex = Excel(excel)
    uni = Uniprot(uniprot)
    uni_IDs = uni.getAnnotated()
    with open(output, "wt") as out:
        for n in range(1, len(ex.lines)):
            line = Excel.Line(ex.lines[n])
            if line.hitID not in uni_IDs:
                out.write(">%s\n%s\n" % (line.cdsID, line.seq))


def startbymethionine(fasta, output):
    ft = Fasta(fasta)
    count = 0
    total_count = 0
    with open(output, "wt") as out:
        for key in ft.dict:
            total_count += 1
            if ft.dict[key].startswith("M"):
                count += 1
                out.write(">%s\n%s\n" % (key, ft.dict[key]))
    print("%s of a total of %s CDS start with methionine." % (count, total_count))


def mergegoterms(goterm1, goterm2):
    merged = goterm1.split("; ") + goterm2.split("; ")
    dd = deduplicate_list(merged)
    return dd


def createCompleteExcel(excel, filtered_fasta, cds_fasta, cds_pep, matrix_FPKM, uniprot, blast2go, interpro, signalp, phobius, tmhmm_mature, tmhmm, output, GOobo_file, list_of_DEfiles):
    #Loading ildes
    ex = Excel(excel)
    ft = Fasta(filtered_fasta)
    pep = Fasta(cds_pep)
    cds = Fasta(cds_fasta)
    fpkm = Ematrix(matrix_FPKM)
    obo = GOobo(GOobo_file)
    uni = Uniprot(uniprot, obo.dict)
    b2go = Blast2GO(blast2go, obo.dict)
    ip = Interpro(interpro, obo.dict)
    tm_mature = TMHMM(tmhmm_mature)
    sp = SignalP4(signalp, tm_mature)
    tm = TMHMM(tmhmm)
    ph = Phobius(phobius)
    DE_dict, FC_dict, pvalue_dict = DEsummarize(list_of_DEfiles)
    uni_count = 0
    b2go_count = 0

    # Header
    with open(output, "wt") as out:
        out.write("ID\tTranscript Sequence\tCDS Sequence\tPeptide Sequence\tFamily\tBestHitID\tBestHitIdentity\tBestHitCoverage\tBestHitBitScore\tBestHitDescription\tBestHitDatabase\tUniprotGO\tUniprotGO names\tKeywords\tKeywords names\tInterproID\tInterpro names\tInterproGO\tInterproGO names\tMergedGO\tMergedGO names\tSignalP\tTMHMM in mature (only for SP)\tTMHMM\tClass\tTissue Specifity\tSG_FC\tMG_FC\tTotal FPKMs\tDE analysis\tDE FC\tPosterior probability of being DE (PPDE)\t%s\n" % ("\t".join(fpkm.order)))
        for key in cds.dict:
            cdsID = key
            transcriptID = key.split(".")[0]
            cds_seq = cds.dict[key]
            pep_seq = pep.dict[key]
            transcript_seq = ft.dict[transcriptID]

            # Excel
            if key in ex.dict:
                ex_line = Excel.Line(ex.dict[key])
                hitID = gethitID(ex_line.row)
                identity = ex_line.row[3]
                coverage = ex_line.row[4]
                bitscore = ex_line.row[5]
                description = ex_line.row[6]
                db = ex_line.row[7]
            else:
                hitID = "None"
                identity = "None"
                coverage = "None"
                bitscore = "None"
                description = "None"
                db = "None"

            # Uniprot and Blast2GO
            if hitID in uni.goterm_dict and uni.goterm_dict[hitID] != "None" and uni.goterm_dict[hitID]!= "":
                GOterm = uni.goterm_dict[hitID]
                GOterm_des = uni.goterm_des_dict[hitID]
                uni_count += 1
            elif cdsID in b2go.dict and b2go.dict[cdsID] != "no GO terms" and b2go.dict[cdsID] != "":
                GOterm = b2go.dict[cdsID]
                GOterm_des = b2go.des_dict[cdsID]
                b2go_count += 1
            else:
                GOterm = "None"
                GOterm_des = "None"
            if hitID in uni.keyword_dict and uni.keyword_dict[hitID] != "None" and uni.keyword_dict[hitID]!= "":
                KW = uni.keyword_dict[hitID]
                KW_des = uni.keyword_des_dict[hitID]
            else:
                KW = "None"
                KW_des = "None"

            # Interpro
            if cdsID in ip.ipro_dict:
                ipro_id = ip.ipro_dict[cdsID]
                ipro_des = ip.ipro_des_dict[cdsID]
            else:
                ipro_id = "None"
                ipro_des = "None"
            if cdsID in ip.goterm_dict:
                ipro_goterm = ip.goterm_dict[cdsID]
                ipro_goterm_des = ip.goterm_des_dict[cdsID]
            else:
                ipro_goterm = "None"
                ipro_goterm_des = "None"

            # Merge GO terms
            if ipro_goterm != "None" and GOterm != "None":
                merged = "; ".join(mergegoterms(GOterm, ipro_goterm))
                merged_des = "; ".join(mergegoterms(GOterm_des, ipro_goterm_des))
            elif ipro_goterm == "None" and GOterm != "None":
                merged = GOterm
                merged_des = GOterm_des
            elif ipro_goterm != "None" and GOterm == "None":
                merged = ipro_goterm
                merged_des = ipro_goterm_des
            else:
                merged = "None"
                merged_des = "None"

            #Family
            if db == "Uniref90" or db == "SwissProt" or db == "ArachnidaDB":
                if hitID in uni.family_dict and uni.family_dict[hitID] != "None":
                    family = uni.family_dict[hitID]
                else:
                    family = " ".join(description.split("=")[0].split(" ")[1:-1])
            elif db == "TickSialoFam":
                family = "%s family" %description.split("|")[0].split(", ")[1]
            elif db == "None" and ipro_des != "None":
                family = ipro_des.split("; ")[0]
            elif db == "None" and merged_des != "None":
                family = merged_des.split("; ")[0]
            else:
                family = "None"



            # SignalP and Phobius
            if cdsID in sp.dict:
                sigp = sp.dict[cdsID]
                if cdsID in tm_mature.dict:
                    sigp_mature_tmhmm = tm_mature.dict[cdsID]
                else:
                    sigp_mature_tmhmm = "None"
            else:
                sigp = "None"
                sigp_mature_tmhmm = "None"

            if cdsID in ph.dict:
                if ph.dict[cdsID] == "Yes" and sigp == "None":
                    sigp = ph.dict[cdsID]


            # TMHMM
            if cdsID in tm.dict:
                tm_domain = tm.dict[cdsID]
            else:
                tm_domain = "None"

            # Class
            if sigp != "None" and sigp_mature_tmhmm == "None" and (merged != "None" or ipro_id != "None"):
                classification = "Putatively Secreted, Annotated"
            elif sigp != "None" and sigp_mature_tmhmm == "None" and merged == "None" and ipro_id == "None":
                classification = "Putatively Secreted, not Annotated"
            elif sigp != "None" and sigp_mature_tmhmm != "None" and merged == "None" and ipro_id == "None":
                classification = "Putatively Non-Secreted, not Annotated"
            elif sigp != "None" and sigp_mature_tmhmm != "None" and (merged != "None" or ipro_id != "None"):
                classification = "Putatively Non-Secreted, Annotated"
            elif sigp == "None" and (merged != "None" or ipro_id != "None"):
                classification = "Putatively Non-Secreted, Annotated"
            elif sigp == "None" and merged == "None" and ipro_id == "None":
                classification = "Putatively Non-Secreted, not Annotated"

            # Modificating classification according to the blast hit
            if "Secreted" in description or "secreted" in description or "|S|" in description:
                if merged != "None" or ipro_id != "None":
                    classification = "Putatively Secreted, Annotated"
                else:
                    classification = "Putatively Secreted, not Annotated"

            # Expression
            evalues = fpkm.dict[transcriptID]
            # Expression class
            evalues_float = []
            for item in evalues:
                evalues_float.append(float(item))
            total_FPKM = np.sum(evalues_float)
            SG_evalues = evalues_float[0:55]
            MG_evalues = evalues_float[55:88]
            SG_mean = np.mean(SG_evalues)
            MG_mean = np.mean(MG_evalues)
            SG_ratio = SG_mean / (MG_mean + SG_mean)
            MG_ratio = 1 - SG_ratio
            SG_FC = np.log2(SG_mean/MG_mean)
            MG_FC = np.log2(MG_mean/SG_mean)
            if SG_FC >= 2:
                eclass = "SG specific"
            elif MG_FC >= 2:
                eclass = "MG specific"
            else:
                eclass = "Neutral"

            #Differential Expression
            if transcriptID in DE_dict:
                DE_str = "; ".join(DE_dict[transcriptID])
                FC = "; ".join(FC_dict[transcriptID])
                pvalue = "; ".join(pvalue_dict[transcriptID])
            else:
                DE_str = "Stable"
                FC = "None"
                pvalue = "None"



            # Writing
            out.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (cdsID, transcript_seq, cds_seq, pep_seq, family, hitID, identity, coverage, bitscore, description, db, GOterm, GOterm_des, KW, KW_des, ipro_id, ipro_des, ipro_goterm, ipro_goterm_des, merged, merged_des, sigp, sigp_mature_tmhmm, tm_domain, classification, eclass, SG_FC, MG_FC, total_FPKM, DE_str, FC, pvalue, "\t".join(evalues)))
    print ("Uniprot GO terms: %s\nB2GO GO terms: %s\n" %(uni_count, b2go_count))

def DEsummarize(list_of_files):
    # The file must have two columns, one for the name of the file, and one for the coding it will have in the string.
    DE_dict = {}
    FC_dict = {}
    pvalue_dict = {}

    with open(list_of_files) as filelist:
        filelist = list(filelist)
        for line in filelist:
            file = line.split("\t")[0]
            name = line.split("\t")[1].replace("\n","")
            de = DE(file)
            for key in de.IDsUPfc:
                if key not in DE_dict:
                    DE_dict[key] = ["%sUP" %name]
                    FC_dict[key] = [str(de.IDsUPfc[key])]
                    pvalue_dict[key] = [str(de.IDsUPpvalue[key])]
                else:
                    DE_dict[key] += ["%sUP" %name]
                    FC_dict[key] += [str(de.IDsUPfc[key])]
                    pvalue_dict[key] += [str(de.IDsUPpvalue[key])]
            for key in de.IDsDOWNfc:
                if key not in DE_dict:
                    DE_dict[key] = ["%sDOWN" %name]
                    FC_dict[key] = [str(de.IDsDOWNfc[key])]
                    pvalue_dict[key] = [str(de.IDsDOWNpvalue[key])]
                else:
                    DE_dict[key] += ["%sDOWN" %name]
                    FC_dict[key] += [str(de.IDsDOWNfc[key])]
                    pvalue_dict[key] += [str(de.IDsDOWNpvalue[key])]
    return DE_dict, FC_dict, pvalue_dict

def DEgetlist(list_of_files, excel_transcriptIDs=False):
    all_dict = {}
    with open(list_of_files) as filelist:
        filelist = list(filelist)
        for line in filelist:
            file = line.split("\t")[0]
            name = line.split("\t")[1].replace("\n","").replace("|","-")
            de = DE(file, excel_transcriptIDs)
            UP_name = "%sUP.txt" %name
            DOWN_name = "%sDOWN.txt" %name
            if UP_name not in all_dict:
                all_dict[UP_name] = de.IDsUP
            else:
                all_dict[UP_name] += de.IDsUP
            if DOWN_name not in all_dict:
                all_dict[DOWN_name] = de.IDsDOWN
            else:
                all_dict[DOWN_name] += de.IDsDOWN
    for key in all_dict:
        with open(key,"wt") as out:
            out.write("\n".join(all_dict[key]))
    files = list(all_dict.keys())
    return files

def enrichment(files,output_label, pwd):
    os.system("mkdir %s/GO %s/KW %s/IP %s/family" %(pwd,pwd,pwd,pwd))
    for file in files:
        outname = "%s_goTerm.tsv" %file.split(".")[0]
        os.system("java -classpath /shared/srna_shared/michael/java/ sRNAfuncTerms.MakeEnrichment %s %s_GO_all.txt" %(file, output_label))
        os.system("mv %s GO/" %outname)
        os.system("java -classpath /shared/srna_shared/michael/java/ sRNAfuncTerms.MakeEnrichment %s %s_KW_all.txt" %(file, output_label))
        os.system("mv %s KW/%s" %(outname, outname.replace("_goTerm.tsv","_KW.tsv")))
        os.system("java -classpath /shared/srna_shared/michael/java/ sRNAfuncTerms.MakeEnrichment %s %s_IP_all.txt" %(file, output_label))
        os.system("mv %s IP/%s" %(outname, outname.replace("_goTerm.tsv","_IP.tsv")))
        os.system("java -classpath /shared/srna_shared/michael/java/ sRNAfuncTerms.MakeEnrichment %s %s_family_all.txt" %(file, output_label))
        os.system("mv %s family/%s" %(outname, outname.replace("_goTerm.tsv","_family.tsv")))

def numberofB2GO(excel,b2go,obo):
    ex = CompleteExcel(excel)
    obo = GOobo(obo)
    b2go = Blast2GO(b2go, obo.dict)
    count = 0
    for n in range(1,len(ex.lines)):
        line = CompleteExcel.Line(ex.lines[n])
        if line.cdsID in b2go.dict and b2go.dict[line.cdsID] != "None":
            count +=1
    print("Number of B2GO hits: %s" %count)


