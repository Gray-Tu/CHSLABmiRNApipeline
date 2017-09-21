#!/usr/local/bin/python
#--2 and 3--
__author__ = "gray"
__date__ = "20170920"
__version__ = "1.0.1"
__aim__ = """
_2_QCRun.py #For miRNAseq use (TruSeq adapter content)

#cutadapt -a <TruSeq 3' adapt> -m 18 -0 <outFastq> <input.gz, default .gz>

#fastq_quality_filter -q <Q:20> -p <70> -o <output> -i <inPutfastq> -Q33

#fastq_quality_trimmer -t <Q:20> -o <output> -i <inPutFastq> -Q33

===================
format:
#fastq_to_fasta -o <output> -i <input> -Q33

#collapse
#fastx_collapser -o <output> -i <input>

"""
__Config__ = """
SAMPLESHEET {GZfile, SampleMarker}
STORE_DIR   {RawDir}
QCWORK_DIR  {QCDir}
QC_3ADA {TGGAATTCTCGGGTGCCAAGG}
QC_MINLEN   {18}
QC_QVALUE   {20}
QC_Pcent    {70}
Q33         {Q33 or ""}
CORE        {20}
"""

import subprocess as sup
import sys
import os
import multiprocessing as mul


def cmdCutadapt(SampleMarker=None, ada3="TGGAATTCTCGGGTGCCAAGG",minilength="18", 
                                               StoreDir="./", TarDir="./"):
    InFile = StoreDir+"/"+SampleMarker+".fastq.gz"
    OutFile = TarDir+"/"+SampleMarker+"_de_ada.fastq"
    #out file suffix: _de_ada
    return  ["cutadapt", "-a", ada3, "-m", minilength, "-o", OutFile, InFile]
                                                         


def cmdFastxTookit(SampleMarker, StoreDir, Qvalue="20", Pcent="70", Q33="", 
                                          TarFileDir="./", minilength="18"):

    Input = StoreDir+"/"+SampleMarker+"_de_ada.fastq"
    filterOut = TarFileDir+"/"+SampleMarker+"_filtered.fastq"
    #fastq filter
    comdFilter = ["fastq_quality_filter", "-q", Qvalue, "-p" ,Pcent, 
                  "-o", filterOut, "-i", Input, Q33]

    #trimmer
    trimOut = TarFileDir+"/"+SampleMarker+"_QC.fastq"
    comdTrim = ["fastq_quality_trimmer", "-t", Qvalue, "-l", minilength,
                "-i", filterOut, "-o", trimOut, Q33]

    #fq to fa
    faOut = TarFileDir+"/"+SampleMarker+"_QC.fasta"
    comdFq2Fa= ["fastq_to_fasta", "-i", trimOut, "-o", faOut, Q33]

    #collapse
    collOut = TarFileDir+"/"+SampleMarker+"_QC_collapse.fasta"
    comdColl = ["fastx_collapser", "-i", faOut, "-o", collOut]

    return [comdFilter, comdTrim, comdFq2Fa, comdColl]


def DoQC(cmdTuple):
    cutada = cmdTuple[0]
    filter = cmdTuple[1]
    trimmer= cmdTuple[2]
    Fq2Fa  = cmdTuple[3]
    coll   = cmdTuple[4]
    comd = ""
    for cmd in cmdTuple:
        print(cmd)
        sup.call(cmd)
'''    for cmd in cmdTuple:
        comd += " ".join(cmd)+";"
        print(comd)
    
    sup.call(comd, shell = True)
'''

def poolUse(CmdList, Core):
    pool = mul.Pool(Core)
    pool.map(DoQC, CmdList)


def ConfigSplit(Con_fp):
    reDict = {} #name to value
    with open(Con_fp, "r") as Fr:
        content = Fr.readlines()
    for line in content:
        if line.strip() != "":
            item = line.strip("\n").split("\t")
            reDict[item[0]] = item[1]
    return reDict
    
if __name__ == "__main__":
    #input
    Config = sys.argv[1]
    if os.path.exists(Config):
        pass
    else:
        print(__Config__)
        sys.exit(1)
    ConFigDict = ConfigSplit(Config)
    SampleSheet = ConFigDict["SAMPLESHEET"]    
    OriDir = ConFigDict["STORE_DIR"]
    QCworkDir = ConFigDict["QCWORK_DIR"]
    QC_3ada = ConFigDict["QC_3ADA"]
    QC_minlen = ConFigDict["QC_MINLEN"]
    QC_qval = ConFigDict["QC_QVALUE"]
    QC_pcent = ConFigDict["QC_Pcent"]
    Q33 = ConFigDict["Q33"]
    Core = int(ConFigDict["CORE"])

    #QCwork dir 
    if os.path.exists(QCworkDir):
        pass
    else:
        sup.call(["mkdir", "-p", QCworkDir])


    CmdList = [] #(a,b,c,d,e)
    with open(SampleSheet, "r") as Fsheet:
        content = Fsheet.readlines() #No header
    for line in content:
        item = line.strip("\n").split("\t")
        Marker = item[1]
        cutCmd = cmdCutadapt(Marker, QC_3ada, QC_minlen, OriDir, QCworkDir)
        fastxCmd = cmdFastxTookit(Marker, QCworkDir, QC_qval, QC_pcent, Q33,
                                                       QCworkDir, QC_minlen)
        CmdList.append([cutCmd]+fastxCmd)
        

    #---------------
    poolUse(CmdList, Core)





