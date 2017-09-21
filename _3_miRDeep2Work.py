#!/usr/local/bin/python
# --2 and 3--
__author__ = "gray"
__date__ = "20170920"
__version__ = "1.0.1"
__aim__ = """
_3_miRDeep2Work.py

(apply mirDeep2 pipeline)

#1: rename sampleMarker_QC_collapse.fasta [content fasta header rename]

#>1-43234 => seq_1_x43234

[cd to dir]
#2: mapper.pl <input> -c -r <20 bowtie -m> -i -j -o 10 -p <genome> -s <outFa> 
    -t <outArf> -v 2> report_mapper_<SampleMarker>.log

[cd to dir]
#3: miRDeep2.pl <input: outFa> <genomeFasta> <input: outArf> <mature_fa> 
                <mature_out_fa> <hairpin_fa> -t <species Name: Human>
                -P -u 2> report_miRDeep2_<SampleMarker>.log

"""
__Config__ = """
SAMPLESHEET\t{GZ\tSampleMarker}
QCWORK_DIR\t{Path}
mirDeep2_DIR\t{Path}
miRDeep2_Config\t{miRDeep2Config_fp}
"""


import sys
import os
import subprocess as sup
import multiprocessing as mul

def ConFigSplit(ConFig):
    reDict = {}
    with open(ConFig, "r") as Fr:
        content = Fr.readlines()
    for line in content:
        if line.strip() != "":
            item = line.strip().split("\t")
            reDict[item[0]] = item[1]
    return reDict

#---------------
def createEnv(SampleMarker, WorkDir="./miRDeep2Work/"): #in config file
    print(" ".join(["mkdir", "-p", WorkDir+"/"+SampleMarker+"/Data/"]))
    sup.call(["mkdir", "-p", WorkDir+"/"+SampleMarker+"/Data/"])

def renameFunction(SampleMarker, QCworkDir, WorkDir="./miRDeep2Work/"):
    oriFile = QCworkDir+"/"+SampleMarker+"_QC_collapse.fasta"
    outFile = WorkDir+"/"+SampleMarker+"/Data/"+SampleMarker+"_rename.fasta"
    #do it
    with open(oriFile, "r") as Fr:
        with open(outFile, "w") as Fw:
            lineCount = 0
            for line in Fr:
                lineCount += 1
                if lineCount%2 == 1:#header #>1-423232 => >seq_1_x423232
                    item = line.strip(">\n").split("-")
                    header = ">seq_"+item[0]+"_x"+item[1]+"\n"
                    Fw.write(header)
                else:
                    Fw.write(line)


def DoRename(parameter):
    SampleMarker = parameter["SampleMarker"]
    WorkDir = parameter["miRDeep2WorkDir"]
    QCDir = parameter["QC_Dir"]
    createEnv(SampleMarker, WorkDir)
    print("RenameWork:"+SampleMarker)
    renameFunction(SampleMarker, QCDir, WorkDir)

#---------------------------------
#mirdeep2
def mapper(ori_data, marker, miRDeep2config):
    #Config:
    miRDeep2Dict = ConFigSplit(miRDeep2config)
    #
    comd = ["mapper.pl", ori_data, "-c",  #fasta
     "-r", "20", #bowtie -m (multiple map
     "-i", "-j", #-i convert to dna, -j remove other letters
     "-o", "10", #bowtie threads (10 *N) < 60
     "-p", miRDeep2Dict["genome"], #bowtie genome
     "-s", marker+"_reads.fa",            #out1 outFa
     "-t", marker+"_reads_vs_mapper.arf", #out2 outArf
     "-v"]#, "2>", "report_mapper_"+marker+".log"]
    p = sup.Popen(comd, stderr=sup.PIPE)
    p.wait()
    p_stderr = p.communicate()[1]
    with open("report_mapper_"+marker+".log", "wb") as Fw:
        Fw.write(p_stderr)

def miRDeep2(marker, miRDeep2config):
    #Config:
    miRDeep2Dict = ConFigSplit(miRDeep2config)
    #
    comd =["miRDeep2.pl", marker+"_reads.fa", #outFa
     miRDeep2Dict["genomeFa"], 
     marker+"_reads_vs_mapper.arf",     #outArf
     miRDeep2Dict["mature_fa"],         #miRNA mature
     miRDeep2Dict["mature_out_fa"],     #out spec miRNA
     miRDeep2Dict["hairpin_fa"],        #hairpin
     "-t", miRDeep2Dict["miRDeep2Spec"],#via quantifier.pl
     "-P", "-u"]#, "2>", "report_miRDeep2_"+marker+".log" ]
    p = sup.Popen(comd, stderr=sup.PIPE)
    p.wait()
    p_stderr = p.communicate()[1]
    with open("report_miRDeep2_"+marker+".log", "wb") as Fw:
        Fw.write(p_stderr)


def miRDeep2env(SampleSheet ,miRDeep2Dir, miRDeep2Config):
    paraList = []
    with open(SampleSheet, "r") as Fsheet:
        content = Fsheet.readlines()
    for line in content: #No header
        Marker = line.strip("\n").split("\t")[1]
        #outFa = Marker+"_read.fa"
        #outArf = Marker+"_reads_vs_mapper.arf"

        #need fork?
        para_mapper = (Marker, miRDeep2Dir, miRDeep2Config)#, outFa, outArf)
        paraList.append(para_mapper)
    return paraList

#Multiple use
def DomiRDeep2(parameter):
    marker = parameter[0]
    miRDeep2Dir = parameter[1]
    miRDeep2config = parameter[2]
    #outFa = parameter[3]
    #outArf = parameter[4]
    #chdir
    workDir = miRDeep2Dir+"/"+marker
    ori_data = workDir+"/Data/"+marker+"_rename.fasta"
    os.chdir(workDir)
    print(os.getcwd())
    
    mapper(ori_data, marker, miRDeep2config)
    miRDeep2(marker, miRDeep2config)


if __name__ == "__main__":
    if os.path.exists(sys.argv[1]):
        pass
    else:
        print(___Config__)
        sys.exit(1)

    ConfigDict = ConFigSplit(sys.argv[1])
    SampleSheet = ConfigDict["SAMPLESHEET"]
    QCdir = ConfigDict["QCWORK_DIR"]
    miRDeep2Dir = ConfigDict["miRDeep2_DIR"]
    #-rename
    renameDictList = []
    with open(SampleSheet, "r") as Fsheet:
        content = Fsheet.readlines()
    for line in content: #no header
        item = line.strip("\n").split("\t")
        Marker = item[1]
        renameDictList.append( {"SampleMarker":Marker, "QC_Dir":QCdir, 
                                "miRDeep2WorkDir":miRDeep2Dir})
    if len(renameDictList) < 20:
        pool = mul.Pool(len(renameDictList))
        pool.map(DoRename, renameDictList)
    else:
        pool = mul.Pool(20)
        pool.map(DoRename, renameDictList)

    #-miRDeep2
    miRDeep2Config = ConfigDict["miRDeep2_Config"]    
    miRDeep2Paralist = miRDeep2env(SampleSheet ,miRDeep2Dir, miRDeep2Config)

    pool.map(DomiRDeep2, miRDeep2Paralist)
    

