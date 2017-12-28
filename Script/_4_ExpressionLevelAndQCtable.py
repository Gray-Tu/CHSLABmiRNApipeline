#!/usr/local/bin/python
# --2 and 3--
__author__ = "gray"
__date__ = "20171228"
__version__ = "1.0.2"
__aim__ = """
_4_ExpressionLevelAndQCtable.py

Create a miRNA profiling, make a QC table for this analysis batch.

input: Global Config file ( SampleSheet:    sample marker
                            Raw_dir:        RawReads
                            QC_dir:         QCReads
                            miRDeep2_dir:   mapping to genome
                            )

ouput:  QCtable_fp,
        miRNAexpression profiling
"""
__Config__ ="""

SAMPLESHEET {Filepath}  #marker
STORE_DIR   {Dirpath}   #+/+marker+"fastq.gz"
QCWORK_DIR  {Dirpath}   #+/+marker+"_QC.fasta"
miRDeep2_DIR{Dirpath}   #+/+marker+/+{miRNAs_expressed_all_samples+*.csv}
                        #           +report_mapper_+"marker"+.log
                    

QCTABLE_fp  {Filepath}
miRNA_EXP   {Filepath}

"""
import sys
import os
import subprocess as sup
import glob

def ConfigSplit(_con_fp):
    reDict = {} #name to value
    with open(_con_fp, "r") as Fr:
        content = Fr.readlines()
    for line in content: #no header
        if line.strip() != "":
            item = line.strip().split("\t")
            if "#" != item[0][0]:
                reDict[item[0]] = item[1]
    return reDict

def SampleMarkerParser(_sample_sheet):
    reList = []
    with open(_sample_sheet, "r") as Fsheet:
        content = Fsheet.readlines() 
    reList = [line.strip().split("\t")[1] for line in content]
    return reList
#---------------------------------

def QCtable_Raw(_rawDir, _marker):
    GZfile = _rawDir+"/"+_marker+".fastq.gz"
    #zcat:
    pZcat = sup.Popen(["zcat",GZfile], stdout=sup.PIPE)
    Reads = sup.check_output(("wc","-l"), stdin=pZcat.stdout).decode("UTF-8")
    pZcat.wait()
    return float(str(Reads).split(" ")[0])/4
 
def QCtable_QC(_qcDir, _marker):
    QCfile = _qcDir+"/"+_marker+"_QC.fasta"
    pWc = sup.check_output(["wc", "-l", QCfile]).decode("UTF-8")
    return float(pWc.split(" ")[0])/2

def QCtable_MapGenome(_miRDeep2Dir, _marker):
    #report_mapper.log:
    ReportLog = _miRDeep2Dir+"/"+_marker+"/report_mapper_"+_marker+".log"
    #QC-reads and map to gemone
    with open(ReportLog, "r") as Fmap:
        content = Fmap.readlines()
    logLineItem = content[-1].strip().split("\t")
    Mapread = logLineItem[1]
    return float(Mapread)

def QCtable_MapMiRNA(_miRDeep2Dir, _marker):
    #read_occ
    ReadOcc = glob.glob(_miRDeep2Dir+"/"+_marker+"/expression_analyses/*/read_occ")
    #print(miRDeep2Dir+"/"+SampleMarker+"/expression_analyses/*/read_occ")
    #check file number
    if len(ReadOcc) != 1:
        print("Error: "+_marker+" read_occ file lost!"+str(ReadOcc))
        sys.exit(1)
    else:
        
        with open(ReadOcc[0], "r") as Fr:
            cont = Fr.readlines()        
    return sum([float(line.split("\t")[0].split("_x")[1]) for line in cont])

def QCtableOut(_sample_list, _rawDir, _qcDir, _miRDeep2Dir, QCtable_fp):
    print("Start: make QCtable:"+QCtable_fp)
    #format: Sample, RAW, QC, QC/RAW, MAPtoGenome, mapGenome/QC, MAPtomiRNA, mapMiRNA/QC, mapMiRNA/RAW
    with open(QCtable_fp, "w") as Fqc:
        #firstline
        Fqc.write("\t".join(["Sample", "RAW_reads", "QC_reads", "QC/RAW", 
                     "mapGenome_reads", "mapGenome/QC", "mapMiRNA_reads",
                     "mapMiRNA/QC", "mapMiRNA/RAW"])+"\n")
        for _marker in _sample_list:
            RAW = QCtable_Raw(_rawDir, _marker)
            QC = QCtable_QC(_qcDir, _marker)
            mapGenome = QCtable_MapGenome(_miRDeep2Dir, _marker)
            mapMiRNA = QCtable_MapMiRNA(_miRDeep2Dir, _marker)
            wl = [  _marker,
                    str(int(RAW)),
                    str(int(QC)),
                    str(QC/RAW),
                    str(int(mapGenome)),
                    str(mapGenome/QC),
                    str(int(mapMiRNA)),
                    str(mapMiRNA/QC),
                    str(mapMiRNA/RAW) ]
            Fqc.write("\t".join(wl)+"\n")
    print("End: QCtable")
        
#-----------------------------------
def MapReadprofile(_miRDeep2Dir, _marker):
    reDict = {} # miRNAname: [read, rpm]
    #check "miRNAs_expressed_all_samples*.csv"
    ExpFile = glob.glob(_miRDeep2Dir+"/"+_marker+"/miRNAs_expressed_all_sample*.csv")
    if len(ExpFile) != 1:
        print("Error: "+_marker+" miRNA expressed file lost")
        sys.exit(1)
    else:
        with open(ExpFile[0], "r") as Fexp:
            content = Fexp.readlines()
        for line in content[1:]:
            item = line.strip().split("\t")
            matureName = item[0]
            rawCount = float(item[-2])
            rpm = float(item[-1])
            #Methond 1, Avg reads and rpm
            if matureName in reDict:
                reDict[matureName].append((rawCount, rpm))
            else:
                reDict[matureName] = [(rawCount, rpm)]
    return reDict

def AllSampleProfile(_sample_list,  _miRDeep2Dir):
    reBigDict = {} # marker: ProfileDict
    for marker in _sample_list:
        reBigDict[marker] = MapReadprofile(_miRDeep2Dir, marker)

    return reBigDict

def OutReadprofile(_profile_dict, rawRead_fp, rpm_fp):
    print("Start: make expression profile: "+rpm_fp)
    Samplelist = list(_profile_dict.keys())
    miRNAlist = list(_profile_dict[Samplelist[0]].keys())
    with open(rawRead_fp, "w") as Fread:
        with open(rpm_fp, "w") as Frpm:
            #firstline
            Fread.write("miRNA\t"+"\t".join(Samplelist)+"\n")
            Frpm.write("miRNA\t"+"\t".join(Samplelist)+"\n")
            for miRNA in miRNAlist:
                Fread.write(miRNA)
                Frpm.write(miRNA)
                for sample in Samplelist:
                    Reads = sum([read for read,rpm in _profile_dict[sample][miRNA] ])/max(len(_profile_dict[sample][miRNA]), 1)
                    RPMs = sum([rpm for read,rpm in _profile_dict[sample][miRNA] ])/max(len(_profile_dict[sample][miRNA]), 1)
                    Fread.write("\t{0:.2f}".format(Reads))
                    Frpm.write("\t{0:.2f}".format(RPMs))
                Fread.write("\n")
                Frpm.write("\n")
    print("End: expression profile")

if __name__ == "__main__":
    #--test
    ConFigDict = ConfigSplit(sys.argv[1])
    SampleSheet = ConFigDict["SAMPLESHEET"]
    SAMPLELIST = SampleMarkerParser(SampleSheet)
    MIRDEEP2Dir = ConFigDict["miRDeep2_DIR"]
    RAWDir = ConFigDict["STORE_DIR"]
    QCDir = ConFigDict["QCWORK_DIR"]
    QCtable_fp = ConFigDict["QCTABLE_fp"]
    miRNA_read_Profle = ConFigDict["MiRNA_READ_Profile_fp"]
    miRNA_rpm_Profile = ConFigDict["MiRNA_RPM_Profile_fp"]

    #---------
    #QCtable
    print(SAMPLELIST, RAWDir, QCDir, MIRDEEP2Dir, QCtable_fp)
    #QCtableOut(SAMPLELIST, RAWDir, QCDir, MIRDEEP2Dir, QCtable_fp)
    
    #exp
    ProfileDict = AllSampleProfile(SAMPLELIST, MIRDEEP2Dir)
    #print(ProfileDict)
    OutReadprofile(ProfileDict, miRNA_read_Profle, miRNA_rpm_Profile)

    #Sample1
    '''
    A = SampleList[0]
    raw = QCtable_Raw(rawDir, A)
    QC = QCtable_QC(QCDir, A)
    mapG = QCtable_MapGenome(miRDeep2Dir, A)
    miRNA = QCtable_MapMiRNA(miRDeep2Dir, A)
    print((raw, QC, mapG, miRNA))
    '''







