#!usr/local/bin/python
#--2 and 3--
import multiprocessing as mul
import sys
import subprocess as sup

_author = "gray"
_date = "20180124"
_version = "1.0.0"
_aim = sys.argv[0]+""" #For miRNAseq use (detector tRNA-RF or tiRNA via blastn

#blastdb path
#Report path

"""
_Config = """
miRDeep2_DIR    {mirdeep2 worspace}
BLASTDB {tRNA or other DB name}
BLASTREPORT {Report reads}
"""

def runBlast(fasta_fp, blastdb):
    """
    run blastn -outfmt 6
    
    Args:
       fasta_fp : in miRDeep2Dir file path, add _reads.fa
        blastdb : for blastn -db using

    """
    comd = ["blastn", "-query", fasta_fp, "-db", blastdb, "-outfmt", "6",
            "-out", fasta_fp+"_blastReport.txt"]
    sup.call(comd)
    print(" ".join(comd))
    print("END:"+fasta_fp)


def makefileList(SAMPLESHEET, miRDeep2Dir):
    """
    create a full file list
    
    Args:
        SAMPLESHEET, samplesheet file path
        miRDeep2Dir, collpase path
    """
    SampleNameList = []
    with open(SAMPLESHEET) as Fr:
        for line in Fr:
            item = line.strip().split("\t")
            file_fp = miRDeep2Dir+"/"+item[1]+"/"+item[1]+"_reads.fa"
            SampleNameList.append(file_fp)
    return SampleNameList

def blastReadCount(FastaSampleList, BLASTreport, id=99):
    """merge blast report (count the reads of BLASTresult)
    
    Args:
        FastaSampleList, with mirDeep2 path. 
                         the  subitem need add "_blastReport.txt"
        BLASTreport, output file path
        id, blast identity
    
    """

    with open(BLASTreport, "w") as Freport:
        for f in FastaSampleList:
            fid = f.split("/")[-1].replace("_reads.fa","")
            Readsum = 0
            StoreAlign = []
            with open(f+"_blastReport.txt") as Fr:
                for line in Fr:
                    item = line.strip().split("\t")
                    if float(item[2]) >= id:
                        if item[0] not in StoreAlign:
                            StoreAlign.append(item[0])
                            Readsum += float(item[0].split("_x")[1])
            Freport.write(fid+"\t"+str(Readsum)+"\n")

        


def ConfigSplit(Con_fp):
    reDict = {} #name to value
    with open(Con_fp, "r") as Fr:
        content = Fr.readlines()
        for line in content:
            if line.strip() != "":
                print(line)
                item = line.strip("\n").split("\t")
                reDict[item[0]] = item[1]
    return reDict



if __name__ == "__main__":
    ConfigFile = sys.argv[1]
    configDict = ConfigSplit(ConfigFile)
    sampleSheet = configDict["SAMPLESHEET"]
    blastdb = configDict["BLASTDB"]
    blastout = configDict["BLASTREPORT"]
    miRDeep2dir = configDict["miRDeep2_DIR"]
    
    #get samplefilelist
    SampleList = makefileList(sampleSheet, miRDeep2dir)

    #runblast
    #create using input
    for file_path in SampleList:
       runBlast(file_path, blastdb)

    #get Report
    blastReadCount(SampleList, blastout, 99)

    print(sys.argv[0]+"END")






