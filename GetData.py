#!/usr/bin/python
#--2 and 3--
__author__ = "gray"
__date__ = "20170919"
__version__ = "1.0.1"
__aim__ = """
GetData.py for miseq pipeline CHSLAB used

Copy file,
Rename file,
unzip file > for QC used

input:
    sample sheet
    project Dir (Target Dir)
   
[sample sheet] format
RawSampleName\tNewSampleName[marker]

"""
import sys
import os
import subprocess as sup


def GetData(SampleSheet, TargetDir="./"):
    #check SampleSheet
    if os.path.exists(SampleSheet):
        pass
    else:
        print("No Find:"+SampleSheet)
        sys.exit(1)
    #--------
    with open(SampleSheet,"r") as Fr:
        #no header
        content = Fr.readlines()
    for line in content:
        item = line.strip("\n").split("\t")
        Oripath = item[0]
        Marker = item[1]
        #---cp change name, (with gz file)
        Comd = "cp "+Oripath+" "+TargetDir+"/"+Marker+".fastq.gz"
        print(Comd)
        sup.call(Comd, shell=True)



#
if __name__ == "__main__":
    SampleSheet = sys.argv[1]
    TargetDir = sys.argv[2]
    #check dir 
    if os.path.exists(TargetDir):
        pass
    else:
        sup.call("mkdir -p "+TargetDir,shell=True)
    GetData(SampleSheet, TargetDir)
    






