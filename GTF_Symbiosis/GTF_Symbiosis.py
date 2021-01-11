from gtfparse import read_gtf
import pandas as pd
import numpy as np
import ensembl_rest
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
from pywaffle import Waffle 
from prettytable import PrettyTable
import math
from matplotlib.gridspec import GridSpec
from copy import deepcopy
from collections import OrderedDict
from BCBio import GFF
from Bio import SeqIO
from itertools import combinations
from venn import venn
import networkx as nx
import itertools
from venn import pseudovenn
from upsetplot import from_contents
from upsetplot import plot
import os

##Collapsed Returns (Returns transcripts and genes with a list of exons in those regions)

##Given a GTF Files which contains Transcripts, Genes, Codons, Exons, etc, it finds all exons associated with a transcript 
##and adds them into an exon column

##colFrame --- GTF file dataframe

def collapsedReturns(colFrame2):
    colFrame2['Exons'] = np.empty((len(colFrame2), 0)).tolist()
    counter = 0
    while counter < len(colFrame2.index):
        exons = []
        if colFrame2.at[counter,"feature"] == "gene" and colFrame2.at[counter+1,"feature"] == "transcript":
            itterator=2
            while counter+itterator < len(colFrame2.index) and colFrame2.at[counter+itterator,"feature"] != "transcript":
                if colFrame2.at[counter+itterator,"feature"] == "exon":
                    colFrame2.at[counter+1,"Exons"].append([colFrame2.at[counter+itterator,"start"],colFrame2.at[counter+itterator,"end"]])
                itterator=itterator+1
            counter = counter + itterator - 1
        elif colFrame2.at[counter,"feature"] == "transcript" and counter == 0:
            itterator=1
            while counter+itterator < len(colFrame2.index) and colFrame2.at[counter+itterator,"feature"] != "transcript":
                if colFrame2.at[counter+itterator,"feature"] == "exon":
                    colFrame2.at[counter,"Exons"].append([colFrame2.at[counter+itterator,"start"],colFrame2.at[counter+itterator,"end"]])
                itterator=itterator+1
            counter = counter + itterator - 1
        elif colFrame2.at[counter,"feature"] == "transcript" and colFrame2.at[counter-1,"feature"] != "gene":
            itterator=1
            while counter+itterator < len(colFrame2.index) and colFrame2.at[counter+itterator,"feature"] != "transcript":
                if colFrame2.at[counter+itterator,"feature"] == "exon":
                    colFrame2.at[counter,"Exons"].append([colFrame2.at[counter+itterator,"start"],colFrame2.at[counter+itterator,"end"]])
                itterator=itterator+1
            counter = counter + itterator - 1
        else:
            counter=counter+1
    return(colFrame2)

##Exon/Splice Junction Matching (returns a list of 0s and 1s coresponding to whether the exons/splice junctions in list1 match any exons/splice junctions in list2, takes a 
##threshold value which determines how close they can be)

##list1 --- list of exons stored as list of two values [start, stop]
##list2 --- list of exons stored as list of two values [start, stop]
##threshold --- distance away between start and stop values that are considered matching

def exonMatching(list1, list2, threshold):
    itNum = len(list1)
    itNum2 = len(list2)
    listofzeroes = [0]*itNum
    itterator = 0
    itterator2 = 0
    while itterator < itNum:
        while itterator2 < itNum2 and list1[itterator][0] >= (list2[itterator2][1] + threshold):
            itterator2 = itterator2 + 1
        if itterator2 < itNum2 and abs(list1[itterator][1] - list2[itterator2][1]) <= threshold and abs(list1[itterator][0] - list2[itterator2][0]) <= threshold:
            listofzeroes[itterator] = 1
            itterator = itterator + 1
        else:
            itterator = itterator + 1
            
    return(listofzeroes)

##Start/Stop Matching (Gives the distance between the start and stops of two given exon lists)

##list1 --- list of exons stored as list of two values [start, stop]
##list2 --- list of exons stored as list of two values [start, stop]

def startStopDistance(list1, list2):
    itNum = len(list1)
    itNum2 = len(list2)
    
    if itNum == 1:
        start1 = list1[0][0]
        end1 = list1[0][1]
    else:
        start1 = list1[0][0]
        end1 = list1[itNum-1][1]
        
    if itNum2 == 1:
        start2 = list2[0][0]
        end2 = list2[0][1]
    else:
        start2 = list2[0][0]
        end2 = list2[itNum2-1][1]
            
    return [start1-start2, end1-end2]

##Best Match (removes all fits below a certain amount, and returns the best fit of the others)

##passedFrame --- a GTF file which has already been matched
##threshold --- a value that each match must pass to be able to be counted
##label2 -- a passed label to be removed from the list of locations

def bestMatch(passedFrame, threshold, label2 = ""):
    gtfFrame = passedFrame.copy()
    gtfFrame['Fit_Value'] = np.empty((len(gtfFrame), 0)).tolist()
    oneExonSize = 0
    for index, row in gtfFrame.iterrows():
        biggestValue=0
        bvIndex=-1
        if not gtfFrame.at[index,"Matched_Transcript"]:
            pass
        elif len(gtfFrame.at[index,"Splice_Junctions"]) == 0:
            oneExonSize = oneExonSize + 1
            i=1
            bvIndex = 0
            onExValue = abs(gtfFrame.at[index, "Start_Stop_Distance"][0][0])+abs(gtfFrame.at[index, "Start_Stop_Distance"][0][1])
            while i < len(gtfFrame.at[index,"Matched_Transcript"]):
                matchValue = abs(gtfFrame.at[index, "Start_Stop_Distance"][i][0])+abs(gtfFrame.at[index, "Start_Stop_Distance"][i][1])
                if matchValue < onExValue:
                    bvIndex = i
                    onExValue= matchValue
                i = i + 1
            gtfFrame.at[index,"Matched_Transcript"] = gtfFrame.at[index,"Matched_Transcript"][bvIndex]
            gtfFrame.at[index,"Start_Stop_Distance"] = gtfFrame.at[index,"Start_Stop_Distance"][bvIndex]
            gtfFrame.at[index,"Splice_Junction_Matches"] = gtfFrame.at[index,"Splice_Junction_Matches"][bvIndex]
            gtfFrame.at[index,"Gene_Exon_Lengths"] = gtfFrame.at[index,"Gene_Exon_Lengths"][bvIndex]
            gtfFrame.at[index,"Gene_SJ_Lengths"] = gtfFrame.at[index,"Gene_SJ_Lengths"][bvIndex]
            gtfFrame.at[index,"Matches"] = gtfFrame.at[index,"Matches"][bvIndex]
            gtfFrame.at[index,"Perfect_Matches"] = gtfFrame.at[index,"Perfect_Matches"][bvIndex]
            gtfFrame.at[index,"Fit_Value"] = -1
        else:
            i=0
            while i < len(gtfFrame.at[index,"Matched_Transcript"]):
                matchValue = 1 - (abs(len(gtfFrame.at[index,"Splice_Junctions"])-gtfFrame.at[index,"Matches"][i]) + abs(len(gtfFrame.at[index,"Splice_Junctions"])-gtfFrame.at[index,"Gene_SJ_Lengths"][i])) / len(gtfFrame.at[index,"Splice_Junctions"])
                if matchValue > biggestValue and matchValue >= threshold:
                    bvIndex = i
                    biggestValue= matchValue
                gtfFrame.at[index, "Fit_Value"].append(matchValue)
                i = i + 1
            if bvIndex == -1:
                gtfFrame.at[index,"Matched_Transcript"] = np.NaN
                gtfFrame.at[index,"Start_Stop_Distance"] = np.NaN
                gtfFrame.at[index,"Splice_Junction_Matches"] = np.NaN
                gtfFrame.at[index,"Gene_Exon_Lengths"] = np.NaN
                gtfFrame.at[index,"Gene_SJ_Lengths"] = np.NaN
                gtfFrame.at[index,"Matches"] = np.NaN
                gtfFrame.at[index,"Perfect_Matches"] = np.NaN
                gtfFrame.at[index,"Fit_Value"] = np.NaN
                if label2 != "":
                    gtfFrame.at[index,"Located_In"] = gtfFrame.at[index,"Located_In"].remove(label2)
                    
            else:
                gtfFrame.at[index,"Matched_Transcript"] = gtfFrame.at[index,"Matched_Transcript"][bvIndex]
                gtfFrame.at[index,"Start_Stop_Distance"] = gtfFrame.at[index,"Start_Stop_Distance"][bvIndex]
                gtfFrame.at[index,"Splice_Junction_Matches"] = gtfFrame.at[index,"Splice_Junction_Matches"][bvIndex]
                gtfFrame.at[index,"Gene_Exon_Lengths"] = gtfFrame.at[index,"Gene_Exon_Lengths"][bvIndex]
                gtfFrame.at[index,"Gene_SJ_Lengths"] = gtfFrame.at[index,"Gene_SJ_Lengths"][bvIndex]
                gtfFrame.at[index,"Matches"] = gtfFrame.at[index,"Matches"][bvIndex]
                gtfFrame.at[index,"Perfect_Matches"] = gtfFrame.at[index,"Perfect_Matches"][bvIndex]
                gtfFrame.at[index,"Fit_Value"] = gtfFrame.at[index,"Fit_Value"][bvIndex]
                
    return gtfFrame, oneExonSize

##Unique Matches From Best Matched GTF File (returns a number count of how many transcipts are not duplicates)

##sumFrame --- a matched GTF dataframe which has already been filtered for best matches with bestMatch


def uniqueMatches(sumFrame):
    
    sumFramedropped = sumFrame.dropna()
    sumFramedropped = sumFramedropped[sumFramedropped.astype(str)['Matched_Transcript'] != '[]']
    
    sumFramedropped = sumFramedropped.sort_values(by=['Fit_Value'])

    sumFramedropped = sumFramedropped.reset_index()
    sumFramedropped = sumFramedropped.drop(columns=['index'])
    sumFramedropped = sumFramedropped.drop_duplicates(subset=['Matched_Transcript'])
    
    return sumFramedropped.shape[0]

##Returns Total Number of Perfect Matches in a given GTF dataframe

##sumFrame --- a matched GTF file

def sumPerfectMatches(sumFrame):
    matchSum = 0
    for index, row in sumFrame.iterrows():
        if not sumFrame.at[index,"Perfect_Matches"]:
            pass
        else:
            if sum(sumFrame.at[index,"Perfect_Matches"]) > 0:
                matchSum = matchSum + 1
    return matchSum


##Total Number of Exons

##sumFrame --- GTF dataframe with an exon column

def sumofExons(sumFrame):
    matchSum = 0
    for index, row in sumFrame.iterrows():
        if str(sumFrame.at[index,'Exons']) == '[]':
            pass
        else:
            matchSum = matchSum + len(sumFrame.at[index,"Exons"])
    return matchSum

##Return a dataframe with only perfect matches

##sumFrame --- a matched gtf dataframe


def returnPerfectMatches(sumFrame):
    indices = []
    for index, row in sumFrame.iterrows():
        if not sumFrame.at[index,"Perfect_Matches"]:
            pass
        else:
            i=0
            while i < len(sumFrame.at[index,"Perfect_Matches"]):
                if sumFrame.at[index,"Perfect_Matches"][i] == 1:
                    indices.append(index)
                    i = len(sumFrame.at[index,"Perfect_Matches"])
                else:
                    i = i + 1
    sumFrame2 = sumFrame.iloc[indices, ]
    return sumFrame2

##Return a dataframe with only multiple perfect matches

##sumFrame --- a matched gtf dataframe

def returnMultiPerfectMatches(sumFrame):
    indices = []
    sumFrame['Perfect_Starts'] = np.empty((len(sumFrame), 0)).tolist()
    sumFrame['Perfect_Stops'] = np.empty((len(sumFrame), 0)).tolist()
    for index, row in sumFrame.iterrows():
        if not sumFrame.at[index,"Perfect_Matches"]:
            pass
        else:
            i=0
            k=0
            first = False
            while i < len(sumFrame.at[index,"Perfect_Matches"]):
                if sumFrame.at[index,"Perfect_Matches"][i] == 1:
                    if sumFrame.at[index,"Matched_Transcript"][i] !=  sumFrame.at[index,"transcript_id"]:
                        sumFrame.at[index,'Perfect_Starts'].append(sumFrame.at[index,"Start_Stop_Distance"][i][0])
                        sumFrame.at[index,'Perfect_Stops'].append(sumFrame.at[index,"Start_Stop_Distance"][i][1])
                    if first == True and k == 0:
                        indices.append(index)
                        k=1
                    first=True
                    i = i + 1
                else:
                    i = i + 1
    sumFrame2 = sumFrame.iloc[indices, ]
    return sumFrame2

##Returns the number of matched Splice Junctions (matched can occur across matched transcripts i.e: there are two splice junctions, in one matched transcript only the first is matched, in another matched transcript the second splice junction is matched, a matched SJ value of 2 will be returned since both were matched),
##also returns total number of splice junctions in the given GTF

##sumFrame --- a matched gtf dataframe

def sumTotalSJMatched(sumFrame):
    matchSum = 0
    masterSize = 0
    for index, row in sumFrame.iterrows():
        if not sumFrame.at[index,"Splice_Junction_Matches"]:
            pass
        else:
            list = sumFrame.at[index,"Splice_Junction_Matches"][0]
            masterSize = masterSize + len(sumFrame.at[index,"Splice_Junction_Matches"][0])
            i=0
            while i < len(sumFrame.at[index,"Splice_Junction_Matches"]):
                j=0
                while j < len(sumFrame.at[index,"Splice_Junction_Matches"][i]):
                    if sumFrame.at[index,"Splice_Junction_Matches"][i][j] == 1:
                        if list[j] == 0:
                            list[j] = 1
                    j = j+1
                i=i+1
            
            matchSum = matchSum + sum(list)
    return matchSum, masterSize

##Splice Junction Column Creator: Given an exon column it makes a new splice junction column based on these exons

##sumFrame --- a GTF file that is not yet matched, and has already run through collapsedReturns

def splicJuncColumn(sumFrame):
    sumFrame['Splice_Junctions'] = np.empty((len(sumFrame), 0)).tolist()
    for index, row in sumFrame.iterrows():
        if not sumFrame.at[index,"Exons"]:
            pass
        else:
            if len(sumFrame.at[index,"Exons"]) == 1:
                pass
            elif len(sumFrame.at[index,"Exons"]) > 1:
                i = 0
                while i < len(sumFrame.at[index,"Exons"]) - 1:
                    start = sumFrame.at[index,"Exons"][i][1]
                    end = sumFrame.at[index,"Exons"][i+1][0]
                    sumFrame.at[index,"Splice_Junctions"].append([start,end])
                    i = i + 1
            
    return(sumFrame)

##Splice Junction Missing: Returns 4 numbers which operate as bins. If the first splice junction was unmatched then count1 increases by 1, if any splice junction between (and including) the second splice junction and the middle splice junction are missing count2 increases by 1, and so on where count3 is middle + 1 to second to last, and count4 is the last splice junction. This accounts for negative strand reads by flipping the values.

##sumFrame --- a gtf dataframe which has already gone throuhg best match


def splicJuncMissing(sumFrame):
    indices=[]
    
    #create new index to help itteration
    sumFrame = sumFrame.reset_index()
    sumFrame = sumFrame.drop(columns=['index'])
    
    sumFrame['Splice_Misses'] = np.zeros((len(sumFrame), 4),dtype=int).tolist()
    for index, row in sumFrame.iterrows():
        if not sumFrame.at[index,"Splice_Junction_Matches"]:
            pass
        elif str(sumFrame.at[index,"Splice_Junction_Matches"]) == "[]":
            pass
        elif len(sumFrame.at[index,"Splice_Junction_Matches"]) == 1 or len(sumFrame.at[index,"Splice_Junction_Matches"]) == 0:
            pass
        elif len(sumFrame.at[index,"Splice_Junction_Matches"]) == 2:
            if sumFrame.at[index,"Splice_Junction_Matches"][0] == 0:
                sumFrame.at[index,"Splice_Misses"][0] = 1
            if sumFrame.at[index,"Splice_Junction_Matches"][1] == 0:
                sumFrame.at[index,"Splice_Misses"][3] = 1
            indices.append(index)
        elif len(sumFrame.at[index,"Splice_Junction_Matches"]) % 2 == 0:
            if sumFrame.at[index,"Splice_Junction_Matches"][0] == 0:
                sumFrame.at[index,"Splice_Misses"][0] = 1
            if sumFrame.at[index,"Splice_Junction_Matches"][len(sumFrame.at[index,"Splice_Junction_Matches"])-1] == 0:
                sumFrame.at[index,"Splice_Misses"][3] = 1
            i = 1
            while i < len(sumFrame.at[index,"Splice_Junction_Matches"])/2:
                if sumFrame.at[index,"Splice_Junction_Matches"][i] == 0:
                    sumFrame.at[index,"Splice_Misses"][1] = 1
                if sumFrame.at[index,"Splice_Junction_Matches"][len(sumFrame.at[index,"Splice_Junction_Matches"])-i] == 0:
                    sumFrame.at[index,"Splice_Misses"][2] = 1
                i = i + 1
            indices.append(index)
        elif len(sumFrame.at[index,"Splice_Junction_Matches"]) % 2 != 0:
            if sumFrame.at[index,"Splice_Junction_Matches"][0] == 0:
                sumFrame.at[index,"Splice_Misses"][0] = 1
            if sumFrame.at[index,"Splice_Junction_Matches"][len(sumFrame.at[index,"Splice_Junction_Matches"])-1] == 0:
                sumFrame.at[index,"Splice_Misses"][3] = 1
            if sumFrame.at[index,"Splice_Junction_Matches"][math.floor(len(sumFrame.at[index,"Splice_Junction_Matches"])/2)+1] == 0:
                sumFrame.at[index,"Splice_Misses"][1] = 1
            i = 1
            while i < math.floor(len(sumFrame.at[index,"Splice_Junction_Matches"])/2):
                if sumFrame.at[index,"Splice_Junction_Matches"][i] == 0:
                    sumFrame.at[index,"Splice_Misses"][1] = 1
                if sumFrame.at[index,"Splice_Junction_Matches"][len(sumFrame.at[index,"Splice_Junction_Matches"])-i] == 0:
                    sumFrame.at[index,"Splice_Misses"][2] = 1
                i = i + 1
            indices.append(index)

    sumFrame2 = sumFrame.iloc[indices,]
    count1 = 0
    count2 = 0
    count3 = 0
    count4 = 0
    for index, row in sumFrame2.iterrows():
        if sumFrame2.at[index, "strand"] == "-":
            sumFrame2.at[index,"Splice_Misses"] = sumFrame2.at[index,"Splice_Misses"][::-1]
        if sumFrame2.at[index,"Splice_Misses"][0] == 1:
            count1 = count1 + 1
        if sumFrame2.at[index,"Splice_Misses"][1] == 1:
            count2 = count2 + 1
        if sumFrame2.at[index,"Splice_Misses"][2] == 1:
            count3 = count3 + 1
        if sumFrame2.at[index,"Splice_Misses"][3] == 1:
            count4 = count4 + 1
    
            
    return count1, count2, count3, count4


##GTF Merge: Combines and collapses given GTF files into one final dataframe

##gtfFrames --- GTF Files, requires at least two


def gtfMerge(gtfFrame1, gtfFrame2, *argv):
    columnNames = gtfFrame1.columns.tolist()

    #create new index to help itteration
    gtfFrame1 = gtfFrame1.reset_index()
    gtfFrame1 = gtfFrame1.drop(columns=['index'])

    gtfFrame2 = gtfFrame2.reset_index()
    gtfFrame2 = gtfFrame2.drop(columns=['index'])
    
    if len(gtfFrame1.columns.tolist()) >= len(gtfFrame2.columns.tolist()):
        columnNames = gtfFrame1.columns.tolist()
    else:
        columnNames = gtfFrame2.columns.tolist()
    columnNames.append("Exons")
    
    print("Collapsing GTF #1:")
    gtfFrame1 = collapsedReturns(gtfFrame1)
    
    print("Collapsing GTF #2:")
    gtfFrame2 = collapsedReturns(gtfFrame2)
    
    gtfFrame1 = gtfFrame1[(gtfFrame1["feature"] == "gene") | (gtfFrame1["feature"] == "transcript")]
    gtfFrame2 = gtfFrame2[(gtfFrame2["feature"] == "gene") | (gtfFrame2["feature"] == "transcript")]
    
    resultFrame = gtfFrame1.append(gtfFrame2)
    
    i=1
    for arg in argv:
        
        print("Collapsing GTF #" + str(i+2) + ":")
        gtfFrameArg = arg.reset_index()
        gtfFrameArg = gtfFrameArg.drop(columns=['index'])
        gtfFrameArg = collapsedReturns(gtfFrameArg)
        gtfFrameArg = gtfFrameArg[(gtfFrameArg["feature"] == "gene") | (gtfFrameArg["feature"] == "transcript")]
        resultFrame = resultFrame.append(gtfFrameArg)
        i=i+1
    
    resultFrame = resultFrame[columnNames]
    resultFrame = resultFrame.sort_values(by=['start', 'end'])
    
    resultFrame = resultFrame.reset_index()
    resultFrame = resultFrame.drop(columns=['index'])
    
    return(resultFrame)

##Reshape GTF: Removes old columns so that a GTF file can be used in symbiosis again

##gtfFrame1 -- gtf file which has already gone through Symbiosis

def reshapeGTF(gtfFrame1):
    gtfFrame1 = gtfFrame1.drop(columns=['Matched_Gene', 'Matched_Transcript', 'Start_Stop_Distance', 'Splice_Junction_Matches', 'Gene_Exon_Lengths', 'Gene_SJ_Lengths', 'Matches', 'Perfect_Matches', "Splice_Junctions"], errors = 'ignore')
    return(gtfFrame1)

##Start Stop Crop (return dataframe where only perfect matches under a certain start and stop distance are accepted)

##sumFrame --- a GTF dataframe which has already been matcehd

##cropDistance --- a distance which you desire to cut off start stop distances at

def startStopCrop(sumFrame, cropDistance):
    indices = []
    for index, row in sumFrame.iterrows():
        if not sumFrame.at[index,"Perfect_Matches"]:
            pass
        else:
            i=0
            while i < len(sumFrame.at[index,"Perfect_Matches"]):
                if sumFrame.at[index,"Perfect_Matches"][i] == 1:
                    k=0
                    for item in sumFrame.at[index,"Start_Stop_Distance"]:
                        if k==0 and abs(item[0]) <= cropDistance and abs(item[1]) <= cropDistance:
                            indices.append(index)
                            i = len(sumFrame.at[index,"Perfect_Matches"])
                            k=1
                    i = i + 1
                else:
                    i = i + 1
    sumFrame2 = sumFrame.iloc[indices, ]
    return sumFrame2

##Write GTF: saves a GTF dataframe as a gtf file

## name --- the name you wish to name the file as

## givenGTF --- the gtf file which you want to save

def writeGTF(name, givenGTF):
    correctColumns = ["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "gene_id", "transcirpt_id" ]
    givenColList = givenGTF.columns.tolist()
    extraCol = []
    
    for item in givenColList:
        if item not in correctColumns:
            extraCol.append(item)

    correctColumns.extend(extraCol)
    givenGTF = givenGTF[correctColumns]
    
    with open(name,'w') as fh:
        for index, row in givenGTF.iterrows():
            string = ""
            k = 0
            for column in givenGTF:
                if k <= 7:
                    k = k + 1
                    string = string + str(givenGTF.at[index, column]) + "\t"
                else:
                    k = k + 1
                    string = string + str(column) + " " + "\"" + str(givenGTF.at[index, column]) + "\"" + "; "
            string = string + "\n"
            fh.write(string)
                        
##Symbiosis Master Function: the primary function of Symbiosis which compares two GTF Files            
          
##main gtf --- the GTF dataframe which will be the basis of comparision

##refGTF -- the GTF dataframe which the main GTf will be compared against, either a reference, or another read file (this is marked by the inTissue flag)

##inTissue --- if False, refGTF is a reference GTF like Ensembl, if True then it is another sample, and main and refGTf will be combined in the end results (i.e. transcripts not matched in refGTF will be added to MainGTF to create a full construct of all transcripts)

##premerged1 --- if True, GTF has already been collapsed

##premerged2 --- if True, GTF has already been collapsed
            
def symbiosisMasterFunction(mainGTF,refGTF, threshold, inTissue = False, preMerged1 = False, preMerged2 = False):
    ###Create an array of each GTFs chromosones
    chromosones1 = mainGTF.seqname.unique().tolist()
    chromosones2 = refGTF.seqname.unique().tolist()

    ##Counter for testing certain properties
    counter = 0

    #Create final dataframe
    MasterFrameGeneMatch = pd.DataFrame()
    MasterFrameGeneMatch2 = pd.DataFrame()

    #Check if chromosones don't match
    if sorted(chromosones1) != sorted(chromosones2):
        print("Non-matching chromosones")

    print("Combining Dataframes:")

    for item in chromosones1:

            #Subset our dataframe by the chromosone
            chromeFrame1 = mainGTF[mainGTF["seqname"] == item]
            chromeFrame2 = refGTF[refGTF["seqname"] == item]
            
            ##Create a list of Column Names to resort at end
            columnNames1 = chromeFrame1.columns.tolist()
            columnNames2 = chromeFrame2.columns.tolist()

            #create new index to help itteration
            chromeFrame1 = chromeFrame1.reset_index()
            chromeFrame1 = chromeFrame1.drop(columns=['index'])

            chromeFrame2 = chromeFrame2.reset_index()
            chromeFrame2 = chromeFrame2.drop(columns=['index'])
            
            
            if preMerged1 == False:
                #Create Exon Column to compare
                chromeFrame1 = collapsedReturns(chromeFrame1)
                
                
            if preMerged2 == False:
                #Create Exon Column to compare
                chromeFrame2 = collapsedReturns(chromeFrame2)
          


            #Look at only transcripts from file we are evaluating and unique transcripts and all genes from reference GTF
            chromeFrame1 = chromeFrame1[chromeFrame1["feature"] == "transcript"]
            chromeFrame2 = chromeFrame2[(chromeFrame2["feature"] == "gene") | (chromeFrame2["feature"] == "transcript")]
            for index, row in chromeFrame2.iterrows():
                if chromeFrame2.at[index, "feature"] == "transcript":
                    if 'transcript_name' in chromeFrame2.columns:
                        chromeFrame2.at[index, "gene_id"] = chromeFrame2.at[index, "transcript_name"]
                    elif 'transcript_id' in chromeFrame2.columns:
                        chromeFrame2.at[index, "gene_id"] = chromeFrame2.at[index, "transcript_id"]
            
            if preMerged1 == False:
                #Create Exon Column to sort
                columnNames1.append("Exons")
            if preMerged2 == False:
                #Create Exon Column to compare
                columnNames2.append("Exons")
            chromeFrame1 = chromeFrame1[columnNames1]
            chromeFrame2 = chromeFrame2[columnNames2]


            #Create Empty List Columns
            chromeFrame1['Matched_Gene'] = np.empty((len(chromeFrame1), 0)).tolist()
            chromeFrame1['Matched_Transcript'] = np.empty((len(chromeFrame1), 0)).tolist()
            chromeFrame1['Start_Stop_Distance'] = np.empty((len(chromeFrame1), 0)).tolist()
            chromeFrame1['Splice_Junction_Matches'] = np.empty((len(chromeFrame1), 0)).tolist()
            chromeFrame1['Gene_Exon_Lengths'] = np.empty((len(chromeFrame1), 0)).tolist()
            chromeFrame1['Gene_SJ_Lengths'] = np.empty((len(chromeFrame1), 0)).tolist()
            chromeFrame1['Splice_Junc_Column_2'] = np.empty((len(chromeFrame1), 0)).tolist()
            
            
            chromeFrame2['Matched_Gene'] = np.empty((len(chromeFrame2), 0)).tolist()
            chromeFrame2['Matched_Transcript'] = np.empty((len(chromeFrame2), 0)).tolist()
            chromeFrame2['Start_Stop_Distance'] = np.empty((len(chromeFrame2), 0)).tolist()
            chromeFrame2['Splice_Junction_Matches'] = np.empty((len(chromeFrame2), 0)).tolist()
            chromeFrame2['Gene_Exon_Lengths'] = np.empty((len(chromeFrame2), 0)).tolist()
            chromeFrame2['Gene_SJ_Lengths'] = np.empty((len(chromeFrame2), 0)).tolist()
            chromeFrame2['Splice_Junc_Column_2'] = np.empty((len(chromeFrame2), 0)).tolist()
            chromeFrame2['found'] = False
            
            
            ##chromeFrame1['Overlap Percentage'] = np.empty((len(chromeFrame1), 0)).tolist()

            #Sort values based on their start position
            chromeFrame1 = chromeFrame1.sort_values(by=['start'])
            chromeFrame2 = chromeFrame2.sort_values(by=['start'])
            
            
            #create new index to help itteration
            chromeFrame1 = chromeFrame1.reset_index()
            chromeFrame1 = chromeFrame1.drop(columns=['index'])

            chromeFrame2 = chromeFrame2.reset_index()
            chromeFrame2 = chromeFrame2.drop(columns=['index'])
            
            ##Loop through first dataframe and check for correpsonding genes in the second
            itterator = 0
            for index, row in chromeFrame2.iterrows():
                Exons2 = chromeFrame2.at[index,"Exons"]
                chromeFrame2.at[index,"Exons"] = sorted(Exons2)
            for index, row in chromeFrame1.iterrows():
                Exons1 = chromeFrame1.at[index,"Exons"]
                chromeFrame1.at[index,"Exons"] = sorted(Exons1)
                
            #Add Splice Junction Columns
            chromeFrame1 = splicJuncColumn(chromeFrame1)
            chromeFrame2 = splicJuncColumn(chromeFrame2)
            
            
            for index, row in chromeFrame1.iterrows():
                

                #Move Along Rows in the Reference Until It Is Past Beginning of Current read
                if itterator != len(chromeFrame2.index) and row[3] > chromeFrame2.at[itterator,"end"]:
                    while itterator < len(chromeFrame2.index) and row[3] > chromeFrame2.at[itterator,"end"]:
                        itterator = itterator + 1
                
                ## Pass If we have reached the end of the reference genome or if the current gene is ahead of our current read
                if itterator == len(chromeFrame2.index) or row[4] < chromeFrame2.at[itterator,"start"]:
                    pass
                elif row[3] >= chromeFrame2.at[itterator,"start"] and row[3] <= chromeFrame2.at[itterator,"end"]: ##If the start is contained within the transcript we know they overlap
                    if chromeFrame2.at[itterator,"feature"]=="gene":
                        chromeFrame1.at[index,"Matched_Gene"].append(chromeFrame2.at[itterator,"gene_id"])
                    elif chromeFrame2.at[itterator,"feature"]=="transcript":
                        chromeFrame1.at[index,"Gene_Exon_Lengths"].append(len(chromeFrame2.at[itterator,"Exons"]))
                        chromeFrame1.at[index,"Gene_SJ_Lengths"].append(len(chromeFrame2.at[itterator,"Splice_Junctions"]))
                        chromeFrame1.at[index,"Matched_Transcript"].append(chromeFrame2.at[itterator,"gene_id"])
                        chromeFrame1.at[index,"Start_Stop_Distance"].append(startStopDistance(chromeFrame1.at[index,"Exons"], chromeFrame2.at[itterator,"Exons"]))
                        chromeFrame1.at[index,"Splice_Junction_Matches"].append(exonMatching(chromeFrame1.at[index,"Splice_Junctions"], chromeFrame2.at[itterator,"Splice_Junctions"], threshold))

                        chromeFrame1.at[index,"Splice_Junc_Column_2"].append(chromeFrame2.at[itterator,"Splice_Junctions"])
                        chromeFrame2.at[itterator,"found"] = True
                                
                    
                    if(chromeFrame2.at[itterator,"end"] > row[4]):
                        compareValue = chromeFrame2.at[itterator,"end"]
                    else:
                        compareValue = row[4]
                    looping = 1
                    while itterator+looping < len(chromeFrame2.index) and compareValue > chromeFrame2.at[itterator+looping,"start"]:
                        if itterator+looping == len(chromeFrame2.index) or row[4] < chromeFrame2.at[itterator+looping,"start"]:
                            pass
                        elif row[3] >= chromeFrame2.at[itterator+looping,"start"] and row[3] <= chromeFrame2.at[itterator+looping,"end"]:
                            if chromeFrame2.at[itterator+looping,"feature"]=="gene":
                                chromeFrame1.at[index,"Matched_Gene"].append(chromeFrame2.at[itterator+looping,"gene_id"])
                            elif chromeFrame2.at[itterator+looping,"feature"]=="transcript":
                                
                                chromeFrame1.at[index,"Gene_Exon_Lengths"].append(len(chromeFrame2.at[itterator+looping,"Exons"]))
                                chromeFrame1.at[index,"Gene_SJ_Lengths"].append(len(chromeFrame2.at[itterator+looping,"Splice_Junctions"]))
                                chromeFrame1.at[index,"Matched_Transcript"].append(chromeFrame2.at[itterator+looping,"gene_id"])
                                chromeFrame1.at[index,"Start_Stop_Distance"].append(startStopDistance(chromeFrame1.at[index,"Exons"], chromeFrame2.at[itterator+looping,"Exons"]))
                                chromeFrame1.at[index,"Splice_Junction_Matches"].append(exonMatching(chromeFrame1.at[index,"Splice_Junctions"], chromeFrame2.at[itterator+looping,"Splice_Junctions"], threshold))

                                chromeFrame1.at[index,"Splice_Junc_Column_2"].append(chromeFrame2.at[itterator+looping,"Splice_Junctions"])
                                chromeFrame2.at[itterator+looping,"found"] = True
                            
                           
                        elif row[4] >= chromeFrame2.at[itterator+looping,"start"] and row[4] <= chromeFrame2.at[itterator+looping,"end"]:
                            if chromeFrame2.at[itterator+looping,"feature"]=="gene":
                                chromeFrame1.at[index,"Matched_Gene"].append(chromeFrame2.at[itterator+looping,"gene_id"])
                            elif chromeFrame2.at[itterator+looping,"feature"]=="transcript":
                                chromeFrame1.at[index,"Gene_Exon_Lengths"].append(len(chromeFrame2.at[itterator+looping,"Exons"]))
                                chromeFrame1.at[index,"Gene_SJ_Lengths"].append(len(chromeFrame2.at[itterator+looping,"Splice_Junctions"]))
                                chromeFrame1.at[index,"Matched_Transcript"].append(chromeFrame2.at[itterator+looping,"gene_id"])
                                chromeFrame1.at[index,"Start_Stop_Distance"].append(startStopDistance(chromeFrame1.at[index,"Exons"], chromeFrame2.at[itterator+looping,"Exons"]))
                                chromeFrame1.at[index,"Splice_Junction_Matches"].append(exonMatching(chromeFrame1.at[index,"Splice_Junctions"], chromeFrame2.at[itterator+looping,"Splice_Junctions"], threshold))

                                chromeFrame1.at[index,"Splice_Junc_Column_2"].append(chromeFrame2.at[itterator+looping,"Splice_Junctions"])
                                chromeFrame2.at[itterator+looping,"found"] = True
                                
                            
                        elif row[3] <= chromeFrame2.at[itterator+looping,"start"] and row[4] >= chromeFrame2.at[itterator+looping,"end"]:
                            if chromeFrame2.at[itterator+looping,"feature"]=="gene":
                                chromeFrame1.at[index,"Matched_Gene"].append(chromeFrame2.at[itterator+looping,"gene_id"])
                            elif chromeFrame2.at[itterator+looping,"feature"]=="transcript":
                                chromeFrame1.at[index,"Gene_Exon_Lengths"].append(len(chromeFrame2.at[itterator+looping,"Exons"]))
                                chromeFrame1.at[index,"Gene_SJ_Lengths"].append(len(chromeFrame2.at[itterator+looping,"Splice_Junctions"]))
                                chromeFrame1.at[index,"Matched_Transcript"].append(chromeFrame2.at[itterator+looping,"gene_id"])
                                chromeFrame1.at[index,"Start_Stop_Distance"].append(startStopDistance(chromeFrame1.at[index,"Exons"], chromeFrame2.at[itterator+looping,"Exons"]))
                                chromeFrame1.at[index,"Splice_Junction_Matches"].append(exonMatching(chromeFrame1.at[index,"Splice_Junctions"], chromeFrame2.at[itterator+looping,"Splice_Junctions"], threshold))

                                chromeFrame1.at[index,"Splice_Junc_Column_2"].append(chromeFrame2.at[itterator+looping,"Splice_Junctions"])
                                chromeFrame2.at[itterator+looping,"found"] = True
                            
                        looping = looping+1
                elif row[4] >= chromeFrame2.at[itterator,"start"] and row[4] <= chromeFrame2.at[itterator,"end"]: #If the end is contained withing the transcript we know they overlap
                    if chromeFrame2.at[itterator,"feature"]=="gene":
                        chromeFrame1.at[index,"Matched_Gene"].append(chromeFrame2.at[itterator,"gene_id"])
                    elif chromeFrame2.at[itterator,"feature"]=="transcript":
                        chromeFrame1.at[index,"Gene_Exon_Lengths"].append(len(chromeFrame2.at[itterator,"Exons"]))
                        chromeFrame1.at[index,"Gene_SJ_Lengths"].append(len(chromeFrame2.at[itterator,"Splice_Junctions"]))
                        chromeFrame1.at[index,"Matched_Transcript"].append(chromeFrame2.at[itterator,"gene_id"])
                        chromeFrame1.at[index,"Start_Stop_Distance"].append(startStopDistance(chromeFrame1.at[index,"Exons"], chromeFrame2.at[itterator,"Exons"]))
                        chromeFrame1.at[index,"Splice_Junction_Matches"].append(exonMatching(chromeFrame1.at[index,"Splice_Junctions"], chromeFrame2.at[itterator,"Splice_Junctions"], threshold))

                        chromeFrame1.at[index,"Splice_Junc_Column_2"].append(chromeFrame2.at[itterator,"Splice_Junctions"])
                        chromeFrame2.at[itterator,"found"] = True
                    
                        
                    if(chromeFrame2.at[itterator,"end"] > row[4]):
                        compareValue = chromeFrame2.at[itterator,"end"]
                    else:
                        compareValue = row[4]
                    looping = 1
                    while itterator+looping < len(chromeFrame2.index) and compareValue > chromeFrame2.at[itterator+looping,"start"]:
                        if itterator+looping == len(chromeFrame2.index) or row[4] < chromeFrame2.at[itterator+looping,"start"]:
                            pass
                        elif row[3] >= chromeFrame2.at[itterator+looping,"start"] and row[3] <= chromeFrame2.at[itterator+looping,"end"]:
                            if chromeFrame2.at[itterator+looping,"feature"]=="gene":
                                chromeFrame1.at[index,"Matched_Gene"].append(chromeFrame2.at[itterator+looping,"gene_id"])
                            elif chromeFrame2.at[itterator+looping,"feature"]=="transcript":
                            
                                chromeFrame1.at[index,"Gene_Exon_Lengths"].append(len(chromeFrame2.at[itterator+looping,"Exons"]))
                                chromeFrame1.at[index,"Gene_SJ_Lengths"].append(len(chromeFrame2.at[itterator+looping,"Splice_Junctions"]))
                                chromeFrame1.at[index,"Matched_Transcript"].append(chromeFrame2.at[itterator+looping,"gene_id"])
                                chromeFrame1.at[index,"Start_Stop_Distance"].append(startStopDistance(chromeFrame1.at[index,"Exons"], chromeFrame2.at[itterator+looping,"Exons"]))
                                chromeFrame1.at[index,"Splice_Junction_Matches"].append(exonMatching(chromeFrame1.at[index,"Splice_Junctions"], chromeFrame2.at[itterator+looping,"Splice_Junctions"], threshold))

                                chromeFrame1.at[index,"Splice_Junc_Column_2"].append(chromeFrame2.at[itterator+looping,"Splice_Junctions"])
                                chromeFrame2.at[itterator+looping,"found"] = True
                                
                           
                        elif row[4] >= chromeFrame2.at[itterator+looping,"start"] and row[4] <= chromeFrame2.at[itterator+looping,"end"]:
                            if chromeFrame2.at[itterator+looping,"feature"]=="gene":
                                chromeFrame1.at[index,"Matched_Gene"].append(chromeFrame2.at[itterator+looping,"gene_id"])
                            elif chromeFrame2.at[itterator+looping,"feature"]=="transcript":

                                chromeFrame1.at[index,"Gene_Exon_Lengths"].append(len(chromeFrame2.at[itterator+looping,"Exons"]))
                                chromeFrame1.at[index,"Gene_SJ_Lengths"].append(len(chromeFrame2.at[itterator+looping,"Splice_Junctions"]))
                                chromeFrame1.at[index,"Matched_Transcript"].append(chromeFrame2.at[itterator+looping,"gene_id"])
                                chromeFrame1.at[index,"Start_Stop_Distance"].append(startStopDistance(chromeFrame1.at[index,"Exons"], chromeFrame2.at[itterator+looping,"Exons"]))
                                chromeFrame1.at[index,"Splice_Junction_Matches"].append(exonMatching(chromeFrame1.at[index,"Splice_Junctions"], chromeFrame2.at[itterator+looping,"Splice_Junctions"], threshold))

                                chromeFrame1.at[index,"Splice_Junc_Column_2"].append(chromeFrame2.at[itterator+looping,"Splice_Junctions"])
                                chromeFrame2.at[itterator+looping,"found"] = True
                                
                        elif row[3] <= chromeFrame2.at[itterator+looping,"start"] and row[4] >= chromeFrame2.at[itterator+looping,"end"]:
                            if chromeFrame2.at[itterator+looping,"feature"]=="gene":
                                chromeFrame1.at[index,"Matched_Gene"].append(chromeFrame2.at[itterator+looping,"gene_id"])
                            elif chromeFrame2.at[itterator+looping,"feature"]=="transcript":
                                chromeFrame1.at[index,"Gene_Exon_Lengths"].append(len(chromeFrame2.at[itterator+looping,"Exons"]))
                                chromeFrame1.at[index,"Gene_SJ_Lengths"].append(len(chromeFrame2.at[itterator+looping,"Splice_Junctions"]))
                                chromeFrame1.at[index,"Matched_Transcript"].append(chromeFrame2.at[itterator+looping,"gene_id"])
                                chromeFrame1.at[index,"Start_Stop_Distance"].append(startStopDistance(chromeFrame1.at[index,"Exons"], chromeFrame2.at[itterator+looping,"Exons"]))
                                chromeFrame1.at[index,"Splice_Junction_Matches"].append(exonMatching(chromeFrame1.at[index,"Splice_Junctions"], chromeFrame2.at[itterator+looping,"Splice_Junctions"], threshold))

                                chromeFrame1.at[index,"Splice_Junc_Column_2"].append(chromeFrame2.at[itterator+looping,"Splice_Junctions"])
                                chromeFrame2.at[itterator+looping,"found"] = True
                               
                        looping = looping+1
                elif row[3] <= chromeFrame2.at[itterator,"start"] and row[4] >= chromeFrame2.at[itterator,"end"]: #If the read encapsulates our transcript
                    if chromeFrame2.at[itterator,"feature"]=="gene":
                        chromeFrame1.at[index,"Matched_Gene"].append(chromeFrame2.at[itterator,"gene_id"])
                    elif chromeFrame2.at[itterator,"feature"]=="transcript":
                        chromeFrame1.at[index,"Gene_Exon_Lengths"].append(len(chromeFrame2.at[itterator,"Exons"]))
                        chromeFrame1.at[index,"Gene_SJ_Lengths"].append(len(chromeFrame2.at[itterator,"Splice_Junctions"]))
                        chromeFrame1.at[index,"Matched_Transcript"].append(chromeFrame2.at[itterator,"gene_id"])
                        chromeFrame1.at[index,"Start_Stop_Distance"].append(startStopDistance(chromeFrame1.at[index,"Exons"], chromeFrame2.at[itterator,"Exons"]))
                        chromeFrame1.at[index,"Splice_Junction_Matches"].append(exonMatching(chromeFrame1.at[index,"Splice_Junctions"], chromeFrame2.at[itterator,"Splice_Junctions"], threshold))

                        chromeFrame1.at[index,"Splice_Junc_Column_2"].append(chromeFrame2.at[itterator,"Splice_Junctions"])
                        chromeFrame2.at[itterator,"found"] = True
                        
                    if(chromeFrame2.at[itterator,"end"] > row[4]):
                        compareValue = chromeFrame2.at[itterator,"end"]
                    else:
                        compareValue = row[4]
                    looping = 1
                    while itterator+looping < len(chromeFrame2.index) and compareValue > chromeFrame2.at[itterator+looping,"start"]:
                        if itterator+looping == len(chromeFrame2.index) or row[4] < chromeFrame2.at[itterator+looping,"start"]:
                            pass
                        elif row[3] >= chromeFrame2.at[itterator+looping,"start"] and row[3] <= chromeFrame2.at[itterator+looping,"end"]:
                            if chromeFrame2.at[itterator+looping,"feature"]=="gene":
                                chromeFrame1.at[index,"Matched_Gene"].append(chromeFrame2.at[itterator+looping,"gene_id"])
                            elif chromeFrame2.at[itterator+looping,"feature"]=="transcript":
                                chromeFrame1.at[index,"Gene_Exon_Lengths"].append(len(chromeFrame2.at[itterator+looping,"Exons"]))
                                chromeFrame1.at[index,"Gene_SJ_Lengths"].append(len(chromeFrame2.at[itterator+looping,"Splice_Junctions"]))
                                chromeFrame1.at[index,"Matched_Transcript"].append(chromeFrame2.at[itterator+looping,"gene_id"])
                                chromeFrame1.at[index,"Start_Stop_Distance"].append(startStopDistance(chromeFrame1.at[index,"Exons"], chromeFrame2.at[itterator+looping,"Exons"]))
                                chromeFrame1.at[index,"Splice_Junction_Matches"].append(exonMatching(chromeFrame1.at[index,"Splice_Junctions"], chromeFrame2.at[itterator+looping,"Splice_Junctions"], threshold))

                                chromeFrame1.at[index,"Splice_Junc_Column_2"].append(chromeFrame2.at[itterator+looping,"Splice_Junctions"])
                                chromeFrame2.at[itterator+looping,"found"] = True
                            
                           
                        elif row[4] >= chromeFrame2.at[itterator+looping,"start"] and row[4] <= chromeFrame2.at[itterator+looping,"end"]:
                            if chromeFrame2.at[itterator+looping,"feature"]=="gene":
                                chromeFrame1.at[index,"Matched_Gene"].append(chromeFrame2.at[itterator+looping,"gene_id"])
                            elif chromeFrame2.at[itterator+looping,"feature"]=="transcript":
                                chromeFrame1.at[index,"Gene_Exon_Lengths"].append(len(chromeFrame2.at[itterator+looping,"Exons"]))
                                chromeFrame1.at[index,"Gene_SJ_Lengths"].append(len(chromeFrame2.at[itterator+looping,"Splice_Junctions"]))
                                chromeFrame1.at[index,"Matched_Transcript"].append(chromeFrame2.at[itterator+looping,"gene_id"])
                                chromeFrame1.at[index,"Start_Stop_Distance"].append(startStopDistance(chromeFrame1.at[index,"Exons"], chromeFrame2.at[itterator+looping,"Exons"]))
                                chromeFrame1.at[index,"Splice_Junction_Matches"].append(exonMatching(chromeFrame1.at[index,"Splice_Junctions"], chromeFrame2.at[itterator+looping,"Splice_Junctions"], threshold))

                                chromeFrame1.at[index,"Splice_Junc_Column_2"].append(chromeFrame2.at[itterator+looping,"Splice_Junctions"])
                                chromeFrame2.at[itterator+looping,"found"] = True
                               
                        elif row[3] <= chromeFrame2.at[itterator+looping,"start"] and row[4] >= chromeFrame2.at[itterator+looping,"end"]:
                            if chromeFrame2.at[itterator+looping,"feature"]=="gene":
                                chromeFrame1.at[index,"Matched_Gene"].append(chromeFrame2.at[itterator+looping,"gene_id"])
                            elif chromeFrame2.at[itterator+looping,"feature"]=="transcript":
                                chromeFrame1.at[index,"Gene_Exon_Lengths"].append(len(chromeFrame2.at[itterator+looping,"Exons"]))
                                chromeFrame1.at[index,"Gene_SJ_Lengths"].append(len(chromeFrame2.at[itterator+looping,"Splice_Junctions"]))
                                chromeFrame1.at[index,"Matched_Transcript"].append(chromeFrame2.at[itterator+looping,"gene_id"])
                                chromeFrame1.at[index,"Start_Stop_Distance"].append(startStopDistance(chromeFrame1.at[index,"Exons"], chromeFrame2.at[itterator+looping,"Exons"]))
                                chromeFrame1.at[index,"Splice_Junction_Matches"].append(exonMatching(chromeFrame1.at[index,"Splice_Junctions"], chromeFrame2.at[itterator+looping,"Splice_Junctions"], threshold))

                                chromeFrame1.at[index,"Splice_Junc_Column_2"].append(chromeFrame2.at[itterator+looping,"Splice_Junctions"])
                                chromeFrame2.at[itterator+looping,"found"] = True
                                
                            
                        looping = looping+1

            MasterFrameGeneMatch = MasterFrameGeneMatch.append(chromeFrame1) #add our matched data to the master dataframe
            MasterFrameGeneMatch2 = MasterFrameGeneMatch2.append(chromeFrame2)
            print(item,end = ' ') #Display each chromosone as we finish it
    
    ##Add Non-matched in other sample if it's inTissue
    if inTissue == True:
        nonMatchedMaster2 = MasterFrameGeneMatch2[MasterFrameGeneMatch2["found"] == False]
        nonMatchedMaster2 = nonMatchedMaster2[nonMatchedMaster2["feature"] == "transcript"]
        nonMatchedMaster2.drop(columns=['found'])


        columnSimple = MasterFrameGeneMatch.columns.tolist()
        nonMatchedMaster2 = nonMatchedMaster2[columnSimple]


        MasterFrameGeneMatch = MasterFrameGeneMatch.append(nonMatchedMaster2)
        
    MasterFrameGeneMatch = MasterFrameGeneMatch.reset_index()
    MasterFrameGeneMatch = MasterFrameGeneMatch.drop(columns=['index'])
    
    ##Create A Column For Number Of Matches and Perfect Matches
    MasterFrameGeneMatch['Matches'] = np.empty((len(MasterFrameGeneMatch), 0)).tolist()
    MasterFrameGeneMatch['Perfect_Matches'] = np.empty((len(MasterFrameGeneMatch), 0)).tolist()
    for index, row in MasterFrameGeneMatch.iterrows():
        if not MasterFrameGeneMatch.at[index,"Splice_Junction_Matches"]:
            pass
        else:
            i=0
            while i < len(MasterFrameGeneMatch.at[index,"Splice_Junction_Matches"]):
                MasterFrameGeneMatch.at[index, "Matches"].append(sum(MasterFrameGeneMatch.at[index,"Splice_Junction_Matches"][i]))
                i = i + 1
                
        if not MasterFrameGeneMatch.at[index,"Splice_Junction_Matches"]:
            pass
        else:
            i=0
            while i < len(MasterFrameGeneMatch.at[index,"Splice_Junction_Matches"]):
                if MasterFrameGeneMatch.at[index, "Matches"][i] == MasterFrameGeneMatch.at[index,"Gene_SJ_Lengths"][i]:
                    if len(MasterFrameGeneMatch.at[index, "Splice_Junctions"]) == MasterFrameGeneMatch.at[index, "Gene_SJ_Lengths"][i]:
                        MasterFrameGeneMatch.at[index, "Perfect_Matches"].append(1)
                    else:
                        MasterFrameGeneMatch.at[index, "Perfect_Matches"].append(0)
                else:
                    MasterFrameGeneMatch.at[index, "Perfect_Matches"].append(0)
                i = i + 1
     
    return(MasterFrameGeneMatch)


##Transcriptome Master Function: Constructs a Transcriptome from a list of GTFs

##gtfList --- an array of gtfs which you wish to be combined

##gtfLabels --- an array of the associated name you wish to label each GTF

##threshold --- threshold you consider acceptable for distance between splice junction start and stops to be counted as a match

##distanceThreshold --- threshold you consider acceptable between start and stop sites to be counted as a match


def transcriptomeMasterFunction2(GTFList, GTFLabels, threshold, distanceThreshold):
    
    if len(GTFList) <= 1:
        print("At least two GTFs are required to assemble transcriptome")
        return
    else:
        i = 0
        chromosoneList = []
        while i < len(GTFList):
            chromosoneList.append(len(GTFList[i].seqname.unique().tolist())) 
            i = i + 1
        if len(set(chromosoneList)) > 1:
            print("Non-matching chromosones:")
            print("Transcriptome will be constructed using chromosones of first GTF")
            print(sorted(GTFList[0].seqname.unique().tolist()))
            print("\n")
            
    ###Create an array of each GTFs chromosones
    chromosones1 = GTFList[0].seqname.unique().tolist()

    #Create final dataframe
    MasterFrameGeneMatch = pd.DataFrame()

    print("Combining Dataframes:")

    
    
    for item in chromosones1:
        
        workingDataframe = pd.DataFrame()
        
        i = 0
        while i < len(GTFList):
            
            #Choose GTF to work with
            mainGTF = GTFList[i]
            
            #Select only the current Chromosone
            chromeFrame1 = mainGTF[mainGTF["seqname"] == item]
            
            ##Create a list of Column Names to resort at end
            columnNames1 = chromeFrame1.columns.tolist()

            #create new index to help itteration
            chromeFrame1 = chromeFrame1.reset_index()
            chromeFrame1 = chromeFrame1.drop(columns=['index'])

            #Create Exon Column to compare
            chromeFrame1 = collapsedReturns(chromeFrame1)         

            #Look at only transcripts
            chromeFrame1 = chromeFrame1[chromeFrame1["feature"] == "transcript"]
            for index, row in chromeFrame1.iterrows():
                if 'transcript_name' in chromeFrame1.columns:
                    chromeFrame1.at[index, "gene_id"] = chromeFrame1.at[index, "transcript_name"] + "/" + GTFLabels[i] + "-" + str(item)
                elif 'transcript_id' in chromeFrame1.columns:
                    chromeFrame1.at[index, "gene_id"] = chromeFrame1.at[index, "transcript_id"] + "/" + GTFLabels[i] + "-" + str(item)
            
            columnNames1.append("Exons")
            
            chromeFrame1 = chromeFrame1[columnNames1]

            #Add Neccesary Columns
            chromeFrame1['Matched_Transcript'] = np.empty((len(chromeFrame1), 0)).tolist()
            chromeFrame1['Start_Stop_Distance'] = np.empty((len(chromeFrame1), 0)).tolist()
            chromeFrame1['Located_In'] = [[GTFLabels[i]] for _ in range(len(chromeFrame1))]
            chromeFrame1['found'] = False
            chromeFrame1['Shared_SJ_Transcript'] = np.empty((len(chromeFrame1), 0)).tolist()
            chromeFrame1['Shared_SJ_Weight'] = np.empty((len(chromeFrame1), 0)).tolist()
            chromeFrame1['Previously_Matched'] = False

            #Sort values based on their start position
            chromeFrame1 = chromeFrame1.sort_values(by=['start'])
            
            #create new index to help itteration
            chromeFrame1 = chromeFrame1.reset_index()
            chromeFrame1 = chromeFrame1.drop(columns=['index'])
            
            ##Ensure all exons are sorted
            for index, row in chromeFrame1.iterrows():
                Exons1 = chromeFrame1.at[index,"Exons"]
                chromeFrame1.at[index,"Exons"] = sorted(Exons1)
                
            #Add Splice Junction Column
            chromeFrame1 = splicJuncColumn(chromeFrame1)
            
            #Add to our working frame
            if i == 0:
                workingDataframe = chromeFrame1.copy()
            else:
                workingDataframe = workingDataframe.append(chromeFrame1)
            i = i + 1
        
        #Sort values based on their start position in our working frame
        workingDataframe = workingDataframe.sort_values(by=['start'])

        #create new index to help itteration in our working frame
        workingDataframe = workingDataframe.reset_index()
        workingDataframe = workingDataframe.drop(columns=['index'])
        
        for index, row in workingDataframe.iterrows():

            ##If the transcript has already been perfectly matched than skip it
            if workingDataframe.at[index, "found"] == False:
                
                ##If the end is greater then the start of the next index they overlap
                if index+1 < len(workingDataframe.index) and workingDataframe.at[index, "end"] >= workingDataframe.at[index+1,"start"]:

                    #Only match one exon length genes with one exon length genes
                    if len(workingDataframe.at[index,"Splice_Junctions"]) == 0 and len(workingDataframe.at[index+1,"Splice_Junctions"]) == 0:
                        if abs(workingDataframe.at[index,"start"] - workingDataframe.at[index+1,"start"])<=distanceThreshold and abs(workingDataframe.at[index,"end"] - workingDataframe.at[index+1,"end"])<=distanceThreshold:
                            workingDataframe.at[index,"Matched_Transcript"].append(workingDataframe.at[index+1,"gene_id"])
                            workingDataframe.at[index,"Start_Stop_Distance"].append([workingDataframe.at[index,"start"] - workingDataframe.at[index+1,"start"], workingDataframe.at[index,"end"] - workingDataframe.at[index+1,"end"]])
                            workingDataframe.at[index + 1,"found"] = True
                            if workingDataframe.at[index+1,"Located_In"][0] not in workingDataframe.at[index,"Located_In"]:
                                workingDataframe.at[index,"Located_In"].append(workingDataframe.at[index+1,"Located_In"][0])
                            

                    #Match multi-exon length genes
                    elif len(workingDataframe.at[index,"Splice_Junctions"]) != 0:
                        exonMatches = exonMatching(workingDataframe.at[index,"Splice_Junctions"], workingDataframe.at[index+1,"Splice_Junctions"], threshold)
                        if sum(exonMatches) == len(workingDataframe.at[index,"Splice_Junctions"]) and len(workingDataframe.at[index,"Splice_Junctions"]) == len(workingDataframe.at[index+1,"Splice_Junctions"]):
                            if abs(workingDataframe.at[index,"start"] - workingDataframe.at[index+1,"start"])<=distanceThreshold and abs(workingDataframe.at[index,"end"] - workingDataframe.at[index+1,"end"])<=distanceThreshold:
                                workingDataframe.at[index,"Matched_Transcript"].append(workingDataframe.at[index+1,"gene_id"])
                                workingDataframe.at[index,"Start_Stop_Distance"].append([workingDataframe.at[index,"start"] - workingDataframe.at[index+1,"start"], workingDataframe.at[index,"end"] - workingDataframe.at[index+1,"end"]])
                                workingDataframe.at[index+1,"found"] = True
                                if workingDataframe.at[index+1,"Located_In"][0] not in workingDataframe.at[index,"Located_In"]:
                                    workingDataframe.at[index,"Located_In"].append(workingDataframe.at[index+1,"Located_In"][0])
                
                    looping = 2
                    while index+looping < len(workingDataframe.index) and workingDataframe.at[index, "end"] >= workingDataframe.at[index+looping,"start"]:

                        #Only match one exon length genes with one exon length genes
                        if len(workingDataframe.at[index,"Splice_Junctions"]) == 0 and len(workingDataframe.at[index+looping,"Splice_Junctions"]) == 0:
                            if abs(workingDataframe.at[index,"start"] - workingDataframe.at[index+looping,"start"])<=distanceThreshold and abs(workingDataframe.at[index,"end"] - workingDataframe.at[index+looping,"end"])<=distanceThreshold:
                                workingDataframe.at[index,"Matched_Transcript"].append(workingDataframe.at[index+looping,"gene_id"])
                                workingDataframe.at[index,"Start_Stop_Distance"].append([workingDataframe.at[index,"start"] - workingDataframe.at[index+looping,"start"], workingDataframe.at[index,"end"] - workingDataframe.at[index+looping,"end"]])
                                workingDataframe.at[index+looping,"found"] = True
                                if workingDataframe.at[index+looping,"Located_In"][0] not in workingDataframe.at[index,"Located_In"]:
                                    workingDataframe.at[index,"Located_In"].append(workingDataframe.at[index+looping,"Located_In"][0])

                        #Match multi-exon length genes
                        elif len(workingDataframe.at[index,"Splice_Junctions"]) != 0:
                            exonMatches = exonMatching(workingDataframe.at[index,"Splice_Junctions"], workingDataframe.at[index+looping,"Splice_Junctions"], threshold)
                            if sum(exonMatches) == len(workingDataframe.at[index,"Splice_Junctions"]) and len(workingDataframe.at[index,"Splice_Junctions"]) == len(workingDataframe.at[index+looping,"Splice_Junctions"]):
                                if abs(workingDataframe.at[index,"start"] - workingDataframe.at[index+looping,"start"])<=distanceThreshold and abs(workingDataframe.at[index,"end"] - workingDataframe.at[index+looping,"end"])<=distanceThreshold:
                                    workingDataframe.at[index,"Matched_Transcript"].append(workingDataframe.at[index+looping,"gene_id"])
                                    workingDataframe.at[index,"Start_Stop_Distance"].append([workingDataframe.at[index,"start"] - workingDataframe.at[index+looping,"start"], workingDataframe.at[index,"end"] - workingDataframe.at[index+looping,"end"]])
                                    workingDataframe.at[index+looping,"found"] = True
                                    if workingDataframe.at[index+looping,"Located_In"][0] not in workingDataframe.at[index,"Located_In"]:
                                        workingDataframe.at[index,"Located_In"].append(workingDataframe.at[index+looping,"Located_In"][0])
                        
                        ##Itterate
                        looping = looping + 1
                        
        workingDataframe = workingDataframe[workingDataframe["found"]==False]
        
        #Sort values based on their start position in our working frame
        workingDataframe = workingDataframe.sort_values(by=['start'])
        
        #create new index to help itteration in our working frame
        workingDataframe = workingDataframe.reset_index()
        workingDataframe = workingDataframe.drop(columns=['index'])
        
        for index, row in workingDataframe.iterrows():
            
            ##If the end is greater then the start of the next index they overlap
            if index+1 < len(workingDataframe.index) and workingDataframe.at[index, "end"] >= workingDataframe.at[index+1,"start"]:

                #Match multi-exon length genes
                if len(workingDataframe.at[index,"Splice_Junctions"]) != 0:
                    exonMatches = exonMatching(workingDataframe.at[index,"Splice_Junctions"], workingDataframe.at[index+1,"Splice_Junctions"], threshold)
                    if sum(exonMatches) >= 1:
                        workingDataframe.at[index, 'Shared_SJ_Transcript'].append(workingDataframe.at[index+1,"gene_id"]) 
                        workingDataframe.at[index, 'Shared_SJ_Weight'].append(1)
                        workingDataframe.at[index+1, 'Shared_SJ_Transcript'].append(workingDataframe.at[index,"gene_id"]) 
                        workingDataframe.at[index+1, 'Shared_SJ_Weight'].append(1)

                looping = 2
                while index+looping < len(workingDataframe.index) and workingDataframe.at[index, "end"] >= workingDataframe.at[index+looping,"start"]:

                    #Match multi-exon length genes
                    if len(workingDataframe.at[index,"Splice_Junctions"]) != 0:
                        exonMatches = exonMatching(workingDataframe.at[index,"Splice_Junctions"], workingDataframe.at[index+looping,"Splice_Junctions"], threshold)
                        if sum(exonMatches) >= 1:
                            workingDataframe.at[index, 'Shared_SJ_Transcript'].append(workingDataframe.at[index+looping,"gene_id"]) 
                            workingDataframe.at[index, 'Shared_SJ_Weight'].append(1)
                            workingDataframe.at[index+looping, 'Shared_SJ_Transcript'].append(workingDataframe.at[index,"gene_id"]) 
                            workingDataframe.at[index+looping, 'Shared_SJ_Weight'].append(1)

                    ##Itterate
                    looping = looping + 1
        
        MasterFrameGeneMatch = MasterFrameGeneMatch.append(workingDataframe)
        
        print(item, end=" ")
        
    #create new index to help itteration in our working frame
    MasterFrameGeneMatch = MasterFrameGeneMatch.reset_index()
    MasterFrameGeneMatch = MasterFrameGeneMatch.drop(columns=['index'])
    
    return(MasterFrameGeneMatch)

##Gtf Transcriptome Reference: compares a transcriptome against a reference

##gtfStart --- the transcriptome dataframe

##gtfReference --- the reference GTf dataframe


def gtfTranscriptomeReference(gtfStart, gtfReference, referenceLabel, threshold, distanceThreshold):
    
    ###Create an array of each GTFs chromosones
    chromosones1 = gtfStart.seqname.unique().tolist()
    chromosones2 = gtfReference.seqname.unique().tolist()

    ##Counter for testing certain properties
    counter = 0

    #Create final dataframe
    MasterFrameGeneMatch = pd.DataFrame()
    MasterFrameGeneMatch2 = pd.DataFrame()

    #Check if chromosones don't match
    if sorted(chromosones1) != sorted(chromosones2):
        print("Non-matching chromosones")

    print("\n")
    print("Comparing to Reference:")

    for item in chromosones1:

            #Subset our dataframe by the chromosone
            chromeFrame1 = gtfStart[gtfStart["seqname"] == item]
            chromeFrame2 = gtfReference[gtfReference["seqname"] == item]

            #create new index to help itteration
            chromeFrame1 = chromeFrame1.reset_index()
            chromeFrame1 = chromeFrame1.drop(columns=['index'])
            
            #print(chromeFrame1.head())

            chromeFrame2 = chromeFrame2.reset_index()
            chromeFrame2 = chromeFrame2.drop(columns=['index'])     
                
            chromeFrame2 = collapsedReturns(chromeFrame2)

            #Look at only transcripts from reference GTF
            chromeFrame2 = chromeFrame2[chromeFrame2["feature"] == "transcript"] 

            chromeFrame1['Matched_Reference_Gene'] = [[] for _ in range(len(chromeFrame1))]
            chromeFrame1['Matched_Reference_Transcript'] = [[] for _ in range(len(chromeFrame1))]
            chromeFrame1['Reference_Transcript_Partial'] = [[] for _ in range(len(chromeFrame1))]
            chromeFrame1['Reference_Gene_Partial'] = [[] for _ in range(len(chromeFrame1))]
            chromeFrame1['Matched_Ref_Start_Stop'] = [[] for _ in range(len(chromeFrame1))]

            #Sort values based on their start position
            chromeFrame1 = chromeFrame1.sort_values(by=['start'])
            chromeFrame2 = chromeFrame2.sort_values(by=['start'])
            
            #create new index to help itteration
            chromeFrame1 = chromeFrame1.reset_index()
            chromeFrame1 = chromeFrame1.drop(columns=['index'])

            chromeFrame2 = chromeFrame2.reset_index()
            chromeFrame2 = chromeFrame2.drop(columns=['index'])
            
            ##Loop through first dataframe and check for correpsonding genes in the second
            for index, row in chromeFrame2.iterrows():
                Exons2 = chromeFrame2.at[index,"Exons"]
                chromeFrame2.at[index,"Exons"] = sorted(Exons2)

            #Add splice junc column for reference gtf
            chromeFrame2 = splicJuncColumn(chromeFrame2)

            itterator = 0
            for index, row in chromeFrame1.iterrows():

                #Move Along Rows in the Reference Until It Is Past Beginning of Current read
                while itterator < len(chromeFrame2.index)-1 and chromeFrame1.at[index, 'start'] > chromeFrame2.at[itterator,"end"]:
                    itterator = itterator + 1
                if (chromeFrame2.at[itterator,"start"] <= chromeFrame1.at[index, 'start'] <= chromeFrame2.at[itterator,"end"]) or (chromeFrame2.at[itterator,"start"] <= chromeFrame1.at[index, 'end'] <= chromeFrame2.at[itterator,"end"]) or (chromeFrame1.at[index,"start"] <= chromeFrame2.at[itterator, 'end'] <= chromeFrame1.at[index,"end"]):

                    #Only match one exon length genes with one exon length genes
                    if len(chromeFrame1.at[index,"Splice_Junctions"]) == 0 and len(chromeFrame2.at[itterator,"Splice_Junctions"]) == 0:
                        if abs(chromeFrame1.at[index, 'start']-chromeFrame2.at[itterator,"start"])<=distanceThreshold and abs(chromeFrame1.at[index, 'end']-chromeFrame2.at[itterator,"end"])<=distanceThreshold:
                            if referenceLabel not in chromeFrame1.at[index,"Located_In"]:
                                chromeFrame1.at[index,"Located_In"].append(referenceLabel)
                            chromeFrame1.at[index, 'Matched_Reference_Gene'].append(chromeFrame2.at[itterator,"gene_id"])
                            chromeFrame1.at[index, 'Matched_Reference_Transcript'].append(chromeFrame2.at[itterator,"transcript_name"])

                    #Match multi-exon length genes if they have over our threshold score for matching
                    elif len(chromeFrame1.at[index,"Splice_Junctions"]) != 0:
                        exonMatches = exonMatching(chromeFrame1.at[index,"Splice_Junctions"], chromeFrame2.at[itterator,"Splice_Junctions"], threshold)
                        if sum(exonMatches) == len(chromeFrame1.at[index,"Splice_Junctions"]) and len(chromeFrame1.at[index,"Splice_Junctions"]) == len(chromeFrame2.at[itterator,"Splice_Junctions"]) and abs(chromeFrame1.at[index, 'start']-chromeFrame2.at[itterator,"start"])<=distanceThreshold and abs(chromeFrame1.at[index, 'end']-chromeFrame2.at[itterator,"end"])<=distanceThreshold:
                            if referenceLabel not in chromeFrame1.at[index,"Located_In"]:
                                chromeFrame1.at[index,"Located_In"].append(referenceLabel)
                            chromeFrame1.at[index, 'Matched_Reference_Gene'].append(chromeFrame2.at[itterator,"gene_id"])
                            chromeFrame1.at[index, 'Matched_Reference_Transcript'].append(chromeFrame2.at[itterator,"transcript_name"])
                            chromeFrame1.at[index,"Matched_Ref_Start_Stop"].append([chromeFrame1.at[index,"start"] - chromeFrame2.at[itterator,"start"], chromeFrame1.at[index,"end"] - chromeFrame2.at[itterator,"end"]])
                        elif sum(exonMatches) >= 1:
                            chromeFrame1.at[index, 'Reference_Transcript_Partial'].append(chromeFrame2.at[itterator,"transcript_name"])
                            chromeFrame1.at[index, 'Reference_Gene_Partial'].append(chromeFrame2.at[itterator,"gene_id"])
                            
                    looping = 1
                    while itterator+looping < len(chromeFrame2.index) and chromeFrame1.at[index, 'end'] >= chromeFrame2.at[itterator+looping,"start"]:
                        
                        if (chromeFrame2.at[itterator+looping,"start"] <= chromeFrame1.at[index, 'start'] <= chromeFrame2.at[itterator+looping,"end"]) or (chromeFrame2.at[itterator+looping,"start"] <= chromeFrame1.at[index, 'end'] <= chromeFrame2.at[itterator+looping,"end"]) or (chromeFrame1.at[index,"start"] <= chromeFrame2.at[itterator+looping, 'end'] <= chromeFrame1.at[index,"end"]):
                            
                            if len(chromeFrame1.at[index,"Splice_Junctions"]) == 0 and len(chromeFrame2.at[itterator+looping,"Splice_Junctions"]) == 0:
                                if abs(chromeFrame1.at[index, 'start']-chromeFrame2.at[itterator+looping,"start"])<=distanceThreshold and abs(chromeFrame1.at[index, 'end']-chromeFrame2.at[itterator+looping,"end"])<=distanceThreshold:
                                    if referenceLabel not in chromeFrame1.at[index,"Located_In"]:
                                        chromeFrame1.at[index,"Located_In"].append(referenceLabel)
                                    chromeFrame1.at[index, 'Matched_Reference_Gene'].append(chromeFrame2.at[itterator+looping,"gene_id"])
                                    chromeFrame1.at[index, 'Matched_Reference_Transcript'].append(chromeFrame2.at[itterator+looping,"transcript_name"])

                            #Match multi-exon length genes if they have over our threshold score for matching
                            elif len(chromeFrame1.at[index,"Splice_Junctions"]) != 0:
                                exonMatches = exonMatching(chromeFrame1.at[index,"Splice_Junctions"], chromeFrame2.at[itterator+looping,"Splice_Junctions"], threshold)
                                if sum(exonMatches) == len(chromeFrame1.at[index,"Splice_Junctions"]) and len(chromeFrame1.at[index,"Splice_Junctions"]) == len(chromeFrame2.at[itterator+looping,"Splice_Junctions"]) and abs(chromeFrame1.at[index, 'start']-chromeFrame2.at[itterator+looping,"start"])<=distanceThreshold and abs(chromeFrame1.at[index, 'end']-chromeFrame2.at[itterator+looping,"end"])<=distanceThreshold:
                                    if referenceLabel not in chromeFrame1.at[index,"Located_In"]:
                                        chromeFrame1.at[index,"Located_In"].append(referenceLabel)
                                    chromeFrame1.at[index, 'Matched_Reference_Gene'].append(chromeFrame2.at[itterator+looping,"gene_id"])
                                    chromeFrame1.at[index, 'Matched_Reference_Transcript'].append(chromeFrame2.at[itterator+looping,"transcript_name"])
                                    chromeFrame1.at[index,"Matched_Ref_Start_Stop"].append([chromeFrame1.at[index,"start"] - chromeFrame2.at[itterator+looping,"start"], chromeFrame1.at[index,"end"] - chromeFrame2.at[itterator+looping,"end"]])
                                elif sum(exonMatches) >= 1:
                                    chromeFrame1.at[index, 'Reference_Transcript_Partial'].append(chromeFrame2.at[itterator+looping,"transcript_name"])
                                    chromeFrame1.at[index, 'Reference_Gene_Partial'].append(chromeFrame2.at[itterator+looping,"gene_id"])
                        looping = looping + 1

            MasterFrameGeneMatch = MasterFrameGeneMatch.append(chromeFrame1) #add our matched data to the master dataframe
            print(item,end = ' ') #Display each chromosone as we finish it
    
    MasterFrameGeneMatch = MasterFrameGeneMatch.reset_index()
    MasterFrameGeneMatch = MasterFrameGeneMatch.drop(columns=['index'])
    
    return(MasterFrameGeneMatch)


##gtf Symbiosis: Master function which both does the Symbiosis and also displays the results



def gftSymbiosis(analyzedGTF,referenceGTF,matchThreshold, matchThreshold2, preMerged1 = False, preMerged2 = False, reshape = False, inTissue = False):
    
    #Find Number Of Exons, Transcripts, and Genes in both GTFs
    if preMerged1 == False:
        exonSizeAn = len(analyzedGTF[analyzedGTF["feature"] == "exon"].index)
        tranSizeAn = len(analyzedGTF[analyzedGTF["feature"] == "transcript"].index)
        geneSizeAn = len(analyzedGTF[analyzedGTF["feature"] == "gene"].index)
    elif preMerged1 == True:
        exonSizeAn = sumofExons(analyzedGTF)
        tranSizeAn = len(analyzedGTF[analyzedGTF["feature"] == "transcript"].index)
        geneSizeAn = len(analyzedGTF[analyzedGTF["feature"] == "gene"].index)

    if preMerged2 == False:
        exonSizeRef = len(referenceGTF[referenceGTF["feature"] == "exon"].index)
        tranSizeRef = len(referenceGTF[referenceGTF["feature"] == "transcript"].index)
        geneSizeRef = len(referenceGTF[referenceGTF["feature"] == "gene"].index)
    elif preMerged2 == True:
        exonSizeRef = sumofExons(referenceGTF)
        tranSizeRef = len(referenceGTF[referenceGTF["feature"] == "transcript"].index)
        geneSizeRef = len(referenceGTF[referenceGTF["feature"] == "gene"].index)
    
    if reshape == True:
        analyzedGTF = reshapeGTF(analyzedGTF)
        referenceGTF = reshapeGTF(referenceGTF)
    
    if inTissue == False:
        #Print Out Basic Data of Both GTFs
        t = PrettyTable(['Name', 'Transcripts', 'Exons', 'Genes'])
        t.add_row(['Analyzed', tranSizeAn, exonSizeAn, geneSizeAn])
        t.add_row(['Reference', tranSizeRef, exonSizeRef, geneSizeRef])
        print(t)
    else:
        #Print Out Basic Data of Both GTFs
        t = PrettyTable(['Name', 'Transcripts', 'Exons', 'Genes'])
        t.add_row(['Analyzed 1', tranSizeAn, exonSizeAn, geneSizeAn])
        t.add_row(['Analyzed 2', tranSizeRef, exonSizeRef, geneSizeRef])
        print(t)
        
    #Run the master comparison Fucntion
    exonMatchedGTF = symbiosisMasterFunction(analyzedGTF, referenceGTF, matchThreshold, inTissue, preMerged1, preMerged2)   
    
    if inTissue == True:
        
        #find the best matches
        bestMatchFrame, oneExonSizeAn = bestMatch(exonMatchedGTF, matchThreshold2)
        
        #Find the Number of Perfect Matches
        perfectMatchNum = sumPerfectMatches(exonMatchedGTF)

        #Create Dataframe of only matched transcripts
        fullFrame = bestMatchFrame.dropna(subset=['Matched_Transcript'])
        fullFrame = fullFrame[fullFrame.astype(str)['Matched_Transcript'] != '[]']

        #Find All Genuine No-Matches
        genuineNoMatch = exonMatchedGTF[exonMatchedGTF.astype(str)['Matched_Transcript'] == '[]'].shape[0]
        
        #Find number of matched transcripts from first gtf
        matched1 = exonMatchedGTF[exonMatchedGTF.astype(str)['Matched_Transcript'] != '[]'].shape[0]
        
        #Find Missing Splice Junction Locations
        count1, count2, count3, count4 = splicJuncMissing(fullFrame)
        
        #Find all matched transcripts from second gtf
        matched2List = []
        for index, row in exonMatchedGTF.iterrows():
            if not exonMatchedGTF.at[index,"Matched_Transcript"]:
                pass
            else: 
                matched2List.extend(exonMatchedGTF.at[index,"Matched_Transcript"])
        #mached2List = list(set(matched2List))
        matched2ListUnique = list(OrderedDict.fromkeys(matched2List)) 
        matched2 = len(matched2ListUnique)
        
        print(matched2)
        print(matched1)
        
        # Make our dataset ready for graphing:
        barWidth = 0.9
        height = [genuineNoMatch, matched1+matched2]
        bars = ['No Matching Transcript', 'Matching (All)']
        x = np.arange(len(bars))
        rPlacement = [1,2,]
        height2 = [count1, count2, count3, count4]
        bars2 = ['First', 'First 1/2', 'Second 1/2', 'Last']
        rPlacement2 = [1,2,3,4]
        x2 = np.arange(len(bars2))

        ##Print histograms of Start and Stop changes
        fig, axs = plt.subplots(1, 2, figsize=(10,10))
        axs[0].bar(rPlacement, height, width = barWidth)
        axs[0].set_title('Transcript Types')
        axs[0].set_ylabel('Count')
        axs[0].set_xticks(x+1)
        axs[0].set_xticklabels(bars, rotation=40, ha='right')

        axs[1].bar(rPlacement2, height2, width = barWidth)
        axs[1].set_title('Missing Splice Junction Locations')
        axs[1].set_ylabel('Count')
        axs[1].set_xticks(x2+1)
        axs[1].set_xticklabels(bars2, rotation=40, ha='right')
        
         # Create names on the x-axis
        label = [str(genuineNoMatch), str(matched1+matched2), str(perfectMatchNum)]

        for i in range(len(rPlacement)):
            axs[0].text(x = rPlacement[i]-.25 , y = height[i]+(height[i]/50), s = label[i], size = 18)

        # Create names on the x-axis
        label2 = [str(count1), str(count2), str(count3), str(count4)]

        for i in range(len(rPlacement2)):
            axs[1].text(x = rPlacement2[i]-.25 , y = height2[i]+(height2[i]/50), s = label2[i], size = 18)


        # Show graphic
        plt.show()

        
    
        ##Create multi-perfect match frame
        perfectMatchFrame = returnMultiPerfectMatches(exonMatchedGTF)
        starts = []
        stops = []
        for index, row in perfectMatchFrame.iterrows():
            i = 0
            while i < len(perfectMatchFrame.at[index, "Perfect_Starts"]):
                if perfectMatchFrame.at[index, "strand"] == "-":
                    starts.append(perfectMatchFrame.at[index, "Perfect_Stops"][i])
                    stops.append(perfectMatchFrame.at[index, "Perfect_Starts"][i])
                elif perfectMatchFrame.at[index, "strand"] == "+":
                    starts.append(perfectMatchFrame.at[index, "Perfect_Starts"][i])
                    stops.append(perfectMatchFrame.at[index, "Perfect_Stops"][i])
                i = i + 1

        #Turn them into numpy arrays
        startArray = np.array(starts, dtype=object)
        stopArray = np.array(stops, dtype=object)
        
        #Flatten them
        startSingles = [item for item in startArray if isinstance(item, int)]
        stopSingles = [item for item in stopArray if isinstance(item, int)]
        
        startSingles = np.array(startSingles)
        stopSingles = np.array(stopSingles)
        
        startGroups = [i for i in startArray if i not in startSingles]
        stopGroups = [i for i in stopArray if i not in stopSingles]
        
        startGroups = np.array(startGroups,dtype=object)
        stopGroups = np.array(stopGroups,dtype=object)
        
        startGroups = np.hstack(startGroups)
        stopGroups = np.hstack(stopGroups)
        
        startArray = np.concatenate((startSingles, startGroups))
        stopArray = np.concatenate((stopSingles, stopGroups))
        
        #Create Dataframes from start and stop values
        startFrame = pd.DataFrame(startArray, columns= ['start'])
        stopFrame = pd.DataFrame(stopArray, columns= ['stop'])
        
        
        # Create Quant Frames for each
        quantFrameStart = startFrame.quantile([.1, .25, .5, .75, .9])
        quantFrameStop = stopFrame.quantile([.1, .25, .5, .75, .9])
    
        print("\n")
    
        #Print Basic Table
        t = PrettyTable(['Quantile', 'Start', 'Stop'])
        t.add_row(['.1', quantFrameStart.iloc[0,0], quantFrameStop.iloc[0,0]])
        t.add_row(['.25', quantFrameStart.iloc[1,0], quantFrameStop.iloc[1,0]])
        t.add_row(['.5', quantFrameStart.iloc[2,0], quantFrameStop.iloc[2,0]])
        t.add_row(['.75', quantFrameStart.iloc[3,0], quantFrameStop.iloc[3,0]])
        t.add_row(['.9', quantFrameStart.iloc[4,0], quantFrameStop.iloc[4,0]])                            
        print(t)
        
        matchedMPGTF = returnMultiPerfectMatches(exonMatchedGTF)
        
        t = PrettyTable(['# Multi-Match Starts/Stops', '# of Multi-Matches'])
        t.add_row([len(startArray),matchedMPGTF.shape[0]])                           
        print(t)
                   
         ##Print histograms of Start and Stop changes
        fig, axs = plt.subplots(1, 2, constrained_layout=True)

        axs[0].hist(startFrame['start'],range=(-2000,2000))
        axs[0].set_title('Start Distance')
        axs[0].set_xlabel('Distance (Base Pairs)')
        axs[0].set_ylabel('Count')
        fig.suptitle('Start and Stop Distance of Multi-Perfect Matches', fontsize=16)

        axs[1].hist(stopFrame['stop'],range=(-2000,2000))
        axs[1].set_title('Stop Distance')
        axs[1].set_xlabel('Distance (Base Pairs)')
        axs[1].set_ylabel('Count')
        #Show graphic
        plt.show()

        ##Print histograms of Start and Stop changes
        fig, axs = plt.subplots(1, 2, constrained_layout=True)

        axs[0].hist(startFrame['start'],range=(-200,200))
        axs[0].set_title('Start Distance')
        axs[0].set_xlabel('Distance (Base Pairs)')
        axs[0].set_ylabel('Count')
        fig.suptitle('Start and Stop Distance of Multi-Perfect Matches', fontsize=16)

        axs[1].hist(stopFrame['stop'],range=(-200,200))
        axs[1].set_title('Stop Distance')
        axs[1].set_xlabel('Distance (Base Pairs)')
        axs[1].set_ylabel('Count')
        #Show graphic
        plt.show()

        ##Print histograms of Start and Stop changes
        fig, axs = plt.subplots(1, 2, constrained_layout=True)

        axs[0].hist(startFrame['start'],range=(-10,10))
        axs[0].set_title('Start Distance')
        axs[0].set_xlabel('Distance (Base Pairs)')
        axs[0].set_ylabel('Count')
        fig.suptitle('Start and Stop Distance of Multi-Perfect Matches', fontsize=16)

        axs[1].hist(stopFrame['stop'],range=(-10,10))
        axs[1].set_title('Stop Distance')
        axs[1].set_xlabel('Distance (Base Pairs)')
        axs[1].set_ylabel('Count')
        #Show graphic
        plt.show()
        
        ##Print histograms of Start and Stop changes
        fig, axs = plt.subplots(1, 2, constrained_layout=True)

        axs[0].hist(startFrame['start'],range=(-5,5))
        axs[0].set_title('Start Distance')
        axs[0].set_xlabel('Distance (Base Pairs)')
        axs[0].set_ylabel('Count')
        fig.suptitle('Start and Stop Distance of Multi-Perfect Matches', fontsize=16)

        axs[1].hist(stopFrame['stop'],range=(-10,10))
        axs[1].set_title('Stop Distance')
        axs[1].set_xlabel('Distance (Base Pairs)')
        axs[1].set_ylabel('Count')
        #Show graphic
        plt.show()
    
    if inTissue == False:
        
        #Find the Number of Perfect Matches
        perfectMatchNum = sumPerfectMatches(exonMatchedGTF)

        #Find Number of non-matched transcripts
        genuineNoMatch = exonMatchedGTF[exonMatchedGTF.astype(str)['Matched_Transcript'] == '[]'].shape[0]

        #find the best matches
        bestMatchFrame, oneExonSizeAn = bestMatch(exonMatchedGTF, matchThreshold2)

        #Total number of Transcripts
        totalSize = bestMatchFrame.shape[0]

        #Create Dataframe of only matched transcripts
        fullFrame = bestMatchFrame.dropna(subset=['Matched_Transcript'])
        fullFrame = fullFrame[fullFrame.astype(str)['Matched_Transcript'] != '[]']


        #Find Missing Splice Junction Locations
        count1, count2, count3, count4 = splicJuncMissing(fullFrame)

        #find size of matched transcripts
        matchesNum = fullFrame.shape[0]

        #Find amount of unmatched to exons transcripts
        noMatchesNum = totalSize - matchesNum - genuineNoMatch
        
        #Find total of matched unique exons
        uniqueSize = uniqueMatches(bestMatchFrame)

        # Make our dataset ready for graphing:
        barWidth = 0.9
        height = [matchesNum, genuineNoMatch, noMatchesNum]
        bars = ['Non-Zero Matched Transcripts', 'No Matching Transcript', 'Match Transcript/No Splice Junctions Matches']
        x = np.arange(len(bars))
        rPlacement = [1,2,3]
        height2 = [count1, count2, count3, count4]
        bars2 = ['First', 'First 1/2', 'Second 1/2', 'Last']
        rPlacement2 = [1,2,3,4]
        x2 = np.arange(len(bars2))

        ##Print histograms of Start and Stop changes
        fig, axs = plt.subplots(1, 2, figsize=(10,10))
        axs[0].bar(rPlacement, height, width = barWidth)
        axs[0].set_title('Transcript Types')
        axs[0].set_ylabel('Count')
        axs[0].set_xticks(x+1)
        axs[0].set_xticklabels(bars, rotation=40, ha='right')

        axs[1].bar(rPlacement2, height2, width = barWidth)
        axs[1].set_title('Missing Splice Junction Locations')
        axs[1].set_ylabel('Count')
        axs[1].set_xticks(x2+1)
        axs[1].set_xticklabels(bars2, rotation=40, ha='right')


        # Create names on the x-axis
        label = [str(matchesNum), str(genuineNoMatch), str(noMatchesNum)]

        for i in range(len(rPlacement)):
            axs[0].text(x = rPlacement[i]-.25 , y = height[i]+(height[i]/50), s = label[i], size = 18)

        # Create names on the x-axis
        label2 = [str(count1), str(count2), str(count3), str(count4)]

        for i in range(len(rPlacement2)):
            axs[1].text(x = rPlacement2[i]-.25 , y = height2[i]+(height2[i]/50), s = label2[i], size = 18)


        # Show graphic
        plt.show()
        
        ##Create perfect match frame
        perfectMatchFrame = fullFrame[fullFrame['Perfect_Matches'] == 1]
        perfectMatchFrame[['start','stop']] = pd.DataFrame(perfectMatchFrame.Start_Stop_Distance.tolist(), index= perfectMatchFrame.index)
        startStopFrame = perfectMatchFrame[['start','stop']]


        #Find Qunatiles of the Perfect Match Frame
        quantFrame = startStopFrame.quantile([.1, .25, .5, .75, .9], axis = 0)

        #Print Basic Table
        t = PrettyTable(['Quantile', 'Start', 'Stop'])
        t.add_row(['.1', quantFrame.iloc[0,0], quantFrame.iloc[0,1]])
        t.add_row(['.25', quantFrame.iloc[1,0], quantFrame.iloc[1,1]])
        t.add_row(['.5', quantFrame.iloc[2,0], quantFrame.iloc[2,1]])
        t.add_row(['.75', quantFrame.iloc[3,0], quantFrame.iloc[3,1]])
        t.add_row(['.9', quantFrame.iloc[4,0], quantFrame.iloc[4,1]])                            
        print(t)

        ##Print histograms of Start and Stop changes
        fig, axs = plt.subplots(1, 2, constrained_layout=True)

        axs[0].hist(perfectMatchFrame['start'],range=(-2000,2000))
        axs[0].set_title('Start Distance')
        axs[0].set_xlabel('Distance (Base Pairs)')
        axs[0].set_ylabel('Count')
        fig.suptitle('Start and Stop Distance of Perfect Matches', fontsize=16)

        axs[1].hist(perfectMatchFrame['stop'],range=(-2000,2000))
        axs[1].set_title('Stop Distance')
        axs[1].set_xlabel('Distance (Base Pairs)')
        axs[1].set_ylabel('Count')
        #Show graphic
        plt.show()

        ##Print histograms of Start and Stop changes
        fig, axs = plt.subplots(1, 2, constrained_layout=True)

        axs[0].hist(perfectMatchFrame['start'],range=(-200,200))
        axs[0].set_title('Start Distance')
        axs[0].set_xlabel('Distance (Base Pairs)')
        axs[0].set_ylabel('Count')
        fig.suptitle('Start and Stop Distance of Perfect Matches', fontsize=16)

        axs[1].hist(perfectMatchFrame['stop'],range=(-200,200))
        axs[1].set_title('Stop Distance')
        axs[1].set_xlabel('Distance (Base Pairs)')
        axs[1].set_ylabel('Count')
        #Show graphic
        plt.show()

        ##Print histograms of Start and Stop changes
        fig, axs = plt.subplots(1, 2, constrained_layout=True)

        axs[0].hist(perfectMatchFrame['start'],range=(-10,10))
        axs[0].set_title('Start Distance')
        axs[0].set_xlabel('Distance (Base Pairs)')
        axs[0].set_ylabel('Count')
        fig.suptitle('Start and Stop Distance of Perfect Matches', fontsize=16)

        axs[1].hist(perfectMatchFrame['stop'],range=(-10,10))
        axs[1].set_title('Stop Distance')
        axs[1].set_xlabel('Distance (Base Pairs)')
        axs[1].set_ylabel('Count')
        #Show graphic
        plt.show()
    
    #Print Basic Table
    t = PrettyTable(['Name', 'Amount'])
    t.add_row(['Perfect Matches', perfectMatchNum])
    t.add_row(['One Exon Length Matches', oneExonSizeAn])
    print(t)
    
    if inTissue == False:
        data = {'Matched Transcripts': round((uniqueSize/tranSizeRef)*100), 'Unmatched Transcripts': 100-round((uniqueSize/tranSizeRef)*100)}
        fig = plt.figure(
            FigureClass=Waffle, 
            rows=5, 
            values=data, 
            colors=("#983D3D", "#232066"),
            title={'label': 'Amount of Transcripts Matched From Reference Transcriptome', 'loc': 'left'},
            labels=["{0} ({1}%)".format(k, v) for k, v in data.items()],
            legend={'loc': 'lower left', 'bbox_to_anchor': (0, -0.4), 'ncol': len(data), 'framealpha': 0}
        )
        fig.gca().set_facecolor('#EEEEEE')
        fig.set_facecolor('#EEEEEE')
        plt.show()
    
    

        #Find Total Amount Of Exons Matched
        totalSJMatch, SJMasterSize = sumTotalSJMatched(exonMatchedGTF)

        #Plot Splice Junction Span in Waffle Plot
        data = {'Matched Splice Junctions': round((totalSJMatch/SJMasterSize)*100), 'Unmatched Splice Junctions': 100-round((totalSJMatch/SJMasterSize)*100)}
        fig = plt.figure(
            FigureClass=Waffle, 
            rows=5, 
            values=data, 
            colors=("#983D3D", "#232066"),
            title={'label': 'Amount Of Splice Junctions Matched From Analyzed GTF', 'loc': 'left'},
            labels=["{0} ({1}%)".format(k, v) for k, v in data.items()],
            legend={'loc': 'lower left', 'bbox_to_anchor': (0, -0.4), 'ncol': len(data), 'framealpha': 0}
        )
        fig.gca().set_facecolor('#EEEEEE')
        fig.set_facecolor('#EEEEEE')
        plt.show()
    
    
    return(exonMatchedGTF)
    
##Gtf Transcriptome 2: Master fucntion which both assembles the Transcriptome and displays the results

##analyzedGTFs --- list of GTFs you wish combine into a transcriptome

##analyzedLabels --- list of names you wish to associate with each gtf

##exonMatchThreshold --- distance which you wish to count as a matching splice junction

##distanceThreshold --- distance which you wish to count as a matching transcript

##gtfReference --- a gtf dataframe of a reference like Ensembl which your final transcirptome can be compared against

##referenceLabel --- the name you wish to associate with your reference

    
def gftTranscriptome2(analyzedGTFs,analyzedLabels, exonMatchThreshold, distanceThreshold, gtfReference = "", referenceLabel =""):
    
    print("Combining:")
    
    #Get the total size of all our GTFs combined
    i = 0
    completeSize = 0
    while i < len(analyzedGTFs):
        print(analyzedLabels[i], end='')
        sizingGTF = analyzedGTFs[i].copy()
        completeSize = completeSize + sizingGTF[sizingGTF["feature"]=="transcript"].shape[0]
        if i < len(analyzedGTFs)-1:
            print(", ", end='')
        i = i + 1
    print("\n")
    
    #Combine gtfs
    exonMatchedGTF = transcriptomeMasterFunction2(analyzedGTFs, analyzedLabels, exonMatchThreshold, distanceThreshold)
    
    #Print a table of the starting and final sizes
    t = PrettyTable(['Name', 'Amount'])
    t.add_row(['Starting Size', completeSize])
    t.add_row(['Final Size', exonMatchedGTF.shape[0]])
    print("\n")
    print(t)
    
    #Create the network to mark our genes
    print("Creating Connections:")
    my_network_data = []
    for index, row in exonMatchedGTF.iterrows():
        sizeDF = len(exonMatchedGTF.index)
        numbersInterval = [round(sizeDF/4), round(2*sizeDF/4), round(3*sizeDF/4), sizeDF-1]
        terms = ["25%", "50%", "75%", "100%"]
        for item in exonMatchedGTF.at[index, "Shared_SJ_Transcript"]:
            if item is not None:
                connection = [item, exonMatchedGTF.at[index, "gene_id"]]
                connection.sort()
                my_network_data.append(connection)
        if index in numbersInterval:
            locationIndex = numbersInterval.index(index)
            print(terms[locationIndex], end=" ")
    print("\n")
    g = nx.Graph()
    g.add_edges_from(my_network_data)
    
    #Print Out the details of the network we made
    print(nx.info(g))
    
    
    #Get ready to add enw gene and transcript names
    exonMatchedGTF2 = exonMatchedGTF.copy()
    exonMatchedGTF2["new_gene_id"] = None
    exonMatchedGTF2["new_transcript_id"] = None
    exonMatchedGTF2 = exonMatchedGTF2.set_index('gene_id')
    fullsize = len(exonMatchedGTF2.index)
    fulllength = len(str(fullsize))
    
    #Get the size of our network to print out percentage intervals
    sizeDF = nx.number_connected_components(g)
    print("Number of Clusters: " + str(sizeDF))
    numbersInterval = [round(sizeDF/4), round(2*sizeDF/4), round(3*sizeDF/4), sizeDF-1]
    terms = ["25%", "50%", "75%", "100%"]
    
    #Add gene names and transcript names
    k = 0
    print("Creating Genes:")
    for c in nx.connected_components(g):
        clist = list(c)
        kstr = str(k)
        zero_filled_k = kstr.zfill(fulllength)
        p = 1
        for item in clist:
            exonMatchedGTF2.at[item, "new_gene_id"] = "COISORAT." + zero_filled_k
            exonMatchedGTF2.at[item, "new_transcript_id"] = "COISORAT." + zero_filled_k + "." + str(p)
            p = p + 1
        if k in numbersInterval:
            locationIndex = numbersInterval.index(k)
            print(terms[locationIndex], end = " ")
        k = k + 1
        
    for index, row in exonMatchedGTF2.iterrows():
        if exonMatchedGTF2.at[index, "new_gene_id"] == None:
            kstr = str(k)
            zero_filled_k = kstr.zfill(fulllength)
            exonMatchedGTF2.at[index, "new_gene_id"] = "COISORAT." + zero_filled_k
            exonMatchedGTF2.at[index, "new_transcript_id"] = "COISORAT." + zero_filled_k + "." + str(1)
            k=k+1
        #Add start and stops if it makes the size larger
        for item in exonMatchedGTF2.at[index, "Start_Stop_Distance"]:
            if exonMatchedGTF2.at[index, "start"] + item[0] < exonMatchedGTF2.at[index, "start"]:
                exonMatchedGTF2.at[index, "start"] = exonMatchedGTF2.at[index, "start"] + item[0]
            if exonMatchedGTF2.at[index, "end"] + item[1] < exonMatchedGTF2.at[index, "end"]:
                exonMatchedGTF2.at[index, "end"] = exonMatchedGTF2.at[index, "end"] - item[1]
        exonMatchedGTF2.at[index, "Matched_Transcript"].append(exonMatchedGTF2.at[index, "transcript_id"])
​
    #Clean Up Columns and column names
    exonMatchedGTF2 = exonMatchedGTF2.reset_index()
    exonMatchedGTF2 = exonMatchedGTF2.drop(columns=['found'])
    exonMatchedGTF2 = exonMatchedGTF2.drop(columns=['transcript_id'])
    exonMatchedGTF2 = exonMatchedGTF2.drop(columns=['gene_id'])
    exonMatchedGTF2 = exonMatchedGTF2.drop(columns=['Previously_Matched'])
    exonMatchedGTF2 = exonMatchedGTF2.rename(columns={"new_gene_id": "gene_id", "new_transcript_id": "transcirpt_id"})
    
    if isinstance(gtfReference, pd.DataFrame):
        
        exonMatchedGTF3 = gtfTranscriptomeReference(exonMatchedGTF2, gtfReference, referenceLabel, exonMatchThreshold, distanceThreshold)
        
        fullRefSize = len(exonMatchedGTF3.index)
        
        nonReference = exonMatchedGTF3[exonMatchedGTF3.astype(str)['Matched_Reference_Gene'] == '[]']
        nonReferenceSize = len(nonReference.index)
        
        nonReferencePart = exonMatchedGTF3[exonMatchedGTF3.astype(str)['Reference_Transcript_Partial'] == '[]']
        nonReferencePartSize = len(nonReferencePart.index)
        
        doubleMatch = exonMatchedGTF3[exonMatchedGTF3.astype(str)['Reference_Transcript_Partial'] != '[]']
        doubleMatch = doubleMatch[doubleMatch.astype(str)['Matched_Reference_Gene'] != '[]']
        doubleMatchSize = len(doubleMatch.index)
        
        #Print a table of the starting and final sizes
        t = PrettyTable(['Name', 'Amount'])
        t.add_row(['Matched To Ensemble (Full)', fullRefSize-nonReferenceSize])
        t.add_row(['Non-Matched To Ensemble (Full)', nonReferenceSize])
        t.add_row(['Matched To Ensemble (Partial)', fullRefSize-nonReferencePartSize])
        t.add_row(['Matched To Ensemble (Full & Partial)', doubleMatchSize])
        print("\n")
        print(t)
        
        
        #Create each set for our venn diagram dictionary
        k = 0
        vennDict = {}
        while k < len(analyzedLabels):
            vennDict[analyzedLabels[k]] = {"default"}
            k = k + 1
        
        vennDict[referenceLabel] = {"default"}
        
        for index, row in exonMatchedGTF3.iterrows():
            for item in exonMatchedGTF3.at[index, "Located_In"]:
                if (item in analyzedLabels) or (item == referenceLabel):
                    vennDict[item].add(exonMatchedGTF.at[index, "gene_id"])
        
        r = dict(vennDict)
        del r[referenceLabel]
        
        while k < len(analyzedLabels):
            vennDict[analyzedLabels[k]].remove("default")
            k = k + 1
        vennDict[referenceLabel].remove("default")
        if len(analyzedLabels) == 6:
            pseudovenn(r)
        elif len(analyzedLabels) <= 6:
            venn(r)
        else:
            print("More than 6 GTFs combined, result will not be displayed as a venn diagram.")
            print("\n")
        
        #Create our PYU plot
        pyuSeries = from_contents(vennDict)
        plot(pyuSeries)
        plt.show()
        
        return(exonMatchedGTF3)
    
    
    #Print Venn Diagram for under 6 GTFs and a pyuplot for everything        
    for index, row in exonMatchedGTF.iterrows():
        for item in exonMatchedGTF.at[index, "Located_In"]:
            if item in analyzedLabels:
                vennDict[item].add(exonMatchedGTF.at[index, "gene_id"])
​
    while k < len(analyzedLabels):
        vennDict[analyzedLabels[k]].remove("default")
        k = k + 1
    if len(analyzedLabels) == 6:
        pseudovenn(vennDict)
    elif len(analyzedLabels) <= 6:
        venn(vennDict)
    else:
        print("More than 6 GTFs combined, result will not be displayed as a venn diagram.")
        print("\n")
    
    
    #Create each set for our venn diagram dictionary
    k = 0
    vennDict = {}
    while k < len(analyzedLabels):
        vennDict[analyzedLabels[k]] = {"default"}
        k = k + 1
    
    #Create our PYU plot
    pyuSeries = from_contents(vennDict)
    plot(pyuSeries)
    plt.show()
    
    return(exonMatchedGTF2)
    