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

## The finale result will look like:
    
    ## feature   :   exons
    ##transcript :  [[124,456], [546, 789]]
    ## exon 1
    ## exon 2

def collapsedReturns(colFrame2):
    
    colFrame2['Exons'] = np.empty((len(colFrame2), 0)).tolist()
    counter = 0
    while counter < len(colFrame2.index):
        
        ##Collect the exons for each transcript
        exons = []
        if colFrame2.at[counter,"feature"] == "transcript" :
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


##Returns the number of matched Splice Junctions (matched can occur across matched transcripts i.e: there are two splice junctions, in one matched transcript only the first is matched, in another matched transcript the second splice junction is matched, a matched SJ value of 2 will be returned since both were matched),
##also returns total number of splice junctions in the given GTF

##sumFrame --- a matched gtf dataframe

def sumTotalSJMatched(sumFrame):
    
    matchSum = 0
    masterSize = 0
    for index, row in sumFrame.iterrows():
        if sumFrame.at[index,"Matched_Reference_Gene"]:
            matchSum = matchSum + len(sumFrame.at[index,"Splice_Junctions"])
            masterSize = masterSize + len(sumFrame.at[index,"Splice_Junctions"])
        elif sumFrame.at[index,"Reference_Gene_Partial"]:
            listG = sumFrame.at[index,"Shared_Ref_Weight"][0][0]
            masterSize = masterSize + len(sumFrame.at[index,"Splice_Junctions"])
            i=0
            while i < len(sumFrame.at[index,"Reference_Transcript_Partial"]):
                j=0
                while j < len(sumFrame.at[index,"Splice_Junctions"]):
                    if sumFrame.at[index,"Shared_Ref_Weight"][i][0][j] == 1:
                        if listG[j] == 0:
                            listG[j] = 1
                    j = j+1
                i=i+1
            
            matchSum = matchSum + sum(list)
        else:
            masterSize = masterSize + len(sumFrame.at[index,"Splice_Junctions"])
            
            
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


def splicJuncMissing(sumFrame, ref = False):
    
    indices=[]
    
    #create new index to help itteration
    sumFrame = sumFrame.reset_index()
    sumFrame = sumFrame.drop(columns=['index'])
    
    if ref == False:
        sumFrame['Splice_Misses'] = np.zeros((len(sumFrame), 4),dtype=int).tolist()
        for index, row in sumFrame.iterrows():
            i = 0
            while i < len(sumFrame.at[index,"Shared_SJ_Weight"]):
                if not sumFrame.at[index,"Shared_SJ_Weight"][i][0]:
                    pass
                elif str(sumFrame.at[index,"Shared_SJ_Weight"][i][0]) == "[]":
                    pass
                elif len(sumFrame.at[index,"Shared_SJ_Weight"][i][0]) == 1 or len(sumFrame.at[index,"Shared_SJ_Weight"][i][0]) == 0:
                    pass
                elif len(sumFrame.at[index,"Shared_SJ_Weight"][i][0]) == 2:
                    if sumFrame.at[index,"Shared_SJ_Weight"][i][0][0] == 0:
                        sumFrame.at[index,"Splice_Misses"][0] = 1
                    if sumFrame.at[index,"Shared_SJ_Weight"][i][0][1] == 0:
                        sumFrame.at[index,"Splice_Misses"][3] = 1
                    indices.append(index)
                elif len(sumFrame.at[index,"Shared_SJ_Weight"][i][0]) % 2 == 0:
                    if sumFrame.at[index,"Shared_SJ_Weight"][i][0][0] == 0:
                        sumFrame.at[index,"Splice_Misses"][0] = 1
                    if sumFrame.at[index,"Shared_SJ_Weight"][i][0][len(sumFrame.at[index,"Shared_SJ_Weight"][i][0])-1] == 0:
                        sumFrame.at[index,"Splice_Misses"][3] = 1
                    j = 1
                    while j < len(sumFrame.at[index,"Shared_SJ_Weight"][i][0])/2:
                        if sumFrame.at[index,"Shared_SJ_Weight"][i][0][j] == 0:
                            sumFrame.at[index,"Splice_Misses"][1] = 1
                        if sumFrame.at[index,"Shared_SJ_Weight"][i][0][len(sumFrame.at[index,"Shared_SJ_Weight"][i][0])-j] == 0:
                            sumFrame.at[index,"Splice_Misses"][2] = 1
                        j = j + 1
                    indices.append(index)
                elif len(sumFrame.at[index,"Shared_SJ_Weight"][i][0]) % 2 != 0:
                    if sumFrame.at[index,"Shared_SJ_Weight"][i][0][0] == 0:
                        sumFrame.at[index,"Splice_Misses"][0] = 1
                    if sumFrame.at[index,"Shared_SJ_Weight"][i][0][len(sumFrame.at[index,"Shared_SJ_Weight"][i][0])-1] == 0:
                        sumFrame.at[index,"Splice_Misses"][3] = 1
                    if sumFrame.at[index,"Shared_SJ_Weight"][i][0][math.floor(len(sumFrame.at[index,"Shared_SJ_Weight"][i][0])/2)+1] == 0:
                        sumFrame.at[index,"Splice_Misses"][1] = 1
                    j = 1
                    while j < math.floor(len(sumFrame.at[index,"Shared_SJ_Weight"][i][0])/2):
                        if sumFrame.at[index,"Shared_SJ_Weight"][i][0][j] == 0:
                            sumFrame.at[index,"Splice_Misses"][1] = 1
                        if sumFrame.at[index,"Shared_SJ_Weight"][i][0][len(sumFrame.at[index,"Shared_SJ_Weight"][i][0])-j] == 0:
                            sumFrame.at[index,"Splice_Misses"][2] = 1
                        j = j + 1
                    indices.append(index)
                i = i + 1

        sumFrame2 = sumFrame.iloc[indices,]
        count1 = 0
        count2 = 0
        count3 = 0
        count4 = 0

        sumFrame2 = sumFrame2.reset_index()
        sumFrame2 = sumFrame2.drop(columns=['index'])

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
    
    elif ref == True:
        
        sumFrame['Splice_Misses'] = np.zeros((len(sumFrame), 4),dtype=int).tolist()
        for index, row in sumFrame.iterrows():
            i = 0
            while i < len(sumFrame.at[index,"Shared_Ref_Weight"]):
                if not sumFrame.at[index,"Shared_Ref_Weight"][i][0]:
                    pass
                elif str(sumFrame.at[index,"Shared_Ref_Weight"][i][0]) == "[]":
                    pass
                elif len(sumFrame.at[index,"Shared_Ref_Weight"][i][0]) == 1 or len(sumFrame.at[index,"Shared_Ref_Weight"][i][0]) == 0:
                    pass
                elif len(sumFrame.at[index,"Shared_Ref_Weight"][i][0]) == 2:
                    if sumFrame.at[index,"Shared_Ref_Weight"][i][0][0] == 0:
                        sumFrame.at[index,"Splice_Misses"][0] = 1
                    if sumFrame.at[index,"Shared_Ref_Weight"][i][0][1] == 0:
                        sumFrame.at[index,"Splice_Misses"][3] = 1
                    indices.append(index)
                elif len(sumFrame.at[index,"Shared_Ref_Weight"][i][0]) % 2 == 0:
                    if sumFrame.at[index,"Shared_Ref_Weight"][i][0][0] == 0:
                        sumFrame.at[index,"Splice_Misses"][0] = 1
                    if sumFrame.at[index,"Shared_Ref_Weight"][i][0][len(sumFrame.at[index,"Shared_Ref_Weight"][i][0])-1] == 0:
                        sumFrame.at[index,"Splice_Misses"][3] = 1
                    j = 1
                    while j < len(sumFrame.at[index,"Shared_Ref_Weight"][i][0])/2:
                        if sumFrame.at[index,"Shared_Ref_Weight"][i][0][j] == 0:
                            sumFrame.at[index,"Splice_Misses"][1] = 1
                        if sumFrame.at[index,"Shared_Ref_Weight"][i][0][len(sumFrame.at[index,"Shared_Ref_Weight"][i][0])-j] == 0:
                            sumFrame.at[index,"Splice_Misses"][2] = 1
                        j = j + 1
                    indices.append(index)
                elif len(sumFrame.at[index,"Shared_Ref_Weight"][i][0]) % 2 != 0:
                    if sumFrame.at[index,"Shared_Ref_Weight"][i][0][0] == 0:
                        sumFrame.at[index,"Splice_Misses"][0] = 1
                    if sumFrame.at[index,"Shared_Ref_Weight"][i][0][len(sumFrame.at[index,"Shared_Ref_Weight"][i][0])-1] == 0:
                        sumFrame.at[index,"Splice_Misses"][3] = 1
                    if sumFrame.at[index,"Shared_Ref_Weight"][i][0][math.floor(len(sumFrame.at[index,"Shared_Ref_Weight"][i][0])/2)+1] == 0:
                        sumFrame.at[index,"Splice_Misses"][1] = 1
                    j = 1
                    while j < math.floor(len(sumFrame.at[index,"Shared_Ref_Weight"][i][0])/2):
                        if sumFrame.at[index,"Shared_Ref_Weight"][i][0][j] == 0:
                            sumFrame.at[index,"Splice_Misses"][1] = 1
                        if sumFrame.at[index,"Shared_Ref_Weight"][i][0][len(sumFrame.at[index,"Shared_Ref_Weight"][i][0])-j] == 0:
                            sumFrame.at[index,"Splice_Misses"][2] = 1
                        j = j + 1
                    indices.append(index)
                i = i + 1

        sumFrame2 = sumFrame.iloc[indices,]
        count1 = 0
        count2 = 0
        count3 = 0
        count4 = 0

        sumFrame2 = sumFrame2.reset_index()
        sumFrame2 = sumFrame2.drop(columns=['index'])

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

    else:
        print("ref variable must equal True or False")

##Reformat Exons: re-adds the exon rows to final output to allow for exporting as a usable GTF.
        
##startFrame: the collapsed dataframe which you wish to readd the exon rows to
        
def reformatExons(startFrame):
    
    ##startFrame: the collapsed dataframe which you wish to readd the exon rows to
    
    endFrame = pd.DataFrame()
    
     ###Create an array of each GTFs chromosones
    chromosones1 = startFrame.seqname.unique().tolist()
    
    mtLoc= False
    if 'MT' in chromosones1:
        mtLoc = True
        
    
    chromosones1.remove('X')
    chromosones1.remove('Y')
    if mtLoc == True:
        chromosones1.remove('MT')
    
    chromosones2 = list(map(int, chromosones1))
    
    chromosones2.sort()
    
    chromosones2.append('X')
    chromosones2.append('Y')
    if mtLoc == True:
        chromosones2.append('MT')
    
    chromosones2 = list(map(str, chromosones2))
    
    for item in chromosones2:
        
        chromeFrame = pd.DataFrame()
        runFrame = startFrame[startFrame["seqname"] == item]
        runFrame = runFrame.reset_index()
        runFrame = runFrame.drop(columns=['index'])
        
        if item == "MT":
            for index, row in runFrame.iterrows():
                runFrame.at[index, "start"] = runFrame.at[index, "Exons"][0][0]
                runFrame.at[index, "end"] = runFrame.at[index, "Exons"][len(runFrame.at[index, "Exons"])-1][1]
        
        for index, row in runFrame.iterrows():
            
            row1 = runFrame.iloc[index]
            rowFrame = pd.DataFrame(row1)
            rowFrame = rowFrame.T
            
            rowFrame = rowFrame.reset_index()
            rowFrame = rowFrame.drop(columns=['index'])
        
            tempFrame = pd.DataFrame()
            k=1
            if rowFrame.at[0,"strand"] == "+":
                for itemExon in rowFrame.at[0,"Exons"]:
                    newRow = pd.DataFrame([[rowFrame.at[0,"seqname"],rowFrame.at[0,"source"],"exon", itemExon[0], itemExon[1], rowFrame.at[0,"score"], rowFrame.at[0,"strand"], rowFrame.at[0,"frame"], rowFrame.at[0,"gene_id"], rowFrame.at[0,"transcript_id"] , k]], columns=['seqname','source','feature', 'start', 'end', 'score', 'strand', 'frame', 'gene_id', 'transcript_id', 'exon_id'])
                    tempFrame = pd.concat([tempFrame, newRow])
                    k = k + 1
            elif rowFrame.at[0,"strand"] == "-":
                for itemExon in reversed(rowFrame.at[0,"Exons"]):
                    newRow = pd.DataFrame([[rowFrame.at[0,"seqname"],rowFrame.at[0,"source"],"exon", itemExon[0], itemExon[1], rowFrame.at[0,"score"], rowFrame.at[0,"strand"], rowFrame.at[0,"frame"], rowFrame.at[0,"gene_id"], rowFrame.at[0,"transcript_id"] , k]], columns=['seqname','source','feature', 'start', 'end', 'score', 'strand', 'frame', 'gene_id', 'transcript_id', 'exon_id'])
                    tempFrame = pd.concat([tempFrame, newRow])
                    k = k + 1

            rowFrame = pd.concat([rowFrame, tempFrame])
            chromeFrame = pd.concat([chromeFrame, rowFrame])
        
        endFrame = pd.concat([endFrame, chromeFrame ])
        print(item, end=" ")
        
    endFrame = endFrame.reset_index()
    endFrame = endFrame.drop(columns=['index'])
    
    return(endFrame)
        
##Write GTF: saves a GTF dataframe as a gtf file

## name --- the name you wish to name the file as

## givenGTF --- the gtf file which you want to save

##extraColumns --- a flag to indicate whether you wish to print more than the basic 11 columns (true will print the extra columns)

def writeGTF(name, givenGTF, extraColumns = True):
    
    correctColumns = ["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "gene_id", "transcript_id", "exon_id" ]
    
    if extraColumns == True:
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
                elif (str(column) == "gene_id") or (str(column) == "transcript_id") or (str(column) == "transcript_name"):
                        k = k + 1
                        string = string + str(column) + " " + "\"" + str(givenGTF.at[index, column]) + "\"" + "; "
                elif extraColumns == True:
                    k = k + 1
                    string = string + str(column) + " " + "\"" + str(givenGTF.at[index, column]) + "\"" + "; "

            
            string = string + "\n"
            fh.write(string)
            
## Smaller or Larger Compared To Reference --- (a function to characterize and compare how a transcript compares to a reference, wether transcripts are on average smaller or larger)

##givenGtf - the gtf dataframe which has already been matched to a reference
##exonThreshold -- distance you qualify as being longer or shorter
                        
def smallLarger(givenGTF, exonThreshold):
    
    perfectEnds = []
    perfectStarts = []
    
    for index, row in givenGTF.iterrows():
        for item in givenGTF.at[index, 'Matched_Ref_Start_Stop']:
            
            if givenGTF.at[index, 'strand'] == "-":
                perfectEnds.append(-item[0])
                perfectStarts.append(item[1])
            else:
                perfectStarts.append(-item[0])
                perfectEnds.append(item[1])
                
        for item in givenGTF.at[index, 'Shared_Ref_Weight']:
            matchedLen = len(givenGTF.at[index, 'Exons'])
            refLen = len(item[0])

            if ((matchedLen-refLen)==0) & (sum(item[2])==matchedLen):
                
                if givenGTF.at[index, 'strand'] == "-":
                    perfectEnds.append(item[1][0] - givenGTF.at[index, 'start'])
                    perfectStarts.append(givenGTF.at[index, 'end'] - item[1][1])
                else:
                    perfectStarts.append(item[1][0] - givenGTF.at[index, 'start'])
                    perfectEnds.append(givenGTF.at[index, 'end'] - item[1][1])
                
    
    perfectsiteData = [perfectEnds, perfectStarts]
    perfectsiteTypes = ["Change in End Site", "Change in Start Site"]
    fullPFrame = pd.DataFrame()
    for i in range(len(perfectsiteData)):
        loopFrame = pd.DataFrame(perfectsiteData[i])
        loopFrame["Site Type"] = perfectsiteTypes[i]
        fullPFrame = fullPFrame.append(loopFrame)
    
    
    
    fullPFrame = fullPFrame.rename(columns={0: "Base Pair Change"})
    
    plt.figure(figsize=(20, 10))
    sns.set(font_scale=2)
    ax = sns.violinplot(x= "Site Type", y="Base Pair Change", data=fullPFrame, inner="quartile")
    
    plt.show()
    
    ##givenGTF: a gtf which has been comapred to a reference which you wish to characterize
    ##exonThreshold: how much larger does a read have to be compared to another to be considered a "larger" transcript, and vice versa for "smaller"
    
    #Select a gtf with only partial matches
    matchGTF  = givenGTF[givenGTF.astype(str)['Shared_Ref_Weight'] != '[]']
    
    #select a gtf with only multiple partial matches
    multipleMatchGTF = matchGTF[matchGTF['Shared_Ref_Weight'].str.len() > 1]
    
    #select a gtf with only single partial matches
    singleMatchGTF = matchGTF[matchGTF['Shared_Ref_Weight'].str.len() == 1]
    
    #add a column of how many matches each transcript has
    matchGTF['Length'] = matchGTF['Shared_Ref_Weight'].str.len()
    
    #get an avergae amount of these matches
    averageNumberMatches = matchGTF["Length"].mean()
    
    #reset index to loop through the single match gtf
    singleMatchGTF = singleMatchGTF.reset_index()
    singleMatchGTF = singleMatchGTF.drop(columns=['index'])
    
    #Get the size of all of these gtfs
    readCounts = matchGTF.shape[0]
    multipleMatchCount = multipleMatchGTF.shape[0]
    singleMatchCount = singleMatchGTF.shape[0]
    
    #print table of results up until now
    t = PrettyTable(['Name', 'Amount', 'Percent'])
    t.add_row(['Total Transcripts With Partial Matches To Reference', readCounts, 100])
    t.add_row(['Transcripts With Multiple Partial Matches', multipleMatchCount, round(((multipleMatchCount)/readCounts)*100, 2)])
    t.add_row(['Transcripts With A Single Partial Match', singleMatchCount, round(((singleMatchCount)/readCounts)*100, 2)])
    t.add_row(['Average # of Matches (For Transcripts With Partial Matches to Reference)', round(averageNumberMatches,2), "N/A"])
    print("\n")
    print(t)

    print("From this point on only single partial match transcripts will be used.")
    
    #create values for counts and list to contain differences
    startBigger = 0
    endBigger = 0
    startSmaller = 0
    endSmaller = 0
    endChange = []
    startChange = []
    
    #Size of the single matches
    totalCount = singleMatchCount
    
    #counts and lists for exon differences
    longer = 0
    shorter = 0
    exonChange = []
    
    zeroExon = []
    zeroExon2 = []
    
    
    for index, row in singleMatchGTF.iterrows():
        for item in singleMatchGTF.at[index, 'Shared_Ref_Weight']:
            matchedLen = len(singleMatchGTF.at[index, 'Exons'])
            refLen = len(item[0])
            
            exonChange.append(matchedLen-refLen)
            if matchedLen > refLen:
                longer = longer + 1
            elif matchedLen < refLen:
                shorter = shorter + 1
            
            if singleMatchGTF.at[index, 'strand'] == "-":
                endChange.append(item[1][0] - singleMatchGTF.at[index, 'start'])
                startChange.append(singleMatchGTF.at[index, 'end'] - item[1][1])
            else:
                startChange.append(item[1][0] - singleMatchGTF.at[index, 'start'])
                endChange.append(singleMatchGTF.at[index, 'end'] - item[1][1])
                
            
            if (item[0][0][0] - singleMatchGTF.at[index, 'Exons'][0][1] >= exonThreshold):
                if singleMatchGTF.at[index, 'strand'] == "-":
                    endBigger = endBigger+1
                else:
                    startBigger = startBigger+1
                
            if (singleMatchGTF.at[index, 'Exons'][matchedLen-1][0] - item[0][refLen-1][1]>= exonThreshold):
                if singleMatchGTF.at[index, 'strand'] == "-":
                    startBigger = startBigger+1
                else:
                    endBigger = endBigger+1
                    
            if (singleMatchGTF.at[index, 'Exons'][0][0] - item[0][0][1] >= exonThreshold):
                if singleMatchGTF.at[index, 'strand'] == "-":
                    endSmaller = endSmaller+1
                else:
                    startSmaller = startSmaller+1
                    
            if (item[0][refLen-1][0] - singleMatchGTF.at[index, 'Exons'][matchedLen-1][1] >= exonThreshold):
                if singleMatchGTF.at[index, 'strand'] == "-":
                    startSmaller = startSmaller+1
                else:
                    endSmaller = endSmaller+1
                    
    #Print another table of the results so far                
    t = PrettyTable(['Name', 'Amount', 'Percent'])
    t.add_row(['Total Partial Matches', totalCount, 100])
    t.add_row(['More Exons (Count)', longer, round(((longer)/totalCount)*100, 3)])
    t.add_row(['Less Exons (Count)', shorter, round(((shorter)/totalCount)*100, 3)])
    t.add_row(['Longer End (Exon Past Reference Read End)', endBigger, round(((endBigger)/totalCount)*100, 3)])
    t.add_row(['Shorter End (Refernce Exon Past Read End)', endSmaller, round(((endSmaller)/totalCount)*100, 3)])
    t.add_row(['Longer Start (Exon Past Reference Read Start)', startBigger, round(((startBigger)/totalCount)*100, 3)])
    t.add_row(['Shorter Start (Refernce Exon Past Read Start)', startSmaller, round(((startSmaller)/totalCount)*100, 3)])
    print("\n")
    print(t)
    
    
    print("max exon change:", max(exonChange))
    print("min exon change:", min(exonChange))
    
    exondf = pd.DataFrame(exonChange)
    
    exonQuant = exondf[0].quantile([.1, .25, .5, .75, .9])
    
    print("Exon Quantitles:")
    print(exonQuant)
    
    plt.figure(figsize=(20, 10))
    sns.set(font_scale=2)
    ax = sns.violinplot(data=exonChange, inner="quartile") 
    
    plt.ylim(-40, 40)
    
    plt.show()
    
    print("max start change:", max(startChange))
    print("min start change:", min(startChange))
    
    print("max end change:", max(endChange))
    print("min end change:", min(endChange))
    
    
    siteData = [endChange, startChange]
    siteTypes = ["Change in End Site", "Change in Start Site"]
    fullSFrame = pd.DataFrame()
    for i in range(len(siteData)):
        loopFrame = pd.DataFrame(siteData[i])
        loopFrame["Site Type"] = siteTypes[i]
        fullSFrame = fullSFrame.append(loopFrame)

    fullSFrame = fullSFrame.rename(columns={0: "Base Pair Change"})
    
    plt.figure(figsize=(20, 10))
    sns.set(font_scale=2)
    ax = sns.violinplot(x= "Site Type", y="Base Pair Change", data=fullSFrame, inner="quartile")
    
    plt.show()                     


##Transcriptome Master Function: Constructs a Transcriptome from a list of GTFs

##gtfList --- an array of gtfs which you wish to be combined

##gtfLabels --- an array of the associated name you wish to label each GTF

##threshold --- threshold you consider acceptable for distance between splice junction start and stops to be counted as a match

##distanceThreshold --- threshold you consider acceptable between start and stop sites to be counted as a match


def transcriptomeMasterFunction2(GTFList, GTFLabels, threshold, distanceThreshold):
    
    if len(GTFList) <= 1:
        print("At least two GTFs are required to assemble transcriptome")
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
                    chromeFrame1.at[index, "gene_id"] = chromeFrame1.at[index, "transcript_name"] + "/" + GTFLabels[i] + "-" + str(item) + str(index)
                else:
                    chromeFrame1.at[index, "gene_id"] = chromeFrame1.at[index, "transcript_id"] + "/" + GTFLabels[i] + "-" + str(item) + str(index)
            
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
            chromeFrame1['Shared_Located_In'] = np.empty((len(chromeFrame1), 0)).tolist()

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
                        workingDataframe.at[index, 'Shared_SJ_Weight'].append([exonMatches, len(workingDataframe.at[index+1,"Splice_Junctions"])])
                        workingDataframe.at[index+1, 'Shared_SJ_Transcript'].append(workingDataframe.at[index,"gene_id"]) 
                        workingDataframe.at[index+1, 'Shared_SJ_Weight'].append([exonMatches, len(workingDataframe.at[index,"Splice_Junctions"])])
                        if workingDataframe.at[index+1,"Located_In"][0] not in workingDataframe.at[index,"Shared_Located_In"]:
                            workingDataframe.at[index,"Shared_Located_In"].append(workingDataframe.at[index+1,"Located_In"][0])
                        if workingDataframe.at[index,"Located_In"][0] not in workingDataframe.at[index+1,"Shared_Located_In"]:
                            workingDataframe.at[index+1,"Shared_Located_In"].append(workingDataframe.at[index,"Located_In"][0])

                looping = 2
                while index+looping < len(workingDataframe.index) and workingDataframe.at[index, "end"] >= workingDataframe.at[index+looping,"start"]:

                    #Match multi-exon length genes
                    if len(workingDataframe.at[index,"Splice_Junctions"]) != 0:
                        exonMatches = exonMatching(workingDataframe.at[index,"Splice_Junctions"], workingDataframe.at[index+looping,"Splice_Junctions"], threshold)
                        if sum(exonMatches) >= 1:
                            workingDataframe.at[index, 'Shared_SJ_Transcript'].append(workingDataframe.at[index+looping,"gene_id"]) 
                            workingDataframe.at[index, 'Shared_SJ_Weight'].append([exonMatches, len(workingDataframe.at[index+looping,"Splice_Junctions"])])
                            workingDataframe.at[index+looping, 'Shared_SJ_Transcript'].append(workingDataframe.at[index,"gene_id"]) 
                            workingDataframe.at[index+looping, 'Shared_SJ_Weight'].append([exonMatches, len(workingDataframe.at[index,"Splice_Junctions"])])
                            if workingDataframe.at[index+looping,"Located_In"][0] not in workingDataframe.at[index,"Shared_Located_In"]:
                                workingDataframe.at[index,"Shared_Located_In"].append(workingDataframe.at[index+looping,"Located_In"][0])
                            if workingDataframe.at[index,"Located_In"][0] not in workingDataframe.at[index+looping,"Shared_Located_In"]:
                                workingDataframe.at[index+looping,"Shared_Located_In"].append(workingDataframe.at[index,"Located_In"][0])

                    ##Itterate
                    looping = looping + 1
        
        MasterFrameGeneMatch = MasterFrameGeneMatch.append(workingDataframe)
        
        print(item, end=" ")
        
    #create new index to help itteration in our working frame
    MasterFrameGeneMatch = MasterFrameGeneMatch.reset_index()
    MasterFrameGeneMatch = MasterFrameGeneMatch.drop(columns=['index'])
    
    return(MasterFrameGeneMatch)



##gtfTranscript Reference Matcher -- Match a collapsed gtf, or a merged gtf with a reference gtf

##gtfStart -- the gtf dataframe you want to compare to the reference

##gtfReference -- the gtf dataframe of your reference

##referenceLabel -- the name you want to use with your reference

##threshold --- threshold you consider acceptable for distance between splice junction start and stops to be counted as a match

##distanceThreshold --- threshold you consider acceptable between start and stop sites to be counted as a match


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
            chromeFrame1['Shared_Ref_Weight'] = np.empty((len(chromeFrame1), 0)).tolist()
            chromeFrame1['Overlap_Ref_Transcript'] = np.empty((len(chromeFrame1), 0)).tolist()
            
            if 'transcript_id' in chromeFrame2.columns:
                if 'transcript_name' not in chromeFrame2.columns:
                    chromeFrame2 = chromeFrame2.rename(columns={"transcript_id": "transcript_name"})

            #Sort values based on their start position
            chromeFrame1 = chromeFrame1.sort_values(by=['start'])
            chromeFrame2 = chromeFrame2.sort_values(by=['start'])
            
            #create new index to help itteration
            chromeFrame1 = chromeFrame1.reset_index()
            chromeFrame1 = chromeFrame1.drop(columns=['index'])

            chromeFrame2 = chromeFrame2.reset_index()
            chromeFrame2 = chromeFrame2.drop(columns=['index'])
            
            for index, row in chromeFrame2.iterrows():
                Exons2 = chromeFrame2.at[index,"Exons"]
                chromeFrame2.at[index,"Exons"] = sorted(Exons2)
                
            chromeFrame2 = chromeFrame2.reset_index()
            chromeFrame2 = chromeFrame2.drop(columns=['index'])

            #Add splice junc column for reference gtf
            chromeFrame2 = splicJuncColumn(chromeFrame2)

            itterator = 0
            for index, row in chromeFrame1.iterrows():

                #Move Along Rows in the Reference Until It Is Past Beginning of Current read
                while itterator < len(chromeFrame2.index)-1 and chromeFrame1.at[index, 'start'] > chromeFrame2.at[itterator,"end"]:
                    itterator = itterator + 1
                    
                if (chromeFrame2.at[itterator,"start"] <= chromeFrame1.at[index, 'start'] <= chromeFrame2.at[itterator,"end"]) or (chromeFrame2.at[itterator,"start"] <= chromeFrame1.at[index, 'end'] <= chromeFrame2.at[itterator,"end"]) or (chromeFrame1.at[index,"start"] <= chromeFrame2.at[itterator, 'end'] <= chromeFrame1.at[index,"end"]):
                    
                    #Only match one exon length genes with one exon length genes
                    if len(chromeFrame1.at[index,"Splice_Junctions"]) == 0:
                        if abs(chromeFrame1.at[index, 'start']-chromeFrame2.at[itterator,"start"])<=distanceThreshold and abs(chromeFrame1.at[index, 'end']-chromeFrame2.at[itterator,"end"])<=distanceThreshold and len(chromeFrame2.at[itterator,"Splice_Junctions"]) == 0:
                            if referenceLabel not in chromeFrame1.at[index,"Located_In"]:
                                chromeFrame1.at[index,"Located_In"].append(referenceLabel)
                            chromeFrame1.at[index, 'Matched_Reference_Gene'].append(chromeFrame2.at[itterator,"gene_id"])
                            chromeFrame1.at[index, 'Matched_Reference_Transcript'].append(chromeFrame2.at[itterator,"transcript_name"])
                            chromeFrame1.at[index, "Matched_Ref_Start_Stop"].append([chromeFrame1.at[index,"start"] - chromeFrame2.at[itterator,"start"], chromeFrame1.at[index,"end"] - chromeFrame2.at[itterator,"end"]])
                        else:
                            chromeFrame1.at[index, 'Overlap_Ref_Transcript'].append(chromeFrame2.at[itterator,"transcript_name"])
                            
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
                            chromeFrame1.at[index, 'Shared_Ref_Weight'].append([chromeFrame2.at[itterator,"Exons"], [chromeFrame2.at[itterator,"start"], chromeFrame2.at[itterator,"end"]], exonMatches])
                        else:
                            chromeFrame1.at[index, 'Overlap_Ref_Transcript'].append(chromeFrame2.at[itterator,"transcript_name"])
                            
                    looping = 1
                    while itterator+looping < len(chromeFrame2.index) and chromeFrame1.at[index, 'end'] >= chromeFrame2.at[itterator+looping,"start"]:
                        
                        if (chromeFrame2.at[itterator+looping,"start"] <= chromeFrame1.at[index, 'start'] <= chromeFrame2.at[itterator+looping,"end"]) or (chromeFrame2.at[itterator+looping,"start"] <= chromeFrame1.at[index, 'end'] <= chromeFrame2.at[itterator+looping,"end"]) or (chromeFrame1.at[index,"start"] <= chromeFrame2.at[itterator+looping, 'end'] <= chromeFrame1.at[index,"end"]):
                            
                            if len(chromeFrame1.at[index,"Splice_Junctions"]) == 0 and len(chromeFrame2.at[itterator+looping,"Splice_Junctions"]) == 0:
                                if abs(chromeFrame1.at[index, 'start']-chromeFrame2.at[itterator+looping,"start"])<=distanceThreshold and abs(chromeFrame1.at[index, 'end']-chromeFrame2.at[itterator+looping,"end"])<=distanceThreshold:
                                    if referenceLabel not in chromeFrame1.at[index,"Located_In"]:
                                        chromeFrame1.at[index,"Located_In"].append(referenceLabel)
                                    chromeFrame1.at[index, 'Matched_Reference_Gene'].append(chromeFrame2.at[itterator+looping,"gene_id"])
                                    chromeFrame1.at[index, 'Matched_Reference_Transcript'].append(chromeFrame2.at[itterator+looping,"transcript_name"])
                                    chromeFrame1.at[index,"Matched_Ref_Start_Stop"].append([chromeFrame1.at[index,"start"] - chromeFrame2.at[itterator+looping,"start"], chromeFrame1.at[index,"end"] - chromeFrame2.at[itterator+looping,"end"]])
                                else:
                                    chromeFrame1.at[index, 'Overlap_Ref_Transcript'].append(chromeFrame2.at[itterator+looping,"transcript_name"])

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
                                    chromeFrame1.at[index, 'Shared_Ref_Weight'].append([chromeFrame2.at[itterator+looping,"Exons"], [chromeFrame2.at[itterator+looping,"start"], chromeFrame2.at[itterator+looping,"end"]], exonMatches])
                                else:
                                    chromeFrame1.at[index, 'Overlap_Ref_Transcript'].append(chromeFrame2.at[itterator+looping,"transcript_name"])
                                    
                        looping = looping + 1

            MasterFrameGeneMatch = MasterFrameGeneMatch.append(chromeFrame1) #add our matched data to the master dataframe
            print(item,end = ' ') #Display each chromosone as we finish it
    
    MasterFrameGeneMatch = MasterFrameGeneMatch.reset_index()
    MasterFrameGeneMatch = MasterFrameGeneMatch.drop(columns=['index'])
    
    return(MasterFrameGeneMatch)


##gtf Symbiosis: Master function which both does the Symbiosis and also displays the results

##analyzedGTFs --- an array of gtfs which you wish to be combined

##analyzedGTFLabels --- an array of the associated name you wish to label each GTF

##threshold --- threshold you consider acceptable for distance between splice junction start and stops to be counted as a match

##distanceThreshold --- threshold you consider acceptable between start and stop sites to be counted as a match

##gtfReference -- the gtf dataframe of your reference

##referenceLabel -- the name you want to use with your reference

##refdistanceThreshold --- threshold you consider acceptable between start and stop sites to be counted as a match for your reference


def gtfSymbiosis(analyzedGTFs, analyzedGTFLabels, threshold, distanceThreshold, referenceGTF = "", referenceGTFLabel = "", refdistanceThreshold = 0):
    
    #Print Out Basic Data of Both GTFs
    t = PrettyTable(['Name', 'Transcripts', 'Exons', 'Genes']) 
    i = 0
    while i < len(analyzedGTFs):
        runGTf = analyzedGTFs[i]
        exonSizeAn = len(runGTf[runGTf["feature"] == "exon"].index)
        tranSizeAn = len(runGTf[runGTf["feature"] == "transcript"].index)
        geneSizeAn = len(runGTf[runGTf["feature"] == "gene"].index)
        t.add_row([analyzedGTFLabels[i], tranSizeAn, exonSizeAn, geneSizeAn])
        i = i + 1
    if isinstance(referenceGTF, pd.DataFrame):
        exonSizeRef = len(referenceGTF[referenceGTF["feature"] == "exon"].index)
        tranSizeRef = len(referenceGTF[referenceGTF["feature"] == "transcript"].index)
        geneSizeRef = len(referenceGTF[referenceGTF["feature"] == "gene"].index)
        t.add_row(['Reference', tranSizeRef, exonSizeRef, geneSizeRef])

    print(t)
        
    #Run the master comparison Fucntion
    
    exonMatchedGTF = transcriptomeMasterFunction2(analyzedGTFs, analyzedGTFLabels, threshold, distanceThreshold)
    print("\n")
    t = PrettyTable(['Name', 'Size']) 
    t.add_row(['Final Condensed Size', len(exonMatchedGTF.index)])
    print(t)
    
    if isinstance(referenceGTF, pd.DataFrame):
        
        exonMatchedGTF2 = gtfTranscriptomeReference(exonMatchedGTF, referenceGTF, referenceGTFLabel, threshold, refdistanceThreshold)
        
        #Find the Number of Perfect Matches
        perfectMatches = exonMatchedGTF2[exonMatchedGTF2.astype(str)['Matched_Reference_Gene'] != '[]']
        perfectMatchNum = perfectMatches.shape[0]

        #Find Number of genuine non-matched transcripts
        nonMatched = exonMatchedGTF2[exonMatchedGTF2.astype(str)['Matched_Reference_Gene'] == '[]']
        nonMatched = nonMatched[nonMatched.astype(str)['Reference_Gene_Partial'] == '[]']
        nonMatchedNum = nonMatched.shape[0]

        #Total number of Transcripts
        totalSize = exonMatchedGTF2.shape[0]
        
        #find size of partial matched transcripts
        partialMatches = exonMatchedGTF2[exonMatchedGTF2.astype(str)['Reference_Gene_Partial'] != '[]']
        partialMatches = partialMatches[partialMatches.astype(str)['Matched_Reference_Gene'] == '[]']
        partialMatchesNum = partialMatches.shape[0]

        #Find Missing Splice Junction Locations
        count1, count2, count3, count4 = splicJuncMissing(exonMatchedGTF2, ref = True)

        # Make our dataset ready for graphing:
        barWidth = 0.9
        height = [nonMatchedNum, perfectMatchNum, partialMatchesNum]
        bars = ['Not Matched To Reference (Partial or Full)', 'Perfect Matches', 'Partial Matches (Perfect Matches Exluded)']
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
        axs[1].set_title('Missing Splice Junction Locations (Compared to Refernce)')
        axs[1].set_ylabel('Count')
        axs[1].set_xticks(x2+1)
        axs[1].set_xticklabels(bars2, rotation=40, ha='right')

        # Create names on the x-axis
        label = [str(nonMatchedNum), str(perfectMatchNum), str(partialMatchesNum)]

        for i in range(len(rPlacement)):
            axs[0].text(x = rPlacement[i]-.25 , y = height[i]+(height[i]/50), s = label[i], size = 18)

        # Create names on the x-axis
        label2 = [str(count1), str(count2), str(count3), str(count4)]

        for i in range(len(rPlacement2)):
            axs[1].text(x = rPlacement2[i]-.25 , y = height2[i]+(height2[i]/50), s = label2[i], size = 18)

        # Show graphic
        plt.show()
        
        ##Create perfect match frame
        perfectMatchFrame = perfectMatches
        
        ##Create multi-perfect match frame
        starts = []
        stops = []
        for index, row in perfectMatchFrame.iterrows():
            i = 0
            while i < len(perfectMatchFrame.at[index, "Matched_Ref_Start_Stop"]):
                if str(perfectMatchFrame.at[index, "strand"]) == "-":
                    starts.append(perfectMatchFrame.at[index, "Matched_Ref_Start_Stop"][i][1]*(-1))
                    stops.append(perfectMatchFrame.at[index, "Matched_Ref_Start_Stop"][i][0]*(-1))
                elif str(perfectMatchFrame.at[index, "strand"]) == "+":
                    starts.append(perfectMatchFrame.at[index, "Matched_Ref_Start_Stop"][i][0])
                    stops.append(perfectMatchFrame.at[index, "Matched_Ref_Start_Stop"][i][1])
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
        
        #Print Basic Table
        t = PrettyTable(['Quantile', 'Start', 'Stop'])
        t.add_row(['.1', quantFrameStart.iloc[0,0], quantFrameStop.iloc[0,0]])
        t.add_row(['.25', quantFrameStart.iloc[1,0], quantFrameStop.iloc[1,0]])
        t.add_row(['.5', quantFrameStart.iloc[2,0], quantFrameStop.iloc[2,0]])
        t.add_row(['.75', quantFrameStart.iloc[3,0], quantFrameStop.iloc[3,0]])
        t.add_row(['.9', quantFrameStart.iloc[4,0], quantFrameStop.iloc[4,0]])                            
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
        
        unmatchedPercentage = round((nonMatchedNum/totalSize)*100)
        data = {'Matched Transcripts': 100-unmatchedPercentage, 'Unmatched Transcripts': unmatchedPercentage}
        fig = plt.figure(
            FigureClass=Waffle, 
            rows=5, 
            values=data, 
            colors=("#983D3D", "#232066"),
            title={'label': 'Amount of Transcripts Matched To Reference Transcriptome', 'loc': 'left'},
            labels=["{0} ({1}%)".format(k, v) for k, v in data.items()],
            legend={'loc': 'lower left', 'bbox_to_anchor': (0, -0.4), 'ncol': len(data), 'framealpha': 0}
        )
        fig.gca().set_facecolor('#EEEEEE')
        fig.set_facecolor('#EEEEEE')
        plt.show()

        #Find Total Amount Of Exons Matched
        totalSJMatch, SJMasterSize = sumTotalSJMatched(exonMatchedGTF2)

        #Plot Splice Junction Span in Waffle Plot
        data = {'Matched Splice Junctions': round((totalSJMatch/SJMasterSize)*100), 'Unmatched Splice Junctions': 100-round((totalSJMatch/SJMasterSize)*100)}
        fig = plt.figure(
            FigureClass=Waffle, 
            rows=5, 
            values=data, 
            colors=("#983D3D", "#232066"),
            title={'label': 'Amount Of Splice Junctions Matched To Reference Transcriptome', 'loc': 'left'},
            labels=["{0} ({1}%)".format(k, v) for k, v in data.items()],
            legend={'loc': 'lower left', 'bbox_to_anchor': (0, -0.4), 'ncol': len(data), 'framealpha': 0}
        )
        fig.gca().set_facecolor('#EEEEEE')
        fig.set_facecolor('#EEEEEE')
        plt.show()
        
        
        #Create each set for our venn diagram dictionary
        k = 0
        vennDict = {}
        labelsList = analyzedGTFLabels
        labelsList.append(referenceGTFLabel)
        while k < len(labelsList):
            vennDict[labelsList[k]] = {"default"}
            k = k + 1
            
        for index, row in exonMatchedGTF.iterrows():
            for item in exonMatchedGTF.at[index, "Located_In"]:
                if (item in labelsList):
                    vennDict[item].add(exonMatchedGTF.at[index, "transcript_id"])
        
        r = dict(vennDict)
        
        k=0
        while k < len(labelsList):
            vennDict[labelsList[k]].remove("default")
            k = k + 1
        if len(labelsList) == 6:
            pseudovenn(r)
        elif len(labelsList) <= 6:
            venn(r)
        else:
            print("More than 6 GTFs combined, result will not be displayed as a venn diagram.")
            print("\n")
        
        #Create our PYU plot
        pyuSeries = from_contents(vennDict)
        plot(pyuSeries)
        plt.show()
        
        
        
        return(exonMatchedGTF2)
    
    elif referenceGTF == "":
        
        #Find the Number of Perfect Matches
        perfectGTF = exonMatchedGTF[exonMatchedGTF.astype(str)['Matched_Transcript'] != '[]']
        noSingleExonPerfectGTF = perfectGTF[perfectGTF.astype(str)['Splice_Junctions'] != '[]']
        perfectNum = len(perfectGTF.index)
        noSingleExonPerfectNum = len(noSingleExonPerfectGTF.index)
        
        print("\n")
        
        t = PrettyTable(['# Perfect Matches (Multi Exon)', '# of Single Exon Matches'])
        t.add_row([noSingleExonPerfectNum,perfectNum - noSingleExonPerfectNum])                           
        print(t)
        
        #Find all non-perfect matches & no matches
        noPerfectGTF = exonMatchedGTF[exonMatchedGTF.astype(str)['Matched_Transcript'] == '[]']
        noPerfectNum = len(noPerfectGTF.index)
        
        #Find no-matches
        genuineNoMatch = noPerfectGTF[noPerfectGTF.astype(str)['Shared_SJ_Transcript'] == '[]']
        genuineNoMatchNum = len(genuineNoMatch.index)
        
        #Find Missing Splice Junction Locations
        allmatchFrame = exonMatchedGTF[exonMatchedGTF.astype(str)['Shared_SJ_Transcript'] != '[]']
        count1, count2, count3, count4 = splicJuncMissing(allmatchFrame)
        
        #Find all matched transcripts from second gtf
        matched2List = []
        for index, row in exonMatchedGTF.iterrows():
            i = 0
            while i < len(exonMatchedGTF.at[index,"Shared_SJ_Transcript"]):
                if not exonMatchedGTF.at[index,"Shared_SJ_Transcript"][i]:
                    pass
                else: 
                    matched2List.extend(exonMatchedGTF.at[index,"Shared_SJ_Transcript"][i])
                i = i + 1
        matched2ListUnique = list(OrderedDict.fromkeys(matched2List)) 
        matched2 = len(matched2ListUnique)
        
        # Make our dataset ready for graphing:
        barWidth = 0.9
        height = [genuineNoMatchNum, (perfectNum+noPerfectNum-genuineNoMatchNum)]
        bars = ['No Matching Transcript', 'Matching (All)']
        x = np.arange(len(bars))
        rPlacement = [1,2]
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
        label = [str(genuineNoMatchNum), str(perfectNum+noPerfectNum-genuineNoMatchNum)]

        for i in range(len(rPlacement)):
            axs[0].text(x = rPlacement[i]-.25 , y = height[i]+(height[i]/50), s = label[i], size = 18)

        # Create names on the x-axis
        label2 = [str(count1), str(count2), str(count3), str(count4)]

        for i in range(len(rPlacement2)):
            axs[1].text(x = rPlacement2[i]-.25 , y = height2[i]+(height2[i]/50), s = label2[i], size = 18)

        # Show graphic
        plt.show()
    
        ##Create multi-perfect match frame
        starts = []
        stops = []
        for index, row in noSingleExonPerfectGTF.iterrows():
            i = 0
            while i < len(noSingleExonPerfectGTF.at[index, "Start_Stop_Distance"]):
                if str(noSingleExonPerfectGTF.at[index, "strand"]) == "-":
                    starts.append(noSingleExonPerfectGTF.at[index, "Start_Stop_Distance"][i][1]*(-1))
                    stops.append(noSingleExonPerfectGTF.at[index, "Start_Stop_Distance"][i][0]*(-1))
                elif str(noSingleExonPerfectGTF.at[index, "strand"]) == "+":
                    starts.append(noSingleExonPerfectGTF.at[index, "Start_Stop_Distance"][i][0])
                    stops.append(noSingleExonPerfectGTF.at[index, "Start_Stop_Distance"][i][1])
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
        
        
        t = PrettyTable(['Total # Multi-Match Starts/Stops'])
        t.add_row([len(startArray)])                           
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
        
        #Create each set for our venn diagram dictionary
        k = 0
        vennDict = {}
        while k < len(analyzedGTFLabels):
            vennDict[analyzedGTFLabels[k]] = {"default"}
            k = k + 1
        
        for index, row in exonMatchedGTF.iterrows():
            for item in exonMatchedGTF.at[index, "Located_In"]:
                if (item in analyzedGTFLabels):
                    vennDict[item].add(exonMatchedGTF.at[index, "gene_id"])
        
        r = dict(vennDict)
        
        k=0
        while k < len(analyzedGTFLabels):
            vennDict[analyzedGTFLabels[k]].remove("default")
            k = k + 1
        if len(analyzedGTFLabels) == 6:
            pseudovenn(r)
        elif len(analyzedGTFLabels) <= 6:
            venn(r)
        else:
            print("More than 6 GTFs combined, result will not be displayed as a venn diagram.")
            print("\n")
        
        #Create our PYU plot
        pyuSeries = from_contents(vennDict)
        plot(pyuSeries)
        plt.show()
        

        return(exonMatchedGTF)
    
    
##Gtf Transcriptome 2: Master fucntion which both assembles the Transcriptome and displays the results

##analyzedGTFs --- list of GTFs you wish combine into a transcriptome

##analyzedLabels --- list of names you wish to associate with each gtf

##exonMatchThreshold --- distance which you wish to count as a matching splice junction

##distanceThreshold --- distance which you wish to count as a matching transcript

##gtfReference --- a gtf dataframe of a reference like Ensembl which your final transcirptome can be compared against

##referenceLabel --- the name you wish to associate with your reference

##oneExonSeperate --- a bianry flag that determines whether you want the final analysis seperated by multi-exon and one-exon (default value is False, which corresponds to not splitting the final analysis)

    
def gftTranscriptome2(analyzedGTFs,analyzedLabels, exonMatchThreshold, distanceThreshold, gtfReference = "", referenceLabel ="", oneExonSeperate = False, vennReference = True):
    
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
    
    #If one exon seperate flag is on, seperate our GTF in two
    if oneExonSeperate == True:
        oneExonMatchedGTF = exonMatchedGTF[(exonMatchedGTF.Splice_Junctions.str.len() == 0)]
        multiexonMatchedGTF = exonMatchedGTF[(exonMatchedGTF.Splice_Junctions.str.len() != 0)]
        t.add_row(['Final Size', exonMatchedGTF.shape[0]])
        t.add_row(['One Exon Size', oneExonMatchedGTF.shape[0]])
        t.add_row(['Multi Exon Size', multiexonMatchedGTF.shape[0]])
    else:
        t.add_row(['Final Size', exonMatchedGTF.shape[0]])
    print("\n")
    print(t)
    
    #Create the network to mark our genes
    print("Creating Connections:")
    my_network_data = []
    
    #Get size to print out intervals
    sizeDF = len(exonMatchedGTF.index)
    numbersInterval = [round(sizeDF/4), round(2*sizeDF/4), round(3*sizeDF/4), sizeDF-1]
    terms = ["25%", "50%", "75%", "100%"] 
    
    exonMatchedGTF = exonMatchedGTF.reset_index()
    exonMatchedGTF = exonMatchedGTF.drop(columns=['index'])
    
    #Add Each connection from our data
    for index, row in exonMatchedGTF.iterrows():
        
        for item in exonMatchedGTF.at[index, "Shared_SJ_Transcript"]:
            if item is not None:
                connection = [item, exonMatchedGTF.at[index, "gene_id"]]
                connection.sort()
                my_network_data.append(connection)
        if index in numbersInterval:
            locationIndex = numbersInterval.index(index)
            print(terms[locationIndex], end=" ")

    #Create Networkx Graph    
    print("\n")
    g = nx.Graph()
    g.add_edges_from(my_network_data)
    
    #Get ready to add new gene and transcript names
    exonMatchedGTF2 = exonMatchedGTF.copy()
    exonMatchedGTF2["new_gene_id"] = "None"
    exonMatchedGTF2["new_transcript_id"] = None
    exonMatchedGTF2["solitaryGene"] = False
    exonMatchedGTF2 = exonMatchedGTF2.set_index('gene_id')
    fullsize = len(exonMatchedGTF2.index)
    fulllength = len(str(fullsize))
    
    #Get the size of our network to print out percentage intervals
    sizeDF = nx.number_connected_components(g)
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
    
    exonMatchedGTF2 = exonMatchedGTF2.reset_index()
    
    #Add labels for solitary genes and add start and stops to make largest possible transcript    
    for index, row in exonMatchedGTF2.iterrows():
        if str(exonMatchedGTF2.at[index, "new_gene_id"]) == "None":
            kstr = str(k)
            zero_filled_k = kstr.zfill(fulllength)
            exonMatchedGTF2.at[index, "new_gene_id"] = "COISORAT." + zero_filled_k
            exonMatchedGTF2.at[index, "new_transcript_id"] = "COISORAT." + zero_filled_k + "." + str(1)
            k=k+1 
            exonMatchedGTF2.at[index, "solitaryGene"] = True
        #Add start and stops if it makes the size larger
        for item in exonMatchedGTF2.at[index, "Start_Stop_Distance"]:
            if not item:
                print(exonMatchedGTF2.iloc[index])
            if exonMatchedGTF2.at[index, "start"] + item[0] < exonMatchedGTF2.at[index, "start"]:
                exonMatchedGTF2.at[index, "start"] = exonMatchedGTF2.at[index, "start"] + item[0]
            if exonMatchedGTF2.at[index, "end"] + item[1] < exonMatchedGTF2.at[index, "end"]:
                exonMatchedGTF2.at[index, "end"] = exonMatchedGTF2.at[index, "end"] - item[1]
        exonMatchedGTF2.at[index, "Matched_Transcript"].append(exonMatchedGTF2.at[index, "transcript_id"])
       
    print("\n")
    if oneExonSeperate == True:
        oneExonMatchedGTF = exonMatchedGTF2[(exonMatchedGTF2.Splice_Junctions.str.len() == 0)]
        multiexonMatchedGTF = exonMatchedGTF2[(exonMatchedGTF2.Splice_Junctions.str.len() != 0)]

        multiexonMatchedGTFGN = multiexonMatchedGTF[multiexonMatchedGTF["solitaryGene"]==False]
        multiExonGNNum = len(multiexonMatchedGTFGN.new_gene_id.unique())
        multiExonSGNum = multiexonMatchedGTF.shape[0]-multiexonMatchedGTFGN.shape[0]

        print("Number of Solitary Genes (One Exon): " + str(oneExonMatchedGTF.shape[0]))
        
        print("Number of Gene Networks (Multi Exon): " + str(multiExonGNNum))
        print("Number of Solitary Genes (Multi Exon): " + str(multiExonSGNum))
        
        print(nx.info(g))
        
    else:
        #Print Out the details of the network we made
        print(nx.info(g))
        print("Number of Gene Networks: " + str(sizeDF))
        print("Number of Genes: " + str(k))

    #Clean Up Columns and column names
    exonMatchedGTF2 = exonMatchedGTF2.drop(columns=['found'])
    exonMatchedGTF2 = exonMatchedGTF2.drop(columns=['transcript_id'])
    exonMatchedGTF2 = exonMatchedGTF2.drop(columns=['gene_id'])
    exonMatchedGTF2 = exonMatchedGTF2.drop(columns=['Previously_Matched'])
    exonMatchedGTF2 = exonMatchedGTF2.rename(columns={"new_gene_id": "gene_id", "new_transcript_id": "transcript_id"})
    
    if isinstance(gtfReference, pd.DataFrame):
        
        ##Match Gtf to our refernce
        exonMatchedGTF3 = gtfTranscriptomeReference(exonMatchedGTF2, gtfReference, referenceLabel, exonMatchThreshold, distanceThreshold)
        exonMatchedGTF4 = exonMatchedGTF3.copy()
        
        
        #if we are seperating out the one exons do it at this stage
        if oneExonSeperate == True:
            oneExonMatchedGTF3 = exonMatchedGTF3[(exonMatchedGTF3.Splice_Junctions.str.len() == 0)]
            exonMatchedGTF3 = exonMatchedGTF3[(exonMatchedGTF3.Splice_Junctions.str.len() != 0)]
            print("\n")
            print("Beginning Multi Exon Analysis:")
        
        #Full size of our final GTF
        fullRefSize = len(exonMatchedGTF3.index)
        
        #Get size of gtf that didnt match anything perfectly in the reference
        nonReference = exonMatchedGTF3[exonMatchedGTF3.astype(str)['Matched_Reference_Gene'] == '[]']
        nonReferenceSize = len(nonReference.index)
        
        #get size of gtf that didnt match anything partially in the reference
        nonReferencePart = exonMatchedGTF3[exonMatchedGTF3.astype(str)['Reference_Transcript_Partial'] == '[]']
        nonReferencePartSize = len(nonReferencePart.index)
        
        #get size of the GTf that matched something paritally and perfectly
        doubleMatch = exonMatchedGTF3[exonMatchedGTF3.astype(str)['Reference_Transcript_Partial'] != '[]']
        doubleMatch = doubleMatch[doubleMatch.astype(str)['Matched_Reference_Gene'] != '[]']
        doubleMatchSize = len(doubleMatch.index)
        
        ##Find Gene Coverage
        genesCoveredGN = []
        genesCoveredSG = []
        for index, row in exonMatchedGTF3.iterrows():
            if exonMatchedGTF3.at[index, "Matched_Reference_Gene"]:
                if exonMatchedGTF3.at[index, "solitaryGene"] == True:
                    genesCoveredSG.append(exonMatchedGTF3.at[index, "gene_id"])
                elif exonMatchedGTF3.at[index, "solitaryGene"] == False:
                    genesCoveredGN.append(exonMatchedGTF3.at[index, "gene_id"])
        
        #Get the gene count lengths
        genesCoveredSGFinal = list(set(genesCoveredSG))
        genesCoveredSGNum = len(genesCoveredSGFinal)
        
        genesCoveredGNFinal = list(set(genesCoveredGN))
        genesCoveredGNNum = len(genesCoveredGNFinal)
        
        #Get all solitary genes and gene networks
        genesCoveredGN2 = []
        genesCoveredSG2 = []
        for index, row in exonMatchedGTF3.iterrows():
            if exonMatchedGTF3.at[index, "Reference_Transcript_Partial"]:
                if exonMatchedGTF3.at[index, "solitaryGene"] == True:
                    genesCoveredSG2.append(exonMatchedGTF3.at[index, "gene_id"])
                elif exonMatchedGTF3.at[index, "solitaryGene"] == False:
                    genesCoveredGN2.append(exonMatchedGTF3.at[index, "gene_id"])
        
        #Find partial matches that don't have perfect matches
        genesCoveredGN3 = [x for x in genesCoveredGN2 if x not in genesCoveredGN]
        genesCoveredSG3 = [x for x in genesCoveredSG2 if x not in genesCoveredSG]
        
        #Get the gene count lengths
        genesCoveredSGFinal2 = list(set(genesCoveredSG3))
        genesCoveredSGNumBoth = len(genesCoveredSGFinal2)
        
        genesCoveredGNFinal2 = list(set(genesCoveredGN3))
        genesCoveredGNNumBoth = len(genesCoveredGNFinal2)
        
        #Get amount of gtf that overlapepd a known transcript but did not match anyhting
        transcriptomeOverlap = exonMatchedGTF3[exonMatchedGTF3['Overlap_Ref_Transcript'].map(len) > 0]
        transcriptomeOverlap = transcriptomeOverlap[transcriptomeOverlap['Reference_Transcript_Partial'].map(len) == 0]
        transcriptomeOverlap = transcriptomeOverlap[transcriptomeOverlap['Matched_Reference_Transcript'].map(len) == 0]
        transcriptomeOverlapNum = len(transcriptomeOverlap.index)
        
        #Define numbers for perfect matches, partial matches and no matches
        matchedPerf = fullRefSize-nonReferenceSize
        matchedPart = fullRefSize-nonReferencePartSize-doubleMatchSize
        noMatch = fullRefSize - matchedPerf - matchedPart - transcriptomeOverlapNum
        
        #Print a table of the  respective transcript sizes
        t = PrettyTable(['Name', 'Amount', 'Percent'])
        t.add_row(['Full Size', fullRefSize, 100])
        t.add_row(['Matched To Reference (Perfect)', matchedPerf, round(((matchedPerf)/fullRefSize)*100, 3)])
        t.add_row(['Matched To Reference (Partial (No Perfect))', matchedPart, round(((matchedPart)/fullRefSize)*100, 3)])
        t.add_row(['Overlap With Known Transcript, No SJ Match (Full or Partial)', transcriptomeOverlapNum, round(((transcriptomeOverlapNum)/fullRefSize)*100, 3)])
        t.add_row(['No Match (Perfect or Partial) & No Overlap', noMatch, round(((noMatch)/fullRefSize)*100, 3)])
        print("\n")
        print(t)
        
        ##Get the size of overlapping solitary genes and of non matching solitary genes
        solitaryGTF = exonMatchedGTF3[exonMatchedGTF3["solitaryGene"] == True]
        solitaryGTFNM = solitaryGTF[solitaryGTF.astype(str)['Matched_Reference_Gene'] == '[]']
        solitaryGTFNM = solitaryGTFNM[solitaryGTFNM.astype(str)['Reference_Transcript_Partial'] == '[]']
        solitaryGTFOverlap = solitaryGTFNM[solitaryGTFNM.astype(str)['Overlap_Ref_Transcript'] != '[]']
        solitaryGTFOverlapNUM = solitaryGTFOverlap.shape[0]
        
        solitaryGTFNM = solitaryGTFNM[solitaryGTFNM.astype(str)['Overlap_Ref_Transcript'] == '[]']
        solitaryGTFNMNum = solitaryGTFNM.shape[0]
        
        #Print Solitary Gene sizes
        t = PrettyTable(['Name', 'Amount', "Percent"])
        t.add_row(['Total Solitary Genes', solitaryGTF.shape[0], 100])
        t.add_row(['Solitary Genes With Perfect Match To Reference', genesCoveredSGNum, round(((genesCoveredSGNum)/solitaryGTF.shape[0])*100, 3)])
        t.add_row(['Solitary Genes With Partial Match To Reference (No Perfect)', genesCoveredSGNumBoth, round(((genesCoveredSGNumBoth)/solitaryGTF.shape[0])*100, 3)])
        t.add_row(['Solitary Genes With Overlap, But No Match To Reference', solitaryGTFOverlapNUM, round(((solitaryGTFOverlapNUM)/solitaryGTF.shape[0])*100, 3)])
        t.add_row(['Solitary Genes With No Match To Reference & No Overlap', solitaryGTFNMNum, round(((solitaryGTFNMNum)/solitaryGTF.shape[0])*100, 3)])
        print("\n")
        print(t)
        
        #Get Number of Gene Networks 
        networkGTF = exonMatchedGTF3[exonMatchedGTF3["solitaryGene"] == False]
        networkGTFNum = len(networkGTF.gene_id.unique())
        
        #Get Number of Overlapping Genes
        networkGTFNM = networkGTF[networkGTF.astype(str)['Matched_Reference_Gene'] == '[]']
        networkGTFNM = networkGTFNM[networkGTFNM.astype(str)['Reference_Transcript_Partial'] == '[]']
        networkGTFOverlap = networkGTFNM[networkGTFNM.astype(str)['Overlap_Ref_Transcript'] != '[]']
        networkGTFOverlapList = networkGTFOverlap.gene_id.unique()
        
        #Make Sure overlapping genese werent already included in any matches
        networkGTFOverlapList = [x for x in networkGTFOverlapList if x not in genesCoveredGN]
        networkGTFOverlapList = [x for x in networkGTFOverlapList if x not in genesCoveredGN3]
        
        networkGTFOverlapNum = len(networkGTFOverlapList)
        
        #Get Number of Non-Matches
        networkGTFNM = networkGTFNM[networkGTFNM.astype(str)['Overlap_Ref_Transcript'] == '[]']
        networkGTFNMList = networkGTFNM.gene_id.unique()
        
        networkGTFNMList = [x for x in networkGTFNMList if x not in genesCoveredGN]
        networkGTFNMList = [x for x in networkGTFNMList if x not in genesCoveredGN3]
        networkGTFNMList = [x for x in networkGTFNMList if x not in networkGTFOverlapList]
        
        networkGTFNMNum = len(networkGTFNMList)
        
        #Print Gene Network sizes
        t = PrettyTable(['Name', 'Amount', "Percent"])
        t.add_row(['Total Gene Networks', networkGTFNum, 100])
        t.add_row(['Gene Networks With Perfect Match To Reference', genesCoveredGNNum, round(((genesCoveredGNNum)/networkGTFNum)*100, 3)])
        t.add_row(['Gene Networks With Partial Match To Reference (No Perfect)', genesCoveredGNNumBoth, round(((genesCoveredGNNumBoth)/networkGTFNum)*100, 3)])
        t.add_row(['Gene Networks With Overlap, But No Match To Reference', networkGTFOverlapNum, round(((networkGTFOverlapNum)/networkGTFNum)*100, 3)])
        t.add_row(['Gene Networks With No Match To Reference & No Overlap', networkGTFNMNum, round(((networkGTFNMNum)/networkGTFNum)*100, 3)])
        print("\n")
        print(t)
        
        #Create each set for our venn diagram dictionary
        k = 0
        vennDict = {}
        vennDict2 = {}
        
        analyzedLabels2 = analyzedLabels
        analyzedLabels2.append(referenceLabel)
        
        geneLists = [ [] for _ in range(len(analyzedLabels2)) ]
        
        listItems = []

        for index, row in exonMatchedGTF3.iterrows():
            for item in exonMatchedGTF3.at[index, "Located_In"]:
                if vennDict.get(item) is None:
                    vennDict[item] = {exonMatchedGTF3.at[index, "gene_id"]}
                    vennDict2[item] = {exonMatchedGTF3.at[index, "transcript_id"]}
                    if item not in listItems:
                        listItems.append(item)
                else:
                    vennDict[item].add(exonMatchedGTF3.at[index, "gene_id"])
                    vennDict2[item].add(exonMatchedGTF3.at[index, "transcript_id"])
        
        if oneExonSeperate == True:
            print("Multi Exon Transcript Results:")
            print("\n")
        else:
            print("All Transcript Results:")
            print("\n")
        
        if vennReference == False:
            del vennDict[referenceLabel]
            del vennDict2[referenceLabel]
            
        r = dict(vennDict)
        r2 = dict(vennDict2)
          
        if len(r2) == 6:
            print("Gene Overlap:")
            pseudovenn(r)
            plt.show()
            print("Transcript Overlap:")
            pseudovenn(r2)
            plt.show()
        elif len(r2) <= 6:
            print("Gene Overlap:")
            venn(r)
            plt.show()
            print("Transcript Overlap:")
            venn(r2)
            plt.show()
        else:
            print("More than 6 GTFs combined, result will not be displayed as a venn diagram.")
            print("\n")
        
        #Create our PYU plot
        pyuSeriesM = from_contents(vennDict)
        pyuSeries2M = from_contents(vennDict2)
        
        print("Gene Overlap:")
        plot(pyuSeriesM, sort_by = "cardinality", show_counts = True)
        plt.savefig('plots/multiGene.png', bbox_inches='tight')
        plt.show()
        
        print("Transcript Overlap:")
        plot(pyuSeries2M, sort_by = "cardinality", show_counts = True)
        plt.savefig('plots/multiTranscript.png', bbox_inches='tight')
        plt.show()
        
        if oneExonSeperate == True:
            
            print("Beginning Single Exon Analysis:")
            
             #Full size of our final GTF
            fullRefSize = len(oneExonMatchedGTF3.index)

            #Get size of gtf that didnt match anything perfectly in the reference
            nonReference = oneExonMatchedGTF3[oneExonMatchedGTF3.astype(str)['Matched_Reference_Gene'] == '[]']
            nonReferenceSize = len(nonReference.index)


            #get size of the GTf that matched something paritally and perfectly
            doubleMatch = oneExonMatchedGTF3[oneExonMatchedGTF3.astype(str)['Reference_Transcript_Partial'] != '[]']
            doubleMatch = doubleMatch[doubleMatch.astype(str)['Matched_Reference_Gene'] != '[]']
            doubleMatchSize = len(doubleMatch.index)



            ##Find Gene Coverage
            genesCovered = []
            for index, row in oneExonMatchedGTF3.iterrows():
                if oneExonMatchedGTF3.at[index, "Matched_Reference_Gene"]:
                    genesCovered.append(oneExonMatchedGTF3.at[index, "gene_id"])

            genesCoveredFinal = list(set(genesCovered))
            genesCoveredNum = len(genesCoveredFinal)
            
            genesCovered2 = []
            for index, row in oneExonMatchedGTF3.iterrows():
                if oneExonMatchedGTF3.at[index, "Reference_Transcript_Partial"]:              
                    genesCovered2.append(oneExonMatchedGTF3.at[index, "gene_id"])

            
            genesCovered2 = [x for x in genesCovered2 if x not in genesCovered]
            
            genesCoveredFinal2 = list(set(genesCovered2))
            genesCoveredNumBoth = len(genesCoveredFinal2)
            
            #Get amount of gtf that overlapepd a known transcript but did not match anyhting
            transcriptomeOverlap = oneExonMatchedGTF3[oneExonMatchedGTF3['Overlap_Ref_Transcript'].map(len) > 0]
            transcriptomeOverlap = transcriptomeOverlap[transcriptomeOverlap['Reference_Transcript_Partial'].map(len) == 0]
            transcriptomeOverlap = transcriptomeOverlap[transcriptomeOverlap['Matched_Reference_Transcript'].map(len) == 0]
            transcriptomeOverlapNum = len(transcriptomeOverlap.index)

            
            matchedPerf = fullRefSize-nonReferenceSize
            noMatch = fullRefSize - matchedPerf - transcriptomeOverlapNum

            #Print a table of the  respective transcript sizes
            t = PrettyTable(['Name', 'Amount', 'Percent'])
            t.add_row(['Full Size', fullRefSize, 100])
            t.add_row(['Perfect Match To Reference', matchedPerf, round(((matchedPerf)/fullRefSize)*100, 3)])
            t.add_row(['Overlap With Known Transcript, No Match', transcriptomeOverlapNum, round(((transcriptomeOverlapNum)/fullRefSize)*100, 3)])
            t.add_row(['No Match & No Overlap', noMatch, round(((noMatch)/fullRefSize)*100, 3)])
            print("\n")
            print(t)
            
            
            ##Get the size of overlapping solitary genes and of non matching solitary genes
            solitaryGTF = oneExonMatchedGTF3[oneExonMatchedGTF3["solitaryGene"] == True]
            solitaryGTFNM = solitaryGTF[solitaryGTF.astype(str)['Matched_Reference_Gene'] == '[]']
            solitaryGTFNM = solitaryGTFNM[solitaryGTFNM.astype(str)['Reference_Transcript_Partial'] == '[]']
            solitaryGTFOverlap = solitaryGTFNM[solitaryGTFNM.astype(str)['Overlap_Ref_Transcript'] != '[]']
            solitaryGTFOverlapNUM = solitaryGTFOverlap.shape[0]

            solitaryGTFNM = solitaryGTFNM[solitaryGTFNM.astype(str)['Overlap_Ref_Transcript'] == '[]']
            solitaryGTFNMNum = solitaryGTFNM.shape[0]

            #Print Solitary Gene sizes
            t = PrettyTable(['Name', 'Amount', "Percent"])
            t.add_row(['Total Solitary Genes', solitaryGTF.shape[0], 100])
            t.add_row(['Solitary Genes With Perfect Match To Reference', genesCoveredNum, round(((genesCoveredNum)/solitaryGTF.shape[0])*100, 3)])
            t.add_row(['Solitary Genes With Overlap, But No Match To Reference', solitaryGTFOverlapNUM, round(((solitaryGTFOverlapNUM)/solitaryGTF.shape[0])*100, 3)])
            t.add_row(['Solitary Genes With No Match To Reference & No Overlap', solitaryGTFNMNum, round(((solitaryGTFNMNum)/solitaryGTF.shape[0])*100, 3)])
            print("\n")
            print(t)


            #Create each set for our venn diagram dictionary
            k = 0
            vennDict = {}
            vennDict2 = {}

            analyzedLabels2 = analyzedLabels
            analyzedLabels2.append(referenceLabel)

            geneLists = [ [] for _ in range(len(analyzedLabels2)) ]

            listItems = []

            for index, row in oneExonMatchedGTF3.iterrows():
                for item in oneExonMatchedGTF3.at[index, "Located_In"]:
                    if vennDict.get(item) is None:
                        vennDict[item] = {oneExonMatchedGTF3.at[index, "gene_id"]}
                        vennDict2[item] = {oneExonMatchedGTF3.at[index, "transcript_id"]}
                        if item not in listItems:
                            listItems.append(item)
                    else:
                        vennDict[item].add(oneExonMatchedGTF3.at[index, "gene_id"])
                        vennDict2[item].add(oneExonMatchedGTF3.at[index, "transcript_id"])

            if vennReference == False:
                del vennDict[referenceLabel]
                del vennDict2[referenceLabel]
                
            print("Single Exon Transcript Results:")
            print("\n")

            r = dict(vennDict)
            r2 = dict(vennDict2)

            if len(r2) == 6:
                print("Gene Overlap:")
                pseudovenn(r)
                plt.show()
                print("Transcript Overlap:")
                pseudovenn(r2)
                plt.show()
            elif len(r2) <= 6:
                print("Gene Overlap:")
                venn(r)
                plt.show()
                print("Transcript Overlap:")
                venn(r2)
                plt.show()
            else:
                print("More than 6 GTFs combined, result will not be displayed as a venn diagram.")
                print("\n")

            #Create our PYU plot
            pyuSeriesS = from_contents(vennDict)
            pyuSeries2S = from_contents(vennDict2)

            print("Gene Overlap:")
            plot(pyuSeriesS, sort_by = "cardinality", show_counts = True)
            plt.savefig('plots/singleGene.png', bbox_inches='tight')
            plt.show()

            print("Transcript Overlap:")
            plot(pyuSeries2S, sort_by = "cardinality", show_counts = True)
            plt.savefig('plots/singleTranscript.png', bbox_inches='tight')
            plt.show()
        
        return(exonMatchedGTF4)
    else:
        
        exonMatchedGTF4 = exonMatchedGTF2.copy
        
        if oneExonSeperate == True:
            oneExonMatchedGTF2 = exonMatchedGTF2[(exonMatchedGTF2.Splice_Junctions.str.len() == 0)]
            exonMatchedGTF2 = exonMatchedGTF2[(exonMatchedGTF2.Splice_Junctions.str.len() != 0)]
            print("\n")
            print("Beginning Multi Exon Analysis:")
        
        #Create each set for our venn diagram dictionary
        k = 0
        vennDict = {}
        vennDict2 = {}

        analyzedLabels2 = analyzedLabels
        analyzedLabels2.append(referenceLabel)

        geneLists = [ [] for _ in range(len(analyzedLabels2)) ]

        listItems = []

        for index, row in exonMatchedGTF2.iterrows():
            for item in exonMatchedGTF2.at[index, "Located_In"]:
                if vennDict.get(item) is None:
                    vennDict[item] = {exonMatchedGTF2.at[index, "gene_id"]}
                    vennDict2[item] = {exonMatchedGTF2.at[index, "transcript_id"]}
                    if item not in listItems:
                        listItems.append(item)
                else:
                    vennDict[item].add(exonMatchedGTF2.at[index, "gene_id"])
                    vennDict2[item].add(exonMatchedGTF2.at[index, "transcript_id"])
                    
                    
        if vennReference == False:
            del vennDict[referenceLabel]
            del vennDict2[referenceLabel]

        if len(analyzedLabels) == 6:
            print("Gene Overlap:")
            pseudovenn(vennDict)
            plt.show()
            print("Transcript Overlap:")
            pseudovenn(vennDict2)
            plt.show()
        elif len(analyzedLabels) <= 6:
            print("Gene Overlap:")
            venn(vennDict)
            plt.show()
            print("Transcript Overlap:")
            venn(vennDict2)
            plt.show()
        else:
            print("More than 6 GTFs combined, result will not be displayed as a venn diagram.")
            print("\n")

        #Create our PYU plot
        pyuSeries = from_contents(vennDict2)
        plot(pyuSeries)
        plt.show()            
        
        if oneExonSeperate == True:
            print("Beginning Single Exon Analysis:")
            k = 0
            vennDict = {}
            vennDict2 = {}

            analyzedLabels2 = analyzedLabels
            analyzedLabels2.append(referenceLabel)

            geneLists = [ [] for _ in range(len(analyzedLabels2)) ]

            listItems = []

            for index, row in oneExonMatchedGTF2.iterrows():
                for item in oneExonMatchedGTF2.at[index, "Located_In"]:
                    if vennDict.get(item) is None:
                        vennDict[item] = {oneExonMatchedGTF2.at[index, "gene_id"]}
                        vennDict2[item] = {oneExonMatchedGTF2.at[index, "transcript_id"]}
                        if item not in listItems:
                            listItems.append(item)
                    else:
                        vennDict[item].add(oneExonMatchedGTF2.at[index, "gene_id"])
                        vennDict2[item].add(oneExonMatchedGTF2.at[index, "transcript_id"])
            
            if vennReference == False:
                del vennDict[referenceLabel]
                del vennDict2[referenceLabel]

            if len(analyzedLabels) == 6:
                print("Gene Overlap:")
                pseudovenn(vennDict)
                plt.show()
                print("Transcript Overlap:")
                pseudovenn(vennDict2)
                plt.show()
            elif len(analyzedLabels) <= 6:
                print("Gene Overlap:")
                venn(vennDict)
                plt.show()
                print("Transcript Overlap:")
                venn(vennDict2)
                plt.show()
            else:
                print("More than 6 GTFs combined, result will not be displayed as a venn diagram.")
                print("\n")

            #Create our PYU plot
            pyuSeries = from_contents(vennDict2)
            plot(pyuSeries)
            plt.show()
            
    
    return(exonMatchedGTF4)
    
    
