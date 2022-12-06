#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#__authors__ = ("Elian Strozyk")
#__contact__ = ("elian.strozyk@etu.umontpellier.fr")
#__version__ = "0.0.1"
#__date__ = "12/06/2022"
#__licence__ ="This program is free software: you can redistribute it and/or modify
#        it under the terms of the GNU General Public License as published by
#        the Free Software Foundation, either version 3 of the License, or
#        (at your option) any later version.
#        This program is distributed in the hope that it will be useful,
#        but WITHOUT ANY WARRANTY; without even the implied warranty of
#        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#        GNU General Public License for more details.
#        You should have received a copy of the GNU General Public License
#        along with this program. If not, see <https://www.gnu.org/licenses/>."

import sys, re

# input: .sam file
# output: .txt file summing up number of different reads

######################## File's quality checking ########################

# Thorough checking of the file's quality
# What's checked:
# file ends in .sam
## Header:
### The very first character is an '@'
## Body:
### File contains at least 11 sections split with \t so 10 \t in total

def nameEndsInSAM(file) :
    # Check whether the file ends in '.sam'
    print("Checking file's extension:", end=' ')
    if file.endswith(".sam") :
        print("done")
        return True
    else : 
        print("The file submitted does not have the correct extension.\nMake sure the file ends with '.sam'.")
        return False

def StartsWithAnAtSymbol(file) :
    # Check whether the very first character is an '@'
    print("Checking file's header:", end=' ')
    if file[0] == "@" :
        print("done")
        return True
    else :
        print("The file submitted does not contain a header whose first line starts with an '@'.")
        return False

def tabNumberChecking(file) :
    # Counting number of \t (at least 10) in the very first line of the file's body
    print("Checking column number:", end=' ')
    tab_counter = 0

    for character in file[0] :
        if character == "\t" :
            tab_counter += 1

    if tab_counter > 10 :
        print("done")
        return True
    else :
        print("The file does not contain enough columns. It must have at least 11 columns.")
        return False

######################## Storing up file's body content ########################

# After checking file's quality:
## file's body is extracted
## file's data are stored in a dictionary

def ExtractsFileBody(file) :
    # extract file's body while removing lines starting with '@'
    print("Extracting file's body:", end=' ')
    file = file.split("\n") # listification after each '\n'
    file_body = []

    for line in file :
        x = re.findall("^@", line)
        if not x :
            file_body.append(line)
    print("done")
    return file_body

def dictionaryCreation(file) :
    # Creating dictionary while taking care of duplicate reads
    print("Creating dictionary:", end=' ')
    reads = []
    readData = []
    dictionary = {}
    
    for line in file :
        reads.append(line.split("\t")[0])
        readData.append(line.split("\t")[1:6])

    # Renaming reads if one's name is similar. 
    # /!\ This method does not work if more than two reads with the same name
    ## If so, read names are redundant, still
    readsUnique = []
    for i in range(0, len(reads) - 1) :
        if reads[i] == reads[i + 1] :
            readsUnique.append(reads[i])
        else :
            readsUnique.append(reads[i] + ";2")

    # dictionary creation
    dictionary = dict(zip(readsUnique, readData))
    print("done")
    return dictionary

######################## Analysis of file's content ########################

# Extracts and translates FLAG
# Filters mapping quality of reads
# Filters reads whose CIGAR are matched at 100%

# The following function does work BUT
## It should have returned a list of binaries and not a single binary
## Either way, probably unoptimised
###########def flagToBinary(flag) :
###########    # Translates a FLAG to a decimal
###########    # flag checking
###########    if flag < 1 or flag > 2048 :
###########        return 0
###########    
###########    # creation of an ordered list containing all powers of 2 from 1 to 2048
###########    binary_list = []
###########    for i in range(0, 12) :
###########        binary_list.append(pow(2, i))
###########
###########    # calculation of the number to substract 
###########    # frame the flag (find the maximal_min and minimal_max within the list)
###########    for a in range(0, len(binary_list)) :
###########        if ((binary_list[a] < flag) and (binary_list[a + 1] > flag)) or binary_list[a] == flag :
###########            numberToSubstract = binary_list[a]
###########    
###########    # recursive substraction
###########    listBinary = []
###########    if flag - numberToSubstract != 0 :
###########        return flagToBinary(flag - numberToSubstract)
###########    else :
###########        return flag

def flagAnalysis(dictionary, decimal) :
    # decimal == 1 : reads paired
    # decimal == 2 : reads mapped in proper pair
    ## for decimal == 2, should also check distance between reads < 150
    counter = 0
    for key in dictionary :
        if decimal == 1 :
            if int(dictionary[key][0]) & 1 == 1 :
                counter += 1
        if decimal == 2 :
            if int(dictionary[key][0]) & 2 == 2 :
                        counter += 1
            
            #for i in dictionary[i] :
            #        if abs(dictionary[i][2] - dictionary[i+1][2]) < 150 :
            #           counter += 1
    
    # '& decimal == decimal': never would I have found it without Arnaud Soulier
    # and my function flagToBinary() does not work correctly
    
    if decimal == 1 :
        return f'Number of paired reads: {counter}'
    elif decimal == 2 : # and distance <= 150
        return f'Number of reads mapped in proper pair: {counter}'

def filterMappingQuality(dictionary) :
    # filters mapping quality depending on an input threshold (or not)
    counter = 0
    print("\tFiltering reads", end = ' ')
    qualityThreshold = input("\n\t\tInsert a mapping quality threshold between 0% and 100%. Base threshold is 20%.")
    
    # threshold's checking not monkeyproof that much
    if len(qualityThreshold) == 0 : #or not (int(qualityThreshold) >= 0 and int(qualityThreshold) >= 100) :
        qualityThreshold = 20
    else : 
        qualityThreshold = int(qualityThreshold)
    
    for key in dictionary :
        if int(dictionary[key][3]) >= qualityThreshold :
            counter += 1
    print("done")
    return f'Number of reads with a mapping quality greater than or equals to {qualityThreshold}%: {counter}'

def CIGARanalysis(dictionary) :
    # filters out reads whose CIGAR is not matched at 100%
    counter = 0
    print("\tCounting number of totally aligned reads", end=' ')
    for key in dictionary :
        if re.findall("^[0-9][0-9][0-9]M$", dictionary[key][4]) or re.findall("^[0-9][0-9]M$", dictionary[key][4]) :
            counter += 1
    print("done")
    return f'Number of totally aligned reads: {counter}'


######################## Output ########################

def makeOutput(dictionary) :
    print("Making summary file...")

    # see dictionaryCreation(): does not work when reads are present more than twice
    counter = 0
    for key in dictionary :
        if not re.findall(";2$", key) :
            counter += 1

    # Writing the new file
    with open(f'{sys.argv[1]}_summary.txt', "w") as output :
        output.write(f'Total number of reads: {len(dictionary)}\n')
        output.write(f'Number of different reads: {counter}\n')
        output.write(f'{flagAnalysis(dictionary, 1)}\n')
        output.write(f'{flagAnalysis(dictionary, 2)}\n')
        output.write(f'{filterMappingQuality(dictionary)}\n')
        output.write(f'{CIGARanalysis(dictionary)}')
    print(f"The output is located in ./{sys.argv[1]}_summary.txt")


######################## main() ########################
def main(file) :
    if nameEndsInSAM(file) :
        with open(file, "r") as fd :
            file_sam = fd.read()
            if StartsWithAnAtSymbol(file_sam) :
                file_body = ExtractsFileBody(file_sam)
                if tabNumberChecking(file_body) :
                    dictionary = dictionaryCreation(file_body)
                    makeOutput(dictionary)
                    

main(sys.argv[1])