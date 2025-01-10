#!/usr/bin/env python3

import os
import sys
import argparse
import numpy as np
from datetime import datetime

multiprocessing = None
if 'linux' in sys.platform:
    multiprocessing = 'multiprocessing'
elif 'darwin' in sys.platform:
    multiprocessing = 'multiprocess'
elif 'win' in sys.platform:
    multiprocessing = 'multiprocessing'
mp = __import__(multiprocessing) # Import multiprocessing as mp

def dictreturn(inputMSA):
    global seqDict
    seqDict = {}
    with open (inputMSA, 'r') as f:
        for i, line in enumerate(f):
            newLine = line.strip()
            if newLine == '':
                continue
            if i >=3 and line[0]!=(' '): 
                newLine = newLine.split()
                if newLine[0] in seqDict.keys():
                    seqDict[newLine[0]]+=newLine[1]
                else:
                    seqDict[newLine[0]] = newLine[1]

    seqDict = dict(sorted(seqDict.items(), key=lambda x: x[0]))
    return seqDict

def seqMreturn():
    global seqMatrix
    seqMatrix = np.zeros([len(seqDict),len(seqDict)])

    start_time = datetime.now()
    for i in range(0, len(seqMatrix)):
        for j in range(0, len(seqMatrix)):
            seqMatrix[i,j] = -1
    end_time = datetime.now()

    print(end_time-start_time)
    return seqMatrix

def differ(output):
    row = [-1 for i in range(len(dictionary))]
    for index, key in enumerate(dictionary.keys()):
        if key == output:
            i = index

    for j, value2 in enumerate(dictionary.values()):
        if i<j :
            difference = 0
            for k in range(len(value2)):
                if value2[k] != dictionary[output][k]:
                    difference += 1
        
            row[j] = difference
    return row
        
def worker_init(seqdict, seqmatrix):
    global dictionary
    dictionary = seqdict
    global matrix
    matrix = seqmatrix

def grouping(seqDict, seqMatrix, outfile):
    key_list = []
    for i, keys in enumerate(seqDict.keys()):
        key_list.append(keys)
    original_seqM = seqMatrix.copy()
    location = []
    for i in range(0,len(original_seqM)):
        location.append([i])

    with open(outfile, "w") as file:
        first_line = "## " + ";".join(key_list) + "\n"
        file.write(first_line)

    tree = {
        0 : key_list.copy()
    }

    start_time = datetime.now()
    while(len(seqMatrix)>2):
        row = 0
        m = len(seqMatrix)
        index = np.argmin(seqMatrix[seqMatrix>-1])
        while(index-(m-row-1)>=0): 
            index -= (m-row)-1
            row += 1
        column = 1
        while(index-1>=0):
            index -= 1
            column += 1
        column += row
    
        #make newick format
        key1 = key_list[row] # Name in the (row)th value of key list
        key2 = key_list[column] # Name in the (column)th value of the key list
        length = seqMatrix[row, column]/2
        
        new_key = "{},{}".format(key1, key2)

        del key_list[column]
        del key_list[row]

        key_list.insert(0,new_key)
        '''
        this is for teao
        '''
        teao_list = key_list.copy()
        for i in range(len(teao_list)):
            if ',' in teao_list[i]:
                teao_list[i] = teao_list[i].split(',')
        tree[length] = teao_list

        # Make new array and update it
        arr2 = seqMatrix
        arr2 = np.delete(arr2, column, axis = 0)
        arr2 = np.delete(arr2, column, axis = 1)

        arr2 = np.delete(arr2, row, axis = 0)
        arr2 = np.delete(arr2, row, axis = 1)
        
        arr2 = np.insert(arr2, 0, -1, axis = 1) # (arr, obj, values, axis)
        arr2 = np.insert(arr2, 0, -1, axis = 0)


        new_name = location[row]+location[column]
        location.insert(0,new_name)

        del location[column+1]
        del location[row+1]

        for l in range(1, len(arr2)):
            sum1 = 0
            for m in range(0, len(location[0])): # [abc]
                for k in range(0, len(location[l])): # [de],[f],[g]
                    if (location[0])[m] > (location[l])[k]:
                        sum1 += original_seqM[(location[l])[k], (location[0])[m]]
                    else:
                        sum1 += original_seqM[(location[0])[m],(location[l])[k]]
                    arr2[0,l] = sum1/((len(location[0]))*(len(location[l])))

        seqMatrix = arr2
    end_time = datetime.now()
    print(end_time-start_time)
    
    # Final newick
    row = 0
    column = 1
    tree[seqMatrix[row,column]/2] = [tree[0]]

    # level = 0
    tree_red = {}
    tree_content = {}

    for j, value in enumerate(tree.values()):
        for i in range(len(value)):
            if type(value[i])==list:
                # print(value[i])
                for k in range(len(value[i])):
                    # print(value[i][k])
                    if type(value[i][k])==str:
                        tree_content[value[i][k]] = hex(i)
                    elif type(value[i][k])==list:  # Duplicate sequences
                        print(value[i][k])
                        sys.exit()
                        tree_content[str(value[i][k])] = hex(i)
            else:
                # print(value[i])
                tree_content[value[i]] = hex(i)
        content = tree_content.copy()
        tree_red[j] = content
        grouping_list = list(content.values())
        grouping_line = ";".join(grouping_list)
        with open(outfile, 'a') as file:
            file.write(grouping_line+"\n")

def check_input(inputMSA):
    """Check whether the input MSA is in CLUSTAL format"""
    extension = os.path.splitext(inputMSA)[1]
    if extension != ".clustal":
        print("Error: Input MSA is not in CLUSTAL format")
        sys.exit(1)
    with open(inputMSA, 'r') as f:
        first_line = f.readline()
        if "CLUSTAL" not in first_line:
            print("Error: Input MSA is not in CLUSTAL format")
            sys.exit(1)    

def main(args):

    check_input(args.msa)
    seqDict = dictreturn(args.msa)

    seqMatrix = seqMreturn()

    with mp.Pool(initializer= worker_init, initargs=(seqDict,seqMatrix), processes = mp.cpu_count()) as executor:
        results = executor.map(differ, seqDict)
    matrix = np.array(results)

    # Create output filename
    basename = os.path.splitext(args.msa)[0]
    outname = os.path.join(basename + "_grouping.txt")

    grouping(seqDict, matrix, outname)
    print("Grouping subfamilies based on UPGMA clustering is done. Output file is saved as:", outname)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="""Grouping subfamilies based on UPGMA clustering. Input is a Multiple Sequence Alignment (CLUSTAL format)
                                                    and the output is a text file containing the subfamily grouping. The output file is named as 
                                                    <input_filename>_grouping.txt""")
    parser.add_argument("-m", "--msa", required=True, help="file location of Multiple Sequence Alignment (CLUSTAL format)")
    args = parser.parse_args()

    main(args)