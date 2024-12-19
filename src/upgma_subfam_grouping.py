#!/usr/bin/env python3

import numpy as np
from datetime import datetime
import sys

multiprocessing = None
if 'linux' in sys.platform:
    multiprocessing = 'multiprocessing'
elif 'darwin' in sys.platform:
    multiprocessing = 'multiprocess'
elif 'win' in sys.platform:
    multiprocessing = 'multiprocessing'
mp = __import__(multiprocessing) #import multiprocessing as mp

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
    # print("empty matrix done!")
    end_time = datetime.now()

    print(end_time-start_time)
    return seqMatrix
    # print(seqDict)

def differ(output):
    # input()

    row = [-1 for i in range(len(dictionary))]

    for index, key in enumerate(dictionary.keys()):
        if key == output:
            i = index
            # print(output)
            # print(i)
    # return(i)
    # input()
    
    for j, value2 in enumerate(dictionary.values()):
        # print(value2)
        if i<j :
            difference = 0
            for k in range(len(value2)):
                # print(dictionary[output][k])
                # input()
                if value2[k] != dictionary[output][k]:
                    difference += 1
            # print("this is difference")
            # print(difference)
            # input()
            # print("this is ij")
            # print(i, j)
            # print("this is  ij")
             
            row[j] = difference
            
            # print("this is matrix update")
            # print(row)


            # print(f"{i}and {j} Done!")
    return row
        
def worker_init(seqdict, seqmatrix):
    global dictionary
    dictionary = seqdict
    global matrix
    matrix = seqmatrix
    # print(dictionary)

def grouping(seqDict, seqMatrix):
    key_list = []
    for i, keys in enumerate(seqDict.keys()):
        key_list.append(keys)
    original_seqM = seqMatrix.copy()
    location = []
    for i in range(0,len(original_seqM)):
        location.append([i])

    with open("grouping.txt", "w") as file:
        first_line = "## " + ";".join(key_list) + "\n"
        file.write(first_line)

    #newick dictionary
    # newick = {}

    tree = {
        0 : key_list.copy()
    }

    start_time = datetime.now()
    while(len(seqMatrix)>2):
        # print(len(seqMatrix))
        row = 0
        m = len(seqMatrix)
        index = np.argmin(seqMatrix[seqMatrix>-1])
        while(index-(m-row-1)>=0): 
            index -= (m-row)-1
            # print(index)
            row += 1
        column = 1
        while(index-1>=0):
            index -= 1
            column += 1
        column += row
    
        #make newick format
        key1 = key_list[row] #name in the (row)th value of key list
        key2 = key_list[column] # name in the (column)th value of the key list
        length = seqMatrix[row, column]/2
        
        # new_newick = "{},{}".format(key1, key2) #name of the key
        # newick[new_newick] = {key1: length,
        #                       key2: length}    

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
        # print(teao_list)
        tree[length] = teao_list

#make new array and update it
        arr2 = seqMatrix
        arr2 = np.delete(arr2, column, axis = 0)
        arr2 = np.delete(arr2, column, axis = 1)

        arr2 = np.delete(arr2, row, axis = 0)
        arr2 = np.delete(arr2, row, axis = 1)
        
        arr2 = np.insert(arr2, 0, -1, axis = 1) #(arr, obj, values, axis)
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
                        # sum1 += original_seqM[k,m]
                        sum1 += original_seqM[(location[l])[k], (location[0])[m]]
                    else:
                        # sum1 += original_seqM[m,k]
                        sum1 += original_seqM[(location[0])[m],(location[l])[k]]
                        # print((location[0])[k],(location[l])[m])
                    arr2[0,l] = sum1/((len(location[0]))*(len(location[l])))

        seqMatrix = arr2
        # print(seqMatrix)
    end_time = datetime.now()
    print(end_time-start_time)
    
    # final newick
    row = 0
    column = 1
    tree[seqMatrix[row,column]/2] = [tree[0]]

    # level = 0
    tree_red = {}
    tree_content = {}

    for j, value in enumerate(tree.values()):
        # print(j)
        # print(value)
        for i in range(len(value)):
            # print(value[i])
            # print(type(value[i]))
            if type(value[i])==list:
                for k in range(len(value[i])):
                    tree_content[value[i][k]] = hex(i)
            else:
                tree_content[value[i]] = hex(i)
        # print(j, tree_content)
        content = tree_content.copy()
        tree_red[j] = content
        grouping_list = list(content.values())
        grouping_line = ";".join(grouping_list)
        with open("grouping.txt", 'a') as file:
            file.write(grouping_line+"\n")
    # print(content)
    # print(tree_red)
    

if __name__ == "__main__":
    input_file = input("file name: ")
    seqDict = dictreturn(input_file)
    seqMatrix = seqMreturn()
    with mp.Pool(initializer= worker_init, initargs=(seqDict,seqMatrix), processes = mp.cpu_count()) as executor:
        results = executor.map(differ, seqDict)
    # print("this is result")
    # print(results)
    # print(type(results))
    matrix = np.array(results)
    # print(matrix)
    # print(seqDict)

    grouping(seqDict, matrix)