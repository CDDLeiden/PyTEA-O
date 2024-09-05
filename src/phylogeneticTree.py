import numpy as np
from datetime import datetime
import os
import sys

multiprocessing = None
if 'linux' in sys.platform:
    multiprocessing = 'multiprocessing'
elif 'darwin' in sys.platform:
    multiprocessing = 'multiprocess'
elif 'win' in sys.platform:
    multiprocessing = 'multiprocessing'
mp = __import__(multiprocessing) #import multiprocessing as mp



def phylo(inputMSA):

#make original matrix
    # start_time = datetime.now()

    # for i, value1 in enumerate(seqDict.values()): #row
    #     for j, value2 in enumerate(seqDict.values()): #column
    #         if i<j: 
    #             difference = 0
    #             for k in range(len(value1)):
    #                 if value1[k] != value2[k]:
    #                     difference += 1
                # seqMatrix[i, j] = difference
        # print(f"{i}and {j} Done!")

    # end_time = datetime.now()
    # print(end_time-start_time)

    


#make list of names
    key_list = []
    for i, keys in enumerate(seqDict.keys()):
        key_list.append(keys)
    original_seqM = seqMatrix.copy()
    location = []
    for i in range(0,len(original_seqM)):
        location.append([i])
    
    #newick dictionary
    newick = {}        

    # newick_format = key_list.copy()

    '''
    This is for tea-o
    if you don't need it, you can delete it
    '''
    tree = {
        0 : key_list.copy()
    }
    '''
    tea-o 
    '''

# make new matrix and update till len>2
    start_time = datetime.now()
    while(len(seqMatrix)>2):
        print(len(seqMatrix))
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
        
        new_newick = "{},{}".format(key1, key2) #name of the key
        newick[new_newick] = {key1: length,
                              key2: length}    

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
        
        # distance1 = (newick[new_newick])[key1]
        # distance2 = (newick[new_newick])[key2]
        # if key1 in newick.keys():
        #     distance1 = (newick[new_newick])[key1] - [x for x in newick[key1].values()][0]

        # if key2 in newick.keys():
        #     distance2 = (newick[new_newick])[key2] - [x for x in newick[key2].values()][0]

        # # input()

        # new_new_newick = "({0}:{1},{2}:{3})".format(newick_format[row], distance1, newick_format[column], distance2)
        # del newick_format[column]
        # del newick_format[row]
        # newick_format.insert(0, new_new_newick)
        

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
    # key1 = key_list[row]
    # key2 = key_list[column]
    
    # if key1 in newick.keys():
    #     final_dist1 = [x for x in newick[key1].values()][0]
    # else:
    #     final_dist1 = 0

    # if key2 in newick.keys():
    #     final_dist2 = [x for x in newick[key2].values()][0]
    # else:
    #     final_dist2 = 0

    # final = "({}:{},{}:{});".format(newick_format[0], seqMatrix[row,column]/2 -final_dist1, newick_format[1] ,seqMatrix[row,column]/2 - final_dist2)


    tree[seqMatrix[row,column]/2] = [tree[0]]

    with open("phylo3.tree", 'w') as file:
        file.write(f"\t{"\t".join([f"{x}" for x in sorted(seqDict.keys())])}\n")
        for i,row in enumerate(original_seqM):
            for j,col in enumerate(row):
                if i == j:
                    # print("here")
                    original_seqM[i,j] = 0
            string = ""
            for x in original_seqM[i]:
                if x != -1:
                    string += f"{x}\t"
                else: 
                    string += "\t"
            file.write(f"{[f"{x}" for x in sorted(seqDict.keys())][i]}\t{string}\n")
    
    # with open("newick.txt", 'w') as file:
    #     file.write(final)
    
    with open("test.txt", 'w') as file:
        for key, value in tree.items():
            file.write(f"{key}: {value}\n")


def differ(arg):
    print("this is differ")
    result = arg
    return result
    
seqDict = {}
inputMSA = 'CS.clustal'
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


seqMatrix = np.zeros([len(seqDict),len(seqDict)])

start_time = datetime.now()
for i in range(0, len(seqMatrix)):
    for j in range(0, len(seqMatrix)):
        seqMatrix[i,j] = -1
print("empty matrix done!")
end_time = datetime.now()

print(end_time-start_time)

if __name__ == "__main__":
    with mp.Pool(processes = 1) as executor:
        inputs = list(seqDict.keys())
        print("yes")
        results = executor.map(differ, inputs)

        print(results)

    input()



phylo('CS.clustal')