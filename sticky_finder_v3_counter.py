#sticky finder counter
#written by Ziqing Ye, modified by Siran Tian @ Trcek Lab @ Johns Hopkins University
#last updated 02/10/2022
#usage: .py fasta_file_randonmized_seqs GC%threshold of the complementary sequences
#output will be transcript name,gene name,sequence length,total # of palindromes,total # of inverted repeats,total number of sense-antisense sequences

import sys
from fasta import FASTAReader

#define RNA class
class RNA(object):
    def __init__(self, fbtr, parent_gene_symbol, transcript_length, sequence_type, sequence): #**
        self.fbtr = fbtr
        self.gene = parent_gene_symbol #nos, pgc, etc.
        self.length = length_tscrip
        self.seq_type = sequence_type #cdna or 3'UTR, etc.
        self.seq = sequence #the string sequence
        #self.structure = structure #the string of binary RNAfold output **


min_len = 6 #defined the minimum length of a complementary sequence is 6-nt for palindrome and sense-antisense, or 6-nt for one arm of the inverted repeats

#define arm class
class Arm(object): #this is a class
    def __init__(self, source_transcript, starting_position, ending_position, sequence, length): #this is a constructor **
        #use: Arm(source_transcript, starting_position, ending_position, sequence, length)
        self.transcript = source_transcript #the source transcript as an RNA class object
        self.start_pos = starting_position #starting nucleotide position of the transcript (nucleotide index starts with 1!)
        self.end_pos = ending_position
        self.seq = sequence #the nucleotide sequence of the arm, a string. e.g. "ATCCG"
        self.length = length #the length of the arm
        #self.structure = structure #binary structure from RNAfold **
        #ending position = starting position + length - 1
        
class ArmPair(object):
    def __init__(self, arm1, arm2):
        self.arm1 = arm1 #class arm
        self.arm2 = arm2 #class arm

def main_function():

    find_irs(rna1, rna2, arm_pairs, min_len) #find ics

    list1=[]
    list2=[]
    list3=[]
    list4=[]

    number_pal=0
    number_irs=0
    number_trans=0

    for ics in arm_pairs:
        arm1 = ics.arm1
        arm2 = ics.arm2

        tuple1=(int(arm1.start_pos),int(arm1.end_pos))
        tuple2=(int(arm2.start_pos),int(arm2.end_pos))
        
        str_1=arm1.seq 

        x=range(int(arm1.start_pos),int(arm1.end_pos))
        y=range(int(arm2.start_pos),int(arm2.end_pos))


        if arm1.seq != arm2.seq and int(arm2.start_pos) > int(arm1.start_pos) and int(arm1.start_pos)+len(arm1.seq) > int(arm2.start_pos) and arm1.end_pos not in list4:
            if (str_1.count("G")+str_1.count("C"))*100/len(str_1) >= float(sys.argv[2]):
                number_trans += 1
                list4.append(arm1.end_pos)
            else:
                continue 

        elif arm1.seq != arm2.seq and int(arm2.start_pos) < int(arm1.start_pos) and int(arm2.start_pos)+len(arm2.seq) > int(arm1.start_pos)and arm1.end_pos not in list4:
            if (str_1.count("G")+str_1.count("C"))*100/len(str_1) >= float(sys.argv[2]):
                number_trans += 1
                list4.append(arm1.end_pos)
            else: 
                continue 


        elif range(int(arm1.start_pos),int(arm1.end_pos)) == range(int(arm2.start_pos),int(arm2.end_pos)): #palindrome
            if (str_1.count("G")+str_1.count("C"))*100/len(str_1) >= float(sys.argv[2]):
                number_pal += 1 
            else:
                continue 
        
        elif tuple1 not in list2 and tuple2 not in list1 and arm1.end_pos not in list3 and range(int(arm1.start_pos),int(arm1.end_pos)) != range(int(arm2.start_pos),int(arm2.end_pos)): 
            if arm2.start_pos in x or arm1.start_pos in y:
                continue 
            
            if (arm1.start_pos-1,arm1.end_pos) in list2 and (arm2.start_pos, arm2.end_pos+1) in list1: # remove the crypic redundant IRs
                continue

            if (arm1.start_pos-2,arm1.end_pos) in list2 and (arm2.start_pos, arm2.end_pos+2) in list1: # remove the crypic redundant IRs
                continue
            
            if (arm1.start_pos-3,arm1.end_pos) in list2 and (arm2.start_pos, arm2.end_pos+3) in list1: # remove the crypic redundant IRs
                continue
            
            if (arm1.start_pos-4,arm1.end_pos) in list2 and (arm2.start_pos, arm2.end_pos+4) in list1: # remove the crypic redundant IRs
                continue
            
            if (arm1.start_pos-5,arm1.end_pos) in list2 and (arm2.start_pos, arm2.end_pos+5) in list1: # remove the crypic redundant IRs
                continue
            
            else:
                if (str_1.count("G")+str_1.count("C"))*100/len(str_1) >= float(sys.argv[2]):

                    list1.append(tuple1)
                    list2.append(tuple2)
                    list3.append(arm1.end_pos)
                    number_irs +=1
                else:
                    continue 


        else:
            continue

    print (rna2.fbtr+","+rna2.gene+","+rna2.length+","+str(number_pal)+","+str(number_irs)+","+str(number_trans))



def find_complement(nucleotide, is_U):
    #given a character nucleotide input (ATCG or AUCG), returns the complementary nucleotide (ATCG)
    if (nucleotide == 'A'):
        if is_U == True:
            return 'U'
        else:
            return 'T'
    elif (nucleotide == 'T' or nucleotide == 'U'):
        return 'A'
    elif (nucleotide == 'C'):
        return 'G'
    elif (nucleotide == 'G'):
        return 'C'
    else:
        return 'NULL'

def is_complement(nucleotide1, nucleotide2):
    if(nucleotide1 == find_complement(nucleotide2, True)): #set to false if ATCG instead of AUCG
        return True
    else:
        return False


def check_whether_is_complementary_seq(rna1_start_ind, rna2_start_ind, rna1, rna2, arm_pairs, min_arm_len):
    #moves past the starting positions and checks whether the complementary
    #sequence meets our minimum length threshold, then add eligible ones to arm_pairs
    #for rna1, goes from left to right, for rna2, goes from right to left
    is_comp = is_complement(rna1.seq[rna1_start_ind], rna2.seq[rna2_start_ind]) #initialize
    curr_length = 0 #initialize, (start_index + curr_length - 1) = ending index, snippet = rna1[start_ind:end_ind+1]
    curr_index1 = rna1_start_ind
    curr_index2 = rna2_start_ind
    while (is_comp == True):
        curr_length += 1
        curr_index1 += 1 #the next to be checked
        curr_index2 -= 1 #for rna2 we go in the inverse direction
        try:
            is_comp = is_complement(rna1.seq[curr_index1], rna2.seq[curr_index2]) #move on to check the next nucleotide
        except: #if throws index out of bound exception, means reached the end of the string
            break
    if (curr_length < min_arm_len):
        return
    else:
        arm1 = Arm(rna1, rna1_start_ind+1, curr_index1, rna1.seq[rna1_start_ind:curr_index1], curr_length)
        arm2 = Arm(rna2, curr_index2+2, rna2_start_ind+1, rna2.seq[curr_index2+1:rna2_start_ind+1], curr_length)
        arm_pairs.append(ArmPair(arm1, arm2))
        


def find_irs(rna1, rna2, arm_pairs, minimum_arm_length):
    #finds all inter-complementary sequences (ics) with minimum arm length throughout rna1 and rna2
    #adds the found ics to "arm_pairs" (a list of objects of the class ArmPair)
    curr_index = 0
    max_len = len(rna1.seq) #saves computing time
    while curr_index < max_len:
        curr_nucleotide = rna1.seq[curr_index] #start by looking at the first nucleotide in arm1
        complement_nucleotide = find_complement(curr_nucleotide, True) #returns ATCG if false, AUCG if true
        #now, look for this nucleotide on the other arm
        #find syntax: str.find(str, beg=0, end=len(string))
        rna2_found_index = rna2.seq.find(complement_nucleotide, 0)
        while(rna2_found_index >= 0):
            check_whether_is_complementary_seq(curr_index, rna2_found_index, rna1, rna2, arm_pairs, minimum_arm_length)
            #the above function moves past the found positions and checks whether the complementary
            #sequence meets our minimum length threshold, then add eligible ones to arm_pairs
            rna2_found_index = rna2.seq.find(complement_nucleotide, rna2_found_index + 1) #finds next complement nucleotide
        #while loop terminates when no more starting complementary indices could be found
        #now we can move on to the next starting position, skipping over...
        #do not skip over for now. Consider skipping over repeats in the next version
        curr_index = curr_index + 1
    return


reader1= FASTAReader(open(sys.argv[1])) 

for ident, var1  in reader1:
    ident2=ident.rstrip("/n").split(";")
    fbgn2=ident2[-4]
    fbgn2=fbgn2.rstrip(";").split("=")
    parent=fbgn2[1]

    length_tscrip=ident2[-5]
    length_tscrip=length_tscrip.rstrip(";").split("=")
    length_tscrip=length_tscrip[1]

    fbtr2=ident2[0]
    fbtr1=fbtr2.split(" ")
    fbtr=fbtr1[0]

    dna2rna=str(var1).replace("T","U") #fasta is DNA seq.

    rna2 = RNA(fbtr, parent, length_tscrip, '3UTR', dna2rna)
    rna1=rna2
    arm_pairs = []
    main_function()
   #the row includes fbtr, fbgn, #pals, #irs, #trans  
