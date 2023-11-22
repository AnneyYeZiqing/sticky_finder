#(Inter-RNA-Species Trans Inverted Repeats) finder
#written by Ziqing Ye, modified by Siran Tian @ Trcek Lab @ Johns Hopkins University
#last updated 02/10/22

from expSS_score import advanced_structure_info

#have a transcript class, which contains sequence, source gene, ...
#lines with asterics ** are lines that involve structure/RNAfold
#in case people don't want RNAfold they can change those ** lines to get rida structure

#define RNA class
class RNA(object):
    def __init__(self, fbtr, parent_gene_symbol, sequence_type, sequence, structure): #**
        self.fbtr = fbtr
        self.gene = parent_gene_symbol #nos, pgc, etc.
        self.seq_type = sequence_type #cdna or 3'UTR, etc.
        self.seq = sequence #the string sequence
        self.structure = structure #the string of binary RNAfold output **

#import sequence
rna2 = RNA('nos-005-nonpalindromic', 'CRISPR', '3UTR',
           'ATTTATTTAAGTAAAAAATTTACTAAATGTTTTCGCGAAAGGCAATTGCCATATGCTTGCCTGGCAAGCATAGCGATGAGTGATCGCCACTCAATGCCTTCAGCATTTGAAGATGCAGCACTTTGCATAAGTGCCGATTGCAATAATCGTTACATTTATTAAACATAACAATTAATCTTGCGATACATGTTTCGCTTTGCGGAATGTCAAAATTTAAAATTTTACAGCCGACGACGAAAGTGTTCCTTGCTATTTCCTTTAGCAAGATTTAAATTTAGATTAAATTCTAATGATACGATTGACAGTTCGAAATTCAAAGTGTTCCTGATCATATCGAAATTTTTCGGCCGCAAGCGAACATTTTACGAAAAGTTGCAACTTTCGTATCGGCCACGACGATTGAACAAGTATTACGATATTGTAAGTCTGTTTCGAAACAGACAGACTTGTGATCGTTCGTTGTCTATACTATAAGATCTATAGGCACGGGATAACGCTCT',
           '....(((((.(((((...))))))))))(((.((.((...(((....))).(((((((((...)))))))))(((.((((((.....)))))).))).....(((((....))))).((((((.....))))))..(((((((((((.((((..(((...)))..)))).))))...)))))))..(((((((((...))))))))).((((((...))))))...(((((..((((((((((..(((((((.......)))))))........(((((......))))).((((.((((((((.((.((...)).)).))))...)))).))))..(((((((((...................))))))))))).)))))))).))))).(((((((((..((((((.((((......))))(((((((....)))))))..))))))))))).))))((((((((..........)))))))))).)).))).....')



rna1 = rna2


#other parameters
min_len = 6 #minimum arm length


#initialize list of arm pairs
arm_pairs = []


#define arm class
class Arm(object): #this is a class
    def __init__(self, source_transcript, starting_position, ending_position, sequence, length, structure): #this is a constructor **
        #use: Arm(source_transcript, starting_position, ending_position, sequence, length)
        self.transcript = source_transcript #the source transcript as an RNA class object
        self.start_pos = starting_position #starting nucleotide position of the transcript (nucleotide index starts with 1!)
        self.end_pos = ending_position
        self.seq = sequence #the nucleotide sequence of the arm, a string. e.g. "ATCCG"
        self.length = length #the length of the arm
        self.structure = structure #binary structure from RNAfold **
        #ending position = starting position + length - 1
        
class ArmPair(object):
    def __init__(self, arm1, arm2):
        self.arm1 = arm1 #class arm
        self.arm2 = arm2 #class arm

#start!
def main_function():
    file1=open('method-nodms_nos-palindromic_3UTR_palir_min=6_S1.csv', 'w')
    file1.write('transcript 1,' + rna1.fbtr+','+rna1.gene+','+rna1.seq_type+',transcript 2,' + rna2.fbtr + ',' + rna2.gene+','+rna2.seq_type+'\n')
    file1.write('arm1 source,starting position,ending position,sequence,structure, arm2 source,starting position,ending position,sequence,structure,arm length,net structure,structure score\n')
    
    file2=open('method-nodms_nos-palindromic_3UTR_trans_min=6_S1.csv', 'w')
    file2.write('transcript 1,' + rna1.fbtr+','+rna1.gene+','+rna1.seq_type+',transcript 2,' + rna2.fbtr + ',' + rna2.gene+','+rna2.seq_type+'\n')
    file2.write('arm1 source,starting position,ending position,sequence,structure, arm2 source,starting position,ending position,sequence,structure,arm length,net structure,structure score\n')

    print('rna1=', rna1, '\nrna2=', rna2)
    find_irs(rna1, rna2, arm_pairs, min_len) #find ics

    list1=[]
    list2=[]
    list3=[]
    list4=[]

    for ics in arm_pairs:
        arm1 = ics.arm1
        arm2 = ics.arm2

        sum_structure = advanced_structure_info(arm1.structure, arm2.structure) #xxoo and max length of consecutive exposed nucleotides


        tuple1=(int(arm1.start_pos),int(arm1.end_pos))
        tuple2=(int(arm2.start_pos),int(arm2.end_pos))

        x=range(int(arm1.start_pos),int(arm1.end_pos))
        y=range(int(arm2.start_pos),int(arm2.end_pos))

        if arm1.seq != arm2.seq and int(arm2.start_pos) > int(arm1.start_pos) and int(arm1.start_pos)+len(arm1.seq) > int(arm2.start_pos) and arm1.end_pos not in list4:

            file2.write(arm1.transcript.gene + ',' + str(arm1.start_pos) + ',' + str(arm1.end_pos) + ',' + arm1.seq + ',' + arm1.structure +
                        ',' + arm2.transcript.gene + ',' + str(arm2.start_pos) + ',' + str(arm2.end_pos) + ',' + arm2.seq + ',' + arm2.structure +
                        ',' + str(arm1.length) + ',' + sum_structure[0] + ',' + str(sum_structure[1]) + '\n') 
            list4.append(arm1.end_pos)
        
        elif arm1.seq != arm2.seq and int(arm2.start_pos) < int(arm1.start_pos) and int(arm2.start_pos)+len(arm2.seq) > int(arm1.start_pos) and arm1.end_pos not in list4:
            file2.write(arm1.transcript.gene + ',' + str(arm1.start_pos) + ',' + str(arm1.end_pos) + ',' + arm1.seq + ',' + arm1.structure +
                ',' + arm2.transcript.gene + ',' + str(arm2.start_pos) + ',' + str(arm2.end_pos) + ',' + arm2.seq + ',' + arm2.structure +
                ',' + str(arm1.length) + ',' + sum_structure[0] + ',' + str(sum_structure[1]) + '\n') 
            list4.append(arm1.end_pos)

        elif range(int(arm1.start_pos),int(arm1.end_pos)) == range(int(arm2.start_pos),int(arm2.end_pos)): #palindrome
            file1.write(arm1.transcript.gene + ',' + str(arm1.start_pos) + ',' + str(arm1.end_pos) + ',' + arm1.seq + ',' + arm1.structure +
            ',' + arm2.transcript.gene + ',' + str(arm2.start_pos) + ',' + str(arm2.end_pos) + ',' + arm2.seq + ',' + arm2.structure +
            ',' + str(arm1.length) + ',' + sum_structure[0] + ',' + str(sum_structure[1]) + '\n')

        
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

                file1.write(arm1.transcript.gene + ',' + str(arm1.start_pos) + ',' + str(arm1.end_pos) + ',' + arm1.seq + ',' + arm1.structure +
                            ',' + arm2.transcript.gene + ',' + str(arm2.start_pos) + ',' + str(arm2.end_pos) + ',' + arm2.seq + ',' + arm2.structure +
                            ',' + str(arm1.length) + ',' + sum_structure[0] + ',' + str(sum_structure[1]) + '\n')
            
                list1.append(tuple1)
                list2.append(tuple2)
                list3.append(arm1.end_pos)


        else:
            continue 


    file1.close()
    file2.close()

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
        #use: Arm(source_transcript, starting_position, ending_position, sequence, length)
        #print('arm1 start pos=', rna1_start_ind+1, ', end pos=', curr_index1, ', seq=', rna1.seq[rna1_start_ind:curr_index1], ', len=', curr_length)
        #print('arm2 start pos=', curr_index2+2, ', end pos=', rna2_start_ind+1, ', seq=', rna2.seq[(curr_index2+1):rna2_start_ind+1], ', len=', curr_length)
        #Arm class usage: Arm('rna1', starting_position, ending_position, sequence, length, structure)
        arm1 = Arm(rna1, rna1_start_ind+1, curr_index1, rna1.seq[rna1_start_ind:curr_index1], curr_length, rna1.structure[rna1_start_ind:curr_index1])
        arm2 = Arm(rna2, curr_index2+2, rna2_start_ind+1, rna2.seq[curr_index2+1:rna2_start_ind+1], curr_length, rna2.structure[(curr_index2+1):rna2_start_ind+1])
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
        curr_index = curr_index + 1
    print('scanning finished')
    return


main_function()