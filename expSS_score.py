#given the output of RNAfold, calculates a score based on
#the max number of continuous sticky nucleotides that can actually
#base pair in trans (i.e. both arms are exposed)

#param structure =a string of structure in rnafold binary output format
#returns (1) structure in ooxx form (o = both arm's matching nts exposed)
#e.g. ATCGG-CCGAT(....&.)... -> xooxo (order based off of 1st arm)
#(2) max number of consecutive trans-base-pair-able nts
def advanced_structure_info(arm1, arm2):
    curr_count = 0 #counter for current consecutive sticky nucleotides
    max_count = 0 #keeps track of the max # of consecutive sticky nucleotides
    strout = ""
    if len(arm1) != len(arm2):
        return "irregular sequence", "NA"
    for i in range(len(arm1)):
        #printing is SUPER inefficient, avoid printing next time
        if arm1[i] == '.' and arm2[-1-i] == '.':
            strout = strout + 'o'
            curr_count += 1 #increment consecutive nucleotide count
        else: # if a non-exposed nucleotide happens
            strout = strout + 'x'
            if curr_count > max_count: 
                max_count = curr_count #update max count if needed
            curr_count = 0 #reset curr count
    if curr_count > max_count: 
        max_count = curr_count #in case it is continuous to the end
    return strout, max_count
