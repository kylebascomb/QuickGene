
import math
'''This function returns count of all bases in a sequence as a dictionary
        Ex: ACGGGTAC -> {'A': 2, 'C': 2, 'G': 3, 'T': 1}
        Parameters:
            seq (str): nucleotide sequence
        Returns:
            returns the count of all bases as a dictionary
'''
def base_counter(seq):
    base_count = {'A':0, 'C':0, 'G':0, 'T':0}

    for i in seq:
        if i == 'A':
            base_count['A'] += 1
        elif i == 'G':
            base_count['G'] += 1
        elif i == 'C':
            base_count['C'] += 1
        elif i == 'T':
            base_count['T'] += 1
    
    return base_count
#print(base_counter("AAGGTTCCAGGT"))

'''This function returns the GC content of a sequence.
        Ex: If the sequence is 100 bases long and you have 20 C’s and 5 G’s, your GC content is 25%
        Parameters:
            seq (str): nucleotide sequence
        Returns:
            returns the GC content 
'''
def gc_content(seq):
    base_count = base_counter(seq)
    gc_content = base_count.get('C') / base_count.get('G')

    
    
    return round(gc_content, 2)

#print(gc_content("AAAAGGGTTCCCCCCCCCC"))


def non_nucleotide_counter(seq):
    '''This function parses a sequence and returns a dictionary of the location of each 
        non ACGT base and the length of unknown bases if they are consecutive
        Ex: ACNGGGNNNTAC -> {2: 1, 6, 3}
        Parameters:
            seq (str): nucleotide sequence
        Returns:
            returns dictionary in the form of {position: length}
    '''
    return {}