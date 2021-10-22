



def gc_content(seq):
    '''This function returns the GC content of a sequence.
        Ex: If the sequence is 100 bases long and you have 20 C’s and 5 G’s, your GC content is 25%
        Parameters:
            seq (str): nucleotide sequence
        Returns:
            returns the GC content 
    '''
    return 0


def base_counter(seq):
    '''This function returns count of all bases in a sequence as a dictionary
        Ex: ACGGGTAC -> {'A': 2, 'C': 2, 'G': 3, 'T': 1}
        Parameters:
            seq (str): nucleotide sequence
        Returns:
            returns the count of all bases as a dictionary
    '''
    return {}


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