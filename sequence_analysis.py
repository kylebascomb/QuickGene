
import math
from Bio.SeqIO import write
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sequences import format_seq, get_seq, write_fasta_file


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
    base_dict = base_counter(seq)
    base_count = base_dict.get('C') + base_dict.get('G') +base_dict.get('A') + base_dict.get('T')
    gc_count = base_dict.get('C') + base_dict.get('G')
    gc_content = gc_count / base_count
    
    return round(gc_content, 2)

#print(gc_content("ATTTTGC"))


def non_nucleotide_counter(seq):
    '''This function parses a sequence and returns a dictionary of the location of each 
        non ACGT base and the length of unknown bases if they are consecutive
        Ex: ACNGGGNNNTAC -> {2: 1, 6: 3}
        Parameters:
            seq (str): nucleotide sequence
        Returns:
            returns dictionary in the form of {position: length}
    '''
    return {}


def find_palindromes(seq):
    '''This function parses the sequence and returns a list of palindromic nucleotide sequences
        within the seq. Returns a list of all unique palindromes

    '''


def create_count_chart(base_count,filepath):
    '''
    This function uses seaborn to create a barchart of the dictionary returned in
    base_counter. It also saves this chart to the path specified
    Parameters:
        base_count (dict): A dictionary with the keys as the base, and the values as the count of each base
            in a sequence
        filepath (str): path to the file to save the chart
    Returns:
        Returns the Figure object
    '''
    plot = sns.barplot(x=list(base_count.keys()), y=list(base_count.values()))
    plot.get_figure().savefig(filepath)
    return plot.get_figure()


def compile_analysis_from_id(id):
    '''
    This function compiles all the analysis into a single dictionary to be used by Flask
    Parameters:
        id (str): Sequence ID
    Returns:
        returns a dictionary of all the sequence information
    '''
    seq_info = format_seq(get_seq(id))
    chart_path = './static/charts/' + id + '.png'
    fasta_path = './static/fasta/' + id + '.fasta'
    seq_info['base_counts'] = base_counter(seq_info['sequence'])
    seq_info['gc_content'] = gc_content(seq_info['sequence'])
    seq_info['chart'] = create_count_chart(seq_info['base_counts'], chart_path)
    seq_info['chart_path'] = chart_path
    seq_info['fasta_path'] = fasta_path

    create_count_chart(seq_info['base_counts'], chart_path)
    write_fasta_file(seq_info, fasta_path)
    return seq_info


#save_chart(create_count_chart({'A': 2, 'C': 2, 'G': 3, 'T': 1}), 'testfig')
print(compile_analysis_from_id('AH002560.3'))

