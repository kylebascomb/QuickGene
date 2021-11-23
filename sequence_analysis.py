
import re
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sequences import format_seq, get_seq, write_fasta_file

NON_NUCLEOTIDE_COLOR = 'red'
PALINDROME_EVEN = 'blue'
PALINDROME_ODD = 'green'

def base_counter(seq):
    '''This function returns count of all bases in a sequence as a dictionary
        Ex: ACGGGTAC -> {'A': 2, 'C': 2, 'G': 3, 'T': 1}
        Parameters:
            seq (str): nucleotide sequence
        Returns:
            returns the count of all bases as a dictionary
'''
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


def gc_content(seq):
    '''This function returns the GC content of a sequence.
        Ex: If the sequence is 100 bases long and you have 20 C’s and 5 G’s, your GC content is 25%
        Parameters:
            seq (str): nucleotide sequence
        Returns:
            returns the GC content 
'''
    base_dict = base_counter(seq)
    base_count = base_dict.get('C') + base_dict.get('G') +base_dict.get('A') + base_dict.get('T')
    gc_count = base_dict.get('C') + base_dict.get('G')
    gc_content = gc_count / base_count
    
    return round(gc_content, 2)


def non_nucleotide_counter(seq):
    '''This function parses a sequence and returns a dictionary of the location of each 
        non ACGT base and the length of unknown bases if they are consecutive
        Ex: ACNGGGNNNTAC -> {2: 1, 6: 3}
        Parameters:
            seq (str): nucleotide sequence
        Returns:
            returns dictionary in the form of {position: length}
    '''
    non_nucleotides = []

    match = re.finditer(r'[^ATCG]{1,}', seq)
    for i in match:
        non_nucleotides.append([i.start(), i.end() -1])
    return non_nucleotides


def compile_html_colors(non_nucleotide, palindrome):
    '''
    This function creates a color code based on the location of non_nucleotides
    and palindromes
    Parameters:
        non_nucleotide (list): List of lists of the location of each non nucleotide in the sequence
        palindrome (list): List of lists of the location of each palindrome in the sequence
    Returns:
        Returns a list of dictionaries. Each dictionary represents the location and color of a subsequence
        example:{
            'color': 'red',
            'start': 100,
            'end': 105
        }

    '''
    colors = []
    for entry in non_nucleotide:
        colors.append({
            'color': NON_NUCLEOTIDE_COLOR,
            'start': entry[0],
            'end': entry[1]
        })
    index = 0
    for entry in palindrome:
        color = PALINDROME_ODD
        if index % 2 == 0:
            color = PALINDROME_EVEN
        colors.append({
            'color': color,
            'start': entry[0],
            'end': entry[1]
        })
        index +=1
    return colors


def find_palindromes(seq):
    '''This function parses the sequence and returns a list of palindromic nucleotide sequences
        within the seq. 
        Parameters:
            seq (str): sequence 
        Returns:
            Returns a list of all the palindromes. Each entry in the list is another with the start position at index 0
            and the end position at index 1

    '''
    palindromes = []
 
    for i in range(len(seq)):
 
        # odd length strings
        expand(seq, i, i, palindromes)
        # even length strings
        expand(seq, i, i + 1, palindromes)
 
    # print all unique palindromic substrings
    return remove_nested_palindromes(palindromes)

# recursive function called in find_palindromes()
def expand(seq, low, high, palindromes):
    '''
    Helper function for finding palindromes
    '''
 
    # run till `s[low.high]` is not a palindrome
    while low >= 0 and high < len(seq) and compare_complimentary(seq[low], seq[high]):
        # update pointers
        low = low - 1
        high = high + 1

    # seq must be length of 3 or longer
    if (high - low  > 3):
        palindromes.append([low +1, high -1])

def remove_nested_palindromes(palindromes):
    '''
    Helper function to remove palindromes that overlap.
    '''
    current_max = palindromes[0][1]
    remove_list = []
    for i in range(len(palindromes)):
        if palindromes[i][0] <= current_max:
            remove_list.append(palindromes[i])
        else:
            current_max = palindromes[i][1]
    for item in remove_list:
        palindromes.remove(item)
    return palindromes



def compare_complimentary(a, b):
    '''
    Helper function to compare if two bases are complementary
    Parameters:
        a (str): first base
        b (str): second base
    Returns:
        Returns True if a and b are complimentary
    '''
    if a == 'A':
        return b == 'T'
    elif a == 'T':
        return b == 'A'
    elif a =='C':
        return b =='G'
    elif a == 'G':
        return b == 'C'
    else:
        return False


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


def format_seq_as_lines(sequence, characters_per_line):
    '''
    This is a helper function that breaks the entire sequence into a list of lines 
    based on the characters_per_line. 
    Parameters:
        sequence (str): Long sequence of characters
        characters_per_line (int): number of characters that each line should contain.
    Returns:
        Returns a list of strings, each with the length specified in characters_per_line. 
    '''
    lines = []
    for i in range(0, len(sequence), characters_per_line):
        lines.append([i, sequence[i:i+characters_per_line]])
    return lines



def format_seq_with_html(sequence, color_code, characters_per_line):
    '''
    This function formats the the sequence with HTML spans according the the color code
    Parameters:
        sequence (str): The sequence of base pairs
        color_code (list): list of dictionaries that represents the colors and location that need to be formatted
        with HTML spans
        characters_per_line (int): number of characters per line 
    Returns:
        Returns a list of all the lines formatted with HTML. Each entry in the list represents one line.
    '''
    lines = format_seq_as_lines(sequence, characters_per_line )
    html_lines=[]

    def wrap_in_p_tags(start, line):
        start_html = '<p>' + str(start) +'   '
        end_html = '</p>'
        return start_html + line + end_html
    
    def wrap_in_grid(start, line):
        start_html = '<div class="row"><div class="col">'+str(start)+'</div>'
        end_html='<div class="col-md-auto"> <p class="text-monospace text-left">'+str(line)+'</p></div></div>'
        return start_html+end_html
    
    def wrap_in_span(color, subseq):
        start_html = '<span style="color: ' + color +'">'
        end_html = '</span>'
        return start_html + subseq + end_html

    def format_line(line_splits):
        line = ''
        for split in line_splits:
            seq = split[0]
            color = split[1]
            if seq != '':
                if color != '':
                    line += wrap_in_span(color, seq)
                else:
                    line+= seq
        return line
    
    def find_color_codes(color_code, start, end):
        codes_in_line = []
        for code in color_code:
            if code['start'] >= start and code['start'] <= end:
                codes_in_line.append(code)      
            elif code['start'] < start and code['end'] <= end and code['end'] > start:
                codes_in_line.append(code)
            elif code['start'] < start and code['end'] > end:
                codes_in_line.append(code)
        return codes_in_line

    
    for line in lines:
        start = line[0]
        subseq = line[1]
        end = start + len(line[1])
        line_splits = []
        codes_in_line = find_color_codes(color_code, start, end)
        
        cur = 0 #cursor

        if len(codes_in_line) == 0:
            line_splits.append([subseq, ''])
        else:
            codes_in_line = sorted(codes_in_line, key=lambda k: k['start'])
            for i in range(len(codes_in_line)): 
                code = codes_in_line[i]
                #if there is no color to start
                if code['start'] > start:
                    line_splits.append([subseq[cur:code['start'] - start ], '']) 
                   #if no overflow
                    if code['end'] <= end:  
                        line_splits.append([subseq[code['start'] - start :code['end']- start + 1], code['color']])
                    else:
                       line_splits.append([subseq[code['start']- start:], code['color']])
                else:
                    #color from start
                    #if no overflow
                    if code['end'] <= end:
                        line_splits.append([subseq[cur:code['end']-start + 1], code['color']])
                    else:
                        line_splits.append([subseq, code['color']])
                if i == (len(codes_in_line) - 1) and code['end'] <= end:
                    line_splits.append([subseq[code['end']-start + 1:], ''])
                cur = code['end']-start + 1
        html_lines.append(wrap_in_grid(start, format_line(line_splits )))

    return html_lines



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
    #html
    colors = compile_html_colors(non_nucleotide_counter(seq_info['sequence']), find_palindromes(seq_info['sequence']))
    seq_info['seq_html'] = format_seq_with_html(seq_info['sequence'], colors, 100)
    #create files
    create_count_chart(seq_info['base_counts'], chart_path)
    write_fasta_file(seq_info, fasta_path)
    return seq_info


#save_chart(create_count_chart({'A': 2, 'C': 2, 'G': 3, 'T': 1}), 'testfig')
#print(compile_analysis_from_id('AH002560.3'))

