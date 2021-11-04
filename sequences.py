from Bio import Entrez, SeqIO
import sys

email = 'kylebascomb@gmail.com'     # should be moved to an environment variable or argument

def get_seq(id, email = email):
    '''This function calls the Entrez efetch function to get the FASTA 
        sequence related to the id.
        Parameters:
            id (str): id of the protein sequence
            email (email) : user email necessary for Entrez usage
        Returns:
            returns a SeqIO seq record of the FASTA sequence
    '''
    Entrez.email = email
    try:
        request = Entrez.efetch(db="nucleotide", id=id, rettype="fasta")
        seq_record = SeqIO.read(request, "fasta")
        return seq_record
    except Exception as e:
        print("ERROR: ID not found:", id)
        print("Error Message:", e)

    return None


def format_seq(seq_record):
    '''
    This function formats the information in a Seq record into a dictionary.
    Parameters:
        seq_record (Bio.Seq): the seq record to be formatted
    Returns:
        returns a dictionary with the following keys: id, name, description, and sequence
    '''
    return {
        'id': seq_record.id,
        'name': seq_record.name,
        'description': seq_record.description,
        'sequence': str(seq_record.seq.upper())
    
    }


print(format_seq(get_seq('AH002560.3')))