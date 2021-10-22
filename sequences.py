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


print(get_seq('AH002560.3').seq)