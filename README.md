# QuickGene
Flask Project for basic bioinformatics analysis using data from NCBI
QuickGene is a Flask based web page that will be used to get a quick overview of a Nucleotide sequence using an Accession ID from NCBI. 
The user will fill in a form on the webpage with the Accession ID. Then, the sequence will be retrieved in the backend using the BioPython Entrez API. 
Once we have the sequence, we will do a few different analyses on the sequence using the concepts we have learned in class. Some ideas for the analysis so far include: 
counting the number of A’s, C’s, G’s, T’s in the sequence and reporting them in a barchart, finding any characters that are not A/C/G/T and reporting their position,
reporting the length of the sequence, and more (we are still brainstorming what would be useful). The analysis would be sent back to the frontend so the user can see the data.

One of the main motivators for this project is to create a user interface that is a bit easier to take in than the default NCBI view. For many beginner Bioinformaticians, the amount of information available on an NCBI page can be overwhelming. The other motivations for the project are mostly for our own educational purposes. We would like to learn how to apply the new Python concepts we have been learning to bioinformatics topics. We are also using Flask for this project to learn more about backend frameworks in Python. 
