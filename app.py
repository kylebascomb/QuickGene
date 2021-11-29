from flask import Flask, render_template, request, url_for, flash, redirect
from sequence_analysis import compile_analysis_from_id
import re

app = Flask(__name__)
app.config['SECRET_KEY'] = 'HJskjAn92nSA2SDA239SDAN1312DSA8'

'''
TODO @Michael
There is a button with the text "Download FASTA File" in index.html.
I need that button to allow the user to download the file with the path stored in seq_info['fasta_path']

TODO @Adam
Create some simple error handling to notify the user if the id is length 0 or if 
the sequence id is invalid. To achieve this you may have to bring the error from get_seq
sequences.py up into the index route to be able to check it. You can find more information 
on this issue in sequences.py in the format_seq defintion

I think flash can be used for this, but any other method will work fine.
'''
@app.route("/", methods=('GET', 'POST'))
def index():
    seq_info={}
    if request.method =='POST':
        id = request.form['id']
        matched = re.match("[A-Z]{2}[0-9]{6}.[0-9]", id)
        is_match = bool(matched)
        if(len(id) == 0):
            flash('Invalid Sequence ID')
        elif(not is_match):
            flash('Invalid Sequence ID')
        else:
            seq_info = compile_analysis_from_id(id)
        
        
    return render_template('index.html', seq_info=seq_info)




if __name__ == '__main__':
    app.run()