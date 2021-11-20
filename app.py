from flask import Flask, render_template, request, url_for, flash, redirect
from sequence_analysis import compile_analysis_from_id


app = Flask(__name__)
app.config['SECRET_KEY'] = 'HJskjAn92nSA2SDA239SDAN1312DSA8'

'''
TODO @Michael
Primary Task
Implement basic HTML for viewing some of the seq_info in 
templates/index.html.

Display the following information:
seq_info['name']
seq_info['description']
Bar Chart (Done)
seq_info['gc_content']
seq_info['sequence']

I have already implemented an <img> for the bar chart, 
but you can make changes to that if you need to fix the formatting.

If you need to test and see actual output you can use the Gene Id: AH002560.3

Additional Task
Create a button to download the file with the path stored in seq_info['fasta_path']
'''
@app.route("/", methods=('GET', 'POST'))
def index():
    seq_info={}
    if request.method =='POST':
        id = request.form['id']
        seq_info = compile_analysis_from_id(id)
        # flash('Errors') TODO
        
    return render_template('index.html', seq_info=seq_info)




if __name__ == '__main__':
    app.run()