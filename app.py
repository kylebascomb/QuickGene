from flask import Flask, render_template, request, send_file, flash
from sequence_analysis import compile_analysis_from_id
import re
import os

app = Flask(__name__)
app.config['SECRET_KEY'] = 'HJskjAn92nSA2SDA239SDAN1312DSA8'


@app.route("/", methods=('GET', 'POST'))
def index():
    seq_info={}
    if request.method =='POST':
        id = request.form['id']
        if(len(id) == 0):
            flash('Invalid Sequence ID')
        else:
            seq_info = compile_analysis_from_id(id)
            if seq_info['errors'] is not None:
                flash('Invalid Sequence ID')   
        
        
    return render_template('index.html', seq_info=seq_info)


@app.route('/static/<filename>', methods=['GET', 'POST'])
def download_fasta(filename):
    return send_file('/static/'+filename, as_attachment=True)


if __name__ == '__main__':
    app.run(port=os.environ.get('PORT', 5000))