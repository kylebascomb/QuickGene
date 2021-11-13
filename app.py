from flask import Flask, render_template, request, url_for, flash, redirect
from sequence_analysis import compile_analysis_from_id


app = Flask(__name__)
app.config['SECRET_KEY'] = 'HJskjAn92nSA2SDA239SDAN1312DSA8'


@app.route("/", methods=('GET', 'POST'))
def index():
    seq_info={}
    if request.method =='POST':
        id = request.form['id']
        seq_info = compile_analysis_from_id(id)
        print(seq_info)
        # flash('Errors') TODO
        
    return render_template('index.html', seq_info=seq_info)

@app.route("/analysis/<id>")
def weekly_pictures(id):
    analysis = compile_analysis_from_id(id)
    return render_template("analysis.html", analysis=analysis)



if __name__ == '__main__':
    app.run()