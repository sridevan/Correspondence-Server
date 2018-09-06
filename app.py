from flask import Flask, render_template
from flask import session
from flask_sqlalchemy import SQLAlchemy

app = Flask(__name__)
app.config['SQLALCHEMY_DATABASE_URI'] = 'mysql+pymysql://root:root@127.0.0.1/rna3dhub-prod?unix_socket=/Applications/MAMP/tmp/mysql/mysql.sock'
db = SQLAlchemy(app)

from models import *

@app.route('/')
def home():
	test_query = UnitInfo.query.filter_by(pdb_id='1ffk')
	return test_query

@app.route('/correspondence/<unitid>')
def correspondence(unitid):
	return unitid
	
@app.route('/welcome')
def welcome():
	return render_template("welcome.html")

if __name__ == '__main__':
	app.run(debug=True)