from flask import Flask, render_template, request, session
from flask_sqlalchemy import SQLAlchemy

app = Flask(__name__)
app.config['SQLALCHEMY_DATABASE_URI'] = 'mysql+pymysql://root:root@127.0.0.1/rna3dhub-prod?unix_socket=/Applications/MAMP/tmp/mysql/mysql.sock'
db = SQLAlchemy(app)

from models import *

@app.route('/')
def home():
	return render_template("welcome.html")

@app.route('/correspondence/<unitid>', methods=['GET', 'POST'])
def correspondence(unitid):
	test_query = UnitCorrespondence.query.filter_by(unit_id_1 = unitid).limit(10).all()
	test_query = str(test_query)
	return test_query

if __name__ == '__main__':
	app.run(debug=True)