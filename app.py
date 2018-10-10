from flask import Flask, json, jsonify, render_template, render_template_string, request, session
from flask_sqlalchemy import SQLAlchemy
#from models import UnitCorrespondence
#from flask_marshmallow import Marshmallow


app = Flask(__name__)
app.config['SQLALCHEMY_DATABASE_URI'] = 'mysql+pymysql://root:root@127.0.0.1/rna3dhub-prod?unix_socket=/Applications/MAMP/tmp/mysql/mysql.sock'

db = SQLAlchemy(app)
#ma = Marshmallow(app)

from models import *

@app.route('/')
def home():
	return render_template("welcome.html")


@app.route('/correspondence')
def correspondence():
	
	#chain_info = '|'.join(unitid.split('|')[:3])
	#print chain_info

	unitid1 = request.args['unitid1']
	unitid2 = request.args['unitid2']

	res1_correspondence = []
	res1_correspondence_query = UnitCorrespondence.query.filter_by(unit_id_1 = unitid1).limit(20).all()
	for row in res1_correspondence_query:
		res1_correspondence.append(row.unit_id_2)

	res2_correspondence = []
	res2_correspondence_query = UnitCorrespondence.query.filter_by(unit_id_1 = unitid2).limit(20).all()
	for row in res2_correspondence_query:
		res2_correspondence.append(row.unit_id_2)

	res1_center_x = []
	res1_center_y = []
	res1_center_z = []

	res1_center_query = UnitCenters.query.filter(UnitCenters.unit_id.in_(res1_correspondence), UnitCenters.name == 'base').limit(20).all()
	for row in res1_center_query:
		res1_center_x.append(row.x)
		res1_center_y.append(row.y)
		res1_center_z.append(row.z)

	res2_center_x = []
	res2_center_y = []
	res2_center_z = []

	res2_center_query = UnitCenters.query.filter(UnitCenters.unit_id.in_(res2_correspondence), UnitCenters.name == 'base').limit(20).all()
	for row in res2_center_query:
		res2_center_x.append(row.x)
		res2_center_y.append(row.y)
		res2_center_z.append(row.z)

	res1_cell_0_0 = []
	res1_cell_0_1 = []
	res1_cell_0_2 = []
	
	res1_cell_1_0 = []
	res1_cell_1_1 = []
	res1_cell_1_2 = []

	res1_cell_2_0 = []
	res1_cell_2_1 = []
	res1_cell_2_2 = []

	res1_rotation_query = UnitRotations.query.filter(UnitRotations.unit_id.in_(res1_correspondence)).limit(20).all()
	for row in res1_rotation_query:
		res1_cell_0_0.append(row.cell_0_0)
		res1_cell_0_1.append(row.cell_0_1)
		res1_cell_0_2.append(row.cell_0_2)
	
		res1_cell_1_0.append(row.cell_1_0)
		res1_cell_1_1.append(row.cell_1_1)
		res1_cell_1_2.append(row.cell_1_2)

		res1_cell_2_0.append(row.cell_2_0)
		res1_cell_2_1.append(row.cell_2_1)
		res1_cell_2_2.append(row.cell_2_2)

	res2_cell_0_0 = []
	res2_cell_0_1 = []
	res2_cell_0_2 = []
	
	res2_cell_1_0 = []
	res2_cell_1_1 = []
	res2_cell_1_2 = []

	res2_cell_2_0 = []
	res2_cell_2_1 = []
	res2_cell_2_2 = []

	res2_rotation_query = UnitRotations.query.filter(UnitRotations.unit_id.in_(res2_correspondence)).limit(20).all()
	for row in res2_rotation_query:
		res2_cell_0_0.append(row.cell_0_0)
		res2_cell_0_1.append(row.cell_0_1)
		res2_cell_0_2.append(row.cell_0_2)
	
		res2_cell_1_0.append(row.cell_1_0)
		res2_cell_1_1.append(row.cell_1_1)
		res2_cell_1_2.append(row.cell_1_2)

		res2_cell_2_0.append(row.cell_2_0)
		res2_cell_2_1.append(row.cell_2_1)
		res2_cell_2_2.append(row.cell_2_2)

	table_rows = zip(res1_correspondence, res1_center_x, res1_center_y, res1_center_z, res2_correspondence, res2_center_x, res2_center_y, res2_center_z)

	table_rot1_rows = zip(res1_correspondence, res1_cell_0_0, res1_cell_0_1, res1_cell_0_2, res1_cell_1_0, res1_cell_1_1, res1_cell_1_2, res1_cell_2_0, res1_cell_2_1, res1_cell_2_2)
	
	table_rot2_rows = zip(res2_correspondence, res2_cell_0_0, res2_cell_0_1, res2_cell_0_2, res2_cell_1_0, res2_cell_1_1, res2_cell_1_2, res2_cell_2_0, res2_cell_2_1, res2_cell_2_2)

	return render_template("correspondence_table.html", query_res1 = unitid1, query_res2 = unitid2, table_items = table_rows, table_rot1_items = table_rot1_rows, table_rot2_items = table_rot2_rows) 

	#return unitid1 + unitid2

if __name__ == '__main__':
	app.run(debug=True)