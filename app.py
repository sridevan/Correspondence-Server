import collections as coll
from flask import Flask, json, jsonify, render_template, render_template_string, request, session
from flask_bootstrap import Bootstrap
from flask_sqlalchemy import SQLAlchemy
import numpy as np
from numpy import nan
np.set_printoptions(threshold=np.nan)
from sqlalchemy.sql.expression import case
#from models import UnitCorrespondence
#from flask_marshmallow import Marshmallow


app = Flask(__name__)
Bootstrap(app)
app.config['SQLALCHEMY_DATABASE_URI'] = 'mysql+pymysql://root:root@127.0.0.1/rna3dhub-prod?unix_socket=/Applications/MAMP/tmp/mysql/mysql.sock'

db = SQLAlchemy(app)

from models import *
from discrepancy import *
from greedyInsertion import *

@app.route('/')
def home():
	return render_template("welcome.html")


@app.route('/correspondence')
def correspondence():
	
	#chain_info = '|'.join(unitid.split('|')[:3])
	#print chain_info

	unitid1 = request.args['unitid1']
	unitid2 = request.args['unitid2']
	unitid3 = request.args['unitid3']

	# Need to automate this part

	res1_correspondence = []
	res1_correspondence_query = UnitCorrespondence.query.filter_by(unit_id_1 = unitid1).limit(10).all()
	for row in res1_correspondence_query:
		res1_correspondence.append(row.unit_id_2)

	res2_correspondence = []
	res2_correspondence_query = UnitCorrespondence.query.filter_by(unit_id_1 = unitid2).limit(10).all()
	for row in res2_correspondence_query:
		res2_correspondence.append(row.unit_id_2)

	res3_correspondence = []
	res3_correspondence_query = UnitCorrespondence.query.filter_by(unit_id_1 = unitid3).limit(10).all()
	for row in res3_correspondence_query:
		res3_correspondence.append(row.unit_id_2)

	combined_res = zip(res1_correspondence, res2_correspondence, res3_correspondence)

	ifes = []
	for elem in  combined_res:
		ifes.append(elem[0][:8])

	coord = []
	for x in combined_res:
		x = ','.join(x)
		coord.append(x)

	ife_full = dict(zip(ifes, coord))

	c1 = []
	c2 = []
	c3 = []

	unitid_c1 = []
	unitid_c2 = []
	unitid_c3 = []

	ordering = case(
    {unit_id: index for index, unit_id in enumerate(res1_correspondence)},
    value=UnitCenters.unit_id
 	)

	res1_center_query = UnitCenters.query.filter(UnitCenters.unit_id.in_(res1_correspondence), UnitCenters.name == 'base').order_by(ordering).limit(10).all()
	for row in res1_center_query:
		c1.append(np.array([row.x, row.y, row.z]))
		unitid_c1.append(row.unit_id)

	ordering = case(
    {unit_id: index for index, unit_id in enumerate(res2_correspondence)},
    value=UnitCenters.unit_id
 	)

	res2_center_query = UnitCenters.query.filter(UnitCenters.unit_id.in_(res2_correspondence), UnitCenters.name == 'base').order_by(ordering).limit(10).all()
	for row in res2_center_query:
		c2.append(np.array([row.x, row.y, row.z]))
		unitid_c2.append(row.unit_id)

	ordering = case(
    {unit_id: index for index, unit_id in enumerate(res3_correspondence)},
    value=UnitCenters.unit_id
 	)

	res3_center_query = UnitCenters.query.filter(UnitCenters.unit_id.in_(res3_correspondence), UnitCenters.name == 'base').order_by(ordering).limit(10).all()
	for row in res3_center_query:
		c3.append(np.array([row.x, row.y, row.z]))
		unitid_c3.append(row.unit_id)

	newcenter = zip(c1, c2, c3)

	zipped_unitid_c = zip(unitid_c1, unitid_c2, unitid_c3)

	r1 = []
	r2 = []
	r3 = []

	unitid_r1 = []
	unitid_r2 = []
	unitid_r3 = []

	ordering = case(
    {unit_id: index for index, unit_id in enumerate(res1_correspondence)},
    value=UnitRotations.unit_id
 	)

	res1_rotation_query = UnitRotations.query.filter(UnitRotations.unit_id.in_(res1_correspondence)).order_by(ordering).limit(10).all()
	for row in res1_rotation_query:
		r1.append(np.array([[row.cell_0_0, row.cell_0_1, row.cell_0_2], 
							[row.cell_1_0, row.cell_1_1, row.cell_1_2], 
							[row.cell_2_0, row.cell_2_1, row.cell_2_2]]))
		unitid_r1.append(row.unit_id)

	ordering = case(
    {unit_id: index for index, unit_id in enumerate(res2_correspondence)},
    value=UnitRotations.unit_id
 	)

 	res2_rotation_query = UnitRotations.query.filter(UnitRotations.unit_id.in_(res2_correspondence)).order_by(ordering).limit(10).all()
	for row in res2_rotation_query:
		r2.append(np.array([[row.cell_0_0, row.cell_0_1, row.cell_0_2], 
							[row.cell_1_0, row.cell_1_1, row.cell_1_2], 
							[row.cell_2_0, row.cell_2_1, row.cell_2_2]]))
		unitid_r2.append(row.unit_id)

	ordering = case(
    {unit_id: index for index, unit_id in enumerate(res3_correspondence)},
    value=UnitRotations.unit_id
 	)

 	res3_rotation_query = UnitRotations.query.filter(UnitRotations.unit_id.in_(res3_correspondence)).order_by(ordering).limit(10).all()
	for row in res3_rotation_query:
		r3.append(np.array([[row.cell_0_0, row.cell_0_1, row.cell_0_2], 
							[row.cell_1_0, row.cell_1_1, row.cell_1_2], 
							[row.cell_2_0, row.cell_2_1, row.cell_2_2]]))
		unitid_r3.append(row.unit_id)
		

	newrotation = zip(r1, r2, r3)

	zipped_unitid_r = zip(unitid_r1, unitid_r2, unitid_r3)

	
	# This part deals with calculating the discepancy for the corresponding instances
	distances = coll.defaultdict(lambda: coll.defaultdict(int))
	corr1 = []
	corr2 = []

	for a in range(0, len(combined_res)):
		for b in range(a+1, len(combined_res)):
			corr1.append(combined_res[a])
			corr2.append(combined_res[b])
			disc = matrix_discrepancy(newcenter[a], newrotation[a], newcenter[b], newrotation[b])
			distances[combined_res[a][0][:8]][combined_res[b][0][:8]] = disc

	dist = np.zeros((len(ifes), len(ifes)))
        for index1, member1 in enumerate(ifes):

            curr = distances.get(member1, {})
            for index2, member2 in enumerate(ifes):
                val = curr.get(member2, None)
                if member2 not in curr:
                    val = None
                dist[index1, index2] = val
	
	ordering, _ , _ = orderWithPathLengthFromDistanceMatrix(dist, 10, scanForNan=True)
                                                            

	neworder = [x for x in sorted(zip(ordering, ifes))]

	coord_ordered = []
	
	#append the coordinates based on new ordering
	for index in neworder:
		for key, val in ife_full.iteritems():
			if index[1] == key:
				coord_ordered.append(val)
	

	#function to get the discrepancy based on the new ordering	
	def get(d, first, second):
		return d.get(second, {}).get(first, 0.0)

	index1 = []
	index2 = []
	ife1 = []
	ife2 = []

	for x in neworder:
		for y in neworder:
			index1.append(x[0])
			ife1.append(x[1])
			index2.append(y[0])
			ife2.append(y[1])

	ife_pairs = zip(ife1, ife2)

	disc_ordered = [get(distances, first, second) or get(distances, second, first) for first, second in ife_pairs]

	heatmap_data = [{"ife1": ife1, "ife1_index": ife1_index, "ife2": ife2, "ife2_index": ife2_index, "discrepancy": discrepancy} for ife1, ife1_index, ife2, ife2_index, discrepancy in zip(ife1, index1, ife2, index2, disc_ordered)]
	
	heatmap_data = json.dumps(heatmap_data, ensure_ascii=False)

	#return json.dumps(unitid_r3) 

	return render_template("correspondence_disc.html", query_res1 = unitid1, query_res2 = unitid2, query_res3 = unitid3, ife = neworder, coord = coord_ordered, data = heatmap_data)

if __name__ == '__main__':
	app.run(debug=True)