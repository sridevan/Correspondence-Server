import collections as coll
from flask import Flask, json, jsonify, render_template, render_template_string, request, session
from sqlalchemy import tuple_
from flask_bootstrap import Bootstrap
from flask_sqlalchemy import SQLAlchemy
import numpy as np
import itertools
import time
from numpy import nan

np.set_printoptions(threshold=np.nan)
from sqlalchemy.sql.expression import case

# from models import UnitCorrespondence
# from flask_marshmallow import Marshmallow


app = Flask(__name__)
Bootstrap(app)
app.config[
    'SQLALCHEMY_DATABASE_URI'] = 'mysql+pymysql://root:root@127.0.0.1/rna3dhub?unix_socket=/Applications/MAMP/tmp/mysql/mysql.sock'

db = SQLAlchemy(app)

from models import *
from discrepancy import *
from greedyInsertion import *



@app.route('/')
def home():
    return render_template("welcome.html")


@app.route('/correspondence')
def correspondence():
    # chain_info = '|'.join(unitid.split('|')[:3])
    # print chain_info

    unitid1 = request.args['start']
    unitid2 = request.args['end']
    # unitid3 = request.args['unitid3']

    pdb_test = (('4YBB', 'AA'), ('5JC9', 'AA'), ('4V9D', 'BA'), ('4V9P', 'FA'), ('3JA1', 'SA'))

    pdb_test2= (('4YBB', 'AA'), ('5JC9', 'AA'), ('4WOI', 'AA'), ('4V9P', 'FA'), ('5J7L', 'BA'), ('5KCS', '1a'), ('5KPS', '27'), ('5UYN', 'A'), ('5KPX', '26'), ('3JBU', 'A'), ('5NP6', 'D'), ('3JA1', 'SA'), ('5U4J', 'a'))

    pdb_small = [('4YBB', 'AA'), ('5JC9', 'AA'), ('4WOI', 'AA'), ('4V9P', 'FA'), ('5J7L', 'BA'), ('4V9P', 'DA'), ('5J8A', 'BA'), ('5J91', 'AA'), ('5JC9', 'BA'), ('4V9O', 'DA'), ('4V57', 'CA'), ('4V52', 'AA'), ('4U27', 'AA'), ('4V7V', 'AA'), ('4WF1', 'AA'), ('4U1U', 'AA'), ('4U26', 'AA'), ('4V9D', 'AA'), ('4V7S', 'AA'), ('5IT8', 'BA'), ('4V9D', 'BA'), ('4U27', 'CA'), ('4V53', 'AA'), ('4V7U', 'AA'), ('4V50', 'CA'), ('4WWW', 'XA'), ('4U20', 'AA'), ('4U1U', 'CA'), ('4U1V', 'AA'), ('4V6C', 'CA'), ('4V4Q', 'AA'), ('4V85', 'AA'), ('4V7T', 'CA'), ('4U24', 'CA'), ('4V6D', 'AA'), ('4V89', 'AA'), ('4V6D', 'CA'), ('4V6E', 'AA'), ('4V9C', 'AA'), ('4V6E', 'CA'), ('4V5B', 'DA'), ('4V4H', 'CA'), ('5J5B', 'BA'), ('5NWY', '0'), ('3R8N', 'A'), ('5WDT', 'a'), ('5MDV', '2'), ('6ENU', 'a'), ('6ENF', 'a'), ('5U9F', 'A'), ('5UYM', 'A'), ('5MDZ', '2'), ('4V80', 'AA'), ('5O2R', 'a'), ('5LZD', 'a'), ('6BU8', 'A'), ('5U4I', 'a'), ('5KCR', '1a'), ('3JBV', 'A'), ('5LZA', 'a'), ('5JU8', 'AA'), ('6ENJ', 'a'), ('3JCD', 'a'), ('6DNC', 'A'), ('5UYP', 'A'), ('5KPW', '26'), ('5KCS', '1a'), ('5UYN', 'A'), ('3JBU', 'A'), ('3JA1', 'SA')]

    pdb_large = [('4YBB', 'AA'), ('5J8A', 'AA'), ('5JC9', 'AA'), ('4WOI', 'DA'), ('4WOI', 'AA'), ('5IT8', 'AA'), ('4V9P', 'FA'), ('4V9O', 'BA'), ('5J7L', 'BA'), ('4V9P', 'HA'), ('4V9P', 'DA'), ('4V9P', 'BA'), ('5J8A', 'BA'), ('4V9O', 'FA'), ('5J91', 'AA'), ('4V6C', 'AA'), ('5JC9', 'BA'), ('4V52', 'CA'), ('4V9O', 'DA'), ('4V7T', 'AA'), ('4V57', 'CA'), ('4WWW', 'QA'), ('4V52', 'AA'), ('4V57', 'AA'), ('4U27', 'AA'), ('5J88', 'AA'), ('4V7V', 'AA'), ('4V54', 'CA'), ('4WF1', 'AA'), ('4YBB', 'BA'), ('4U1U', 'AA'), ('4V64', 'CA'), ('4U26', 'AA'), ('4V53', 'CA'), ('4V9D', 'AA'), ('4U24', 'AA'), ('4V7S', 'AA'), ('4U25', 'AA'), ('5IT8', 'BA'), ('4V56', 'CA'), ('4V9D', 'BA'), ('4V54', 'AA'), ('4U27', 'CA'), ('4V64', 'AA'), ('4V53', 'AA'), ('4V56', 'AA'), ('4V7U', 'AA'), ('4V55', 'CA'), ('4V50', 'CA'), ('4V9O', 'HA'), ('4WWW', 'XA'), ('4V4Q', 'CA'), ('4U20', 'AA'), ('4V55', 'AA'), ('4U1U', 'CA'), ('4V50', 'AA'), ('4U1V', 'AA'), ('4U25', 'CA'), ('4V6C', 'CA'), ('5J88', 'BA'), ('4V4Q', 'AA'), ('4U26', 'CA'), ('4V85', 'AA'), ('5J91', 'BA'), ('4V7T', 'CA'), ('4WF1', 'CA'), ('4U24', 'CA'), ('4U1V', 'CA'), ('4V6D', 'AA'), ('4U20', 'CA'), ('4V89', 'AA'), ('4V7S', 'CA'), ('4V6D', 'CA'), ('4V7V', 'CA'), ('4V6E', 'AA'), ('4V9C', 'CA'), ('4V9C', 'AA'), ('4V7U', 'CA'), ('4V6E', 'CA'), ('4V5B', 'BA'), ('4V5B', 'DA'), ('5J5B', 'AA'), ('4V4H', 'CA'), ('4V4H', 'AA'), ('5J5B', 'BA'), ('5AFI', 'a'), ('5NWY', '0'), ('3R8O', 'A'), ('3R8N', 'A'), ('5H5U', 'h'), ('5WDT', 'a'), ('5WFS', 'a'), ('5MDV', '2'), ('5MGP', 'a'), ('6ENU', 'a'), ('5WE4', 'a'), ('6ENF', 'a'), ('5U9G', 'A'), ('5U9F', 'A'), ('3JCE', 'a'), ('5UYM', 'A'), ('5MDW', '2'), ('5MDZ', '2'), ('6C4I', 'a'), ('4V80', 'AA'), ('4V80', 'CA'), ('5O2R', 'a'), ('5WFK', 'a'), ('5LZD', 'a'), ('5WE6', 'a'), ('6BU8', 'A'), ('5MDY', '2'), ('5U4I', 'a'), ('5UYL', 'A'), ('5KCR', '1a'), ('5IQR', '2'), ('3JBV', 'A'), ('5WF0', 'a'), ('5LZA', 'a'), ('5JTE', 'AA'), ('5JU8', 'AA'), ('3JCJ', 'g'), ('6ENJ', 'a'), ('3J9Z', 'SA'), ('3JCD', 'a'), ('5L3P', 'a'), ('6DNC', 'A'), ('5UYQ', 'A'), ('5UYP', 'A'), ('5UYK', 'A'), ('5KPW', '26'), ('3J9Y', 'a'), ('5KCS', '1a'), ('5KPS', '27'), ('5UYN', 'A'), ('5KPX', '26'), ('3JBU', 'A'), ('5NP6', 'D'), ('3JA1', 'SA'), ('5U4J', 'a')]

    # ife1 = '|'.join(unitid1.split('|')[:3])
    pdb = unitid1.split('|')[0]
    chain_info = unitid1.split('|')[2]
    chain_index1 = unitid1.split('|')[-1]
    chain_index2 = unitid2.split('|')[-1]

    # get list of unit ids as specified by the range
    units_list = []
    units_query = UnitInfo.query.filter_by(pdb_id = pdb, chain = chain_info).\
                                      filter(UnitInfo.chain_index.between(chain_index1, chain_index2))\
                                      .order_by(UnitInfo.chain_index).all()
    for row in units_query:
        units_list.append(row.unit_id)

    # query nts as a string
    query_nts = ', '.join(units_list)

    query_len = len(units_list)

    ordering = case(
        {id: index for index, id in enumerate(units_list)},
        value=UnitCorrespondence.unit_id_1
    )
    '''
    # Database query for getting the corresponding units from the query nts
    Need to write a function to get the members from an EC given the query
    The query can be a
    i) A string of unit_ids
    ii) Loop_ids
    iii) Or range of nucleotides
    '''
    correspondence_query = UnitCorrespondence.query.filter(UnitCorrespondence.unit_id_1.in_(units_list))\
                                                   .order_by(ordering)\
                                                   .filter(tuple_(UnitCorrespondence.pdb_id_2, UnitCorrespondence.chain_name_2)\
                                                   .in_(pdb_test2))\
                                                   #.limit(100)


    result = [[unit.unit_id_2 for unit in units] for unit_id_1, units in itertools.groupby(correspondence_query, lambda x: x.unit_id_1)]

    newresult = zip(*result)

    # Create lists for residue type and number
    unit_list = []
    res_num = []
    res_type = []
    for units in newresult:
        unit_list.append(units[0])
        for unit in units:
            res_num.append(unit.split('|')[-1])
            res_type.append(unit.split('|')[-2])
        # ife = '|'.join(units[0].split('|')[:3])
        # unit_list.append(ife)

    res_num_list = [res_num[i:i + query_len] for i in range(0, len(res_num), query_len)]
    res_type_list = [res_type[i:i + query_len] for i in range(0, len(res_type), query_len)]

    # res_data = [list(i) for i in zip(res_type_list, res_num_list)]

    res_info_list = []
    for base in res_type_list:
        for num in res_num_list:
            res_info = ["{}{}".format(k, v) for k, v in zip(base, num)]
            res_info_list.append(res_info)

    # Create list of IFES
    ife_list = []
    for elem in unit_list:
        ife = '|'.join(elem.split('|')[:3])
        ife_list.append(ife)

    # Create list of coordinates as strings
    coord_unordered = []
    for x in newresult:
        x = ','.join(x)
        coord_unordered.append(x)

    # Create a dictionary of ifes with coordinate data
    ife_coord = dict(zip(ife_list, coord_unordered))

    # Create list of res number as strings
    res_num_ordered = []
    for res in res_num_list:
        res = ','.join(res)
        res_num_ordered.append(res)

    # Create list to store the centers np array
    units_center = []
    units_num_center = []

    # This section of the code deals with the database query to get the centers data
    for units in newresult:

        ordering = case(
            {id: index for index, id in enumerate(units)},
            value=UnitCenters.unit_id
        )

        centers_query = UnitCenters.query.filter(UnitCenters.unit_id.in_(units),
                                                 UnitCenters.name == 'base').order_by(ordering)
        for row in centers_query:
            units_center.append(np.array([row.x, row.y, row.z]))
            units_num_center.append(row.unit_id)

    units_center_list = [units_center[i:i+query_len] for i in range(0, len(units_center), query_len)]

    # Create list to store the rotation np array
    units_rotation = []
    units_num_rotation = []

    # This section of the code deals with the database query to get the rotation data
    for units in newresult:

        ordering = case(
            {id: index for index, id in enumerate(units)},
            value=UnitRotations.unit_id
        )

        rotation_query = UnitRotations.query.filter(UnitRotations.unit_id.in_(units)).order_by(ordering)

        for row in rotation_query:
            units_rotation.append(np.array([[row.cell_0_0, row.cell_0_1, row.cell_0_2],
                                            [row.cell_1_0, row.cell_1_1, row.cell_1_2],
                                            [row.cell_2_0, row.cell_2_1, row.cell_2_2]]))
            units_num_rotation.append(row.unit_id)

    units_rotation_list = [units_rotation[i:i+query_len] for i in range(0, len(units_rotation), query_len)]

    start_time = time.time()

    # This section of the code deals with calculating the discrepancy for the corresponding instances
    distances = coll.defaultdict(lambda: coll.defaultdict(int))

    for a in range(0, len(ife_list)):
        for b in range(0, len(ife_list)):
        #for b in range(a + 1, len(ife_list)):
            disc = matrix_discrepancy(units_center_list[a], units_rotation_list[a], units_center_list[b], units_rotation_list[b])
            ife_a = '|'.join(ife_list[a].split('|')[:3])
            ife_b = '|'.join(ife_list[b].split('|')[:3])
            distances [ife_a][ife_b] = disc

    end_time = time.time()

    disc_time = end_time - start_time


    dist = np.zeros((len(ife_list), len(ife_list)))
    for index1, member1 in enumerate(ife_list):

        curr = distances.get(member1, {})
        for index2, member2 in enumerate(ife_list):
            val = curr.get(member2, None)
            if member2 not in curr:
                val = None
            dist[index1, index2] = val

    ordering, _, _ = orderWithPathLengthFromDistanceMatrix(dist, 10, scanForNan=True)

    # Order the list of ifes based on the new ordering
    ifes_ordered = [x for x in sorted(zip(ordering, ife_list))]

    coord_ordered = []
    # append the coordinates based on new ordering
    for index in ifes_ordered:
        for key, val in ife_coord.iteritems():
            if index[1] == key:
                coord_ordered.append(val)

    # function to get the discrepancy based on the new ordering
    def get(d, first, second):
        return d.get(second, {}).get(first, 0.0)

    index1 = []
    index2 = []
    ife1 = []
    ife2 = []

    for member1 in ifes_ordered:
        for member2 in ifes_ordered:
            index1.append(member1[0])
            ife1.append(member1[1])
            index2.append(member2[0])
            ife2.append(member2[1])

    ife_pairs = zip(ife1, ife2)

    disc_ordered = [get(distances, first, second) or get(distances, second, first) for first, second in ife_pairs]

    heatmap_data = [
        {"ife1": ife1, "ife1_index": ife1_index, "ife2": ife2, "ife2_index": ife2_index, "discrepancy": discrepancy}
        for ife1, ife1_index, ife2, ife2_index, discrepancy in zip(ife1, index1, ife2, index2, disc_ordered)
        ]

    heatmap_data = json.dumps(heatmap_data, ensure_ascii=False)

    return render_template("correspondence_disc.html", query_pdb=pdb, disc_time=disc_time, query_nts=query_nts, coord=coord_ordered, ifes=ifes_ordered, res_list=coord_ordered, data=heatmap_data)
    # return json.dumps(disc_time)


if __name__ == '__main__':
    app.run(debug=True)
