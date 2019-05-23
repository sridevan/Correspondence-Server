import collections as coll
from flask import Flask, json, jsonify, render_template, render_template_string, request, session
from sqlalchemy import tuple_
from flask_bootstrap import Bootstrap
from flask_sqlalchemy import SQLAlchemy
import numpy as np
import itertools
import random
import time
import math
from numpy import nan
#import scipy

np.set_printoptions(threshold=np.nan)
#from sqlalchemy.sql.expression import case
from sqlalchemy import case
#from sqlalchemy.orm import relationship
# from models import UnitCorrespondence
# from flask_marshmallow import Marshmallow

app = Flask(__name__)
app.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = False
#app.config['SQLALCHEMY_DATABASE_URI'] = 'mysql://webfr3d:nrw0FhuKwY2CUYa2TDPU@localhost/rna3dhub-prod'
app.config['SQLALCHEMY_DATABASE_URI'] = 'mysql+pymysql://root:root@127.0.0.1/rna3dhub?unix_socket=/Applications/MAMP/tmp/mysql/mysql.sock'
Bootstrap(app)
db = SQLAlchemy(app)

from models import *
from discrepancy import *
from greedyInsertion import *
from process_input import *
from ordering import *
from queries import *


@app.route('/')
def home():
    # Debug statement
    return render_template("home.html")

@app.route('/correspondence')
def correspondence():
    # chain_info = '|'.join(unitid.split('|')[:3])
    # print chain_info

    pdb_test = [('3JCD', 'a'), ('3R8N', 'A'), ('3JCJ', 'g'), ('5NP6', 'D'), ('6H4N', 'a')]

    data = request.args['units']

    query_list = input_type(data)

    query_ife = '|'.join(query_list[0][0].split('|')[:3])
    query_pdb = query_list[0][0].split('|')[0]
    query_chain = query_list[0][0].split('|')[2]

    reject_list = []

#######################################################################################################
    
    #This section of the code deals with getting the members of Equivalence Class from the query chain
    ife_list = NrChains.query.join(NrClasses, NrReleases)\
        .filter(NrChains.ife_id == query_ife).filter(NrClasses.resolution == '4.0')\
        .order_by(NrReleases.date.desc()).limit(1) 
       
    for row in ife_list:
        class_id = row.nr_class_id

    ec_query = NrClasses.query.filter_by(nr_class_id=class_id)

    for row in ec_query:
        equivalence_class = row.name
        nr_release = row.nr_release_id

    members_query = NrChains.query.filter_by(nr_class_id=class_id)

    ife_members = []
    for row in members_query:
        ife_members.append(row.ife_id)

    rejected_ife = []

    for i, v in enumerate(ife_members):
        if any(c in '+' for c in v):
            rejected_ife.append(ife_members[i])
            del ife_members[i]
        else:
            pass

    members_pdb = []
    members_chain = []

    for ife in ife_members:
        members_pdb.append(ife.split('|')[0])
        members_chain.append(ife.split('|')[-1])

    members_info = zip(members_pdb, members_chain)

#####################################################################################################
    units_complete_list = []

    for range_num in query_list:
        start_range = range_num[0].split('|')[-1]
        stop_range = range_num[1].split('|')[-1]
        units_query = UnitInfo.query.filter_by(pdb_id=query_pdb, chain=query_chain). \
            filter(UnitInfo.chain_index.between(start_range, stop_range)) \
            .order_by(UnitInfo.chain_index).all()

        for row in units_query:
            units_complete_list.append(row.unit_id)

    # query nts as a string
    query_nts = ', '.join(units_complete_list)

    query_complete_len = len(units_complete_list)

#####################################################################################################

    #### This section deals with getting the units of unmodified nucleotides

    standard_nts = ('A', 'C', 'G', 'U')

    units_std_list = []

    for unit in units_complete_list:
        k = unit.split('|')[-2]
        if k in standard_nts:
            units_std_list.append(unit)

    query_std_len= len(units_std_list)

#####################################################################################################

#### This section of the code deals with getting the complete corresponding unit_ids 

    correspondence_complete = UnitCorrespondence.query.filter(UnitCorrespondence.unit_id_1.in_(units_complete_list)) \
        .filter(tuple_(UnitCorrespondence.pdb_id_2, UnitCorrespondence.chain_name_2) \
        .in_(members_info[:100])) 

    result_complete = [[unit.unit_id_2 for unit in units] for unit_id_1, units in
              itertools.groupby(correspondence_complete, lambda x: x.unit_id_1)]

    corr_complete = zip(*result_complete)

    # Create lists for residue type and number
    unit_list = []
    res_num = []
    res_type = []
    for units in corr_complete:
        unit_list.append(units[0])
        #for unit in units:
            #res_num.append(unit.split('|')[-1])
            #res_type.append(unit.split('|')[-2])
        # ife = '|'.join(units[0].split('|')[:3])
        # unit_list.append(ife)

    # res_num_list = [res_num[i:i + query_len] for i in range(0, len(res_num), query_len)]
    # res_type_list = [res_type[i:i + query_len] for i in range(0, len(res_type), query_len)]
    # res_list = [res_num[i:i + query_complete_len] for i in xrange(0, len(res_num), query_complete_len)]

    # Create list of IFES
    ife_list = []
    for elem in unit_list:
        ife = '|'.join(elem.split('|')[:3])
        ife_list.append(ife)

    # Create list of coordinates as strings
    coord_unordered = []
    for x in corr_complete:
        x = ','.join(x)
        coord_unordered.append(x)

    # Create a dictionary of ifes with coordinate data
    ife_coord = dict(zip(ife_list, coord_unordered))

##################################################################################

#### This section of the code deals with getting the corresponding unit_ids without modified nucleotides 

    #ordering = case({id: index for index, id in enumerate(units_std)}, value=UnitCorrespondence.unit_id_1)

    correspondence_std = UnitCorrespondence.query.filter(UnitCorrespondence.unit_id_1.in_(units_std_list)) \
        .filter(tuple_(UnitCorrespondence.pdb_id_2, UnitCorrespondence.chain_name_2) \
        .in_(members_info[:100]))

    result_std = [[unit.unit_id_2 for unit in units] for unit_id_1, units in
              itertools.groupby(correspondence_std, lambda x: x.unit_id_1)]

    corr_std = zip(*result_std)

##################################################################################

#### Queries to obtain center and rotation data for calculating discrepacy

    # Create list to store the centers np array
    units_center = []
    units_num_center = []

    # This section of the code deals with the database query to get the centers data
    for units in corr_std:

        #ordering = case(
            #{id: index for index, id in enumerate(units)},
            #value=UnitCenters.unit_id
        #)

        centers_query = UnitCenters.query.filter(UnitCenters.unit_id.in_(units),
                                                 UnitCenters.name == 'base')
                                                #.order_by(ordering)
        for row in centers_query:
            units_center.append(np.array([row.x, row.y, row.z]))
            units_num_center.append(row.unit_id)

    units_center_list = [units_center[i:i + query_std_len] for i in xrange(0, len(units_center), query_std_len)]

    # Create list to store the rotation np array
    units_rotation = []
    units_num_rotation = []

    # This section of the code deals with the database query to get the rotation data
    for units in corr_std:

        #ordering = case(
            #{id: index for index, id in enumerate(units)},
            #value=UnitRotations.unit_id
        #)

        rotation_query = UnitRotations.query.filter(UnitRotations.unit_id.in_(units))

        for row in rotation_query:
            units_rotation.append(np.array([[row.cell_0_0, row.cell_0_1, row.cell_0_2],
                                            [row.cell_1_0, row.cell_1_1, row.cell_1_2],
                                            [row.cell_2_0, row.cell_2_1, row.cell_2_2]]))
            units_num_rotation.append(row.unit_id)

    units_rotation_list = [units_rotation[i:i + query_std_len] for i in xrange(0, len(units_rotation), query_std_len)]

    rotation_size = len(units_rotation_list)

####################################################################################

#### Build a distance matrix using the queries from above and order them based on
#### similarity

    # This section of the code deals with calculating the discrepancy for the corresponding instances
    distances = coll.defaultdict(lambda: coll.defaultdict(int))

    for a in range(0, len(ife_list)):
        for b in range(a+1, len(ife_list)):
            disc = matrix_discrepancy(units_center_list[a], units_rotation_list[a], units_center_list[b],
                                      units_rotation_list[b])
            distances[ife_list[a]][ife_list[b]] = disc

    # Empty list to append pairs of IFE with NaN discrepancy
    ife_nan = []

    for k, v in distances.items():
        for a, b in v.items():
            if math.isnan(b):
                ife_nan.append((k, a))
                v[a] = 1

    dist = np.zeros((len(ife_list), len(ife_list)))
    for index1, member1 in enumerate(ife_list):
        curr = distances.get(member1, {})
        for index2, member2 in enumerate(ife_list):
            dist[index1, index2] = curr.get(member2, 0)       

    dist = (dist + np.swapaxes(dist, 0, 1))

    #ordering, _, _ = orderWithPathLengthFromDistanceMatrix(dist, 10, scanForNan=True)
    disc_order = optimalLeafOrder(dist)

    new_order = []
    idx_new = []

    for idx, order in enumerate(disc_order):
        new_order.append(ife_list[order])
        idx_new.append(idx)

    ifes_ordered = zip(idx_new, new_order)

    # Order the list of ifes based on the new ordering
    #ifes_ordered = [x for x in sorted(zip(disc_order, ife_list))]

    coord_ordered = []
    # append the coordinates based on new ordering
    for index in ifes_ordered:
        for key, val in ife_coord.iteritems():
            if index[1] == key:
                coord_ordered.append(val)

#########################################################################################

#### Logic to order and build the heatmap data

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
    
###########################################################################################     

    return render_template("correspondence_display.html", query_pdb=query_pdb, query_nts=query_nts,
                          coord=coord_ordered, ifes=ifes_ordered, res_list=coord_ordered, 
                          ec=equivalence_class, release=nr_release, data=heatmap_data)

if __name__ == '__main__':
    app.run(debug=True)