import collections as coll
import csv
from collections import OrderedDict, defaultdict
from copy import deepcopy
from flask import Flask, json, jsonify, render_template, render_template_string, request, session
from sqlalchemy import tuple_
from flask_bootstrap import Bootstrap
from flask_sqlalchemy import SQLAlchemy
import numpy as np
import itertools
from itertools import combinations, groupby
import random
import time
import math
from numpy import nan
import StringIO
from definitions import ribosome_subunits
# import scipy

np.set_printoptions(threshold=np.nan)
# from sqlalchemy.sql.expression import case
from sqlalchemy import case

# from sqlalchemy.orm import relationship
# from models import UnitCorrespondence
# from flask_marshmallow import Marshmallow

app = Flask(__name__)
app.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = False
# app.config['SQLALCHEMY_DATABASE_URI'] = 'mysql://webfr3d:nrw0FhuKwY2CUYa2TDPU@localhost/rna3dhub-prod'
app.config[
    'SQLALCHEMY_DATABASE_URI'] = 'mysql+pymysql://root:root@127.0.0.1/rna3dhub?unix_socket=/Applications/MAMP/tmp/mysql/mysql.sock'
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

    annotation = [('4V4Q|1|AA', 'Empty', 'NO', '-', '-'), ('4V4Q|1|CA', 'Empty', 'NO', '-', '-'),
                  ('4V50|1|AA', 'P-site of SSU only', 'NO', '-', '-'),
                  ('4V50|1|CA', 'P-site of SSU only', 'NO', '-', '-'), ('4V5B|1|BA', 'Empty', 'NO', '-', '-'),
                  ('4V5B|1|DA', 'Empty', 'NO', '-', '-'),
                  ('4V9D|1|AA', 'P/E', 'Ribosome Recycling (early intermediate)', 'RRF (AY)', '-'),
                  ('4V9D|1|BA', 'P/P', 'Posttermination', 'none', '-'),
                  ('4V9O|1|BA', 'Empty', 'Elongation - pretranslocation', 'EF-G (BV)', 'Viomycin (BW)'),
                  ('4V9O|1|DA', 'Empty', 'Elongation - pretranslocation', 'EF-G (DV)', 'Viomycin (DW)'),
                  ('4V9O|1|FA', 'Empty', 'Elongation - pretranslocation', 'EF-G (FV)', 'Viomycin (FW)'),
                  ('4V9O|1|HA', 'Empty', 'Elongation - pretranslocation', 'EF-G (HV)', 'Viomycin (HW)'),
                  ('4V9P|1|BA', 'Empty', 'Elongation - pretranslocation', 'EF-G (BV)', 'Viomycin (BW)'),
                  ('4V9P|1|DA', 'Empty', 'Elongation - pretranslocation', 'EF-G (DV)', 'Viomycin (DW)'),
                  ('4V9P|1|FA', 'Empty', 'Elongation - pretranslocation', 'EF-G (FV)', 'Viomycin (FW)'),
                  ('4V9P|1|HA', 'Empty', 'Elongation - pretranslocation', 'EF-G (HV)', '-'),
                  ('4YBB|1|AA', 'Empty', 'NO', '-', '-'), ('4YBB|1|BA', 'Empty', 'NO', '-', '-'),
                  ('5IT8|1|AA', 'Empty', 'NO', '-', '-'), ('5IT8|1|BA', 'Empty', 'NO', '-', '-'),
                  ('5J5B|1|AA', 'Empty', 'Bound to antibiotic', '-', 'Tetracycline (AA)'),
                  ('5J5B|1|BA', 'Empty', 'Bound to antibiotic', '-', 'Tetracycline (BA)'),
                  ('5J7L|1|AA', 'Empty', 'Bound to antibiotic', '-', 'Tetracycline (AA)'),
                  ('5J7L|1|BA', 'Empty', 'Bound to antibiotic', '-', 'Tetracycline (BA)'),
                  ('5J88|1|AA', 'Empty', 'NO', '-', '-'), ('5J88|1|BA', 'Empty', 'NO', '-', '-'),
                  ('5J8A|1|AA', 'Empty', 'Bound to antibiotic', '-', 'Tigecycline (AA)'),
                  ('5J8A|1|BA', 'Empty', 'Bound to antibiotic', '-', 'Tigecycline (BA)'),
                  ('5J91|1|AA', 'Empty', 'Bound to antibiotic', '-', 'Tigecycline (AA)'),
                  ('5J91|1|BA', 'Empty', 'Bound to antibiotic', '-', 'Tigecycline (BA)'),
                  ('5JC9|1|AA', 'Empty', 'NO', '-', '-'), ('5JC9|1|BA', 'Empty', 'NO', '-', '-'),
                  ('5MDZ|1|2', 'P-site ', 'NO', '-', '-'), ('6BU8|1|A', 'Appears to be A/A, P/P, E/E ', 'NO', '-', '-'),
                  ('6GWT|1|a', 'P/RF1 ', 'Termination', 'RF1 (v), RF3 (w)', 'Apidaecin (z)'),
                  ('6GXM|1|a', 'P/RF1 ', 'Termination', 'RF1 (v), RF3 (w)', 'Apidaecin (z)'),
                  ('6GXN|1|a', 'P/RF1', 'Termination', 'RF1 (v), RF3 (w)', 'Apidaecin (z)'),
                  ('6GXO|1|a', 'P/E', 'Termination', 'RF1 (v), RF3 (w)', 'Apidaecin (z)'),
                  ('6I7V|1|AA', 'Empty', 'Ribosome heterogeneity ', '-', '-'),
                  ('6I7V|1|BA', 'Empty', 'Ribosome heterogeneity ', '-', '-'),
                  ('4U1U|1|AA', 'Empty', 'Bound to antibiotic', '-', 'Quinupristin (B6)'),
                  ('4U1U|1|CA', 'Empty', 'Bound to antibiotic', '-', 'Quinupristin (D6)'),
                  ('4U1V|1|AA', 'Empty', 'Bound to antibiotic', '-', 'Linopristin (B6)'),
                  ('4U1V|1|CA', 'Empty', 'Bound to antibiotic', '-', 'Linopristin (D6)'),
                  ('4U20|1|AA', 'Empty', 'Bound to antibiotic', '-', 'Flopristin (BA)'),
                  ('4U20|1|CA', 'Empty', 'Bound to antibiotic', '-', 'Flopristin (DA)'),
                  ('4U24|1|AA', 'Empty', 'Bound to antibiotic', '-', 'Dalfopristin (BA)'),
                  ('4U24|1|CA', 'Empty', 'Bound to antibiotic', '-', 'Dalfopristin (DA)'),
                  ('4U25|1|AA', 'Empty', 'Bound to antibiotic', '-', 'Virginiamycin (BA)'),
                  ('4U25|1|CA', 'Empty', 'Bound to antibiotic', '-', 'Virginiamycin (DA)'),
                  ('4U26|1|AA', 'Empty', 'Bound to antibiotic', '-', 'Dalfopristin (B4) and Quinupristin (B6)'),
                  ('4U26|1|CA', 'Empty', 'Bound to antibiotic', '-', 'Dalfopristin (D4) and Quinupristin (D6)'),
                  ('4U27|1|AA', 'Empty', 'Bound to antibiotic', '-', 'Flopristin (BA) and Linopristin (B6)'),
                  ('4U27|1|CA', 'Empty', 'Bound to antibiotic', '-', 'Flopristin (DA) and Linopristin (D6)'),
                  ('4V4H|1|AA', 'Empty', 'Bound to antibiotic', '-', 'Kasugamycin (AA)'),
                  ('4V4H|1|CA', 'Empty', 'Bound to antibiotic', '-', 'Kasugamycin (CA)'),
                  ('4V52|1|AA', 'Empty', 'Bound to antibiotic', '-', 'Neomycin (AA & BB)'),
                  ('4V52|1|CA', 'Empty', 'Bound to antibiotic', '-', 'Neomycin (CA & DB)'),
                  ('4V53|1|AA', 'Empty', 'Bound to antibiotic', '-', 'Gentamycin (AA & BB)'),
                  ('4V53|1|CA', 'Empty', 'Bound to antibiotic', '-', 'Gentamycin (CA & DB)'),
                  ('4V54|1|AA', 'Empty', 'Ribosome recycling', 'RRF (B6)', '-'),
                  ('4V54|1|CA', 'Empty', 'Ribosome recycling', 'RRF (D6)', '-'),
                  ('4V55|1|AA', 'Empty', 'Ribosome recycling', 'RRF (B6)', 'Gentamycin (AA & BB)'),
                  ('4V55|1|CA', 'Empty', 'Ribosome recycling', 'RRF (D6)', 'Gentamycin (CA & DB)'),
                  ('4V56|1|AA', 'Empty', 'Bound to antibiotic', '-', 'Spectinomycin (AA)'),
                  ('4V56|1|CA', 'Empty', 'Bound to antibiotic', '-', 'Spectinomycin (CA)'),
                  ('4V57|1|AA', 'Empty', 'Bound to antibiotic', '-', 'Spectinomycin (AA) and Neomycin (AA & BB)'),
                  ('4V57|1|CA', 'Empty', 'Bound to antibiotic', '-', 'Spectinomycin (CA) and Neomycin (CA & DB)'),
                  ('4V64|1|AA', 'Empty', 'Bound to antibiotic', '-', 'Hygromycin B (AA)'),
                  ('4V64|1|CA', 'Empty', 'Bound to antibiotic', '-', 'Hygromycin B (CA)'),
                  ('4V7S|1|AA', 'Empty', 'Bound to antibiotic', '-', 'Telithromycin (BA)'),
                  ('4V7S|1|CA', 'Empty', 'NB', '-', '-'),
                  ('4V7T|1|AA', 'Empty', 'Bound to antibiotic', '-', 'Chloramphenicol (BA)'),
                  ('4V7T|1|CA', 'Empty', 'NB', '-', '-'),
                  ('4V7U|1|AA', 'Empty', 'Bound to antibiotic', '-', 'Erythromycin A (BA)'),
                  ('4V7U|1|CA', 'Empty', 'NB', '-', '-'),
                  ('4V7V|1|AA', 'Empty', 'Bound to antibiotic', '-', 'Clindamycin'),
                  ('4V7V|1|CA', 'Empty', 'NB', '-', '-'), ('4V80|1|AA', '?', '?', '?', '?'),
                  ('4V80|1|CA', '?', '?', '?', '?'), ('4V9C|1|AA', 'P/P', 'NO', '-', 'neomycin (AA & BA)'),
                  ('4V9C|1|CA', 'P/E', 'NO', 'RRF(CY)', 'neomycin (CA & DA)'), ('4WF1|1|AA', '-', 'NB', '-', '-'),
                  ('4WF1|1|CA', '-', 'Bound to antibiotic', '-', 'Negamycin (CA)'),
                  ('4WOI|1|AA', 'P/E', 'Bound to antibiotic', '-', 'Paromomycin (AA & BA)'),
                  ('4WOI|1|DA', 'P/P', 'Bound to antibiotic', '-', 'Paromomycin (CA & DA)'),
                  ('4WWW|1|QA', 'Empty', 'Bound to antibiotic', '-', 'CEM-101 (RA)'),
                  ('4WWW|1|XA', 'Empty', 'NO', '-', '-'),
                  ('5KCR|1|1a', 'P/P ', 'Bound to antibiotic', '-', 'Avilamycin C (1A)'),
                  ('5KCS|1|1a', 'P/P ', 'Bound to antibiotic', 'TetM (1w)', 'Evernimycin (1A)'),
                  ('3J9Y|1|a', 'P/P', 'NO', 'TetM (w)', '-'),
                  ('3JCD|1|a', 'P/4 (8),  E/E (9)', 'Elongation - Backtranslocation (Post EF4) ', 'EF4 (x)', '-'),
                  ('3JCE|1|a', 'A/4 (6), P/4 (8), E/E (9)', 'Elongation - Backtranslocation (Pre EF4)', 'EF4 (x)', '-'),
                  ('5AFI|1|a', 'A/T (y), P/P (v), E/E (w)', 'Elongation - Decoding', 'EF-Tu (z)', 'Kirromycin (z)'),
                  ('5UYK|1|A', 'T tRNA (Y), P/P (W), E/E (X)', 'Elongation - Decoding', 'EF-Tu (Z)', '-'),
                  ('5UYL|1|A', 'A*/T (Y), P/P (W), E/E (X)', 'Elongation - Decoding', 'EF-Tu (Z)', '-'),
                  ('5UYM|1|A', 'A/T (Y), P/P (W), E/E (X)', 'Elongation - Decoding', 'EF-Tu (Z)', '-'),
                  ('5UYN|1|A', 'T tRNA (Y), P/P (W), E/E (X)', 'Elongation - Decoding', 'EF-Tu (Z)', '-'),
                  ('5UYP|1|A', 'A*/T (Y), P/P (W), E/E (X)', 'Elongation - Decoding', 'EF-Tu (Z)', '-'),
                  ('5UYQ|1|A', 'A/T (Y), P/P (W), E/E (X)', 'Elongation - Decoding', 'EF-Tu (Z)', '-'),
                  ('5WDT|1|a', 'A/T (y), P/P (v), E/E (w)', 'Elongation - Decoding', 'EF-Tu H84A (z)', '-'),
                  ('5WE4|1|a', 'A/T (y), P/P (v), E/E (w)', 'Elongation - Decoding', 'EF-Tu (z)', '-'),
                  ('5WE6|1|a', 'A/T (y), P/P (v), E/E (w)', 'Elongation - Decoding', 'EF-Tu H84A (z)', '-'),
                  ('5WF0|1|a', 'A/T (y), P/P (v), E/E (w)', 'Elongation - Decoding', 'EF-Tu H84A (z)', '-'),
                  ('5WFK|1|a', 'A/T (y), P/P (v), E/E (w)', 'Elongation - Decoding', 'EF-Tu H84A (z)', '-'),
                  ('5WFS|1|a', 'A/T (y), P/P (v), E/E (w)', 'Elongation - Decoding', 'EF-Tu H84A (z)', '-'),
                  ('3J9Z|1|SA', 'P/P (S6), E/E (S7)', 'Elongation - Translocation', 'EF-G (S1)', '-'),
                  ('3JA1|1|SA', 'P/E (S2)', 'Elongation - Translocation', 'EF-G (S3)', '-'),
                  ('3R8N|1|A', '?', 'Elongation - Translocation', '?', '-'),
                  ('3R8O|1|A', '?', 'Elongation - Translocation', '?', '-'),
                  ('3JCJ|1|g', 'P/I* ', 'Initiation', 'IF-2 (f)', '-'),
                  ('6DNC|1|A', 'E/E (D), P/P (LA)', 'Premature Termination', 'RF1 (MA)', '-'),
                  ('5NP6|1|D', 'P/P (B)', 'Recoding - Bypassing', '-', '-'), (
                      '6H4N|1|a', 'E/E (w)', 'Ribosome hibernation',
                      'Ribosome modulation factor (v), ribosome hibernation promoting factor (x)', '-'),
                  ('5H5U|1|h', 'P/P (5)', 'Ribosome Rescue', 'ArfA (3), RF2 (4)', '-'),
                  ('5MDV|1|2', 'P/P (5)', 'Ribosome Rescue', 'ArfA (6), RF2 (7)', '-'),
                  ('5MDW|1|2', 'P/P (5)', 'Ribosome Rescue', 'ArfA (6), RF2 (7)', '-'),
                  ('5MDY|1|2', 'P/P (5)', 'Ribosome Rescue', 'ArfA (6), RF2 (7)', '-'),
                  ('5MGP|1|a', 'P/P (x)', 'Ribosome Rescue', 'ArfA (w), RF2 (z)', '-'),
                  ('5U4I|1|a', 'P/P (x), E*/E (y)', 'Ribosome Rescue', 'ArfA (w), RF2 (v)', '-'),
                  ('5U4J|1|a', '?', 'Ribosome Rescue', 'ArfA (w), RF2 (v)', '-'),
                  ('5U9F|1|A', 'P/P(W), E*/E (X)', 'Ribosome Rescue', 'ArfA (Y), RF2 (Z)', '-'),
                  ('5U9G|1|A', 'P/P(W), E*/E (X)', 'Ribosome Rescue', 'ArfA (Y), RF2 (Z)', '-'),
                  ('6ENF|1|a', 'P/P (x)', 'Ribosome Rescue', '-', '-'),
                  ('6ENJ|1|a', 'A/A (9), P/P (x)', 'Ribosome Rescue', 'EF-P (w)', '-'),
                  ('6ENU|1|a', 'P/P (x)', 'Ribosome Rescue', 'EF-P (w)', '-'),
                  ('6C4I|1|a', 'P/P (x), E*/E (y)', 'Ribosome rescue', 'ArfA (w), RF-2 (v)', '-'),
                  ('3JBU|1|A', 'P/P (v)', 'Ribosome Stalling', 'SecM (z)', '-'), (
                      '3JBV|1|A', 'A/P* (V), P/E intermediate (W)', 'Ribosome Stalling', 'SecM (z)',
                      'Chloroamphenicol (b)'), ('5JTE|1|AA', 'A/A (AW), P/P (AX), E/E* (AY)', 'Ribosome Stalling', '-',
                                                'Erythromycin A (BA), ErmBL (B5)'),
                  ('5JU8|1|AA', 'P/P (AX), E*/E* (AY)', 'Ribosome Stalling', '-', 'Erythromycin A (BA), ErmBL (B5)'),
                  ('5NWY|1|0', 'P/P-VemP (M)', 'Ribosome Stalling', '-', 'VemP (s)'),
                  ('5O2R|1|a', 'P/P (x)', 'Ribosome Stalling', 'RF-1 (v)', 'Apidaecin (z)'),
                  ('5LZA|1|a', 'P/P (v)', 'SelB activation', '-', '-'),
                  ('5LZD|1|a', 'A/SelB (y), P/P (v)', 'SelB activation', 'SelB (z)', '-'),
                  ('5LZE|1|a+ 5LZE|1|y+ 5LZE|1|v+ 5LZE|1|x', 'tst', 'SelB activation', 'tst', 'tst'),
                  ('5IQR|1|2', 'A/R (6), P/P (5), E/E (4)', 'Stringent Control', 'RelA (8)', 'Paromomycine (2)'),
                  ('5KPS|1|27', 'P/P (31), E/E (32)', 'Stringent Control', 'RelA (A)', '-'),
                  ('5KPW|1|26', 'A/R (30), P/P (31), E/E (32)', 'Stringent Control', 'RelA (33)', '-'),
                  ('5KPX|1|26', 'A/R (30), P/P (31), E/E (32)', 'Stringent Control', 'RelA (33)', '-'),
                  ('5L3P|1|a', 'A/R( y), P/P (x)', 'Stringent Control', 'RelA (z)', '-'),
                  ('4V85|1|AA', '-', 'Termination', 'RF3 (AW)', 'Viomycin (AY)'),
                  ('4V89|1|AA', '-', 'Termination', 'RF3 (AW)', '-'),
                  ('4V6C|1|AA', '-', 'Translation - elongation', '-', '-'),
                  ('4V6C|1|CA', '-', 'Translation - elongation', '-', '-'),
                  ('4V6D|1|AA', 'P-site ASL fragment', 'Translation - elongation', '-', '-'),
                  ('4V6D|1|CA', 'P-site ASL fragment', 'Translation - elongation', '-', '-'), (
                      '4V6E|1|AA', 'A-site ASL fragment (AX), P-site ASL fragment (AV)', 'Translation - elongation',
                      '-',
                      '-'), (
                      '4V6E|1|CA', 'A-site ASL fragment (CX), P-site ASL fragment (CV)', 'Translation - elongation',
                      '-',
                      '-')]

    annotation_LSU = [('4V4Q|1|BB', 'Empty', 'NO', '-', '-'), ('4V4Q|1|DB', 'Empty', 'NO', '-', '-'),
                      ('4V50|1|BB', 'P-site of SSU only', 'NO', '-', '-'),
                      ('4V50|1|DB', 'P-site of SSU only', 'NO', '-', '-'), ('4V5B|1|AB', 'Empty', 'NO', '-', '-'),
                      ('4V5B|1|CB', 'Empty', 'NO', '-', '-'),
                      ('4V9D|1|CA', 'P/E', 'Ribosome Recycling (early intermediate)', 'RRF (AY)', '-'),
                      ('4V9D|1|DA', 'P/P', 'Posttermination', 'none', '-'),
                      ('4V9O|1|AA', 'Empty', 'Elongation - pretranslocation', 'EF-G (BV)', 'Viomycin (BW)'),
                      ('4V9O|1|CA', 'Empty', 'Elongation - pretranslocation', 'EF-G (DV)', 'Viomycin (DW)'),
                      ('4V9O|1|EA', 'Empty', 'Elongation - pretranslocation', 'EF-G (FV)', 'Viomycin (FW)'),
                      ('4V9O|1|GA', 'Empty', 'Elongation - pretranslocation', 'EF-G (HV)', 'Viomycin (HW)'),
                      ('4V9P|1|AA', 'Empty', 'Elongation - pretranslocation', 'EF-G (BV)', 'Viomycin (BW)'),
                      ('4V9P|1|CA', 'Empty', 'Elongation - pretranslocation', 'EF-G (DV)', 'Viomycin (DW)'),
                      ('4V9P|1|EA', 'Empty', 'Elongation - pretranslocation', 'EF-G (FV)', 'Viomycin (FW)'),
                      ('4V9P|1|GA', 'Empty', 'Elongation - pretranslocation', 'EF-G (HV)', '-'),
                      ('4YBB|1|DA', 'Empty', 'NO', '-', '-'), ('4YBB|1|CA', 'Empty', 'NO', '-', '-'),
                      ('5IT8|1|DA', 'Empty', 'NO', '-', '-'), ('5IT8|1|CA', 'Empty', 'NO', '-', '-'),
                      ('5J5B|1|DA', 'Empty', 'Bound to antibiotic', '-', 'Tetracycline (AA)'),
                      ('5J5B|1|CA', 'Empty', 'Bound to antibiotic', '-', 'Tetracycline (BA)'),
                      ('5J7L|1|DA', 'Empty', 'Bound to antibiotic', '-', 'Tetracycline (AA)'),
                      ('5J7L|1|CA', 'Empty', 'Bound to antibiotic', '-', 'Tetracycline (BA)'),
                      ('5J88|1|DA', 'Empty', 'NO', '-', '-'), ('5J88|1|CA', 'Empty', 'NO', '-', '-'),
                      ('5J8A|1|DA', 'Empty', 'Bound to antibiotic', '-', 'Tigecycline (AA)'),
                      ('5J8A|1|CA', 'Empty', 'Bound to antibiotic', '-', 'Tigecycline (BA)'),
                      ('5J91|1|DA', 'Empty', 'Bound to antibiotic', '-', 'Tigecycline (AA)'),
                      ('5J91|1|CA', 'Empty', 'Bound to antibiotic', '-', 'Tigecycline (BA)'),
                      ('5JC9|1|DA', 'Empty', 'NO', '-', '-'), ('5JC9|1|CA', 'Empty', 'NO', '-', '-'),
                      ('5MDZ|1|1', 'P-site ', 'NO', '-', '-'),
                      ('6BU8|1|01', 'Appears to be A/A, P/P, E/E ', 'NO', '-', '-'),
                      ('6GWT|1|A', 'P/RF1 ', 'Termination', 'RF1 (v), RF3 (w)', 'Apidaecin (z)'),
                      ('6GXM|1|A', 'P/RF1 ', 'Termination', 'RF1 (v), RF3 (w)', 'Apidaecin (z)'),
                      ('6GXN|1|A', 'P/RF1', 'Termination', 'RF1 (v), RF3 (w)', 'Apidaecin (z)'),
                      ('6GXO|1|A', 'P/E', 'Termination', 'RF1 (v), RF3 (w)', 'Apidaecin (z)'),
                      ('6I7V|1|CA*', 'Empty', 'Ribosome heterogeneity ', '-', '-'),
                      ('6I7V|1|DA*', 'Empty', 'Ribosome heterogeneity ', '-', '-'),
                      ('4U1U|1|BA', 'Empty', 'Bound to antibiotic', '-', 'Quinupristin (B6)'),
                      ('4U1U|1|DA', 'Empty', 'Bound to antibiotic', '-', 'Quinupristin (D6)'),
                      ('4U1V|1|BA', 'Empty', 'Bound to antibiotic', '-', 'Linopristin (B6)'),
                      ('4U1V|1|DA', 'Empty', 'Bound to antibiotic', '-', 'Linopristin (D6)'),
                      ('4U20|1|BA', 'Empty', 'Bound to antibiotic', '-', 'Flopristin (BA)'),
                      ('4U20|1|DA', 'Empty', 'Bound to antibiotic', '-', 'Flopristin (DA)'),
                      ('4U24|1|BA', 'Empty', 'Bound to antibiotic', '-', 'Dalfopristin (BA)'),
                      ('4U24|1|DA', 'Empty', 'Bound to antibiotic', '-', 'Dalfopristin (DA)'),
                      ('4U25|1|BA', 'Empty', 'Bound to antibiotic', '-', 'Virginiamycin (BA)'),
                      ('4U25|1|DA', 'Empty', 'Bound to antibiotic', '-', 'Virginiamycin (DA)'),
                      ('4U26|1|BA', 'Empty', 'Bound to antibiotic', '-', 'Dalfopristin (B4) and Quinupristin (B6)'),
                      ('4U26|1|DA', 'Empty', 'Bound to antibiotic', '-', 'Dalfopristin (D4) and Quinupristin (D6)'),
                      ('4U27|1|BA', 'Empty', 'Bound to antibiotic', '-', 'Flopristin (BA) and Linopristin (B6)'),
                      ('4U27|1|DA', 'Empty', 'Bound to antibiotic', '-', 'Flopristin (DA) and Linopristin (D6)'),
                      ('4V4H|1|BB', 'Empty', 'Bound to antibiotic', '-', 'Kasugamycin (AA)'),
                      ('4V4H|1|DB', 'Empty', 'Bound to antibiotic', '-', 'Kasugamycin (CA)'),
                      ('4V52|1|BB', 'Empty', 'Bound to antibiotic', '-', 'Neomycin (AA & BB)'),
                      ('4V52|1|DB', 'Empty', 'Bound to antibiotic', '-', 'Neomycin (CA & DB)'),
                      ('4V53|1|BB', 'Empty', 'Bound to antibiotic', '-', 'Gentamycin (AA & BB)'),
                      ('4V53|1|DB', 'Empty', 'Bound to antibiotic', '-', 'Gentamycin (CA & DB)'),
                      ('4V54|1|BB', 'Empty', 'Ribosome recycling', 'RRF (B6)', '-'),
                      ('4V54|1|DB', 'Empty', 'Ribosome recycling', 'RRF (D6)', '-'),
                      ('4V55|1|BB', 'Empty', 'Ribosome recycling', 'RRF (B6)', 'Gentamycin (AA & BB)'),
                      ('4V55|1|DB', 'Empty', 'Ribosome recycling', 'RRF (D6)', 'Gentamycin (CA & DB)'),
                      ('4V56|1|BB', 'Empty', 'Bound to antibiotic', '-', 'Spectinomycin (AA)'),
                      ('4V56|1|DB', 'Empty', 'Bound to antibiotic', '-', 'Spectinomycin (CA)'),
                      ('4V57|1|BB', 'Empty', 'Bound to antibiotic', '-', 'Spectinomycin (AA) and Neomycin (AA & BB)'),
                      ('4V57|1|DB', 'Empty', 'Bound to antibiotic', '-', 'Spectinomycin (CA) and Neomycin (CA & DB)'),
                      ('4V64|1|BB', 'Empty', 'Bound to antibiotic', '-', 'Hygromycin B (AA)'),
                      ('4V64|1|DB', 'Empty', 'Bound to antibiotic', '-', 'Hygromycin B (CA)'),
                      ('4V7S|1|BA', 'Empty', 'Bound to antibiotic', '-', 'Telithromycin (BA)'),
                      ('4V7S|1|DA', 'Empty', 'NB', '-', '-'),
                      ('4V7T|1|BA', 'Empty', 'Bound to antibiotic', '-', 'Chloramphenicol (BA)'),
                      ('4V7T|1|DA', 'Empty', 'NB', '-', '-'),
                      ('4V7U|1|BA', 'Empty', 'Bound to antibiotic', '-', 'Erythromycin A (BA)'),
                      ('4V7U|1|DA', 'Empty', 'NB', '-', '-'),
                      ('4V7V|1|BA', 'Empty', 'Bound to antibiotic', '-', 'Clindamycin'),
                      ('4V7V|1|DA', 'Empty', 'NB', '-', '-'), ('4V9C|1|BA', 'P/P', 'NO', '-', 'neomycin (AA & BA)'),
                      ('4V9C|1|DA', 'P/E', 'NO', 'RRF(CY)', 'neomycin (CA & DA)'), ('4WF1|1|BA', '-', 'NB', '-', '-'),
                      ('4WF1|1|DA', '-', 'Bound to antibiotic', '-', 'Negamycin (CA)'),
                      ('4WOI|1|BA', 'P/E', 'Bound to antibiotic', '-', 'Paromomycin (AA & BA)'),
                      ('4WOI|1|CA', 'P/P', 'Bound to antibiotic', '-', 'Paromomycin (CA & DA)'),
                      ('4WWW|1|RA', 'Empty', 'Bound to antibiotic', '-', 'CEM-101 (RA)'),
                      ('4WWW|1|YA', 'Empty', 'NO', '-', '-'),
                      ('5KCR|1|1A', 'P/P ', 'Bound to antibiotic', '-', 'Avilamycin C (1A)'),
                      ('5KCS|1|1A', 'P/P ', 'Bound to antibiotic', 'TetM (1w)', 'Evernimycin (1A)'),
                      ('3J9Y|1|A', 'P/P', 'NO', 'TetM (w)', '-'),
                      ('3JCD|1|A', 'P/4 (8),  E/E (9)', 'Elongation - Backtranslocation (Post EF4) ', 'EF4 (x)', '-'), (
                          '3JCE|1|A', 'A/4 (6), P/4 (8), E/E (9)', 'Elongation - Backtranslocation (Pre EF4)',
                          'EF4 (x)',
                          '-'),
                      ('5AFI|1|A', 'A/T (y), P/P (v), E/E (w)', 'Elongation - Decoding', 'EF-Tu (z)', 'Kirromycin (z)'),
                      ('5UYK|1|01', 'T tRNA (Y), P/P (W), E/E (X)', 'Elongation - Decoding', 'EF-Tu (Z)', '-'),
                      ('5UYL|1|01', 'A*/T (Y), P/P (W), E/E (X)', 'Elongation - Decoding', 'EF-Tu (Z)', '-'),
                      ('5UYM|1|01', 'A/T (Y), P/P (W), E/E (X)', 'Elongation - Decoding', 'EF-Tu (Z)', '-'),
                      ('5UYN|1|01', 'T tRNA (Y), P/P (W), E/E (X)', 'Elongation - Decoding', 'EF-Tu (Z)', '-'),
                      ('5UYP|1|01', 'A*/T (Y), P/P (W), E/E (X)', 'Elongation - Decoding', 'EF-Tu (Z)', '-'),
                      ('5UYQ|1|01', 'A/T (Y), P/P (W), E/E (X)', 'Elongation - Decoding', 'EF-Tu (Z)', '-'),
                      ('5WDT|1|A', 'A/T (y), P/P (v), E/E (w)', 'Elongation - Decoding', 'EF-Tu H84A (z)', '-'),
                      ('5WE4|1|A', 'A/T (y), P/P (v), E/E (w)', 'Elongation - Decoding', 'EF-Tu (z)', '-'),
                      ('5WE6|1|A', 'A/T (y), P/P (v), E/E (w)', 'Elongation - Decoding', 'EF-Tu H84A (z)', '-'),
                      ('5WF0|1|A', 'A/T (y), P/P (v), E/E (w)', 'Elongation - Decoding', 'EF-Tu H84A (z)', '-'),
                      ('5WFK|1|A', 'A/T (y), P/P (v), E/E (w)', 'Elongation - Decoding', 'EF-Tu H84A (z)', '-'),
                      ('5WFS|1|A', 'A/T (y), P/P (v), E/E (w)', 'Elongation - Decoding', 'EF-Tu H84A (z)', '-'),
                      ('3J9Z|1|LA', 'P/P (S6), E/E (S7)', 'Elongation - Translocation', 'EF-G (S1)', '-'),
                      ('3JA1|1|LA', 'P/E (S2)', 'Elongation - Translocation', 'EF-G (S3)', '-'),
                      ('3JCJ|1|A', 'P/I* ', 'Initiation', 'IF-2 (f)', '-'),
                      ('6DNC|1|B', 'E/E (D), P/P (LA)', 'Premature Termination', 'RF1 (MA)', '-'),
                      ('5NP6|1|Y', 'P/P (B)', 'Recoding - Bypassing', '-', '-'), (
                          '6H4N|1|A', 'E/E (w)', 'Ribosome hibernation',
                          'Ribosome modulation factor (v), ribosome hibernation promoting factor (x)', '-'),
                      ('5H5U|1|A', 'P/P (5)', 'Ribosome Rescue', 'ArfA (3), RF2 (4)', '-'),
                      ('5MDV|1|1', 'P/P (5)', 'Ribosome Rescue', 'ArfA (6), RF2 (7)', '-'),
                      ('5MDW|1|1', 'P/P (5)', 'Ribosome Rescue', 'ArfA (6), RF2 (7)', '-'),
                      ('5MDY|1|1', 'P/P (5)', 'Ribosome Rescue', 'ArfA (6), RF2 (7)', '-'),
                      ('5MGP|1|A', 'P/P (x)', 'Ribosome Rescue', 'ArfA (w), RF2 (z)', '-'),
                      ('5U4I|1|A', 'P/P (x), E*/E (y)', 'Ribosome Rescue', 'ArfA (w), RF2 (v)', '-'),
                      ('5U4J|1|A', '?', 'Ribosome Rescue', 'ArfA (w), RF2 (v)', '-'),
                      ('5U9F|1|01', 'P/P(W), E*/E (X)', 'Ribosome Rescue', 'ArfA (Y), RF2 (Z)', '-'),
                      ('5U9G|1|01', 'P/P(W), E*/E (X)', 'Ribosome Rescue', 'ArfA (Y), RF2 (Z)', '-'),
                      ('6ENF|1|A', 'P/P (x)', 'Ribosome Rescue', '-', '-'),
                      ('6ENJ|1|A', 'A/A (9), P/P (x)', 'Ribosome Rescue', 'EF-P (w)', '-'),
                      ('6ENU|1|A', 'P/P (x)', 'Ribosome Rescue', 'EF-P (w)', '-'),
                      ('6C4I|1|A', 'P/P (x), E*/E (y)', 'Ribosome rescue', 'ArfA (w), RF-2 (v)', '-'),
                      ('3JBU|1|b', 'P/P (v)', 'Ribosome Stalling', 'SecM (z)', '-'), (
                          '3JBV|1|b', 'A/P* (V), P/E intermediate (W)', 'Ribosome Stalling', 'SecM (z)',
                          'Chloroamphenicol (b)'),
                      ('5JTE|1|BA', 'A/A (AW), P/P (AX), E/E* (AY)', 'Ribosome Stalling', '-',
                       'Erythromycin A (BA), ErmBL (B5)'), (
                          '5JU8|1|BA', 'P/P (AX), E*/E* (AY)', 'Ribosome Stalling', '-',
                          'Erythromycin A (BA), ErmBL (B5)'),
                      ('5NWY|1|N', 'P/P-VemP (M)', 'Ribosome Stalling', '-', 'VemP (s)'),
                      ('5O2R|1|A', 'P/P (x)', 'Ribosome Stalling', 'RF-1 (v)', 'Apidaecin (z)'),
                      ('5LZA|1|A', 'P/P (v)', 'SelB activation', '-', '-'),
                      ('5LZD|1|A', 'A/SelB (y), P/P (v)', 'SelB activation', 'SelB (z)', '-'),
                      ('5IQR|1|1', 'A/R (6), P/P (5), E/E (4)', 'Stringent Control', 'RelA (8)', 'Paromomycine (2)'),
                      ('5KPS|1|28', 'P/P (31), E/E (32)', 'Stringent Control', 'RelA (A)', '-'),
                      ('5KPW|1|27', 'A/R (30), P/P (31), E/E (32)', 'Stringent Control', 'RelA (33)', '-'),
                      ('5KPX|1|27', 'A/R (30), P/P (31), E/E (32)', 'Stringent Control', 'RelA (33)', '-'),
                      ('5L3P|1|A', 'A/R( y), P/P (x)', 'Stringent Control', 'RelA (z)', '-'),
                      ('4V85|1|BA', '-', 'Termination', 'RF3 (AW)', 'Viomycin (AY)'),
                      ('4V89|1|BA', '-', 'Termination', 'RF3 (AW)', '-'),
                      ('4V6C|1|BA', '-', 'Translation - elongation', '-', '-'),
                      ('4V6C|1|DA', '-', 'Translation - elongation', '-', '-'),
                      ('4V6D|1|BA', 'P-site ASL fragment', 'Translation - elongation', '-', '-'),
                      ('4V6D|1|DA', 'P-site ASL fragment', 'Translation - elongation', '-', '-'), (
                          '4V6E|1|BA', 'A-site ASL fragment (AX), P-site ASL fragment (AV)', 'Translation - elongation',
                          '-', '-'), (
                          '4V6E|1|DA', 'A-site ASL fragment (CX), P-site ASL fragment (CV)', 'Translation - elongation',
                          '-', '-'), ('3J7Z|1|A', 'Empty', 'Ribosome Stalling', 'ErmCL (a)', 'Erythromycin A (A)'),
                      ('6GC0|1|A', 'Empty', 'Ribosome Assembly', '-', '-'),
                      ('4UY8|1|A', 'Empty', 'Ribosome Stalling', 'TnaC (7)', '-'),
                      ('6GC8|1|A', 'Empty', 'Ribosome Assembly', '-', '-'), ('4V80|1|BA', '?', '?', '?', '?'),
                      ('4V80|1|DA', '?', '?', '?', '?'), ('6GBZ|1|A', 'Empty', 'Ribosome Assembly', '-', '-'),
                      ('6C4H|1|A', 'P-site ', 'Termination', 'RF2 (v)', '-'),
                      ('5GAE|1|A', '?', 'Co-translational protein targeting', '-', '-'), (
                          '5GAD|1|A', 'ESRP (1)', 'Co-translational protein targeting',
                          'SRP protein (i), SRP receptor (l)',
                          '-'), ('5GAH|1|A', 'ESRP (1)', 'Co-translational protein targeting', 'SRP protein (i)', '-'),
                      ('5GAG|1|A', 'ESRP (1)', 'Co-translational protein targeting', 'SRP protein (i)', '-'),
                      ('5LZE|1|A', '?', 'Pre-translocation SelB', '?', '?')]

    data = request.args['units']

    query_list = input_type(data)

    reject_list = ['5AFI|1|A', '5LZE|1|A', '4WRO|1|3L', '4WSD|1|1K', '4WSD|1|1L', '4WSD|1|3L', '4WT1|1|1K', '4WT1|1|1L', '4WT1|1|3K', '4WT1|1|3L', '4WZO|1|3K', '4WZO|1|3L', '4Y4P|1|1w']

    #######################################################################################################
    def getKey(item):
        return item[0]

    #######################################################################################################
    def get_chain_idx(query):

        range_selection = []
        for elem in query:
            range_selection.append(elem)

        chain_idx = []
        for sublist in query:
            units_query = UnitInfo.query.filter(UnitInfo.unit_id.in_(sublist))

            for rows in units_query:
                chain_idx.append(rows.chain_index)

        return chain_idx

    query_type = check_query(query_list)

    if query_type != 'loop_id':
        query_ife = '|'.join(query_list[0][0].split('|')[:3])
        query_pdb = query_list[0][0].split('|')[0]
        query_chain = query_list[0][0].split('|')[2]
    else:
        pass

    #######################################################################################################
    def get_sorted_units(units):
        unsorted_units = units.split(',')
        sorted_units = sorted(unsorted_units, key=lambda x: int(x.split('|')[4]))
        return sorted_units
    ########################################################################################################
    def check_insertion(corr_units):
        for unit in corr_units:
            ife_num = unit.split('|')[-1]
            try:
                num = int(ife_num)
            except ValueError:
                corr_units.remove(unit)

        return corr_units
    #########################################################################################################
    def custom_order(dct, spec):
        res = OrderedDict()
        for key in spec:
            if key in dct:
                res[key] = dct.pop(key)
        res.update(dct.items())

        return res
    #########################################################################################################
    def group_corr(corr_list):
        corr_list.sort()
        keyf = lambda x: '|'.join(x.split('|')[:3])
        corr_grouped = [list(items) for gr, items in groupby(corr_list, key=keyf)]

        return corr_grouped
    #########################################################################################################
    def order_num(corr_list):

        rej_sub = []
        for unit in corr_list:
            try:
                unit.sort(key=lambda x: int(x.split('|')[-1]))
            except ValueError:
                rej_sub.append(unit)
                corr_list.remove(unit)

        return corr_list
    ##########################################################################################################
    
    units_complete_list = []

    if query_type == 'single_range':

        chain_idx = get_chain_idx(query_list)
        chain_idx.sort()

        units_query = UnitInfo.query.filter_by(pdb_id=query_pdb, chain=query_chain). \
            filter(UnitInfo.chain_index.between(chain_idx[0], chain_idx[1])) \
            .order_by(UnitInfo.chain_index).all()

        for row in units_query:
            units_complete_list.append(row.unit_id)

    elif query_type == 'multiple_ranges':
        chain_idx = get_chain_idx(query_list)
        chain_idx.sort()
        partition_size = len(query_list)
        # Partition the list into a list of lists containing the start and end units of each range
        chain_idx = [chain_idx[i:i + 2] for i in range(0, len(chain_idx), 2)]

        for i in chain_idx:
            units_query = UnitInfo.query.filter_by(pdb_id=query_pdb, chain=query_chain). \
                filter(UnitInfo.chain_index.between(i[0], i[1])) \
                .order_by(UnitInfo.chain_index).all()
            for row in units_query:
                units_complete_list.append(row.unit_id)

        units_complete_list = list(OrderedDict.fromkeys(units_complete_list))

    elif query_type == 'units_str':

        for unit in query_list:
            units_complete_list.append(unit[0])

    # todo work to do
    elif query_type == 'loop_id':

        loop_id = query_list[0][0]
        units_query = LoopInfo.query.filter_by(loop_id=loop_id)

        for row in units_query:
            unsorted_units = row.unit_ids
            loop_position = row.loop_name

        units_complete_list = get_sorted_units(unsorted_units)
        query_ife = '|'.join(units_complete_list[0].split('|')[:3])
        query_pdb = units_complete_list[0].split('|')[0]

    ##########################################################################################################

    # This section of the code deals with getting the members of Equivalence Class from the query chain
    ife_list = NrChains.query.join(NrClasses, NrReleases) \
        .filter(NrChains.ife_id == query_ife).filter(NrClasses.resolution == '4.0') \
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

    # remove ifes that are joined (+)
    rejected_ife = []
    for i, v in enumerate(ife_members):
        if any(c in '+' for c in v):
            rejected_ife.append(ife_members[i])
            del ife_members[i]

    for elem in reject_list:
        for i, v in enumerate(ife_members):
            if elem == v:
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
    #######################################################################################################
    # query nts as a string
    query_nts = ', '.join(units_complete_list)

    query_complete_len = len(units_complete_list)
    #####################################################################################################

    # This section deals with getting the units of unmodified nucleotides
    '''
    standard_nts = ('A', 'C', 'G', 'U')

    units_std_list = []

    for unit in units_complete_list:
        k = unit.split('|')[-2]
        if k in standard_nts:
            units_std_list.append(unit)

    query_std_len = len(units_std_list)
    '''

    #####################################################################################################

    # This section of the code deals with getting the complete corresponding unit_ids
    ordering = case(
        {unit: index for index, unit in enumerate(units_complete_list)},
        value=UnitCorrespondence.unit_id_1
    )

    correspondence_complete = UnitCorrespondence.query.filter(UnitCorrespondence.unit_id_1.in_(units_complete_list)) \
        .order_by(ordering) \
        .filter(tuple_(UnitCorrespondence.pdb_id_2, UnitCorrespondence.chain_name_2) \
                .in_(members_info))

    # result_complete = [[unit.unit_id_2 for unit in units] for unit_id_1, units in
    # itertools.groupby(correspondence_complete, lambda x: x.unit_id_1)]

    # corr_complete = zip(*result_complete)
    # Append the units of the query motif
    # corr_complete.append(units_complete_list)

    corr_units = []

    for row in correspondence_complete:
        corr_units.append(row.unit_id_2)

    corr_filtered = check_insertion(corr_units)
    corr_grouped = group_corr(corr_filtered)
    corr_grouped = [x for x in corr_grouped if len(x) == query_complete_len]
    corr_complete = order_num(corr_grouped)
    corr_complete.append(units_complete_list)
    corr_std = deepcopy(corr_complete)

    accepted_seq = ['A', 'C', 'G', 'U']
    mod_idx = []
    for elem1 in corr_complete:
      for elem2 in elem1:
        seq = elem2.split('|')[3]
        if seq not in accepted_seq:
          mod_idx.append(elem1.index(elem2))

    mod_unique = set(mod_idx)
    mod_idx = list(mod_unique)

    for elem1 in corr_std:
      for ele in sorted(mod_idx, reverse = True): 
        del elem1[ele]

    query_std_len = len(corr_std[0])

    return json.dumps(corr_complete)

    '''

    # first_elem = []

    # for sublist in corr_ordered:
    # first_elem.append(sublist[0].split('|')[-1])

    # Create lists for residue type and number
    unit_list = []
    res_num = []
    res_type = []
    for units in corr_complete:
        unit_list.append(units[0])
        for unit in units:
            res_num.append(unit.split('|')[-1])
            res_type.append(unit.split('|')[-2])
        # ife = '|'.join(units[0].split('|')[:3])
        # unit_list.append(ife)

    res_num_list = [res_num[i:i + query_complete_len] for i in range(0, len(res_num), query_complete_len)]
    res_type_list = [res_type[i:i + query_complete_len] for i in range(0, len(res_type), query_complete_len)]
    # res_list = [res_num[i:i + query_complete_len] for i in xrange(0, len(res_num), query_complete_len)]

    residue_num_unordered = []
    for a in range(0, len(res_num_list)):
        residue_num_unordered.append(["{}{}".format(x, y) for x, y in zip(res_type_list[a], res_num_list[a])])

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

    # Create a dictionary of ifes with residue list
    ife_res = dict(zip(ife_list, residue_num_unordered))

    pdb_updated = []
    chain_updated = []

    for ife in ife_list:
        pdb_updated.append(ife.split('|')[0])
        chain_updated.append(ife.split('|')[-1])

    members_info_updated = zip(pdb_updated, chain_updated)
    ##################################################################################

    # Comment out
    # Get the list of corresponding unit-ids without modified nucleotides
    
    ordering = case(
        {unit: index for index, unit in enumerate(units_std_list)},
        value=UnitCorrespondence.unit_id_1
    )

    correspondence_std = UnitCorrespondence.query.filter(UnitCorrespondence.unit_id_1.in_(units_std_list)) \
        .order_by(ordering) \
        .filter(tuple_(UnitCorrespondence.pdb_id_2, UnitCorrespondence.chain_name_2) \
                .in_(members_info_updated))

    corr_std_units = []
    for row in correspondence_std:
        corr_std_units.append(row.unit_id_2)

    corr_std_grouped = group_corr(corr_std_units)
    corr_std_grouped = [x for x in corr_std_grouped if len(x) == query_std_len]
    corr_std = order_num(corr_std_grouped)

    # Append the standard units of the query motif
    corr_std.append(units_std_list)
    

    ##################################################################################

    # Logic to get and display pairwise interactions from the database
    bps_comb = []
    for a in range(0, len(corr_complete)):
        bps_comb.append([(map(str, comb)) for comb in combinations(corr_complete[a], 2)])

    unit1 = []
    unit2 = []
    bpair = []
    bstack = []
    bphosphate = []
    bribose = []
    for a in range(0, len(corr_complete)):
        bps_list = UnitPairInteractions.query.filter(
            tuple_(UnitPairInteractions.unit_id_1, UnitPairInteractions.unit_id_2) \
                .in_(bps_comb[a]))

        for row in bps_list:
            unit1.append(row.unit_id_1)
            unit2.append(row.unit_id_2)
            bpair.append(row.f_lwbp)
            bstack.append(row.f_stacks)
            bphosphate.append(row.f_bphs)
            bribose.append(row.f_brbs)

        pairwise_info = zip(unit1, unit2, bpair, bstack, bphosphate, bribose)

    filtered_pw_info = []
    for elem in pairwise_info:
        a = list(filter(lambda a: a != None, elem))
        filtered_pw_info.append(a)

    # return json.dumps(filtered_pw_info)

    units_order = {}
    for idx, unit in enumerate(units_complete_list):
        unit = unit.split('|')[-1]
        units_order[unit] = idx + 1

    n1 = []
    n2 = []

    for n in filtered_pw_info:
        n_1 = n[0].split('|')[-1]
        n_2 = n[1].split('|')[-1]
        try:
            n1.append(int(n_1))
            n2.append(int(n_2))
        except:
            pass

    possible_pw = zip(n1, n2)
    unique_pw = list(set(possible_pw))
    pw_sorted = sorted(unique_pw, key=getKey)

    pw_info = {k: OrderedDict({t: '-' for t in pw_sorted}) for k in ife_list}

    for sub_lst in filtered_pw_info:
        k0, k1 = '|'.join(sub_lst[0].split('|')[:3]), '|'.join(sub_lst[1].split('|')[:3])
        if k0 == k1 and k0 in pw_info:
            try:
                sub_key = (int(sub_lst[0][sub_lst[0].rfind('|') + 1:]), int(sub_lst[1][sub_lst[1].rfind('|') + 1:]))
                pw_info[k0][sub_key] = sub_lst[2] if len(sub_lst) == 3 else ';'.join(sub_lst[2:])
            except:
                pass

    ######################################################################################

    # Get center and rotation data for calculating discrepancy

    # Create list to store the centers np array
    units_center = []
    units_num_center = []

    # This section of the code deals with the database query to get the centers data
    for units in corr_std:

        ordering = case(
            {id: index for index, id in enumerate(units)},
            value=UnitCenters.unit_id
        )

        centers_query = UnitCenters.query.filter(UnitCenters.unit_id.in_(units),
                                                 UnitCenters.name == 'base').order_by(ordering)
        for row in centers_query:
            units_center.append(np.array([row.x, row.y, row.z]))
            units_num_center.append(row.unit_id)

    units_center_list = [units_center[i:i + query_std_len] for i in xrange(0, len(units_center), query_std_len)]

    # Create list to store the rotation np array
    units_rotation = []
    units_num_rotation = []

    # This section of the code deals with the database query to get the rotation data
    for units in corr_std:

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

    units_rotation_list = [units_rotation[i:i + query_std_len] for i in xrange(0, len(units_rotation), query_std_len)]

    rotation_size = len(units_rotation_list)

    ####################################################################################

    # This section of the code deals with calculating the discrepancy for the corresponding instances
    distances = coll.defaultdict(lambda: coll.defaultdict(int))

    for a in range(0, len(ife_list)):
        for b in range(a + 1, len(ife_list)):
            disc = matrix_discrepancy(units_center_list[a], units_rotation_list[a], units_center_list[b],
                                      units_rotation_list[b])
            distances[ife_list[a]][ife_list[b]] = disc

    # Empty list to append pairs of IFE with NaN discrepancy
    ife_nan = []

    for k, v in distances.items():
        for a, b in v.items():
            if math.isnan(b):
                ife_nan.append((k, a))
                v[a] = -0.1

    dist = np.zeros((len(ife_list), len(ife_list)))
    for index1, member1 in enumerate(ife_list):
        curr = distances.get(member1, {})
        for index2, member2 in enumerate(ife_list):
            dist[index1, index2] = curr.get(member2, 0)

    dist = (dist + np.swapaxes(dist, 0, 1))

    # ordering, _, _ = orderWithPathLengthFromDistanceMatrix(dist, 10, scanForNan=True)
    disc_order = optimalLeafOrder(dist)

    new_ordering = []
    idx_ordering = []

    for idx, order in enumerate(disc_order):
        new_ordering.append(ife_list[order])
        idx_ordering.append(idx)

    ifes_ordered = zip(idx_ordering, new_ordering)

    coord_ordered = []
    # append the coordinates based on new ordering
    for index in ifes_ordered:
        for key, val in ife_coord.iteritems():
            if index[1] == key:
                coord_ordered.append(val)

    res_list_ordered = []
    for index in ifes_ordered:
        for key, val in ife_res.iteritems():
            if index[1] == key:
                res_list_ordered.append(val)

    #########################################################################################

    # Logic to order and build the heatmap data

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

    disc_formatted = []
    for disc in disc_ordered:
        disc = '%.4f' % disc
        disc_formatted.append(disc)

    a = np.array(disc_formatted)
    a = a.astype(np.float)
    p1 = np.percentile(a, 90)
    p2 = np.percentile(a, 95)
    p3 = np.percentile(a, 99)
    maxDisc = np.amax(a)

    
    # Need to comment out
    rows = zip(ife1, ife2, disc_formatted)

    with open('/Applications/mamp/htdocs/corr-server/Disc/SSU-J/' + loop_id + '_disc.csv', "w") as f:
        writer = csv.writer(f)
        writer.writerow(["ID1", "ID2", "Disc"])
        for row in rows:
            writer.writerow(row)
    

    # return 'Max value for discrepancy is: {} and 95th percentile is {}'.format(maxDisc, p)

    heatmap_data = [
        {"ife1": ife1, "ife1_index": ife1_index, "ife2": ife2, "ife2_index": ife2_index, "discrepancy": discrepancy}
        for ife1, ife1_index, ife2, ife2_index, discrepancy in zip(ife1, index1, ife2, index2, disc_formatted)
    ]

    trna_occupancy = []
    functional_state = []
    factors_bound = []
    antibiotic_bound = []
    for order in ifes_ordered:
        for state in annotation:
            if order[1] == state[0]:
                trna_occupancy.append(state[1])
                functional_state.append(state[2])
                factors_bound.append(state[3])
                antibiotic_bound.append(state[4])

    new_order = []
    for elem in ifes_ordered:
        new_order.append(elem[1])

    pw_info_ordered = custom_order(pw_info, new_order)

    ###########################################################################################
    return render_template("correspondence_display.html", query_pdb=query_pdb, query_nts=query_nts,
                           coord=coord_ordered, ifes=ifes_ordered, maxDisc=maxDisc, p1=p1, p2=p2, p3=p3,
                           # loop_position=loop_position,
                           ec=equivalence_class, release=nr_release, residue_info=res_list_ordered, data=heatmap_data,
                           trna_occupancy=trna_occupancy, functional_state=functional_state,
                           factors_bound=factors_bound,
                           antibiotic_bound=antibiotic_bound, pw_info=pw_info_ordered, pw_list=pw_sorted)


@app.route('/bridges')
def bridges():

    def getKey(item):
        return item[0]

    ife_collection = ["3J9Y|1|a", "3J9Z|1|SA", "3JA1|1|SA", "3JBU|1|A", "3JBV|1|A", "3JCD|1|a", "3JCE|1|a", "3JCJ|1|g",
                      "4U1U|1|AA", "4U1U|1|CA", "4U1V|1|AA", "4U1V|1|CA", "4U20|1|AA",
                      "4U20|1|CA", "4U24|1|AA", "4U24|1|CA", "4U25|1|AA", "4U25|1|CA", "4U26|1|AA", "4U26|1|CA",
                      "4U27|1|AA", "4U27|1|CA", "4V4H|1|AA", "4V4H|1|CA", "4V4Q|1|AA", "4V4Q|1|CA", "4V50|1|AA",
                      "4V50|1|CA", "4V52|1|AA", "4V52|1|CA", "4V53|1|AA", "4V53|1|CA", "4V54|1|AA", "4V54|1|CA",
                      "4V55|1|AA", "4V55|1|CA", "4V56|1|AA", "4V56|1|CA", "4V57|1|AA", "4V57|1|CA", "4V5B|1|BA",
                      "4V5B|1|DA", "4V64|1|AA", "4V64|1|CA", "4V6C|1|AA", "4V6C|1|CA", "4V6D|1|AA", "4V6D|1|CA",
                      "4V6E|1|AA", "4V6E|1|CA", "4V7S|1|AA", "4V7S|1|CA", "4V7T|1|AA", "4V7T|1|CA", "4V7U|1|AA",
                      "4V7U|1|CA", "4V7V|1|AA", "4V7V|1|CA", "4V85|1|AA", "4V89|1|AA",
                      "4V9C|1|AA", "4V9C|1|CA", "4V9D|1|AA", "4V9D|1|BA", "4V9O|1|BA", "4V9O|1|DA", "4V9O|1|FA",
                      "4V9O|1|HA", "4V9P|1|BA", "4V9P|1|DA", "4V9P|1|FA", "4V9P|1|HA", "4WF1|1|AA", "4WF1|1|CA",
                      "4WOI|1|AA", "4WOI|1|DA", "4WWW|1|QA", "4WWW|1|XA", "4YBB|1|AA", "4YBB|1|BA", "5H5U|1|h",
                      "5IQR|1|2", "5IT8|1|AA", "5IT8|1|BA", "5J5B|1|AA", "5J5B|1|BA", "5J7L|1|BA", "5J88|1|AA",
                      "5J88|1|BA", "5J8A|1|AA", "5J8A|1|BA", "5J91|1|AA", "5J91|1|BA", "5JC9|1|AA", "5JC9|1|BA",
                      "5JTE|1|AA", "5JU8|1|AA", "5KCR|1|1a", "5KCS|1|1a", "5KPS|1|27", "5KPW|1|26", "5KPX|1|26",
                      "5L3P|1|a", "5MDV|1|2", "5MDW|1|2", "5MDY|1|2", "5MDZ|1|2", "5MGP|1|a", "5NP6|1|D", "5NWY|1|0",
                      "5O2R|1|a", "5U4I|1|a", "5U4J|1|a", "5U9F|1|A", "5U9G|1|A", "5UYK|1|A", "5UYL|1|A", "5UYM|1|A",
                      "5UYN|1|A", "5UYP|1|A", "5UYQ|1|A", "5WDT|1|a", "5WE4|1|a", "5WE6|1|a", "5WF0|1|a", "5WFS|1|a",
                      "6BU8|1|A", "6C4I|1|a", "6DNC|1|A", "6ENF|1|a", "6ENJ|1|a", "6ENU|1|a", "5J7L|1|AA"]

    unit1 = []
    unit2 = []
    bpair = []
    bstack = []
    bphosphate = []
    bribose = []
    fcrossing = []

    for elem in range(0, len(ribosome_subunits)):
        bridge_list = UnitPairInteractions.query.filter(
            UnitPairInteractions.unit_id_1.like(ribosome_subunits[elem][0] + '%') &
            UnitPairInteractions.unit_id_2.like(ribosome_subunits[elem][1] + '%'))
        for row in bridge_list:
            unit1.append(row.unit_id_1)
            unit2.append(row.unit_id_2)
            bpair.append(row.f_lwbp)
            bstack.append(row.f_stacks)
            bphosphate.append(row.f_bphs)
            bribose.append(row.f_brbs)
            fcrossing.append(row.f_crossing)

        pairwise_info = zip(unit1, unit2, bpair, bstack, bphosphate, bribose)

    bridging_interactions = []
    for elem in pairwise_info:
        a = list(filter(lambda a: a != None, elem))
        bridging_interactions.append(a)

    unit1 = []
    unit2 = []
    for i in bridging_interactions:
        n_1 = i[0].split('|')[-1]
        n_2 = i[1].split('|')[-1]
        try:
            unit1.append(int(n_1))
            unit2.append(int(n_2))
        except:
            pass

    possible_pw = zip(unit1, unit2)
    unique_pw = list(set(possible_pw))
    pw_sorted = sorted(unique_pw, key=lambda element: (element[0], element[1]))

    pw_info = {k: OrderedDict({t: '-' for t in pw_sorted}) for k in ife_collection}

    for sub_lst in bridging_interactions:
        k0, k1 = '|'.join(sub_lst[0].split('|')[:3]), '|'.join(sub_lst[1].split('|')[:3])
        if k0 in pw_info:
            v1, v2 = sub_lst[0], sub_lst[1]
            # `sub_key` is aimed to be a key for inner dict of the predefined `pw_info` dict
            # thus it's composed as a tuple of trailing numbers of the first 2 items
            # in sub_list (ex. `(262, 263)`)
            sub_key = (int(v1[v1.rfind('|') + 1:]), int(v2[v2.rfind('|') + 1:]))
            pw_info[k0][sub_key] = sub_lst[2] if len(sub_lst) == 3 else ';'.join(sub_lst[2:])

    SSU_chain = pw_info.keys()
    LSU_chain = []
    for elem1 in SSU_chain:
        for elem2 in ribosome_subunits:
            if elem1 == elem2[0]:
                LSU_chain.append(elem2[1])

    return render_template("bridge_table.html", SSU_chain=SSU_chain, LSU_chain=LSU_chain, pw_info=pw_info, pw_list=pw_sorted)
    # return json.dumps(str(LSU_chain))

@app.route('/trna_interactions')
def trna_interactions():

    # ife_pairs = [('6BU8|1|A', '6BU8|1|Y'), ('3JCE|1|a', '3JCE|1|6'), ('5UYK|1|A', '5UYK|1|Y'), ('5UYL|1|A', '5UYL|1|Y'), ('5UYM|1|A', '5UYM|1|Y'), ('5UYN|1|A', '5UYN|1|Y'), ('5UYP|1|A', '5UYP|1|Y'), ('5UYQ|1|A', '5UYQ|1|Y'), ('5WDT|1|a', '5WDT|1|y'), ('5WE4|1|a', '5WE4|1|y'), ('5WE6|1|a', '5WE6|1|y'), ('5WF0|1|a', '5WF0|1|y'), ('5WFK|1|a', '5WFK|1|y'), ('5WFS|1|a', '5WFS|1|y'), ('6ENJ|1|a', '6ENJ|1|9'), ('3JBV|1|A', '3JBV|1|V'), ('5JTE|1|AA', '5JTE|1|AW'), ('5LZD|1|a', '5LZD|1|y'), ('5IQR|1|2', '5IQR|1|6'), ('5KPW|1|26', '5KPW|1|30'), ('5KPX|1|26', '5KPX|1|30'), ('5L3P|1|a', '5L3P|1|y'), ('4V6E|1|AA', '4V6E|1|AX'), ('4V6E|1|CA', '4V6E|1|CX')]
    # ife_pairs = [('6BU8|1|A', '6BU8|1|W'), ('3JCE|1|a', '3JCE|1|8'), ('5UYK|1|A', '5UYK|1|W'), ('5UYL|1|A', '5UYL|1|W'), ('5UYM|1|A', '5UYM|1|W'), ('5UYN|1|A', '5UYN|1|W'), ('5UYP|1|A', '5UYP|1|W'), ('5UYQ|1|A', '5UYQ|1|W'), ('5WDT|1|a', '5WDT|1|v'), ('5WE4|1|a', '5WE4|1|v'), ('5WE6|1|a', '5WE6|1|v'), ('5WF0|1|a', '5WF0|1|v'), ('5WFK|1|a', '5WFK|1|v'), ('5WFS|1|a', '5WFS|1|v'), ('6ENJ|1|a', '6ENJ|1|x'), ('3JBV|1|A', '3JBV|1|W'), ('5JTE|1|AA', '5JTE|1|AX'), ('5LZD|1|a', '5LZD|1|v'), ('5IQR|1|2', '5IQR|1|5'), ('5KPW|1|26', '5KPW|1|31'), ('5KPX|1|26', '5KPX|1|31'), ('5L3P|1|a', '5L3P|1|x'), ('4V6E|1|AA', '4V6E|1|AV'), ('4V6E|1|CA', '4V6E|1|CV')]
    # ife_pairs = [('6BU8|1|A', '6BU8|1|X'), ('3JCD|1|a', '3JCD|1|9'), ('3JCE|1|a', '3JCE|1|9'), ('5UYK|1|A', '5UYK|1|X'), ('5UYL|1|A', '5UYL|1|X'), ('5UYM|1|A', '5UYM|1|X'), ('5UYN|1|A', '5UYN|1|X'), ('5UYP|1|A', '5UYP|1|X'), ('5UYQ|1|A', '5UYQ|1|X'), ('5WDT|1|a', '5WDT|1|w'), ('5WE4|1|a', '5WE4|1|w'), ('5WE6|1|a', '5WE6|1|w'), ('5WF0|1|a', '5WF0|1|w'), ('5WFK|1|a', '5WFK|1|w'), ('5WFS|1|a', '5WFS|1|w'), ('3J9Z|1|SA', '3J9Z|1|S7'), ('6DNC|1|A', '6DNC|1|D'), ('6H4N|1|a', '6H4N|1|w'), ('5U4I|1|a', '5U4I|1|y'), ('5U9F|1|A', '5U9F|1|X'), ('5U9G|1|A', '5U9G|1|X'), ('6C4I|1|a', '6C4I|1|y'), ('5JTE|1|AA', '5JTE|1|AY'), ('5JU8|1|AA', '5JU8|1|AY'), ('5IQR|1|2', '5IQR|1|4'), ('5KPS|1|27', '5KPS|1|32'), ('5KPW|1|26', '5KPW|1|32'), ('5KPX|1|26', '5KPX|1|32')]
    ife_pairs = [('6BU8|1|01', '6BU8|1|Y'), ('3JCE|1|A', '3JCE|1|6'), ('5UYK|1|01', '5UYK|1|Y'), ('5UYL|1|01', '5UYL|1|Y'), ('5UYM|1|01', '5UYM|1|Y'), ('5UYN|1|01', '5UYN|1|Y'), ('5UYP|1|01', '5UYP|1|Y'), ('5UYQ|1|01', '5UYQ|1|Y'), ('5WDT|1|A', '5WDT|1|y'), ('5WE4|1|A', '5WE4|1|y'), ('5WE6|1|A', '5WE6|1|y'), ('5WF0|1|A', '5WF0|1|y'), ('5WFK|1|A', '5WFK|1|y'), ('5WFS|1|A', '5WFS|1|y'), ('6ENJ|1|A', '6ENJ|1|9'), ('3JBV|1|b', '3JBV|1|V'), ('5JTE|1|BA', '5JTE|1|AW'), ('5LZD|1|A', '5LZD|1|y'), ('5IQR|1|1', '5IQR|1|6'), ('5KPW|1|27', '5KPW|1|30'), ('5KPX|1|27', '5KPX|1|30'), ('5L3P|1|A', '5L3P|1|y'), ('4V6E|1|BA', '4V6E|1|AX'), ('4V6E|1|DA', '4V6E|1|CX')]

    ife_key = []
    for i in ife_pairs:
        ife_key.append(i[1])

    unit1 = []
    unit2 = []
    bpair = []
    bstack = []
    bphosphate = []
    bribose = []
    fcrossing = []

    for elem in range(0, len(ife_pairs)):
        bridge_list = UnitPairInteractions.query.filter(
            UnitPairInteractions.unit_id_1.like(ife_pairs[elem][0] + '%') &
            UnitPairInteractions.unit_id_2.like(ife_pairs[elem][1] + '%'))
        for row in bridge_list:
            unit1.append(row.unit_id_1)
            unit2.append(row.unit_id_2)
            bpair.append(row.f_lwbp)
            bstack.append(row.f_stacks)
            bphosphate.append(row.f_bphs)
            bribose.append(row.f_brbs)
            fcrossing.append(row.f_crossing)

        pairwise_info = zip(unit1, unit2, bpair, bstack, bphosphate, bribose)

    pw_interactions = []
    for elem in pairwise_info:
        a = list(filter(lambda a: a != None, elem))
        pw_interactions.append(a)

    unit1 = []
    unit2 = []
    for i in pw_interactions:
        n_1 = i[0].split('|')[-1]
        n_2 = i[1].split('|')[-1]
        try:
            unit1.append(int(n_1))
            unit2.append(int(n_2))
        except:
            pass

    possible_pw = zip(unit1, unit2)
    unique_pw = list(set(possible_pw))
    pw_sorted = sorted(unique_pw, key=lambda element: (element[0], element[1]))

    pw_info = {k: OrderedDict({t: '-' for t in pw_sorted}) for k in ife_key}


    for sub_lst in pw_interactions:
        k0, k1 = '|'.join(sub_lst[0].split('|')[:3]), '|'.join(sub_lst[1].split('|')[:3])
        if k1 in pw_info:
            v1, v2 = sub_lst[0], sub_lst[1]
            # `sub_key` is aimed to be a key for inner dict of the predefined `pw_info` dict
            # thus it's composed as a tuple of trailing numbers of the first 2 items
            # in sub_list (ex. `(262, 263)`)
            sub_key = (int(v1[v1.rfind('|') + 1:]), int(v2[v2.rfind('|') + 1:]))
            pw_info[k1][sub_key] = sub_lst[2] if len(sub_lst) == 3 else ';'.join(sub_lst[2:])

    ife_trna = []
    for i in ife_pairs:
        ife_trna.append(i[1])

    index_map = {v: i for i, v in enumerate(ife_trna)}
    pw_info = OrderedDict(sorted(pw_info.items(), key=lambda pair: index_map[pair[0]]))

    return render_template("trna_interaction.html", ife_pairs=ife_pairs, pw_info=pw_info, pw_list=pw_sorted)
   
    '''
if __name__ == '__main__':
    app.run(debug=True)
