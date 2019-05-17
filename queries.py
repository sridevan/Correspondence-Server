"""
This file will contain the logic for getting and processing
database queries
"""
from models import *

def get_class_id(query_ife):

    ife_list = NrChains.query.join(NrClasses, NrReleases)\
        .filter(NrChains.ife_id == query_ife).filter(NrClasses.resolution == '4.0')\
        .order_by(NrReleases.date.desc()).limit(1) 
       
    for row in ife_list:
        class_id = row.nr_class_id

    return class_id