"""
This file will contain the logic for getting and processing
database queries
"""
from models import *

def get_chain_idx(query):

        range_selection = []
        for elem in query:
            range_selection.append(elem)

        chain_idx = []
        for elem in range_selection:
            units_query = UnitInfo.query.filter_by(pdb_id=query_pdb, chain=query_chain).filter(UnitInfo.unit_id.in_(elem))

            for rows in units_query:
                chain_idx.append(rows.chain_index)

        return chain_idx

def return_units(query_type, query_pdb, query_chain, query_list):

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

        # Partition the list into a list of lists containing the start and end units of each range
        chain_idx = [chain_idx[i:i + len(query_list)] for i in range(0, len(chain_idx), len(query_list))]

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

    elif query_type == 'loop_id':
        
        loop_id = query_list[0][0]
    
    return units_complete_list