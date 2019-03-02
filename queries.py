"""
This file will contain the logic for getting and processing
database queries
"""

def ranges(data, pdb, chain):

    units_list = []

    if len(data) > 2:

        for range in data:
            units_list.append(range[0])

        '''
        for range in data:
            units_query = UnitInfo.query.filter_by(pdb_id=pdb, chain=chain). \
                filter(UnitInfo.chain_index.between(range[0], range[1])) \
                .order_by(UnitInfo.chain_index).all()

            for row in units_query:
                units_list.append(row.unit_id)
        '''

    return units_list