"""
This file will contain the logic for determining the type of input
Generally, input can be divided into three types: string of unit_ids,
ranges and loop_id
"""

import logging
import collections as coll

from werkzeug.exceptions import BadRequest

PARTS_LIMIT = 15

def input_type(data):

    parts = data.split(',')
    if len(parts) == 0 or len(parts) > PARTS_LIMIT:
        raise BadRequest("Must give 1 to %s parts in a collection" %
                         PARTS_LIMIT)

    processed = []
    for part in parts:
        if not part:
            raise BadRequest("Cannot give empty part of a list")
        units = part.split(':')
        if len(units) == 1:
            if not units[0]:
                raise BadRequest("Cannot give empty unit")
            processed.append(tuple([units[0], units[0]]))
        elif len(units) == 2:
            if not units[0] or not units[1]:
                raise BadRequest("Both units in range must exist")
            processed.append(tuple(units))
        else:
            raise BadRequest("Range should must have 1 or 2 endpoints")

    return processed
