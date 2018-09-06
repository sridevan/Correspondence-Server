
from app import db

class UnitInfo(db.Model):

	__tablename__ = "unit_info"

	unit_id = db.Column(db.String, primary_key=True)
	pdb_id = db.Column(db.String, nullable=True)
	unit_type_id = db.Column(db.String, nullable=True)

	def __init__(self, pdb_id, unit_type_id):
		self.pdb_id = pdb_id
		self.unit_type_id = unit_type_id

	def __repr__(self):
		return "pdb id {} molecule type {}".format(self.pdb_id, self.unit_type_id)
    
