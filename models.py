from app import db

class UnitCorrespondence(db.Model):

	__tablename__ = "correspondence_units"

	correspondence_id = db.Column(db.String, primary_key=True)
	unit_id_1 = db.Column(db.String, primary_key=True)
	unit_id_2 = db.Column(db.String, primary_key=True)

	def __init__(self, unit_id_1, unit_id_2):
		self.unit_id_1 = unit_id_1
		self.unit_id_2 = unit_id_2

	def __repr__(self):
		return "{}:{}".format(self.unit_id_1, self.unit_id_2)
		#equivalent_residues = {self.unit_id_1: self.unit_id_2}
		#return equivalent_residues
		


    
