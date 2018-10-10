from app import db

class UnitCorrespondence(db.Model):

	__tablename__ = "correspondence_units"

	#correspondence_id = db.Column(db.String, primary_key=True)
	unit_id_1 = db.Column(db.String, primary_key=True)
	unit_id_2 = db.Column(db.String, primary_key=True)
	pdb_id_1 = db.Column(db.String, primary_key=True)
	pdb_id_2 = db.Column(db.String, primary_key=True)

	#def __init__(self, unit_id_1, unit_id_2, pdb_id_1, pdb_id_2):
		#self.unit_id_1 = unit_id_1
		#self.unit_id_2 = unit_id_2
		#self.pdb_id_1 = pdb_id_1
		#self.pdb_id_2 = pdb_id_2

class UnitCenters(db.Model):

	__tablename__ = "unit_centers"

	#correspondence_id = db.Column(db.String, primary_key=True)
	unit_center_id = db.Column(db.Integer, primary_key=True)
	unit_id = db.Column(db.String)
	pdb_id = db.Column(db.String)
	name = db.Column(db.String)
	x = db.Column(db.Float)
	y = db.Column(db.Float)
	z = db.Column(db.Float)

class UnitRotations(db.Model):

	__tablename__ = "unit_rotations"

	#correspondence_id = db.Column(db.String, primary_key=True)
	unit_id = db.Column(db.String, primary_key=True)
	pdb_id = db.Column(db.String)
	cell_0_0 = db.Column(db.Float)
	cell_0_1 = db.Column(db.Float)
	cell_0_2 = db.Column(db.Float)
	cell_1_0 = db.Column(db.Float)
	cell_1_1 = db.Column(db.Float)
	cell_1_2 = db.Column(db.Float)
	cell_2_0 = db.Column(db.Float)
	cell_2_1 = db.Column(db.Float)
	cell_2_2 = db.Column(db.Float)

    
