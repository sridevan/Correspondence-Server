from app import db


class UnitInfo(db.Model):
	__tablename__ = "unit_info"

	unit_id = db.Column(db.String, primary_key=True)
	pdb_id = db.Column(db.String)
	chain = db.Column(db.String)
	chain_index = db.Column(db.Integer)


class UnitCorrespondence(db.Model):
    __tablename__ = "correspondence_units"

    # correspondence_id = db.Column(db.String, primary_key=True)
    unit_id_1 = db.Column(db.String, primary_key=True)
    unit_id_2 = db.Column(db.String, primary_key=True)
    pdb_id_1 = db.Column(db.String, primary_key=True)
    pdb_id_2 = db.Column(db.String, primary_key=True)
    chain_name_2 = db.Column(db.String, primary_key=True)

    # def __init__(self, unit_id_1, unit_id_2, pdb_id_1, pdb_id_2):
    # self.unit_id_1 = unit_id_1
    # self.unit_id_2 = unit_id_2
    # self.pdb_id_1 = pdb_id_1
    # self.pdb_id_2 = pdb_id_2


class UnitCenters(db.Model):
    __tablename__ = "unit_centers"

    # correspondence_id = db.Column(db.String, primary_key=True)
    unit_center_id = db.Column(db.Integer, primary_key=True)
    unit_id = db.Column(db.String)
    pdb_id = db.Column(db.String)
    name = db.Column(db.String)
    x = db.Column(db.Float)
    y = db.Column(db.Float)
    z = db.Column(db.Float)


class UnitRotations(db.Model):
    __tablename__ = "unit_rotations"

    # correspondence_id = db.Column(db.String, primary_key=True)
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

class LoopInfo(db.Model):
    __tablename__ = "loop_info"

    loop_id = db.Column(db.String, primary_key=True)
    unit_ids = db.Column(db.Text)
    loop_name = db.Column(db.Text)

class NrReleases(db.Model):
    __tablename__ = "nr_releases"

    nr_release_id = db.Column(db.String, primary_key=True)
    date = db.Column(db.Date)
    classes = db.relationship("NrClasses", backref='nr_releases', lazy=True)

class NrClasses(db.Model):
    __tablename__ = "nr_classes"

    nr_class_id = db.Column(db.Integer, primary_key=True)
    name = db.Column(db.String)
    nr_release_id = db.Column(db.String, db.ForeignKey("nr_releases.nr_release_id"))
    resolution = db.Column(db.String)
    chains = db.relationship("NrChains", backref='nr_classes', lazy=True)

class NrChains(db.Model):
    __tablename__ = "nr_chains"

    nr_chain_id = db.Column(db.Integer, primary_key=True)
    ife_id = db.Column(db.String)
    nr_class_id = db.Column(db.Integer, db.ForeignKey("nr_classes.nr_class_id"))
    nr_release_id = db.Column(db.String, db.ForeignKey("nr_releases.nr_release_id"))

class UnitPairInteractions(db.Model):
    __tablename__ = "unit_pairs_interactions"

    unit_pairs_interactions_id = db.Column(db.Integer, primary_key=True)
    unit_id_1 = db.Column(db.String)
    unit_id_2 = db.Column(db.String)
    pdb_id = db.Column(db.String)
    f_lwbp = db.Column(db.String)
    f_stacks = db.Column(db.String)
    f_bphs = db.Column(db.String)
    f_brbs = db.Column(db.String)





