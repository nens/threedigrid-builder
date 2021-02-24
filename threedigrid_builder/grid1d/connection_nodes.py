__all__ = ["ConnectionNodes"]


class ConnectionNodes:
    def __init__(self, the_geom, id, code, storage_area):
        self.the_geom = the_geom
        self.id = id
        self.code = code
        self.storage_area = storage_area

    def __repr__(self):
        return "<ConnectionNodes object (len:{})>".format(len(self.the_geom))
