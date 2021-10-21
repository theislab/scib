class RootCellError(Exception):
    def __init__(self, message):
        self.message = message


class NeighborsError(Exception):
    def __init__(self, message):
        self.message = message
