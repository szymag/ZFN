from src.io.DataWriter import DataWriter

class StdOutDataWriter(DataWriter):
    def __init__(self):
        super().__init__()

    def write(self, data):
        print(data)