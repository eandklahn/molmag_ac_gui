# https://www.programiz.com/python-programming/user-defined-exception

class FileFormatError(Exception):
    """Exception raised for errors in the reading of files

    Attributes:
        filename -- name of file that did not conform to the implied standard
        message -- explanation of the error
    """

    def __init__(self, filename):
        self.filename = filename
        self.message = '{} is not formatted to be read according to the implied standard'.format(self.filename)
        super().__init__(self.message)
        
