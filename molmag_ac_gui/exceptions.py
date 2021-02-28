# https://www.programiz.com/python-programming/user-defined-exception

class FileFormatError(Exception):
    """Exception raised for errors in the reading of files

    Attributes:
        salary -- input salary which caused the error
        message -- explanation of the error
    """

    def __init__(self, filename):
        self.filename = filename
        self.message = '{} is not formatted to be read like this'.format(self.filename)
        super().__init__(self.message)
        
