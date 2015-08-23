# encoding: utf-8

class PyseidonError(Exception):
    """
    Custom error for PySeidon library
    """
    def __init__(self, arg):
        # Call the base class constructor with the parameters it needs
        super(PyseidonError, self).__init__(arg)
