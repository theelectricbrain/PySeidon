# encoding: utf-8

class PyseidonError(Exception):
    def __init__(self, arg):
        # Set some exception infomation
        self.msg = arg
