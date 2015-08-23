class ObjectFromDict(object):
    """
    Turns any given class object into a dictionnary
    """
    def __init__(self, d):
        self.__dict__ = d
