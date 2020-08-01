



class SubdivisionError(Exception):
    def __init__(self,msg):
        if msg:
            self.message = msg
        else:
            self.message = None

    def __str__(self):
        return 'SubdivisionError: {0}'.format(self.message)

class IntervalError(Exception):
    def __init__(self,msg):
        if msg:
            self.message = msg
        else:
            self.message = None

    def __str__(self):
        return 'IntervalError: {0}'.format(self.message)

class IterationError(Exception):
    def __init__(self,msg):
        if msg:
            self.message = msg
        else:
            self.message = None

    def __str__(self):
        return 'IterationError: {0}'.format(self.message)

