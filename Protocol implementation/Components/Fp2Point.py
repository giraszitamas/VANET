class Fp2Point(object):
    """Egy projektív térbeli pont az elliptikus görbén (x,y)"""

    def __init__(self, x, y):
        """
        Projektív térbeli pont inicializálása
        :param x: x tengelybeli érték;
        :param y: y tengelybeli érték;
        """
        self.x = x
        self.y = y
        
    def isInfinity(self):
        """
        :return: végtelen pont-e vagy sem;
        """
        return self.x is None and self.y is None
    
    def toString(self):
        """
        :return: a pont string reprezentációja;
        """
        return str(self.x) + " " + str(self.y) 
