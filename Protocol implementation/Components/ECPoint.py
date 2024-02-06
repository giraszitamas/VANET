# A point on an elliptic curve: (x,y)
class ECPoint:
    """Egy affin térbeli pont az elliptikus görbén (x,y)"""

    def __init__(self, curve, x, y):
        """
        Elliptikus görbén lévő pont inicializálása
        :param curve: az alkalmazott elliptikus görbe;
        :param x: x tengelybeli értéke;
        :param y: y tengelybeli értéke;
        """
        self.curve = curve
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
    
    def getCompact(self):
        """
        :return: a pont string reprezentációja;
        """
        return self.x
