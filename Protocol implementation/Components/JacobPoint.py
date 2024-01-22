class JacobPoint(object):
    """Egy Jacobian  pont az elliptikus görbén (x,y,z)"""

    def __init__(self, x, y, z):
        """
        Jacobian térbeli pont inicializálása
        :param x: x tengelybeli érték;
        :param y: y tengelybeli érték;
        :param z: z tengelybeli érték;
        """
        self.x = x
        self.y = y
        self.z = z

    def isInfinity(self):
        """
        :return: végtelen pont-e vagy sem;
        """
        return self.z == 0
