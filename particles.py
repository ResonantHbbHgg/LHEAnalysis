from math import log, tan, acos, pi, copysign

class particle(object):
    def __init__(self, index, pid, status, mother1, mother2, color1, color2, px, py, pz, e, mass, f1, f2):
        self.index = index
        self.pid = pid
        self.status = status
        self.mother1 = mother1
        self.mother2 = mother2
        self.color1 = color1
        self.color2 = color2
        self.px = px
        self.py = py
        self.pz = pz
        self.e = e
        self.mass = mass
        self.f1 = f1
        self.f2 = f2
        # Computed on the fly
        self.p2 = self.px**2 + self.py**2 + self.pz**2
        self.p = self.p2 ** 0.5
        self.pt = (self.px**2 + self.py**2)**0.5
        self.theta = acos(self.pz / self.p)
        try:
            self.phi = acos(self.px / self.pt)
        except ZeroDivisionError:
            self.phi = 0.
        try:
            self.eta = -log(tan(self.theta / 2.))
        except ValueError:
            self.eta = copysign(9999.0, self.pz / self.p)

        if self.py < 0:
            self.phi = -self.phi

def deltaR(a, b):
    deta = a.eta - b.eta
    dphi = a.phi - b.phi
    if dphi > pi:
        dphi -= pi
    if dphi < -pi:
        dphi += pi
    return (deta**2. + dphi**2.)**.5


