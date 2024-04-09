import ROOT
import shipunit as u 

class Particle:
    
    def __init__(self):
        self.pdg = ROOT.TDatabasePDG.Instance()
        
    def PDGname(self, particle):
        """
        Change particle name for use with the PDG database
        """
        if '1' in particle:
            particle = particle.replace('1',"'")
        if particle == 'nu': return 'nu_e' # simple workaround to handle 3nu channel
        return particle

    def mass(self, particle):
        """
        Read particle mass from PDG database
        """
        particle = self.PDGname(particle)
        tPart = self.pdg.GetParticle(particle)
        print(particle)
        return tPart.Mass()

    def lifetime(self, particle):
        """
        Read particle lifetime from PDG database
        """
        particle = self.PDGname(particle)
        tPart = self.pdg.GetParticle(particle)
        return tPart.Lifetime()

class CKMmatrix():
    """
    CKM matrix, from http://pdg.lbl.gov/2017/reviews/rpp2016-rev-ckm-matrix.pdf
    """
    def __init__(self):
        self.Vud = 0.9743
        self.Vus = 0.2251
        self.Vub = 3.6e-03
        self.Vcd = 0.2249
        self.Vcs = 0.9735
        self.Vcb = 4.11e-02
        self.Vtd = 8.6e-03
        self.Vts = 4.0e-02
        self.Vtb = 0.999
        
        


class constants():
    """
    Store some constants useful for HNL physics
    """
    def __init__(self):
        self.decayConstant = {'pi+':0.130*u.GeV,
                            'pi0':0.130*u.GeV,
                            'eta':1.2*0.130*u.GeV,
                            'eta1':0.152*u.GeV,
                            'eta_c':0.335*u.GeV,
                            'rho+':0.209*u.GeV, 
                            'rho0':0.209*u.GeV,
                            'omega':0.195*u.GeV,
                            'phi':0.229*u.GeV,
                            'D_s+':0.249*u.GeV,
                            'D*_s+':0.315*u.GeV} # decay constants f_h of pseudoscalar and vector mesons
        self.GF = 1.166379e-05/(u.GeV*u.GeV) # Fermi's constant (GeV^-2)
        self.s2thetaw = 0.23126 # square sine of the Weinberg angle
        self.heV = 6.58211928*pow(10.,-16) # no units or it messes up!!
        self.hGeV = self.heV * pow(10.,-9) # no units or it messes up!!
