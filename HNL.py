import os
import sys
import pandas as pd

# Chemin vers l'installation de ROOT
root_install_path = 'root_install'

# Configuration des variables d'environnement pour ROOT
os.environ['ROOTSYS'] = root_install_path
os.environ['PATH'] = os.path.join(root_install_path, 'bin') + ':' + os.environ['PATH']
os.environ['LD_LIBRARY_PATH'] = os.path.join(root_install_path, 'lib') + ':' + os.environ.get('LD_LIBRARY_PATH', '')

# Ajout du chemin des bibliothèques ROOT à sys.path
sys.path.append(os.path.join(root_install_path, 'lib'))

import ROOT
from const_utils import Particle, CKMmatrix, constants
from particle import ParticleConfig
import math

import shipunit as u
import yaml
# Load some useful constants
c = constants()


# Configuration pour HNL
class HNLConfig(ParticleConfig):
    def __init__(self, params):
        self.params = params
        self.config_lines =  []
        self.particle = Particle()
        
    def generate_config(self):
        mass, couplings, process_selection = self.params['mass'], self.params['couplings'], self.params['process_selection']
        # self.config_lines.append("9900015:new = N2 N2 2 0 0 {:.12} 0.0 0.0 0.0 {:.12}  0   1   0   1   0\n".format(mass, self.compute_ctau()))
        may_decay = self.params['HNL_decay']
        datafile = 'hnl_production.yaml'
        with open(datafile, 'rU') as f:
            data = yaml.load(f, Loader=yaml.FullLoader)
        all_channels  = data['channels']
        
        self.add_hnl(mass, couplings, may_decay)
        
        
        histograms = self.make_interpolators('branchingratios.dat')
        self.config_lines.append("Next:numberCount    =  0\n")

        # Inclusive
        # =========

        if process_selection=='inclusive':
            self.setup_pythia_inclusive()
        
        if process_selection=='c':

            selection = data['selections']['c']
            print(selection)
            for cmd in selection['parameters']:
                self.config_lines.append(cmd+'\n')

            # Add new charmed particles
            # -------------------------

            # Select all charmed particles (+ tau lepton)
            c_particles = selection['particles']
            tau_id = 15 # tau- Monte-Carlo ID
            self.add_particles(c_particles + [tau_id], data)

            # Add HNL production channels from charmed particles
            # --------------------------------------------------

            # Find charm and tau decays to HNLs
            c_channels = [ch for ch in all_channels if ch['id'] in c_particles]
            tau_channels = [ch for ch in all_channels if ch['id'] == tau_id]
            # Standard model process: tau+ production from D_s+ decay
            ds_id = 431 # D_s+ Monte-Carlo ID
            ds_tau_br = 0.0548 # SM branching ratio Br(D_s+ -> tau+ nu_tau) (source: PDG 2018)

            # Compute the branching ratio scaling factor, taking into account
            # secondary production from tau
            # Decay chains are encoded as follows:
            #     [(top level id A, [br A -> B, br B -> C, ...]), ...]

            # Most charm particles directly decay to HNLs
            primary_decays = [(ch['id'], [self.get_br(histograms, ch, mass, couplings)])
                            for ch in c_channels]
            # The D_s+ can indirectly produce a HNL by first producing a tau+
            secondary_decays = [(ds_id, [ds_tau_br, self.get_br(histograms, ch, mass, couplings)])
                                for ch in tau_channels]
            all_decays = primary_decays + secondary_decays

            # Compute maximum total branching ratio (to rescale all BRs)
            max_total_br = self.compute_max_total_br(all_decays)
            
            self.exit_if_zero_br(max_total_br, process_selection, mass)
            self.print_scale_factor(1/max_total_br)

            # Add charm decays
            for ch in c_channels:
                self.add_channel(ch, histograms, mass, couplings, 1/max_total_br)
            # Add tau production from D_s+
            # We can freely rescale Br(Ds -> tau) and Br(tau -> N X...) as long as
            # Br(Ds -> tau -> N X...) remains the same.
            # Here, we set Br(tau -> N) to unity to make event generation more efficient.
            # The implicit assumption here is that we will disregard the tau during the analysis.
            total_tau_br = sum(branching_ratios[1] for (_, branching_ratios) in secondary_decays)
            assert(ds_tau_br*total_tau_br <= max_total_br + 1e-12)
            self.config_lines.append("431:addChannel      1  {:.12}    0      -15       16\n"\
                                .format(ds_tau_br*total_tau_br/max_total_br))
            # Add secondary HNL production from tau
            for ch in tau_channels:
                # Rescale branching ratios only if some are non-zero. Otherwise leave them at zero.
                self.add_tau_channel(ch, histograms, mass, couplings, 1/(total_tau_br or 1))

            # Add dummy channels in place of SM processes
            self.fill_missing_channels(max_total_br, all_decays)

            # # List channels to confirm that Pythia has been properly set up
            # P8gen.List(9900015)
        
        if process_selection in ['b', 'bc']:
            selection = data['selections'][process_selection]
            for cmd in selection['parameters']:
                self.config_lines.append(cmd+'\n')
            # Add particles
            particles = selection['particles']
            self.add_particles(particles, data)

            # Find all decay channels
            channels = [ch for ch in all_channels if ch['id'] in particles]
            decays = [(ch['id'], [self.get_br(histograms, ch, mass, couplings)]) for ch in channels]

            # Compute scaling factor
            max_total_br = self.compute_max_total_br(decays)
            self.exit_if_zero_br(max_total_br, process_selection, mass)
            self.print_scale_factor(1/max_total_br)

            # Add beauty decays
            for ch in channels:
                self.add_channel(ch, histograms, mass, couplings, 1/max_total_br)

            # Add dummy channels in place of SM processes
            self.fill_missing_channels(max_total_br, decays)

        if process_selection == 'Z':
            selection = data['selections']['Z']
            print(selection)
            
            tau_id = 15 # tau- Monte-Carlo ID
            self.add_particles([tau_id], data)
            
            for cmd in selection['parameters']:
                self.config_lines.append(cmd+'\n')
                
            Z_particles = selection['particles']

            # self.add_particles(Z_particles, data)

            # Add HNL production channels from charmed particles
            # --------------------------------------------------

            # Find charm and tau decays to HNLs
            Z_channels = [ch for ch in all_channels if ch['id'] in Z_particles]
            
            all_decays = [(ch['id'], [self.get_br(histograms, ch, mass, couplings)])
                            for ch in Z_channels]
            
            max_total_br = self.compute_max_total_br(all_decays)
            
            self.exit_if_zero_br(max_total_br, process_selection, mass)
            self.print_scale_factor(1/max_total_br)

            # Add charm decays
            for ch in Z_channels:
                self.add_channel(ch, histograms, mass, couplings, 1/max_total_br)
                
            self.fill_missing_channels(max_total_br, all_decays)
            
            
        if process_selection == 'W':
            selection = data['selections']['W']
            print(selection)
            
            tau_id = 15 # tau- Monte-Carlo ID
            self.add_particles([tau_id], data)
            
            for cmd in selection['parameters']:
                self.config_lines.append(cmd+'\n')
                
            W_particles = selection['particles']
            
            df_production = self.get_csv_br("W_dec_leps.csv", "24")
   
            
            df_production = df_production.loc[mass]
            df_production = pd.DataFrame(df_production)
            
            i = 0
            for index,elem in df_production.iterrows():
                self.config_lines.append(f"24:addChannel      1  {elem.iloc[0]}    101      9900015       -{index.split()[-1]}\n")
                i+=1
            # --------------------------------------------------

            # # Find charm and tau decays to HNLs
            # W_channels = [ch for ch in all_channels if ch['id'] in W_particles]
            
            # all_decays = [(ch['id'], [self.get_br(histograms, ch, mass, couplings)])
            #                 for ch in Z_channels]
            
            # max_total_br = self.compute_max_total_br(all_decays)
            
            # self.exit_if_zero_br(max_total_br, process_selection, mass)
            # self.print_scale_factor(1/max_total_br)

            # # Add charm decays
            # for ch in W_channels:
            #     self.add_channel(ch, histograms, mass, couplings, 1/max_total_br)
                
            # self.fill_missing_channels(max_total_br, all_decays)
        # ... Autres configurations
        
        return self.config_lines

    
    
    def add_hnl(self, mass, couplings, may_decay):
        hnl_instance = HNL(mass, couplings, debug=True)
        ctau = hnl_instance.computeNLifetime(system="FairShip") * u.c_light * u.cm
    
        self.config_lines.append("9900015:new = N2 N2 2 0 0 {:.12} 0.0 0.0 0.0 {:.12}  0   1   0   1   0\n".format(mass, ctau/u.mm))
        self.config_lines.append("9900015:isResonance = false\n")
        # Configuring decay modes...
        if may_decay:
            self.addHNLdecayChannels(hnl_instance, conffile=os.path.expandvars('DecaySelection.conf'), verbose=False)
        # Finish HNL setup...
        if may_decay:
            self.config_lines.append(f"9900015:mayDecay = on\n")
        else:
            self.config_lines.append(f"9900015:mayDecay = off\n")
        
    def addHNLdecayChannels(self, hnl, conffile=os.path.expandvars('DecaySelection.conf'), verbose=True):
        """
        Configures the HNL decay table in Pythia8
        Inputs:
        - P8Gen: an instance of ROOT.HNLPythia8Generator()
        - hnl: an instance of hnl.HNL()
        - conffile: a file listing the channels one wishes to activate
        """
        # First fetch the list of kinematically allowed decays
        allowed = hnl.allowedChannels()
        # Then fetch the list of desired channels to activate
        wanted = self.load(conffile=conffile, verbose=verbose)
        # Add decay channels
        for dec in allowed:
            if dec not in wanted:
                print('addHNLdecayChannels ERROR: channel not configured!\t', dec)
                quit()
            if allowed[dec] == 'yes' and wanted[dec] == 'yes':
                particles = [p for p in dec.replace('->',' ').split()]
                children = particles[1:]
                childrenCodes = [self.PDGcode(p) for p in children]
                BR = hnl.findBranchingRatio(dec)
                # Take care of Majorana modes
                BR = BR/2.
                codes = ' '.join([str(code) for code in childrenCodes])
                self.config_lines.append('9900015:addChannel =  1 {:.12} 0 {}\n'.format(BR, codes))
                # Charge conjugate modes
                codes = ' '.join([(str(-1*code) if self.particle.pdg.GetParticle(-code)!=None else str(code)) for code in childrenCodes])
                self.config_lines.append('9900015:addChannel =  1 {:.12} 0 {}\n'.format(BR, codes))
                # print "debug readdecay table",particles,children,BR
    
    def add_channel(self, ch, histograms, mass, couplings, scale_factor):
        "Add to PYTHIA a leptonic or semileptonic decay channel to HNL."
        if 'idlepton' in ch:
            br = self.get_br(histograms, ch, mass, couplings)
            if br <= 0: # Ignore kinematically closed channels
                return
            if 'idhadron' in ch: # Semileptonic decay
                self.config_lines.append('{}:addChannel      1  {:.17}   22      {}       9900015   {}\n'.format(ch['id'], br*scale_factor, ch['idlepton'], ch['idhadron']))
            else: # Leptonic decay
                self.config_lines.append('{}:addChannel      1  {:.17}    0       9900015      {}\n'.format(ch['id'], br*scale_factor, ch['idlepton']))
        else: # Wrong decay
            raise ValueError("Missing key 'idlepton' in channel {0}".format(ch))

            
    
    
    def add_tau_channel(self, ch, histograms, mass, couplings, scale_factor):
        "Add to PYTHIA a tau decay channel to HNL."
        if 'idhadron' in ch:
            br = self.get_br(histograms, ch, mass, couplings)
            if br <= 0: # Ignore kinematically closed channels
                return
            if 'idlepton' in ch: # 3-body leptonic decay
                self.config_lines.append('{}:addChannel      1  {:.16}    1531       9900015      {} {}\n'.format(ch['id'], br*scale_factor, ch['idlepton'], ch['idhadron']))
            else: # 2-body semileptonic decay
                self.config_lines.append('{}:addChannel      1  {:.16}    1521       9900015      {}\n'.format(ch['id'], br*scale_factor, ch['idhadron']))
        else:
            raise ValueError("Missing key 'idhadron' in channel {0}".format(ch))
   
   
class HNLbranchings(Particle):
    """
    Lifetime and total and partial decay widths of an HNL
    """
    def __init__(self, mass, couplings, debug=False):
        """
        Initialize with mass and couplings of the HNL

        Inputs:
        mass (GeV)
        couplings (list of [U2e, U2mu, U2tau])
        """
        print('a')
        super().__init__()
        
        self.U2 = couplings
        self.U = [math.sqrt(ui) for ui in self.U2]
        self.MN = mass*u.GeV
        self.CKM = CKMmatrix()
        self.CKMelemSq = {'pi+':self.CKM.Vud**2.,
                        'rho+':self.CKM.Vud**2.,
                        'D_s+':self.CKM.Vcs**2.,
                        'D*_s+':self.CKM.Vcs**2.,
                        'K+':self.CKM.Vus**2.,
                        # Quarks (up,down)
                        # up: 1=u, 2=c, 3=t
                        # down: 1=d, 2=s, 3=b
                        (1,1):self.CKM.Vud**2., (1,2):self.CKM.Vus**2., (1,3):self.CKM.Vub**2.,
                        (2,1):self.CKM.Vcd**2., (2,2):self.CKM.Vcs**2., (2,3):self.CKM.Vcb**2.,
                        (3,1):self.CKM.Vtd**2., (3,2):self.CKM.Vts**2., (3,3):self.CKM.Vtb**2.}
        self.decays = [ 'N -> nu nu nu',
                        'N -> e- e+ nu_e',
                        'N -> e- e+ nu_mu',
                        'N -> e- e+ nu_tau',
                        'N -> e- mu+ nu_mu',
                        'N -> mu- e+ nu_e',
                        'N -> pi0 nu_e',
                        'N -> pi0 nu_mu',
                        'N -> pi0 nu_tau',
                        'N -> pi+ e-',
                        'N -> mu- mu+ nu_e',
                        'N -> mu- mu+ nu_mu',
                        'N -> mu- mu+ nu_tau',
                        'N -> pi+ mu-',
                        'N -> eta nu_e',
                        'N -> eta nu_mu',
                        'N -> eta nu_tau',
                        'N -> rho0 nu_e',
                        'N -> rho0 nu_mu',
                        'N -> rho0 nu_tau',
                        'N -> rho+ e-',
                        'N -> omega nu_e',
                        'N -> omega nu_mu',
                        'N -> omega nu_tau',
                        'N -> rho+ mu-',
                        'N -> eta1 nu_e',
                        'N -> eta1 nu_mu',
                        'N -> eta1 nu_tau',
                        'N -> phi nu_e',
                        'N -> phi nu_mu',
                        'N -> phi nu_tau',
                        'N -> e- tau+ nu_tau',
                        'N -> tau- e+ nu_e',
                        'N -> mu- tau+ nu_tau',
                        'N -> tau- mu+ nu_mu',
                        'N -> D_s+ e-',
                        'N -> D_s+ mu-',
                        'N -> D*_s+ e-',
                        'N -> D*_s+ mu-',
                        'N -> eta_c nu_e',
                        'N -> eta_c nu_mu',
                        'N -> eta_c nu_tau' ]
        if self.MN>=1.:
            self.QCD_corr = self.QCD_correction()
            
        if debug:
            print("HNLbranchings instance initialized with couplings:")
            print("\tU2e   = %s"%self.U2[0])
            print("\tU2mu  = %s"%self.U2[1])
            print("\tU2tau = %s"%self.U2[2])
            print("and mass:")
            print("\tm = %s GeV"%(self.MN))
    
    def sqrt_lambda(self,a,b,c):
        """
        Useful function for decay kinematics. Returns 0 for kinematically forbidden region
        """
        l = a**2 + b**2 + c**2 - 2*a*b - 2*b*c - 2*a*c
        if l<0: 
            return 0
        else:
            return math.sqrt(l)
        
    def QCD_correction(self):
        """
        Returns 3-loops QCD correction to HNL decay width into quarks
        """
        alpha_s = ROOT.TGraph( os.path.expandvars('alpha_s.dat') )
        a_s = alpha_s.Eval(self.MN)
        qcd_corr = a_s / math.pi
        qcd_corr += 5.2 * (a_s / math.pi)**2.
        qcd_corr += 26.4 * (a_s / math.pi)**3.
        return qcd_corr
    
    def Width_3nu(self):
        """
        Returns the HNL decay width into three neutrinos
        """
        width = (c.GF**2.)*(self.MN**5.)*sum(self.U2)/(192.*(u.pi**3.))
        width = 2.*width # Majorana case (charge conjugate channels)
        return width
    
    def Width_nu_f_fbar(self, alpha, beta):
        """
        Returns the HNL tree level decay width into the fermion-antifermion pair and a neutrino (decay through Z boson or Z and W)

        Inputs:
        - alpha is a flavour of neutrino (1,2 or 3), determine the mixing angle U_alpha
        - beta is a flavour of the fermions: (1, 2 or 3) for charge leptons or (4, 5, 6, 7, 8, 9) for quarks
        """
        if (alpha not in [1,2,3]) or (beta not in [1,2,3,4,5,6,7,8,9]):
            print('Width_nu_f_fbar ERROR: unknown channel alpha =',alpha,' beta =',beta)
            quit()
        l = [None,'e-','mu-','tau-','u','d','s','c','b','t']
        x = self.mass(l[beta]) / self.MN
        if x > 0.5: # the decay is kinematically forbidden
            return 0
        width = (c.GF**2.)*(self.MN**5.)*self.U2[alpha-1]/(192.*(math.pi**3.))
        L = 0.
        if x>0.01:
            logContent = (1. - 3.*x**2. - (1.-x**2.)*math.sqrt(1. - 4.*x**2.) ) / ( (x**2.)*(1 + math.sqrt(1. - 4.*x**2.)) )
        else:
            logContent = x**4 + 4.*x**6 + 14.*x**8
        if logContent > 0:
            L = math.log( logContent )
        if beta < 4: # lepton decay
            NZ = 1
            if alpha == beta: # interference case
                C1 = 0.25*(1. + 4.*c.s2thetaw + 8.*c.s2thetaw**2.)
                C2 = 0.5*c.s2thetaw*(2.*c.s2thetaw +1.)
            else: # no interference
                C1 = 0.25*(1. - 4.*c.s2thetaw + 8.*c.s2thetaw**2.)
                C2 = 0.5*c.s2thetaw*(2.*c.s2thetaw -1.)
        else: # decay into quarks
            NZ = 3
            if (beta in [4,7,9]): # up quarks
                C1 = 0.25*(1. - 8./3.*c.s2thetaw + 32./9.*c.s2thetaw**2.)
                C2 = 1./3.*c.s2thetaw*(4./3.*c.s2thetaw - 1.)
            if (beta in [5,6,8]): # down quarks
                C1 = 0.25*(1. - 4./3.*c.s2thetaw + 8./9.*c.s2thetaw**2.)
                C2 = 1./6.*c.s2thetaw*(2./3.*c.s2thetaw - 1.)
        width = width*NZ*( C1*( (1.-14.*x**2. -2.*x**4. -12.*x**6.)*math.sqrt(1-4.*x**2) +12.*x**4. *(-1.+x**4.)*L )
                        + 4.*C2*( x**2. *(2.+10.*x**2. -12.*x**4.) * math.sqrt(1.-4.*x**2) + 6.*x**4. *(1.-2.*x**2+2.*x**4)*L ) )
        width = 2.*width # Majorana case (charge conjugate channels)
        return width
    
    def Integrand(self,xx,xi):
        """
        Function to integrate needed for numerical integration using ROOT. Needed for 3-body decays trough W boson.
        First argument is a list of 1 number, argument. Second is the list of 3 numbers xi = mi/MN
        """
        x = xx[0]
        xl2 = xi[0]**2
        xu2 = xi[1]**2
        xd2 = xi[2]**2
        if x==0: # Workaround for division by zero
            return 0
        res = 1./x
        res *= (x - xl2 - xd2)
        res *= (1. + xu2 - x)
        res *= self.sqrt_lambda(x,xl2,xd2)
        res *= self.sqrt_lambda(1.,x,xu2)
        return res
    
    def I(self,x1,x2,x3):
        """
        Numerical integral needed for 3-body decays trough W boson.
        xi = mi/MN
        """ 
        theFunction = self.Integrand
        func = ROOT.TF1('func',theFunction,0,1,3) # Setting function
        func.SetParameters(x1,x2,x3)
        xmin = (x1 + x3)**2
        xmax = (1. - x2)**2
        wf1 = ROOT.Math.WrappedTF1(func) # Setting simple Gauss integration method
        ig = ROOT.Math.GaussIntegrator()
        ig.SetFunction(wf1)
        ig.SetRelTolerance(0.001)
        res = 12. * ig.Integral(xmin,xmax)
        return res
    
    def Width_l1_l2_nu2(self, alpha, beta):
        """
        Returns the HNL decay width into two different flavour leptons and a neutrino (decay through W boson)

        Inputs:
        - alpha is a flavour of the first lepton (1,2 or 3), determine the mixing angle U_alpha
        - beta is a flavour of the second lepton and neutrino (1, 2 or 3)
        """
        if (alpha not in [1,2,3]) or (beta not in [1,2,3]):
            print('Width_l1_l2_nu2 ERROR: unknown channel alpha =',alpha,' beta =',beta)
            quit()
        if alpha==beta: # The interference case is handled by Width_nu_f_fbar function, workaround for a total width calculation
            return 0
        l = [None,'e-','mu-','tau-']
        x1 = self.mass(l[alpha])/self.MN
        x2 = self.mass(l[beta])/self.MN
        if x1+x2>1: # The decay is kinematically forbidden
            return 0
        width = (c.GF**2.)*(self.MN**5.)*self.U2[alpha-1]/(192.*(math.pi**3.))
        width = width*self.I(x1,x2,0)
        width = 2.*width # Majorana case (charge conjugate channels)
        return width
    
    def Width_l_u_d(self, alpha, beta, gamma):
        """
        Returns the HNL tree level decay width into a charged lepton, up quark and down quark (decay through W boson)

        Inputs:
        - alpha is a flavour of the charged lepton (1,2 or 3), determine the mixing angle U_alpha
        - beta is the up quark generation (1, 2, 3 for u, c, t)
        - gamma is the down quark generation (1, 2, 3 for d, s, b)
        """
        if (alpha not in [1,2,3]) or (beta not in [1,2,3]) or (gamma not in [1,2,3]):
            print('Width_l_u_d ERROR: unknown channel alpha =',alpha,' beta =',beta,' gamma =',gamma)
            quit()
        l = [None,'e-','mu-','tau-']
        u = [None,'u','c','t']
        d = [None,'d','s','b']
        xl = self.mass(l[alpha])/self.MN
        xu = self.mass(u[beta])/self.MN
        xd = self.mass(d[gamma])/self.MN
        if xl+xu+xd>1: # The decay is kinematically forbidden
            return 0
        width = (c.GF**2.)*(self.MN**5.)*self.U2[alpha-1]/(192.*(math.pi**3.))
        width *= 3*self.CKMelemSq[(beta, gamma)]
        width *= self.I(xl,xu,xd)
        width = 2.*width # Majorana case (charge conjugate channels)
        return width

    def Width_H0_nu(self, H, alpha):
        """
        Returns the HNL decay width into neutral meson and neutrino

        Inputs:
        - H is a string (name of the meson)
        - alpha is a flavour of the nu (1, 2 or 3), determine the mixing angle U_alpha
        """
        if (H not in ['pi0','eta','rho0','omega','eta1','phi','eta_c']) or (alpha not in [1,2,3]):
            print('Width_H0_nu ERROR: unknown channel H0 =',H,' alpha =',alpha)
            quit()
        x = self.mass(H)/self.MN
        if x > 1: # The decay is kinematically forbidden
            return 0.
        width = (c.GF**2.)*(c.decayConstant[H]**2.)*(self.MN**3.)*self.U2[alpha-1]/(32.*math.pi)
        if H in ['pi0','eta','eta1','eta_c']: # pseudoscalar mesons
            width *= (1 - x**2.)**2.
        if H in ['rho0','omega','phi']: # vector mesons
            width *= (1. + 2*x**2.)*(1. - x**2.)**2.
            if H=='rho0':
                width *= (1. - 2.*c.s2thetaw)**2.
            if H=='omega':
                width *= (16.*c.s2thetaw**2.)/9.
            if H=='phi':
                width *= (1. - 4./3.*c.s2thetaw)**2.
        width = 2.*width # Majorana case (charge conjugate channels)
        return width

    def Width_H_l(self, H, alpha):
        """
        Returns the HNL decay width into charged meson and lepton

        Inputs:
        - H is a string (name of the positively charged meson)
        - alpha is the lepton flavour (1, 2 or 3), determine the mixing angle U_alpha
        """
        if (H not in ['pi+','rho+','D_s+','D*_s+']) or (alpha not in [1,2,3]):
            print('Width_H_l ERROR: unknown channel H =',H,' alpha =',alpha)
            quit()
        l = [None,'e-','mu-','tau-']
        xl = self.mass(l[alpha])/self.MN
        xh = self.mass(H)/self.MN
        if xl+xh > 1: # The decay is kinematically forbidden
            return 0.
        width = (c.GF**2.)*(c.decayConstant[H]**2.)*self.CKMelemSq[H]*(self.MN**3.)*self.U2[alpha-1]/(16.*math.pi)
        width *= self.sqrt_lambda(1., xl**2, xh**2)
        if H in ['pi+','D_s+']: # pseudoscalar mesons
            width *= ( (1. - xl**2.)**2. - xh**2. *(1. + xl**2.) )
        if H in ['rho+','D*_s+']: # vector mesons
            width *= ( (1. - xl**2.)**2. + xh**2.*(1. + xl**2.) - 2.*xh**4. )
        width = 2.*width # Majorana case (charge conjugate channels)
        return width
    
    def Width_charged_leptons(self):
        """
        Returns the total HNL leptonic decay width with charged leptons in final state
        """
        width = 0.
        for l1 in [1,2,3]:
            for l2 in [1,2,3]:
                width += self.Width_nu_f_fbar(l1,l2)
                width += self.Width_l1_l2_nu2(l1,l2)
        return width
    
    def Width_neutral_mesons(self):
        """
        Returns the total HNL decay width into a neutral meson and a neutrino
        """
        mesons = ['pi0','eta','rho0','omega','eta1','phi','eta_c']
        width = 0.
        for m in mesons:
            for l in [1,2,3]:
                width += self.Width_H0_nu(m,l)
        return width
    
    def Width_charged_mesons(self):
        """
        Returns the total HNL decay width into a charged meson and a charged lepton
        """
        mesons = ['pi+','rho+','D_s+','D*_s+']
        width = 0.
        for m in mesons:
            for l in [1,2,3]:
                width += self.Width_H_l(m,l)
        return width

    def Width_quarks_neutrino(self):
        """
        Returns the total HNL decay width with quarks and neutrino in final state. Uses 3-loops alpha_s correction
        """
        if self.MN<1.: # decay into quarks is not applicable at low HNL mass
            return 0
        width = 0.
        for l in [1,2,3]:
            for q in [4,5,6,7,8,9]:
                width += self.Width_nu_f_fbar(l,q)
        width *= (1. + self.QCD_corr)
        return width

    def Width_quarks_lepton(self):
        """
        Returns the total HNL decay width with quarks and charged lepton in final state. Uses 3-loops alpha_s correction
        """
        if self.MN<1.: # decay into quarks is not applicable at low HNL mass
            return 0
        width = 0.
        for l in [1,2,3]:
            for u in [1,2,3]:
                for d in [1,2,3]:
                    width += self.Width_l_u_d(l,u,d)
        width *= (1. + self.QCD_corr)
        return width

    def NDecayWidth(self):
        """
        Returns the total HNL decay width
        """
        leptonWidth = self.Width_3nu() + self.Width_charged_leptons()
        mesonWidth = self.Width_neutral_mesons() + self.Width_charged_mesons()
        quarkWidth = self.Width_quarks_neutrino() + self.Width_quarks_lepton()
        totalWidth = leptonWidth + max([mesonWidth, quarkWidth])
        return totalWidth

    def findBranchingRatio(self, decayString):
        """
        Returns the branching ratio of the selected HNL decay channel

        Inputs:
        - decayString is a string describing the decay, in the form 'N -> stuff1 ... stuffN'
        """
        br = 0.
        totalWidth = self.NDecayWidth()
        
        if (decayString not in self.decays) and (decayString not in ['N -> hadrons','N -> charged hadrons']):
            print('findBranchingRatio ERROR: unknown decay %s'%decayString)
            quit()
        
        if decayString == 'N -> nu nu nu' or decayString == 'N -> 3nu': br = self.Width_3nu() / totalWidth # inclusive
        if decayString == 'N -> e- e+ nu_e': br = self.Width_nu_f_fbar(1,1) / totalWidth
        if decayString == 'N -> e- e+ nu_mu': br = self.Width_nu_f_fbar(2,1) / totalWidth
        if decayString == 'N -> e- e+ nu_tau': br = self.Width_nu_f_fbar(3,1) / totalWidth
        if decayString == 'N -> e- mu+ nu_mu': br = self.Width_l1_l2_nu2(1,2) / totalWidth
        if decayString == 'N -> mu- e+ nu_e': br = self.Width_l1_l2_nu2(2,1) / totalWidth
        if decayString == 'N -> pi0 nu_e': br = self.Width_H0_nu('pi0',1) / totalWidth
        if decayString == 'N -> pi0 nu_mu': br = self.Width_H0_nu('pi0',2) / totalWidth
        if decayString == 'N -> pi0 nu_tau': br = self.Width_H0_nu('pi0',3) / totalWidth
        if decayString == 'N -> pi+ e-': br = self.Width_H_l('pi+',1) / totalWidth
        if decayString == 'N -> mu- mu+ nu_e': br = self.Width_nu_f_fbar(1,2) / totalWidth
        if decayString == 'N -> mu- mu+ nu_mu': br = self.Width_nu_f_fbar(2,2) / totalWidth
        if decayString == 'N -> mu- mu+ nu_tau': br = self.Width_nu_f_fbar(3,2) / totalWidth
        if decayString == 'N -> pi+ mu-': br = self.Width_H_l('pi+',2) / totalWidth
        if decayString == 'N -> eta nu_e': br = self.Width_H0_nu('eta',1) / totalWidth
        if decayString == 'N -> eta nu_mu': br = self.Width_H0_nu('eta',2) / totalWidth
        if decayString == 'N -> eta nu_tau': br = self.Width_H0_nu('eta',3) / totalWidth
        if decayString == 'N -> rho0 nu_e': br = self.Width_H0_nu('rho0',1) / totalWidth
        if decayString == 'N -> rho0 nu_mu': br = self.Width_H0_nu('rho0',2) / totalWidth
        if decayString == 'N -> rho0 nu_tau': br = self.Width_H0_nu('rho0',3) / totalWidth
        if decayString == 'N -> rho+ e-': br = self.Width_H_l('rho+',1) / totalWidth
        if decayString == 'N -> omega nu_e': br = self.Width_H0_nu('omega',1) / totalWidth
        if decayString == 'N -> omega nu_mu': br = self.Width_H0_nu('omega',2) / totalWidth
        if decayString == 'N -> omega nu_tau': br = self.Width_H0_nu('omega',3) / totalWidth
        if decayString == 'N -> rho+ mu-': br = self.Width_H_l('rho+',2) / totalWidth
        if decayString == 'N -> eta1 nu_e': br = self.Width_H0_nu('eta1',1) / totalWidth
        if decayString == 'N -> eta1 nu_mu': br = self.Width_H0_nu('eta1',2) / totalWidth
        if decayString == 'N -> eta1 nu_tau': br = self.Width_H0_nu('eta1',3) / totalWidth
        if decayString == 'N -> phi nu_e': br = self.Width_H0_nu('phi',1) / totalWidth
        if decayString == 'N -> phi nu_mu': br = self.Width_H0_nu('phi',2) / totalWidth
        if decayString == 'N -> phi nu_tau': br = self.Width_H0_nu('phi',3) / totalWidth
        if decayString == 'N -> e- tau+ nu_tau': br = self.Width_l1_l2_nu2(1,3) / totalWidth
        if decayString == 'N -> tau- e+ nu_e': br = self.Width_l1_l2_nu2(3,1) / totalWidth
        if decayString == 'N -> mu- tau+ nu_tau': br = self.Width_l1_l2_nu2(2,3) / totalWidth
        if decayString == 'N -> tau- mu+ nu_mu': br = self.Width_l1_l2_nu2(3,2) / totalWidth
        if decayString == 'N -> D_s+ e-': br = self.Width_H_l('D_s+',1) / totalWidth
        if decayString == 'N -> D_s+ mu-': br = self.Width_H_l('D_s+',2) / totalWidth
        if decayString == 'N -> D*_s+ e-': br = self.Width_H_l('D*_s+',1) / totalWidth
        if decayString == 'N -> D*_s+ mu-': br = self.Width_H_l('D*_s+',2) / totalWidth
        if decayString == 'N -> eta_c nu_e': br = self.Width_H0_nu('eta_c',1) / totalWidth
        if decayString == 'N -> eta_c nu_mu': br = self.Width_H0_nu('eta_c',2) / totalWidth
        if decayString == 'N -> eta_c nu_tau': br = self.Width_H0_nu('eta_c',3) / totalWidth
        
        if decayString == 'N -> hadrons':
            mesonWidth = self.Width_neutral_mesons() + self.Width_charged_mesons()
            quarkWidth = self.Width_quarks_neutrino() + self.Width_quarks_lepton()
            br = max([mesonWidth, quarkWidth]) / totalWidth
        if decayString == 'N -> charged hadrons':
            br = max([self.Width_charged_mesons(), self.Width_quarks_lepton()]) / totalWidth
        return br

    def allowedChannels(self):
        """
        Returns a dictionary of kinematically allowed decay channels

        Inputs:
        - decayString is a string describing the decay, in the form 'N -> stuff1 ... stuffN'
        """
        m = self.MN
        allowedDecays = {'N -> nu nu nu':'yes'}
        if m > 2.*self.mass('e-'):
            allowedDecays.update({'N -> e- e+ nu_e':'yes'})
            allowedDecays.update({'N -> e- e+ nu_mu':'yes'})
            allowedDecays.update({'N -> e- e+ nu_tau':'yes'})
        if m > self.mass('e-') + self.mass('mu-'):
            allowedDecays.update({'N -> e- mu+ nu_mu':'yes'})
            allowedDecays.update({'N -> mu- e+ nu_e':'yes'})
        if m > self.mass('pi0'):
            allowedDecays.update({'N -> pi0 nu_e':'yes'})
            allowedDecays.update({'N -> pi0 nu_mu':'yes'})
            allowedDecays.update({'N -> pi0 nu_tau':'yes'})
        if m > self.mass('pi+') + self.mass('e-'):
            allowedDecays.update({'N -> pi+ e-':'yes'})
        if m > 2.*self.mass('mu-'):
            allowedDecays.update({'N -> mu- mu+ nu_e':'yes'})
            allowedDecays.update({'N -> mu- mu+ nu_mu':'yes'})
            allowedDecays.update({'N -> mu- mu+ nu_tau':'yes'})
        if m > self.mass('pi+') + self.mass('mu-'):
            allowedDecays.update({'N -> pi+ mu-':'yes'})
        if m > self.mass('eta'):
            allowedDecays.update({'N -> eta nu_e':'yes'})
            allowedDecays.update({'N -> eta nu_mu':'yes'})
            allowedDecays.update({'N -> eta nu_tau':'yes'})
        if m > self.mass('rho0'):
            allowedDecays.update({'N -> rho0 nu_e':'yes'})
            allowedDecays.update({'N -> rho0 nu_mu':'yes'})
            allowedDecays.update({'N -> rho0 nu_tau':'yes'})
        if m > self.mass('rho+') + self.mass('e-'):
            allowedDecays.update({'N -> rho+ e-':'yes'})
        if m > self.mass('omega'):
            allowedDecays.update({'N -> omega nu_e':'yes'})
            allowedDecays.update({'N -> omega nu_mu':'yes'})
            allowedDecays.update({'N -> omega nu_tau':'yes'})
        if m > self.mass('rho+') + self.mass('mu-'):
            allowedDecays.update({'N -> rho+ mu-':'yes'})
        if m > self.mass('eta1'):
            allowedDecays.update({'N -> eta1 nu_e':'yes'})
            allowedDecays.update({'N -> eta1 nu_mu':'yes'})
            allowedDecays.update({'N -> eta1 nu_tau':'yes'})
        if m > self.mass('phi'):
            allowedDecays.update({'N -> phi nu_e':'yes'})
            allowedDecays.update({'N -> phi nu_mu':'yes'})
            allowedDecays.update({'N -> phi nu_tau':'yes'})
        if m > self.mass('e-') + self.mass('tau-'):
            allowedDecays.update({'N -> e- tau+ nu_tau':'yes'})
            allowedDecays.update({'N -> tau- e+ nu_e':'yes'})
        if m > self.mass('mu-') + self.mass('tau-'):
            allowedDecays.update({'N -> mu- tau+ nu_tau':'yes'})
            allowedDecays.update({'N -> tau- mu+ nu_mu':'yes'})
        if m > self.mass('D_s+') + self.mass('e-'):
            allowedDecays.update({'N -> D_s+ e-':'yes'})
        if m > self.mass('D_s+') + self.mass('mu-'):
            allowedDecays.update({'N -> D_s+ mu-':'yes'})
        if m > self.mass('D*_s+') + self.mass('e-'):
            allowedDecays.update({'N -> D*_s+ e-':'yes'})
        if m > self.mass('D*_s+') + self.mass('mu-'):
            allowedDecays.update({'N -> D*_s+ mu-':'yes'})
        if m > self.mass('eta_c'):
            allowedDecays.update({'N -> eta_c nu_e':'yes'})
            allowedDecays.update({'N -> eta_c nu_mu':'yes'})
            allowedDecays.update({'N -> eta_c nu_tau':'yes'})
       
        for decay in self.decays:
            if decay not in allowedDecays:
                allowedDecays.update({decay:'no'})
        return allowedDecays


class HNL(HNLbranchings):
    """
    HNL physics according to the nuMSM
    """
    def __init__(self, mass, couplings, debug=False):
        """
        Initialize with mass and couplings of the HNL

        Inputs:
        mass (GeV)
        couplings (list of [U2e, U2mu, U2tau])
        """
        HNLbranchings.__init__(self, mass, couplings, debug)
    def computeNLifetime(self,system="SI"):
        """
        Compute the HNL lifetime

        Inputs:
        - system: choose between default (i.e. SI, result in s) or FairShip (result in ns)
        """
        self.NLifetime = c.hGeV / self.NDecayWidth()
        if system == "FairShip": self.NLifetime *= 1.e9
        return self.NLifetime