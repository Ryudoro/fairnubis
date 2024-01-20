from abc import ABC, abstractmethod
import yaml
import shipunit as u 
import ROOT
import os
import math
import csv
import re
import numpy as np
import scipy
import six
from HNL import HNLbranchings
from HNL import HNL
from const_utils import Particle
import sys


# Interface de Configuration de Particule
class ParticleConfig(ABC):
    @abstractmethod
    def generate_config(self):
        pass

# Configuration pour HNL
class HNLConfig(ParticleConfig):
    def __init__(self, params):
        self.params = params
        self.config_lines =  []
        self.pdg = ROOT.TDatabasePDG.Instance()
        self.particle = Particle()
        
    def generate_config(self):
        mass, couplings, process_selection = self.params['mass'], self.params['couplings'], self.params['process_selection']
        self.config_lines.append("9900015:new = N2 N2 2 0 0 {:.12} 0.0 0.0 0.0 {:.12}  0   1   0   1   0\n".format(mass, self.compute_ctau()))
        
        datafile = 'hnl_production.yaml'
        with open(datafile, 'rU') as f:
            data = yaml.load(f, Loader=yaml.FullLoader)
        all_channels  = data['channels']
        
        self.add_hnl(mass, couplings)
        
        if process_selection=='inclusive':
            self.setup_pythia_inclusive()
        
        histograms = self.make_interpolators('branchingratios.dat')
        pdg = ROOT.TDatabasePDG.Instance()
        self.config_lines.append("Next:numberCount    =  0\n")
    
        datafile = 'hnl_production.yaml'
        with open(datafile, 'rU') as f:
            data = yaml.load(f, Loader=yaml.FullLoader)
        all_channels  = data['channels']

        # Inclusive
        # =========

        if process_selection=='inclusive':
            self.setup_pythia_inclusive()
        
        if process_selection=='c':

            selection = data['selections']['c']
            for cmd in selection['parameters']:
                self.config_lines.append(cmd+'\n')
            self.add_hnl(mass, couplings)

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
        
        # ... Autres configurations
        return self.config_lines

    def compute_ctau(self):
        return 0.1  # Exemple simplifié
    
    def setup_pythia_inclusive(self):
        self.config_lines.append("SoftQCD:inelastic = on\n")
        self.config_lines.append("PhotonCollision:gmgm2mumu = on\n")
        self.config_lines.append("PromptPhoton:all = on\n")
        self.config_lines.append("WeakBosonExchange:all = on\n")
    
    def add_hnl(self, mass, couplings ):
        hnl_instance = HNL(mass, couplings, debug=True)
        ctau = hnl_instance.computeNLifetime(system="FairShip") * u.c_light * u.cm
    
        self.config_lines.append("9900015:new = N2 N2 2 0 0 {:.12} 0.0 0.0 0.0 {:.12}  0   1   0   1   0\n".format(mass, ctau/u.mm))
        self.config_lines.append("9900015:isResonance = false\n")
        # Configuring decay modes...
        self.addHNLdecayChannels(hnl_instance, conffile=os.path.expandvars('DecaySelection.conf'), verbose=False)
        # Finish HNL setup...
        self.config_lines.append("9900015:mayDecay = on\n")
        
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

            
    def add_particles(self,particles, data):
        """
        Adds the corresponding particles to PYTHIA.

        `particles` must be a list containing either the particles PDG IDs, or
        their PYTHIA names. The commands needed to add the particles are queried
        from `data`.

        If the particle is not self-conjugate, the antiparticle is automatically
        added by PYTHIA.
        """
        for particle_id in particles:
            # Find particle in database (None: particle not found)
            particle = next((p for p in data['particles']
                            if particle_id in [p['id'], p['name']]), None)
            if particle is None:
                raise ValueError("Could not find particle ID {0} in file {1}"
                                .format(particle, data))
            # Add the particle
            self.config_lines.append(particle['cmd']+'\n')
    
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
    
    def get_br(self,histograms, channel, mass, couplings):
        """
        Utility function used to reliably query the branching ratio for a given
        channel at a given mass, taking into account the correct coupling.
        """
        hist = histograms[channel['decay']]
        coupling = couplings[channel['coupling']]
        normalized_br = hist(mass)
        return normalized_br * coupling
    
    def exit_if_zero_br(self,max_total_br, selection, mass, particle='HNL'):
        if max_total_br <= 0:
            print("No phase space for {0} from {1} at this mass: {2}. Quitting."
                .format(particle, selection, mass))
            sys.exit()
        
    def compute_max_total_br(self,decay_chains):
        """
        This function computes the maximum total branching ratio for all decay chains.

        In order to make the event generation as efficient as possible when
        studying very rare processes, it is necessary to rescale the branching
        ratios, while enforcing the invariant that any total branching ratio must
        remain lower that unity.

        This is accomplished by computing, for each particle, the total branching
        ratio to processes of interest, and then dividing all branching ratios by
        the highest of those.

        Note: the decay chain length must be at most 2.
        """
        # For each top-level charmed particle, sum BR over all its decay chains
        top_level_particles = self.get_top_level_particles(decay_chains)
        total_branching_ratios = [self.compute_total_br(particle, decay_chains)
                                for particle in top_level_particles]
        # Find the maximum total branching ratio
        return max(total_branching_ratios)


    def compute_total_br(self,particle, decay_chains):
        """
        Returns the total branching ratio to HNLs for a given particle.
        """
        return sum(np.prod(branching_ratios)
                for (top, branching_ratios) in decay_chains
                if top == particle)
    
    
    def get_top_level_particles(self,decay_chains):
        """
        Returns the set of particles which are at the top of a decay chain.
        """
        return set(top for (top, branching_ratios) in decay_chains)


    def fill_missing_channels(self, max_total_br, decay_chains, epsilon=1e-6):
        """
        Add dummy channels for correct rejection sampling.

        Even after rescaling the branching ratios, they do not sum up to unity
        for most particles since we are ignoring SM processes.

        This function adds a "filler" channel for each particle, in order to
        preserve the ratios between different branching ratios.
        """
        top_level_particles = self.get_top_level_particles(decay_chains)
        for particle in top_level_particles:
            my_total_br = self.compute_total_br(particle, decay_chains)
            remainder = 1 - my_total_br / max_total_br
            assert(remainder > -epsilon)
            assert(remainder < 1 + epsilon)
            if remainder > epsilon:
                print(particle)
                self.add_dummy_channel(particle, remainder)
    
    def print_scale_factor(self,scaling_factor):
        "Prints the scale factor used to make event generation more efficient."
        print("One simulated event per {0:.4g} meson decays".format(scaling_factor))
            
    def add_dummy_channel(self, particle, remainder):
        """
        Add a dummy channel to PYTHIA, with branching ratio equal to `remainder.`

        The purpose of this function is to compensate for the absence of SM
        channels, which are ignored when investigating rare processes. A dummy
        decay channel is instead added to each particle in order to maintain the
        correct ratios between the branching ratios of each particle to rare
        processes. This is usually combined with a global reweighting of the
        branching ratios.

        In order to keep PYTHIA from complaining about charge conservation, a
        suitable process is added which conserves electric charge.

        All dummy channels can be identified by the presence of a photon among the
        decay products.
        """
        # pdg = self.pdg
        # charge = pdg.charge(particle)
        if particle > 0:
            self.config_lines.append('{}:addChannel      1   {:.16}    0       22      -11\n'.format(particle, remainder))
        elif particle < 0:
            self.config_lines.append('{}:addChannel      1   {:.16}    0       22       11\n'.format(particle, remainder))
        else:
            self.config_lines.append('{}:addChannel      1   {:.16}    0       22      22\n'.format(particle, remainder))


           
    def PDGcode(self, particle):
        """
        Read particle ID from PDG database
        """
        particle = self.particle.PDGname(particle)
        tPart = self.particle.pdg.GetParticle(particle)
        return int(tPart.PdgCode())
    
    def load(self, conffile = os.path.expandvars('DecaySelection.conf'), verbose=True):
        f = open(conffile,'r')
        reader = csv.reader(f, delimiter=':')
        configuredDecays = {}
        for row in reader:
            if not row: continue # skip empty lines
            if str(row[0]).strip().startswith('#'):
                continue # skip comment / header lines
            channel = str(row[0]).strip()
            flag = str(row[-1]).partition('#')[0].strip() # skip line comments
            configuredDecays[channel] = flag
        if verbose:
            print('Activated decay channels (plus charge conjugates): ')
            for channel in configuredDecays.keys():
                if configuredDecays[channel] == 'yes':
                    print('\t'+channel)        
        return configuredDecays

    def make_particles_stable(self, P8gen, above_lifetime):
        # FIXME: find time unit and add it to the docstring
        """
        Make the particles with a lifetime above the specified one stable, to allow
        them to decay in Geant4 instead.
        """
        p8 = P8gen.getPythiaInstance()
        n=1
        while n!=0:
            n = p8.particleData.nextId(n)
            p = p8.particleData.particleDataEntryPtr(n)
            if p.tau0() > above_lifetime:
                command = "{}:mayDecay = false".format(n)
                p8.readString(command)
                print("Pythia8 configuration: Made {} stable for Pythia, should decay in Geant4".format(p.name()))
                
    def make_interpolators(self, filepath, kind='linear'):
        """
        This function reads a file containing branching ratio histograms, and
        returns a dictionary of interpolators of the branching ratios, indexed by
        the decay string.
        """
        histogram_data = self.parse_histograms(filepath)
        histograms = {}
        for (hist_string, (masses, br)) in six.iteritems(histogram_data):
            histograms[hist_string] = scipy.interpolate.interp1d(
                masses, br, kind=kind, bounds_error=False, fill_value=0, assume_sorted=True)
        return histograms

    def parse_histograms(self,filepath):
        """
        This function parses a file containing histograms of branching ratios.

        It places them in a dictionary indexed by the decay string (e.g. 'd_K0_e'),
        as a pair ([masses...], [branching ratios...]), where the mass is expressed
        in GeV.
        """
        with open(filepath, 'r') as f:
            lines = f.readlines()
        # Define regular expressions matching (sub-)headers and data lines
        th1f_exp      = re.compile(r'^TH1F\|.+')
        header_exp    = re.compile(r'^TH1F\|(.+?)\|B(?:R|F)/U2(.+?)\|.+? mass \(GeV\)\|?')
        subheader_exp = re.compile(r'^\s*?(\d+?),\s*(\d+?\.\d+?),\s*(\d+\.\d+)\s*$')
        data_exp      = re.compile(r'^\s*(\d+?)\s*,\s*(\d+\.\d+)\s*$')
        # Locate beginning of each histogram
        header_line_idx = [i for i in range(len(lines)) if th1f_exp.match(lines[i]) is not None]
        # Iterate over histograms
        histograms = {}
        for offset in header_line_idx:
            # Parse header
            mh = header_exp.match(lines[offset])
            if mh is None or len(mh.groups()) != 2:
                raise ValueError("Malformed header encountered: {0}".format(lines[offset]))
            decay_code = mh.group(1)
            # Parse sub-header (min/max mass and number of points)
            ms = subheader_exp.match(lines[offset+1])
            if ms is None or len(ms.groups()) != 3:
                raise ValueError("Malformed sub-header encountered: {0}".format(lines[offset+1]))
            npoints  = int(ms.group(1))
            min_mass = float(ms.group(2))
            max_mass = float(ms.group(3))
            masses = np.linspace(min_mass, max_mass, npoints, endpoint=False)
            branching_ratios = np.zeros(npoints)
            # Now read the data lines (skipping the two header lines)
            for line in lines[offset+2:offset+npoints+1]:
                md = data_exp.match(line)
                if md is None or len(md.groups()) != 2:
                    raise ValueError("Malformed data row encountered: {0}".format(line))
                idx = int(md.group(1))
                br  = float(md.group(2))
                branching_ratios[idx] = br
            histograms[decay_code] = (masses, branching_ratios)
        return histograms
    
    
class ParticleConfigFactory:
    particle_config_map = {"HNL" : HNLConfig}

    @classmethod
    def register_particle_config(cls, particle_type, config_class):
        cls.particle_config_map[particle_type] = config_class

    @staticmethod
    def get_particle_config(particle_type, params):
        config_class = ParticleConfigFactory.particle_config_map.get(particle_type)
        if config_class:
            return config_class(params)
        raise ValueError(f"Particle type '{particle_type}' not supported")

# Décorateur pour enregistrer une configuration
def register_config(particle_type):
    def decorator(config_class):
        ParticleConfigFactory.register_particle_config(particle_type, config_class)
        return config_class
    return decorator

# Classe Principale de Simulation
class PythiaSimulation:
    def __init__(self, particle_config):
        self.particle_config = particle_config
        self.base_config = self.read_base_config("template.cmnd")

    def read_base_config(self, filepath):
        with open(filepath, 'r') as file:
            #print(file.readlines()[0])
            return file.readlines()

    def setup_simulation(self):
        particle_config_lines = self.particle_config.generate_config()
        final_config = self.base_config + particle_config_lines
        self.write_config(final_config)

    def write_config(self, config_lines):
        with open("pythia_config.cmnd", "w") as f:
            for line in config_lines:
                f.write(line)

# Utilisation
hnl_params = {"mass": 1.0, "couplings": [0.447e-9,7.15e-9,1.88e-9], "process_selection": "c"}
hnl_config = ParticleConfigFactory.get_particle_config("HNL", hnl_params)
simulation = PythiaSimulation(hnl_config)
simulation.setup_simulation()