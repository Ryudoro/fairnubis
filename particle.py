from abc import ABC, abstractmethod
import ROOT
import os
import math
import csv
import re
import numpy as np
import scipy
import six
import yaml
import shipunit as u 
import sys

# Interface de Configuration de Particule
class ParticleConfig(ABC):
    @abstractmethod
    def generate_config(self):
        pass
    def compute_ctau(self):
        return 0.1  # Exemple simplifié
    
    def setup_pythia_inclusive(self):
        self.config_lines.append("SoftQCD:inelastic = on\n")
        self.config_lines.append("PhotonCollision:gmgm2mumu = on\n")
        self.config_lines.append("PromptPhoton:all = on\n")
        self.config_lines.append("WeakBosonExchange:all = on\n")

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