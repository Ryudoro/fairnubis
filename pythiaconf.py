#!/usr/bin/env python3
from HNL import HNLConfig
from DarkPhoton import DarkPhotonConfig
import argparse
    
class ParticleConfigFactory:
    particle_config_map = {"HNL" : HNLConfig, "DarkPhoton" : DarkPhotonConfig}

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
            return file.readlines()

    def setup_simulation(self):
        particle_config_lines = self.particle_config.generate_config()
        final_config = self.base_config + particle_config_lines
        self.write_config(final_config)

    def write_config(self, config_lines):
        with open("pythia_config.cmnd", "w") as f:
            for line in config_lines:
                f.write(line)


if __name__ == '__main__':
    # Utilisation
    
    parser = argparse.ArgumentParser(description="Pythia Simulation for HNL Particles")
    parser.add_argument("--type", type=str, default="HNL", help="particle")
    parser.add_argument("--mass", type=float, default=1., help="Mass of the HNL particle")
    parser.add_argument("--coupling", nargs=3, type=float, default=[0.447e-9, 7.15e-9, 1.88e-9], help="Three couplings for the HNL particle")
    parser.add_argument("--process", default="W", help="Process selection for the simulation")
    parser.add_argument("--HNL_decay", default=False, help="True or False, are we interested in particule decays")
    parser.add_argument("--epsilon", default = 0.00000008, help="epsilon mixing value for DarkPhoton")
    parser.add_argument("--MesonMother",  help="Choose DP production meson source", required=False,  default=True)
    
    0.00000008
    args = parser.parse_args()
    
    hnl_params = {
        "type" : args.type,
        "mass": args.mass,
        "couplings": args.coupling,
        "process_selection": args.process,
        "HNL_decay": args.HNL_decay,
        "epsilon" : args.epsilon,
        "mothermode" : args.MesonMother
    }

    hnl_params['mothermode'] = "eta11"
    hnl_config = ParticleConfigFactory.get_particle_config(hnl_params["type"], hnl_params)
    simulation = PythiaSimulation(hnl_config)
    simulation.setup_simulation()
