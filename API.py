from flask import Flask, request, jsonify, send_file
from pythiaconf import PythiaSimulation 
from pythiaconf import ParticleConfigFactory
import matplotlib.pyplot as plt
import pandas as pd
from io import BytesIO
import os,sys
import math
import numpy as np
from collections import Counter
import subprocess

print(os.getcwd())
cfg = open("pythia8/Makefile.inc")
lib = "../lib"
for line in cfg:
    if line.startswith("PREFIX_LIB="): lib = line[11:-1]; break
sys.path.insert(0, lib)
print(lib)
import pythia8

def plot_to_image(data, title, xlabel, ylabel, plot_type='histogram', labels=None, rotation=0, figsize=(10, 6)):
    """
    Generates a plot (histogram or bar plot) and saves it to a BytesIO object.

    Parameters:
        data (list or tuple): Data to be plotted. For histograms, a list of values; for bar plots, a tuple (positions, counts).
        title (str): The title of the plot.
        xlabel (str): The x-axis label.
        ylabel (str): The y-axis label.
        plot_type (str): 'histogram' or 'bar' to specify the type of plot.
        labels (list): Labels for the bars in a bar plot, if applicable.
        rotation (int): Rotation angle for x-axis labels.
        figsize (tuple): Size of the plot, e.g., (width, height).
    """
    plt.figure(figsize=figsize)
    
    if plot_type == 'histogram':
        bins = np.histogram_bin_edges(data, bins='auto')
        plt.hist(data, bins=bins, alpha=0.7, edgecolor='black', linewidth=1.2, color='royalblue')
    elif plot_type == 'bar':
        positions, counts = data
        plt.bar(positions, counts, align='center', alpha=0.7, color='royalblue', edgecolor='black')
        if labels is not None:
            plt.xticks(positions, labels, rotation=rotation)
    
    plt.title(title, fontsize=16)
    plt.xlabel(xlabel, fontsize=14)
    plt.ylabel(ylabel, fontsize=14)
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    
    img = BytesIO()
    plt.savefig(img, format='png', bbox_inches="tight")
    img.seek(0)
    plt.close()
    return img


def run_pythia_simulation(config_file, n_events):
    pythia = pythia8.Pythia()
    pythia.readFile(config_file)
    pythia.readString(f"Main:numberOfEvents = {n_events}")
    pythia.init()

    hnl_mother_ids = []
    hnl_energies = []
    hnl_etas = []
    hnl_phis = []
    hnl_distances_R = []

    for event in range(n_events):
        if not pythia.next():
            continue

        for i in range(pythia.event.size()):
            particle = pythia.event[i]
            if particle.id() == 4900023 or particle.id() == 9900015:  # Supposons que 9900015 est l'ID de HNL
                hnl_energies.append(particle.e())
                hnl_etas.append(particle.eta())
                hnl_phis.append(particle.phi())
                distance_R = math.sqrt(particle.px()**2 + particle.py()**2)
                hnl_distances_R.append(distance_R)

                if particle.mother1() != 0:
                    mother_id = pythia.event[particle.mother1()].id()
                    hnl_mother_ids.append(mother_id)
                if particle.mother2() != 0 and particle.mother2() != particle.mother1():
                    mother_id = pythia.event[particle.mother2()].id()
                    hnl_mother_ids.append(mother_id)

    mother_counts = Counter(hnl_mother_ids)
    mothers = list(mother_counts.keys())
    counts = list(mother_counts.values())

    df_mothers = pd.DataFrame({
        'MotherID': mothers,
        'Counts': counts
    })
    
    if not os.path.exists("Results"):
        os.makedirs("Results")
    df_mothers.to_csv(os.path.join("Results", f"result_{config_file[:-5]}_mother.csv"), index=False)

    df = pd.DataFrame({
        'Energy': hnl_energies,
        'Eta': hnl_etas,
        'Phi': hnl_phis,
        'DistanceR': hnl_distances_R
    })
    df.to_csv(os.path.join("Results", f"result_{config_file[:-5]}_hnl.csv"), index=False)


    return "Simulation completed and results saved."


app = Flask(__name__)

@app.route('/simulate', methods=['POST'])
def simulate():
    """
    Endpoint to initiate a simulation. Expects JSON data specifying the particle type and parameters.

    Returns:
        JSON: A confirmation message indicating whether the simulation was successful.
    """
    data = request.json

    convert = {
        "HNL" : 9900015,
        "DarkPhoton" : 4900023
    }
    hnl_params = {
        "type" : data.get("type", 'HNL'),
        "mass": data.get("mass", 1.0),
        "couplings": data.get("coupling", [0.447e-9, 7.15e-9, 1.88e-9]),
        "process_selection": data.get("process", "qcd"),
        "HNL_decay": data.get("HNL_decay", False),
        "epsilon": data.get("epsilon", 0.00000008),
        "mothermode": data.get("MesonMother", True)
    }
    path = data.get("path", os.path.join("cmnd", "pythia_config.cmnd"))
    
    if not os.path.exists("cmnd"):
        os.makedirs("cmnd")
    hnl_config = ParticleConfigFactory.get_particle_config(hnl_params['type'], hnl_params)
    simulation = PythiaSimulation(hnl_config, path)
    simulation.setup_simulation()


    return jsonify({"message": "Simulation completed"})


@app.route('/run_simulation', methods=['POST'])
def run_simulation():
    """
    Configures and runs the PYTHIA simulation based on the provided config file and number of events.

    Parameters:
        config_file (str): Path to the PYTHIA configuration file.
        n_events (int): Number of events to simulate.

    Returns:
        str: Confirmation message that the simulation has been completed and results saved.
    """
    data = request.json
    config_files = data.get('config_file')
    n_events = data.get('n_events')
    
    for config_file in config_files:
        result = run_pythia_simulation(config_file, n_events)

    return jsonify(result)



@app.route('/run_simulation_hep_mc', methods=['POST'])
def run_simulation_hep_mc():
    """
    Endpoint to run a HepMC-format simulation, compiling and executing simulation code.

    Returns:
        JSON: Simulation output or error messages.
    """
    simulation_params = request.get_json()

    if not os.path.exists("hepmc"):
        os.makedirs("hepmc")
    if not os.path.exists("lhe"):
        os.makedirs("lhe")
        
    
    infiles = simulation_params.get('config_file', [os.path.join("cmnd",'pythia_config.cmnd')])
    
    suffix = simulation_params.get('suffix', '')
    mode = simulation_params.get('mode', 'hnl')
    nevents = simulation_params.get('n_events', '100')

    print(nevents)
    project_directory = ''
    build_directory = os.path.join(project_directory, 'build')
    
    make_command = ['make']
    make_result = subprocess.run(make_command, capture_output=True, text=True, cwd='.')

    if make_result.returncode != 0:
        print(make_result)
        return {'message': 'Compilation failed', 'error': make_result.stderr}, 500
    
    
    for infile in infiles:
        cmnd_file = str(os.path.join("cmnd", infile))
        # outfileLHE = simulation_params.get('outfileLHE', '')
        # outfileHepMC = simulation_params.get('outfileHepMC', 'output.hepmc')
        outfile_LHE = str(os.path.join("lhe", f"lhe_result_{infile[:-5]}.lhe"))
        outfile_HepMC = str(os.path.join("hepmc", f"hepmc_result_{infile[:-5]}.hepmc"))
        print(cmnd_file, outfile_LHE)
        command = [
            './build/run',
            '--infile', cmnd_file,
            '--outfileLHE', outfile_LHE,
            '--outfileHepMC', outfile_HepMC,
            '--suffix', suffix,
            '--mode', mode,
            '--nevents', str(nevents)
        ]
        
        result = subprocess.run(command, capture_output=True, text=True)

        if result.returncode == 0:
            continue
        else:
            print(result)
            return {'message': 'Simulation failed', 'error': result.stderr}, 500
    return {'message': 'Simulation completed successfully', 'output': result.stdout}, 200
    

@app.route('/get_plot/<plot_type>')
def get_plot(plot_type):
    """
    Endpoint to retrieve a plot based on the specified type, reading data from previously generated CSV files.

    Returns:
        PNG Image: Generated plot image or error message for invalid types.
    """
    if plot_type == "mothers":
        df_mothers = pd.read_csv('simulation_mothers.csv')
        positions = range(len(df_mothers))
        img = plot_to_image((positions, df_mothers['Counts'].tolist()), "Bar plot des particules mères des HNLs", "ID des particules mères", "Nombre de HNLs", plot_type='bar', labels=df_mothers['MotherID'].tolist(), rotation=45)
    else:
        df = pd.read_csv('simulation_results.csv')
        if plot_type == "energy":
            img = plot_to_image(df['Energy'], 'Energy of HNL', 'Energy [GeV]', 'Particles count', plot_type='histogram')
        elif plot_type == "eta":
            img = plot_to_image(df['Eta'], 'Eta Distribution of HNL', 'Eta', 'Number of Particles', plot_type='histogram')
        elif plot_type == "phi":
            img = plot_to_image(df['Phi'], 'Phi Distribution of HNL', 'Phi', 'Number of Particles', plot_type='histogram')
        elif plot_type == "distanceR":
            img = plot_to_image(df['DistanceR'], 'Distance R Distribution from IP for HNL', 'Distance R [m]', 'Number of Particles', plot_type='histogram')
        else:
            return "Invalid plot type", 400

    return send_file(img, mimetype='image/png')


if __name__ == '__main__':
    app.run(debug=True)
