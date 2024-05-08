import streamlit as st
import requests
import os
from PIL import Image
from io import BytesIO


col1, col2 = st.columns([6, 4])  # Adjust proportions based on your content and logo size

# Use the second column for the logo
with col2:
    st.image("Untitled.jpeg", width=300)
with col1:
    st.image("pythia-logo-b.png", width = 100)
st.title('Pythia Simulation Interface')
# Sidebar for mode selection
mode = st.sidebar.radio("Mode", ['Configure', 'Run', 'MadGraph'])

# Function to display file content
def display_file_content(filename):
    """
    Displays the content of a file in a text area within the app.
    
    Args:
        filename (str): The path to the file to be displayed.
    """
    try:
        with open(filename, "r") as file:
            content = file.read()
            st.text_area("File Content", content, height=300)
    except IOError:
        st.error(f"Error: The file {filename} could not be read.")

def get_and_display_plot(plot_type):
    """
    Fetches and displays a plot image from a given URL based on plot type.
    
    Args:
        plot_type (str): The type of plot to retrieve and display.
    """
    response = requests.get(f'http://127.0.0.1:5000/get_plot/{plot_type}')
    if response.status_code == 200:
        img = Image.open(BytesIO(response.content))
        st.image(img)
    else:
        st.error("Erreur lors de la récupération du graphique.")


def check_hnl_decay_status(file_path):
    """
    Searches a Pythia configuration file to determine if HNL decay is enabled.
    
    Args:
        file_path (str): Path to the Pythia configuration file.
    
    Returns:
        bool: True if HNL decay is enabled, False otherwise.
    """
    try:
        with open(file_path, 'r') as file:
            for line in file:
                if line.strip().startswith('9900015:mayDecay'):
                    return 'on' in line.split('=')[1].strip()
        return False  # Retourne False si la ligne spécifique n'est pas trouvée
    except FileNotFoundError:
        print(f"Le fichier {file_path} n'a pas été trouvé.")
        return False
    
# Configure mode
if mode == 'Configure':
    st.subheader('Simulation Configuration')
    
    if 'type_' not in st.session_state:
        st.session_state['type_'] = 'HNL'

    # Particle type selection outside form to dynamically update other options
    type_ = st.selectbox("Type of particle", ['DarkPhoton', 'HNL'], key='type_selector', on_change=lambda: update_selection())

    def update_selection():
        st.session_state['type_'] = st.session_state['type_selector']
    
    ranges = st.radio("Range of masses or single input", ('True', 'False'), index=0)
    if not ranges == 'True':
        with st.form(key='simulation_form'):
            if st.session_state['type_'] == 'DarkPhoton':
                process_selection = st.selectbox('Process Selection for DarkPhoton', ['qcd', 'meson'], key='process_selection_dark')
            else:
                process_selection = st.selectbox('Process Selection for HNL', ['c', 'b', 'W'], key='process_selection_hnl')
        
            if process_selection == "c":
                mass = st.number_input('Mass (GeV)', value=1.0, format="%.2f", max_value=2.)
            elif process_selection == "b":
                mass = st.number_input('Mass (GeV)', value=1.0, format="%.2f", max_value=7.)
            elif process_selection == "meson":
                mass = st.number_input('Mass (GeV)', value=1.0, format="%.2f", max_value=7.)
            elif process_selection == "qcd":
                mass = st.number_input('Mass (GeV)', value=100.0, format="%.2f", min_value=10.)
            elif process_selection == "W":
                mass = st.number_input('Mass (GeV)', value=1.0, format="%.2f", max_value=50.)
            couplings = st.text_input('Couplings (comma-separated)', value='0.447e-9, 7.15e-9, 1.88e-9')
            path = st.text_input("Path to cmnd",value = "pythia_config.cmnd" )
            HNL_decay = st.checkbox('Enable HNL Decay')
            if type_ == "DarkPhoton":
                epsilon = st.number_input('Epsilon Mixing Value for Dark Photon', value=0.00000008, format="%.8f")
                MesonMother = st.radio("Choose DP production meson source", ('True', 'False'), index=0)
            else:
                epsilon = None
                MesonMother = None
            submit_button = st.form_submit_button('Save Configuration')

            if submit_button:
                data = {
                    "type" : st.session_state['type_'],
                    "mass": mass,
                    "coupling": [float(x.strip()) for x in couplings.split(',')],
                    "process": process_selection,
                    "HNL_decay": HNL_decay,
                    "epsilon": epsilon,
                    "MesonMother": MesonMother == 'True',
                    "path" : os.path.join("cmnd", path)
                }

                response = requests.post('http://127.0.0.1:5000/simulate', json=data)
                
                if response.status_code == 200:
                    st.success("Simulation configured successfully!")
                else:
                    st.error(f"Failed to configure simulation. Error: {response.text}")
    else:
        with st.form(key='simulation_form'):
            if st.session_state['type_'] == 'DarkPhoton':
                process_selection = st.selectbox('Process Selection for DarkPhoton', ['qcd', 'meson'], key='process_selection_dark')
            else:
                process_selection = st.selectbox('Process Selection for HNL', ['c', 'b', 'W'], key='process_selection_hnl')
        
            if process_selection == "c":
                min_mass = st.number_input('Min Mass (GeV)', value=1.0, format="%.2f", max_value=2.)
                max_mass = st.number_input('Max Mass (GeV)', value=1.0, format="%.2f", max_value=2.)
                step_mass = st.number_input('Step Mass (GeV)', value=1.0, format="%.2f", max_value=2.)
            elif process_selection == "b":
                min_mass = st.number_input('Min Mass (GeV)', value=1.0, format="%.2f", max_value=7.)
                max_mass = st.number_input('Max Mass (GeV)', value=1.0, format="%.2f", max_value=7.)
                step_mass = st.number_input('Step Mass (GeV)', value=1.0, format="%.2f", max_value=7.)
            elif process_selection == "meson":
                min_mass = st.number_input('Min Mass (GeV)', value=1.0, format="%.2f", max_value=7.)
                max_mass = st.number_input('Max Mass (GeV)', value=1.0, format="%.2f", max_value=7.)
                step_mass = st.number_input('Step Mass (GeV)', value=1.0, format="%.2f", max_value=7.)
            elif process_selection == "qcd":
                min_mass = st.number_input('Min Mass (GeV)', value=100.0, format="%.2f", min_value=10.)
                max_mass = st.number_input('Max Mass (GeV)', value=100.0, format="%.2f", min_value=10.)
                step_mass = st.number_input('Step Mass (GeV)', value=100.0, format="%.2f", min_value=10.)
            elif process_selection == "W":
                min_mass = st.number_input('Min Mass (GeV)', value=1.0, format="%.2f", max_value=50.)
                max_mass = st.number_input('Max Mass (GeV)', value=1.0, format="%.2f", max_value=50.)
                step_mass = st.number_input('Step Mass (GeV)', value=1.0, format="%.2f", max_value=50.)
            couplings = st.text_input('Couplings (comma-separated)', value='0.447e-9, 7.15e-9, 1.88e-9')
            path = st.text_input("Path to cmnd",value = "pythia_config.cmnd" )
            HNL_decay = st.checkbox('Enable HNL Decay')
            if type_ == "DarkPhoton":
                epsilon = st.number_input('Epsilon Mixing Value for Dark Photon', value=0.00000008, format="%.8f")
                MesonMother = st.radio("Choose DP production meson source", ('True', 'False'), index=0)
            else:
                epsilon = None
                MesonMother = None
            submit_button = st.form_submit_button('Save Configuration')
            
            values_list = list(range(int(min_mass*1000), int(max_mass*1000 + 1), int(step_mass*1000)))
            if submit_button:
                for elem in values_list:
                    data = {
                        "type" : st.session_state['type_'],
                        "mass": elem/1000.,
                        "coupling": [float(x.strip()) for x in couplings.split(',')],
                        "process": process_selection,
                        "HNL_decay": HNL_decay,
                        "epsilon": epsilon,
                        "MesonMother": MesonMother == 'True',
                        "path" : os.path.join("cmnd", f"pythia_config_{elem}.cmnd")
                    }

                    response = requests.post('http://127.0.0.1:5000/simulate', json=data)
                    
                if response.status_code == 200:
                    st.success("Simulation configured successfully!")
                else:
                    st.error(f"Failed to configure simulation. Error: {response.text}")
# Run mode
# Mode Run
elif mode == "Run":
    st.subheader('Run Pythia Simulation')

    simulation_started = False
    folder = "cmnd"
    files = os.listdir(folder)
    files_selected = st.multiselect("Choose the files you want to use for simulation", files)
    
    simulation_params = {
            'config_file': files_selected,
            'n_events': 10000
        }
    
    simulation_params['n_events'] = st.number_input('Mass (GeV)', value=10000, step = 1, min_value=10, max_value=1000000)

    # simulation_params["paths"] = files_selected
    
    # Button to start the simulation
    if st.button('Start Simulation for plotting'):
        st.write('Running the simulation...')

        
        # Replace with the appropriate URL of your API
        response = requests.post('http://127.0.0.1:5000/run_simulation', json=simulation_params)
        
        if response.status_code == 200:
            st.success('Simulation completed successfully!')
            simulation_started = True
        else:
            st.error('An error occurred during the simulation.')
    
    if st.button('Start Simulation for download hepmc and lhe files'):
        st.write('Running the simulation...')
        
        
        
        # Replace with the appropriate URL of your API
        response = requests.post('http://127.0.0.1:5000/run_simulation_hep_mc', json=simulation_params)
        
        if response.status_code == 200:
            st.success('Simulation completed successfully!')
            simulation_started = True
        else:
            st.error('An error occurred during the simulation.')
            
    # Check if the simulation results file exists
    
    if not os.path.exists("Results"):
        os.makedirs("Results")
        
    liste_res = os.listdir("Results")
    
    option_res = st.selectbox("Choose an option:", liste_res)
    
    if option_res is not None:
        st.write("Simulation results are available.")

        
        plot_type_options = {
            "Select": None,
            "Show HNL Energy": "energy",
            "Show HNL Mothers": "mothers",
            "Show HNL Eta Distribution": "eta",
            "Show HNL Phi Distribution": "phi",
            "Show HNL Distance R": "distanceR"
        }
        
        plot_type_selection = st.selectbox("Select Plot Type", list(plot_type_options.keys()))

        if plot_type_selection != "Select":
            get_and_display_plot(plot_type_options[plot_type_selection])
        
        

# MadGraph mode
elif mode == 'MadGraph':
    st.subheader('MadGraph Tools')
    # MadGraph related functionalities here

# Sidebar for additional options
st.sidebar.subheader('Additional Options')

if not os.path.exists("cmnd"):
    os.makedirs("cmnd")
    
liste_cmnd = list(os.listdir("cmnd"))

option = st.sidebar.selectbox("Choose an option:", liste_cmnd)

if option != "" and option != None:
    cmnd_file = os.path.join("cmnd",str(option))
    if os.path.exists(cmnd_file):
        display_file_content(cmnd_file)
    else:
        st.sidebar.error("The .cmnd file does not exist.")
