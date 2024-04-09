import streamlit as st
import requests
import os
from PIL import Image
from io import BytesIO


col1, col2 = st.columns([6, 4])  # L'ajustement des proportions dépend de votre contenu et de la taille du logo

# Utilisation de la deuxième colonne pour le logo
with col2:
    st.image("Untitled.jpeg", width=300)  # Ajustez le chemin et la taille selon vos besoins
with col1:
    st.image("pythia-logo-b.png", width = 100)
st.title('Pythia Simulation Interface')
# Sidebar for mode selection
mode = st.sidebar.radio("Mode", ['Configure', 'Run', 'MadGraph'])

# Function to display file content
def display_file_content(filename):
    try:
        with open(filename, "r") as file:
            content = file.read()
            st.text_area("File Content", content, height=300)
    except IOError:
        st.error(f"Error: The file {filename} could not be read.")

def get_and_display_plot(plot_type):
    response = requests.get(f'http://127.0.0.1:5000/get_plot/{plot_type}')
    if response.status_code == 200:
        img = Image.open(BytesIO(response.content))
        st.image(img)
    else:
        st.error("Erreur lors de la récupération du graphique.")


def check_hnl_decay_status(file_path):
    """
    Recherche dans le fichier de configuration Pythia pour trouver le statut de désintégration du HNL.
    
    :param file_path: Chemin vers le fichier de configuration Pythia.
    :return: True si HNL peut se désintégrer, False sinon.
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

    # Sélection du type de particule en dehors du formulaire pour mise à jour dynamique
    type_ = st.selectbox("Type of particle", ['DarkPhoton', 'HNL'], key='type_selector', on_change=lambda: update_selection())

    def update_selection():
        st.session_state['type_'] = st.session_state['type_selector']
        
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
                "MesonMother": MesonMother == 'True'
            }

            # Envoi de la requête à l'API
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
    simulation_params = {
            'config_file': 'pythia_config.cmnd',  # Ensure this path is correct
            'n_events': 10000  # This parameter can be made configurable in Streamlit
        }
    
    simulation_params['n_events'] = st.number_input('Mass (GeV)', value=10000, step = 1, min_value=10, max_value=1000000)

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
    if os.path.exists('simulation_results.csv') or simulation_started:
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
option = st.sidebar.selectbox("Choose an option:", ["", "Display .cmnd File"])

if option == "Display .cmnd File":
    cmnd_file = "pythia_config.cmnd"
    if os.path.exists(cmnd_file):
        display_file_content(cmnd_file)
    else:
        st.sidebar.error("The .cmnd file does not exist.")
