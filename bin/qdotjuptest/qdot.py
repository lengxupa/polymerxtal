#import os,sys
#sys.path.insert(0, '/opt/hublib')
import hublib.ui as ui
import ipywidgets as w
import papermill as pm
import qdot_output as qd
from hublib.tool import get_inputs, run_simtool
import uuid
from tqdm.auto import tqdm

with open('intro_page.html') as f:
    intro_text = f.read()
intro = w.HTML(value=intro_text)

nstates = ui.Integer(name='Number of States',
                     desc='Number of eigen-energies to be solved for the system.',
                     min=1,
                     max=20,
                     value=8)
shape = ui.Dropdown(
    name='Shape',
    description='Select the shape of the quantum dot.',
    value='Cuboid',
    options=['Cuboid', 'Cylinder', 'Dome', 'Cone', 'Pyramid'],
    #     cb = update_shape
)
with open("images/cubic_qd_simple.png", "rb") as f:
    image = f.read()
shape_image = w.Image(
    value=image,
    format='png',
    width=200,
    height=200,
)
dotx = ui.Number(name='X Dimensions',
                 desc='Input field for the length of the dot in the x-direction.',
                 min=1,
                 max=50,
                 value=10,
                 units='nm')
doty = ui.Number(name='Y Dimensions',
                 desc='Input field for the length of the dot in the y-direction.',
                 min=1,
                 max=50,
                 value=10.5,
                 units='nm')
dotz = ui.Number(name='Z Dimensions',
                 desc='Input field for the length of the dot in the z-direction.',
                 min=1,
                 max=50,
                 value=5,
                 units='nm')
lattice = ui.Number(
    name='Lattice Constant',
    desc=
    '''In the single-band model considered here, this constant represents the discretization spacing of a simple cubic lattice. The shape of the quantum dot will be composed of integer multiples of this quantity.
        (Users can choose a discretization constant for a different material by clicking on the drop down menu - a small inverted triangle at the end of the input box.)''',
    min=.3,
    max=10,
    value=0.565,
    units='nm')
emass = ui.Number(
    name='Effective Mass',
    desc='''Input field for the effective mass of the quantum dot material.
Typical values are: GaAs=0.067
        		InAs=0.023

A truly 'free' electron would have an effective mass of 1.0
        The host material modifies the forces that are exerted on the valence electrons and these forces cause the electrons to not move like a truly free particle anymore.  Their 'effective mass' is modified.

        Please note that this simple model here is based on a single effective mass model, valence bands are not included and all 'optical' transitions are only considered within one band, so called intraband transitions.

        (Users can choose an effective mass for a different material by clicking on the drop down menu - a small inverted triangle at the end of the input box.)''',
    min=.005,
    max=3,
    value=0.067,
)
Eg = ui.Number(name='Energy Gap',
               desc='''Energy difference between the conduction band minimum and valence band maximum of the material.

        (Users can choose an energy gap for a different material by clicking on the drop down menu - a small inverted triangle at the end of the input box.)''',
               min=0,
               max=20,
               value=1.43,
               units='eV')

single_options = ui.Form([nstates, shape, shape_image, dotx, doty, dotz, lattice, emass, Eg],
                         name='Simple Quantum Dot Options')


def update_qopts(_w, val):
    if val == 'Multi-Layer Quantum Dot':
        structure.children = [qds, multi_options]
    elif val == 'Simple Quantum Dot':
        structure.children = [qds, single_options]
    else:
        structure.children = []


qds = ui.Dropdown(name='Type of Quantum Dot Structure',
                  description="Select between a simple quantum dot structure or a multi-layer quantum dot structure.",
                  value='Simple Quantum Dot',
                  options=['Simple Quantum Dot', 'Multi-Layer Quantum Dot', 'Import Quantum Dot'],
                  cb=update_qopts)

structure = ui.Form([qds, single_options], name='Quantum Dot Structure')

light_pic = w.HTML(value='<img src="images/light.gif">')
theta = ui.Number(
    name='Angle Theta',
    desc=
    'Angle of the light source with respect to the vertical axis (z-axis) of the quantum dot, in degrees.  The value theta=0 is shining straight down on the quantum dot.  The value theta=90 is shining parallel to the x-y plane.',
    min=0,
    max=90,
    value=45,
    units='deg')
phi = ui.Number(
    name='Angle Phi',
    desc=
    'Angle of the light source around the quantum dot, in degrees.  This is the angle of rotation around the vertical axis--i.e., around the z-axis, in the x-y plane.',
    min=0,
    max=360,
    value=0,
    units='deg')
polar = ui.Form([theta, phi], name='Light Polarization')
absolute_Ef = ui.Checkbox(
    'Absolute Fermi Level',
    desc=
    'Setting this to true makes the Fermi level an absolute energy value. If this option is not enabled, then Fermi level value is taken with respect to the lowest eigenstate. For example, if Ef = 0 eV, then it will correspond to the lowest eigenenergy (which could be 0.5eV for example).',
    value=False)
Ef = ui.Number(
    'Electron Fermi Level',
    desc=
    'Input field for the position of the Fermi energy relative to the lowest energy state. States below the Fermi level are occupied with electrons, and can therefore be excited to the unoccupied states above the Fermi level.',
    min=-10,
    max=10,
    value=0,
    units='eV')
temperature = ui.Number('Temperature',
                        desc='Input field for the ambient temperature of the quantum dot.',
                        min=1,
                        max=1000,
                        value=300,
                        units='K')
gamma = ui.Number('State Broadening',
                  desc='Input field for the width of the Lorentzian in the absorption calculation.',
                  min=.0001,
                  max=.1,
                  value=.01)
absorption = ui.Form([absolute_Ef, Ef, temperature, gamma], name='Absorption')

sweep = ui.Dropdown(
    name='Sweep Parameter',
    description='Input field for the sweep parameter.',
    value="Angle theta in units of 'degree'",
    options=[
        "Fermi level in units of 'eV'", "Temperature in units of 'K'", "Angle theta in units of 'degree'",
        "Angle phi in units of 'degree'"
    ],
)

sweep_min = ui.Number('Minimum', desc='Minimum value of the parameter being swept.', value=0)
sweep_max = ui.Number('Maximum', desc='Maximum value of the parameter being swept.', value=90)
npts = ui.Integer('Number of Points', desc='Number of points between min/max in the sweep.', min=2, max=20, value=3)

simulate_button = w.Button(description='Simulate', style={'button_color': 'lightblue'})


# general tab callback to trigger resize event
def tab_change(w):
    w['owner'].layout.width = '99%'
    w['owner'].layout.width = '100%'


out_tab = None


def do_simulate(b):
    global out_tab
    run_folder = 'RUNS'
    simulate_button.description = 'Running'
    simulate_button.style = {'button_color': 'orange'}

    nb = '/apps/simtools/qdot/qdot_simtool.ipynb'
    inputs = get_inputs(nb)

    run_name = str(uuid.uuid4()) + '_run.ipynb'
    #run_name = '29737d23-5099-4e3b-a6bd-622977bfce11_run.ipynb'
    simulate_out.clear_output()
    with simulate_out:
        print("Running Simulation: ", end='')
        run_simtool(nb, run_name, inputs, outdir=run_folder, append=True)
        print('Plotting')
        pbar = tqdm(total=12, bar_format='{n}/|/{n_fmt}/{total_fmt} {postfix}')
        display(pbar)

    notebook = run_folder + '/' + run_name
    execution = pm.read_notebook(notebook)

    children = []
    title = []

    title.append('Parameters')
    pbar.set_postfix_str(title[-1])
    out = w.Output(layout={'width': '100%', 'height': '100%'})
    with out:
        display(inputs)
    children.append(out)
    pbar.update(1)

    title.append('Energy Levels')
    pbar.set_postfix_str(title[-1])
    children.append(qd.fig_Energy_states(execution))
    pbar.update(1)

    title.append('Energy Values')
    pbar.set_postfix_str(title[-1])
    children.append(qd.fig_Energy_values(execution, inputs.nstates.value))
    pbar.update(1)

    title.append('X_polarized')
    pbar.set_postfix_str(title[-1])
    children.append(qd.fig_X_polarized(execution))
    pbar.update(1)

    title.append('Y_polarized')
    pbar.set_postfix_str(title[-1])
    children.append(qd.fig_Y_polarized(execution))
    pbar.update(1)

    title.append('Z_polarized')
    pbar.set_postfix_str(title[-1])
    children.append(qd.fig_Z_polarized(execution))
    pbar.update(1)

    title.append('Angle_polarized')
    pbar.set_postfix_str(title[-1])
    children.append(qd.fig_Angle_polarized(execution, inputs.phi.value, inputs.theta.value))
    pbar.update(1)

    title.append('Absorption')
    pbar.set_postfix_str(title[-1])
    children.append(qd.fig_Absorption(execution, inputs.phi.value, inputs.theta.value))
    pbar.update(1)

    title.append('Absorption sweep')
    pbar.set_postfix_str(title[-1])
    children.append(qd.fig_Absorption_sweep(execution, inputs.sweep.value))
    pbar.update(1)

    title.append('Integrated Absorption')
    pbar.set_postfix_str(title[-1])
    children.append(qd.fig_Integrated_absorption(execution, inputs.sweep.value))
    pbar.update(1)

    title.append('input_deck')
    pbar.set_postfix_str(title[-1])
    children.append(qd.fig_Inputdeck(execution))
    pbar.update(1)

    title.append('Output Log')
    pbar.set_postfix_str(title[-1])
    children.append(qd.fig_Log(execution))
    pbar.update(1)

    titles = {i: w for i, w in enumerate(title)}
    new_accord = w.Tab(children=children, layout=w.Layout(width='100%', height='100%'), _titles=titles)
    new_accord.observe(tab_change, 'selected_index')

    if out_tab is None:
        out_tab = w.Tab([], layout={'width': '100%', 'height': '100%'})
        tab.children = [*tab.children, out_tab]
        tab.set_title(len(tab.children) - 1, 'Output')

    out_tab.children = [*out_tab.children, new_accord]
    out_tab.set_title(len(out_tab.children) - 1, run_name)

    simulate_button.description = 'Simulate'
    simulate_button.style = {'button_color': 'lightblue'}
    with simulate_out:
        print('Done')
    pbar.close()


simulate_button.on_click(do_simulate)
simulate_out = w.Output(layout={'width': '100%'})

Sweep = ui.Form([sweep, sweep_min, sweep_max, npts], name='Sweep')
optical = ui.Form([light_pic, polar, absorption, Sweep], name='Light Source')
simulate = w.VBox([simulate_button, simulate_out], layout={'width': '100%', 'height': '100%'})

tab = w.Tab(layout={'height': '100%', 'width': '100%'})
tab.children = [intro, structure, optical, simulate]
tab.set_title(0, 'Introduction')
tab.set_title(1, 'Structure')
tab.set_title(2, 'Optical')
tab.set_title(3, 'Simulate')
tab
