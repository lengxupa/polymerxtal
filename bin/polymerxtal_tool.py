import os
import sys
import re
import subprocess
import time
import hublib.tool as tool
import math
import uuid
import numpy as np
from plotly import __version__
from plotly.offline import download_plotlyjs, init_notebook_mode, plot, iplot
from plotly.graph_objs import Figure, FigureWidget
import plotly.figure_factory as FF
import plotly.graph_objs as go
from skimage import measure
import matplotlib
from simtool import findInstalledSimToolNotebooks, searchForSimTool, findSimToolNotebook
from simtool import getSimToolInputs, getSimToolOutputs, Run
import papermill as pm
from ipywidgets import (
    HBox,
    VBox,
    Button,
    Layout,
    FloatProgress,
    Output,
    Accordion,
    Tab,
    Textarea,
    interactive,
)
from IPython.display import clear_output
import hublib.ui as ui
import random, string
from string import Template
import imolecule

init_notebook_mode(connected=True)


def PLUGIN_PLOT_ISOSURFACES(X_COLUMN, Y_COLUMN, Z_COLUMN, V_COLUMN, NUM_CONTOURS):
    PLOT_ARRAY = []
    PLOT_LAYOUT = {}
    X_DATA = np.unique(X_COLUMN)
    Y_DATA = np.unique(Y_COLUMN)
    Z_DATA = np.unique(Z_COLUMN)
    X_DIMENSION = len(X_DATA)
    Y_DIMENSION = len(Y_DATA)
    Z_DIMENSION = len(Z_DATA)
    X_DISTANCE = (max(X_DATA) - min(X_DATA)) / X_DIMENSION
    Y_DISTANCE = (max(Y_DATA) - min(Y_DATA)) / Y_DIMENSION
    Z_DISTANCE = (max(Z_DATA) - min(Z_DATA)) / Z_DIMENSION
    min_val = np.amin(V_COLUMN)
    max_val = np.amax(V_COLUMN)

    vol = V_COLUMN.reshape(int(X_DIMENSION), int(Y_DIMENSION), int(Z_DIMENSION))
    data = []
    num_iso = int(NUM_CONTOURS)
    viridis_cmap = matplotlib.cm.get_cmap("viridis")
    norm = matplotlib.colors.Normalize(vmin=min_val, vmax=max_val)
    show_colorbar = True
    colormap = []
    # plotly hack to include custom colorbars
    trace = go.Scatter3d(
        x=[0 for i in range(num_iso + 2)],
        y=[0 for i in range(num_iso + 2)],
        z=[0 for i in range(num_iso + 2)],
        mode="markers",
        marker=dict(
            size=0.1,
            color=np.linspace(min_val, max_val, num=num_iso + 2),
            colorscale="Viridis",
            opacity=1,
            colorbar=dict(title="$\\Psi^2$"),
        ),
    )
    for value in np.linspace(min_val, max_val, num=num_iso + 2):
        color_rgb = matplotlib.colors.colorConverter.to_rgb(viridis_cmap(norm(value)))
        colormap.append(color_rgb)
        if value > min_val and value < max_val:
            vertices, simplices, faces, values = measure.marching_cubes_lewiner(
                vol, value, (float(X_DISTANCE), float(Y_DISTANCE), float(Z_DISTANCE))
            )
            x, y, z = zip(*vertices)
            fig = FF.create_trisurf(
                x=x,
                y=y,
                z=z,
                plot_edges=False,
                colormap=[color_rgb, color_rgb],
                simplices=simplices,
                show_colorbar=True,
                title="Isosurface",
            )
            fig["data"][0]["opacity"] = 0.4
            show_colorbar = False
            data.append(fig["data"][0])
        data.append(trace)
    PLOT_ARRAY = data
    fig["layout"]["showlegend"] = False
    PLOT_LAYOUT = fig["layout"]
    PLOT_FIG = FigureWidget(data=PLOT_ARRAY, layout=PLOT_LAYOUT)
    return PLOT_FIG


def PLUGIN_PLOT_ISOVOLUMES(X_COLUMN, Y_COLUMN, Z_COLUMN, V_COLUMN, NUM_CONTOURS):
    PLOT_ARRAY = []
    PLOT_LAYOUT = {}
    data = []
    num_iso = int(NUM_CONTOURS)
    colormap = []
    min_val = np.amin(V_COLUMN)
    max_val = np.amax(V_COLUMN)
    n_value = np.linspace(min_val, max_val, num=NUM_CONTOURS + 2)
    fig = {
        "type": "volume",
        "x": X_COLUMN,
        "y": Y_COLUMN,
        "z": Z_COLUMN,
        "value": V_COLUMN,
        "isomin": n_value[1],
        "isomax": n_value[len(n_value) - 2],
        "spaceframe": {"show": False, "fill": 0.0},
        "surface": {"show": True, "fill": 1.0, "count": NUM_CONTOURS},
        "caps": {
            "x": {"show": True, "fill": 1.0},
            "y": {"show": True, "fill": 1.0},
            "z": {"show": True, "fill": 1.0},
        },
        "slices": {
            "x": {
                "show": False,
            },
            "y": {
                "show": False,
            },
            "z": {
                "show": False,
            },
        },
        "colorscale": "Viridis",
        "reversescale": False,
        "opacity": 0.8,
        "opacityscale": "max",
        "lighting": {
            "ambient": 1.0,
            "diffuse": 0.0,
            "specular": 0.0,
            "roughness": 0.0,
            "fresnel": 0.0,
        },
        "lightposition": {
            "x": 10000,
            "y": 10000,
            "z": 0,
        },
    }
    data.append(fig)
    PLOT_ARRAY = data
    PLOT_LAYOUT = {
        "template": "plotly_white",
        "margin": {"t": 5, "b": 5, "l": 5, "r": 5},
    }
    # display(PLOT_ARRAY)
    PLOT_FIG = FigureWidget(data=PLOT_ARRAY, layout=PLOT_LAYOUT)
    return PLOT_FIG


def PLUGIN_PLOT_ISOSLICES(X_COLUMN, Y_COLUMN, Z_COLUMN, V_COLUMN, NUM_CONTOURS):
    PLOT_ARRAY = []
    PLOT_LAYOUT = {}
    data = []
    num_iso = int(NUM_CONTOURS)
    colormap = []
    min_val = np.amin(V_COLUMN)
    max_val = np.amax(V_COLUMN)
    n_value = np.linspace(min_val, max_val, num=NUM_CONTOURS + 2)
    fig = {
        "type": "isosurface",
        "colorscale": "Viridis",
        "surface": {"show": False, "fill": 1.0},
        "spaceframe": {"show": False, "fill": 0.0},
        "slices": {
            "x": {"show": True, "fill": 0.1},
            "y": {"show": True, "fill": 0.1},
            "z": {"show": True, "fill": 0.1},
        },
        "caps": {
            "x": {"show": False, "fill": 1.0},
            "y": {"show": False, "fill": 1.0},
            "z": {"show": False, "fill": 1.0},
        },
        "contour": {"show": False, "width": 1},
        "isomin": n_value[1],
        "isomax": n_value[len(n_value) - 2],
        "value": V_COLUMN,
        "x": X_COLUMN,
        "y": Y_COLUMN,
        "z": Z_COLUMN,
        "opacity": 0.8,
    }
    data.append(fig)
    PLOT_ARRAY = data
    PLOT_LAYOUT = {"template": "plotly_white"}
    PLOT_FIG = FigureWidget(data=PLOT_ARRAY, layout=PLOT_LAYOUT)
    return PLOT_FIG


def PLUGIN_PLOT_2D(X_COLUMNS, Y_COLUMNS, layout=None):
    PLOT_ARRAY = []
    PLOT_LAYOUT = {}
    import numpy as np
    import os
    import plotly.graph_objs as go

    for indx in range(len(Y_COLUMNS)):
        trace = go.Scatter(
            x=X_COLUMNS[indx],
            y=Y_COLUMNS[indx],
        )
        PLOT_ARRAY.append(trace)
    PLOT_LAYOUT = {}
    if layout is not None:
        PLOT_LAYOUT = layout
    PLOT_LAYOUT["template"] = "plotly_white"
    PLOT_FIG = FigureWidget(data=PLOT_ARRAY, layout=PLOT_LAYOUT)
    return PLOT_FIG


##nb = os.environ["TOOLDIR"] + "/" + 'qdot_basic.ipynb'
# nb = findSimToolNotebook(os.environ["TOOLDIR"] + "/" + 'Chain_Helix.ipynb')
simToolName = "polymerxtal"
nb = searchForSimTool("polymerxtal")
inputs = getSimToolInputs(nb)


def UI_SET_VALUE(widget, value):
    widget.value = value


QDOT_UI = {}

QDOT_UI["infinite"] = ui.Checkbox(
    name="Periodic Infinite Chain",
    value=True,
)
QDOT_UI["infinite"].dd.observe(
    lambda obj: UI_SET_VALUE(inputs.infinite, obj.new), names="value"
)
QDOT_UI["polymer_type"] = ui.Dropdown(
    name="Polymer Type",
    value=inputs.polymer_type.value,
    options=inputs.polymer_type.options,
)
QDOT_UI["polymer_type"].dd.observe(
    lambda obj: UI_SET_VALUE(inputs.polymer_type, obj.new), names="value"
)
QDOT_UI["helicity"] = ui.Dropdown(
    name="Helicity", value=inputs.helicity.value, options=inputs.helicity.options
)
QDOT_UI["helicity"].dd.observe(
    lambda obj: UI_SET_VALUE(inputs.helicity, obj.new), names="value"
)
QDOT_UI["monomers"] = ui.Integer(
    name="Monomers",
    value=inputs.monomers.value,
    min=inputs.monomers.min,
    max=inputs.monomers.max,
)
QDOT_UI["monomers"].dd.observe(
    lambda obj: UI_SET_VALUE(inputs.monomers, obj.new), names="value"
)
QDOT_UI["tacticity"] = ui.Dropdown(
    name="Tacticity", value=inputs.tacticity.value, options=inputs.tacticity.options
)
QDOT_UI["tacticity"].dd.observe(
    lambda obj: UI_SET_VALUE(inputs.tacticity, obj.new), names="value"
)
QDOT_UI["chiriality"] = ui.Dropdown(
    name="Chiriality", value=inputs.chiriality.value, options=inputs.chiriality.options
)
QDOT_UI["chiriality"].dd.observe(
    lambda obj: UI_SET_VALUE(inputs.chiriality, obj.new), names="value"
)

QDOT_UI["defect_checkbox"] = ui.Checkbox(
    name="Include Connectivity Defects",
    description = 'Include head-to-head & tail-to-tail connectivity defects in crystals',
    value=True,
)

QDOT_UI["connectivity"] = ui.Number(
    name="Defect Ratio",
    description="Defect Ratio for head-to-head & tail-to-tail connections",
    value=inputs.head_tail_defect_ratio.value,
    min=inputs.head_tail_defect_ratio.min,
    max=inputs.head_tail_defect_ratio.max
)
QDOT_UI["connectivity"].dd.observe(
    lambda obj: UI_SET_VALUE(inputs.head_tail_defect_ratio, obj.new), names="value"
)
QDOT_UI["configs"] = ui.Integer(
    name="Configs",
    description="Number of attempts to find a configuration that does not violate excluded region",
    value=inputs.configs.value,
    min=inputs.configs.min,
    max=inputs.configs.max,
)
QDOT_UI["configs"].dd.observe(
    lambda obj: UI_SET_VALUE(inputs.configs, obj.new), names="value"
)

QDOT_UI["defects"] = VBox(
    [
        QDOT_UI["connectivity"],
        QDOT_UI["configs"]
    ]
)

QDOT_UI["defect"] = VBox(
    [
        QDOT_UI["defect_checkbox"],
        QDOT_UI["defects"]
    ]
)

def on_defect_checkbox(change):
    if change['new']:
        QDOT_UI["defects"].layout.display = ''
    else:
        QDOT_UI["defects"].layout.display = 'none'

QDOT_UI["defect_checkbox"].dd.observe(on_defect_checkbox, names="value")

QDOT_UI["create_lmpdata_file"] = ui.Checkbox(
    name="LAMMPS Data File",
    value=True,
)

ffieldB = ui.Dropdown(options = ['Dreiding'],  #, 'PCFF'
                   value = 'Dreiding', 
                   name = 'Force field')

bondB = ui.Number(value = 1.1, name = 'Bondscale',description='Applied to equilibrium bond lengths',min=1,
    max=2)

chargeB = ui.Dropdown(options = ['Gasteiger'],
                   value = 'Gasteiger', 
                   name = 'Charge')
QDOT_UI["create_lmpinput_file"] = ui.Checkbox(
    name="LAMMPS Input File",
    value=True,
)
QDOT_UI["create_lmpinput_file"].dd.observe(
    lambda obj: UI_SET_VALUE(inputs.create_lmpinput_file, obj.new), names="value"
)
box = VBox([ffieldB, bondB, chargeB, QDOT_UI["create_lmpinput_file"]])

def on_lmpdata_checkbox(change):
    if change['new']:
        box.layout.display = ''
    else:
        box.layout.display = 'none'
    UI_SET_VALUE(inputs.create_lmpdata_file, change.new)
        
QDOT_UI["create_lmpdata_file"].dd.observe(on_lmpdata_checkbox, names="value"
)

QDOT_UI["lmpdata_file"] = VBox([QDOT_UI["create_lmpdata_file"],box])


QDOT_UI["button"] = Button(description="Build")
QDOT_UI["button"].layout = Layout(width="99%")
QDOT_UI["sim_progress"] = Output()

QDOT_UI["a1"] = VBox(
    [
        QDOT_UI["infinite"],
        QDOT_UI["polymer_type"],
        QDOT_UI["helicity"],
        QDOT_UI["monomers"],
        QDOT_UI["tacticity"],
        QDOT_UI["chiriality"],
        QDOT_UI["defect"]
    ]
)

QDOT_UI["a2"] = VBox(
    [
        QDOT_UI["lmpdata_file"]
    ]
)
QDOT_UI["a"] = Accordion(children=[QDOT_UI["a1"], QDOT_UI["a2"]])
QDOT_UI["a"].set_title(0, "Structure")
QDOT_UI["a"].set_title(1, "Output")

QDOT_UI["l1"] = VBox([QDOT_UI["a"], QDOT_UI["button"], QDOT_UI["sim_progress"]])
QDOT_UI["l2"] = Tab(children=[])
QDOT_UI["bs"] = HBox([QDOT_UI["l1"], QDOT_UI["l2"]])
QDOT_UI["l1"].layout = Layout(width="500px", border="1px")
QDOT_UI["l2"].layout = Layout(width="100%", border="1px")
QDOT_UI["bs"].layout = Layout(width="100%", border="1px")

QDOT_UI["display"] = ui.Form(
    [QDOT_UI["bs"]], name="PolymerXtal - Polymer Helix Chain Builder"
)


def fig_Energy_values(execution, states):
    out_data = np.array(execution.read("Eigenfunctions"))
    list_child = [Output() for i in range(states)]
    new_accord = Tab(children=list_child)
    new_accord.layout = Layout(width="100%")
    out_data_e = np.array(execution.read("Energy states"))
    for i, o in enumerate(out_data_e):
        new_accord.set_title(i, "l" + str(o))
    for i in range(states):
        out_data2 = out_data[i]
        out_fig = PLUGIN_PLOT_ISOVOLUMES(
            out_data2[:, 0], out_data2[:, 1], out_data2[:, 2], out_data2[:, 3], 5
        )
        with new_accord.children[i]:
            display(out_fig)
    return new_accord


def fig_Energy_value(execution, state):
    out_data = np.array(execution.read("Eigenfunctions"))
    out_data2 = out_data[state]
    out_fig = PLUGIN_PLOT_ISOVOLUMES(
        out_data2[:, 0], out_data2[:, 1], out_data2[:, 2], out_data2[:, 3], 5
    )
    # out_fig = PLUGIN_PLOT_ISOSURFACES(out_data2[:,0], out_data2[:,1], out_data2[:,2], out_data2[:,3], 5)
    return out_fig


def fig_Energy_states(execution):
    layout = {
        "xaxis": {"title": ""},
        "yaxis": {
            "title": "Energy (eV)",
        },
        "title": "Energy State",
    }
    out_data = np.array(execution.read("Energy states"))
    out_fig = PLUGIN_PLOT_2D(
        [[0, 1] for o in out_data], [[o, o] for o in out_data], layout
    )
    return out_fig


def fig_X_polarized(execution):
    layout = {
        "xaxis": {
            "title": "Energy (eV)",
        },
        "yaxis": {
            "title": "Light/Dark Transition Strength",
            "type": "log",
            "exponentformat": "e",
            "showexponent": "all",
        },
        "title": "Light and Dark Transitions (X-polarized))",
    }
    out_data = np.array(execution.read("X-polarized"))
    out_fig = PLUGIN_PLOT_2D([out_data[:, 0]], [out_data[:, 1]], layout)
    return out_fig


def fig_X_polarized(execution):
    layout = {
        "xaxis": {
            "title": "Energy (eV)",
        },
        "yaxis": {
            "title": "Light/Dark Transition Strength",
            "type": "log",
            "exponentformat": "e",
            "showexponent": "all",
        },
        "title": "Light and Dark Transitions (X-polarized))",
    }
    out_data = np.array(execution.read("X-polarized"))
    out_fig = PLUGIN_PLOT_2D([out_data[:, 0]], [out_data[:, 1]], layout)
    return out_fig


def fig_Y_polarized(execution):
    layout = {
        "xaxis": {
            "title": "Energy (eV)",
        },
        "yaxis": {
            "title": "Light/Dark Transition Strength",
            "type": "log",
            "exponentformat": "e",
            "showexponent": "all",
        },
        "title": "Light and Dark Transitions (Y-polarized))",
    }
    out_data = np.array(execution.read("Y-polarized"))
    out_fig = PLUGIN_PLOT_2D([out_data[:, 0]], [out_data[:, 1]], layout)
    return out_fig


def fig_Z_polarized(execution):
    layout = {
        "xaxis": {
            "title": "Energy (eV)",
        },
        "yaxis": {
            "title": "Light/Dark Transition Strength",
            "type": "log",
            "exponentformat": "e",
            "showexponent": "all",
        },
        "title": "Light and Dark Transitions (Z-polarized))",
    }
    out_data = np.array(execution.read("Z-polarized"))
    out_fig = PLUGIN_PLOT_2D([out_data[:, 0]], [out_data[:, 1]], layout)
    return out_fig


def fig_Angle_polarized(execution):
    layout = {
        "xaxis": {
            "title": "Energy (eV)",
        },
        "yaxis": {
            "title": "Light/Dark Transition Strength",
            "type": "log",
            "exponentformat": "e",
            "showexponent": "all",
        },
        "title": "Light and Dark Transitions (phi="
        + str(inputs.phi.value)
        + " theta="
        + str(inputs.theta.value)
        + "))",
    }
    out_data = np.array(execution.read("Angle-polarized"))
    out_fig = PLUGIN_PLOT_2D([out_data[:, 0]], [out_data[:, 1]], layout)
    return out_fig


def fig_Absorption(execution):
    layout = {
        "xaxis": {
            "title": "Energy (eV)",
        },
        "yaxis": {
            "title": "Absorption",
            "type": "log",
            "exponentformat": "e",
            "showexponent": "all",
        },
        "title": "Absorption (phi="
        + str(inputs.phi.value)
        + " theta="
        + str(inputs.theta.value)
        + "))",
    }
    out_data = np.array(execution.read("Absorption"))
    out_fig = PLUGIN_PLOT_2D([out_data[:, 0]], [out_data[:, 1]], layout)
    return out_fig


def fig_Absorption_sweep(execution):
    layout = {
        "xaxis": {
            "title": "Energy (eV)",
        },
        "yaxis": {
            "title": "Absorption",
            "type": "log",
            "exponentformat": "e",
            "showexponent": "all",
        },
        "title": "Absorption Sweep of " + str(inputs.sweep.value),
    }
    out_data = np.array(execution.read("Absorption Sweep"))
    out_fig = PLUGIN_PLOT_2D(
        [o[:, 0] for o in out_data], [o[:, 1] for o in out_data], layout
    )
    return out_fig


def fig_Integrated_absorption(execution):
    layout = {
        "xaxis": {
            "title": str(inputs.sweep.value),
        },
        "yaxis": {
            "title": "Absorption",
            "type": "log",
            "exponentformat": "e",
            "showexponent": "all",
        },
        "title": "Integrated Absortion",
    }
    out_data = np.array(execution.read("Integrated Absortion"))
    out_fig = PLUGIN_PLOT_2D([out_data[:, 0]], [out_data[:, 1]], layout)
    return out_fig


def fig_Inputdeck(execution):
    out_data = execution.read("input_deck")
    out_fig = Textarea(value=out_data)
    out_fig.layout = Layout(width="100%", height="400px")
    return out_fig


def fig_Log(execution):
    out_data = execution.read("Output Log")
    out_fig = Textarea(value=out_data)
    out_fig.layout = Layout(width="100%", height="400px")
    return out_fig


def fig_Structure(execution):
    pdbFile = execution.read("PDBview", raw=True)
    return imolecule.draw(pdbFile[7:])
    # helix_name = execution.db.read('Helix Name')
    # return imolecule.draw(execution.outdir + "/" + helix_name + ".pdb")  #, camera_type="orthographic")


def fig_PDB(execution):
    pdbFile = execution.read("PDB", raw=True)
    out_data=''
    with open(pdbFile[7:], 'r') as input_text:
        out_data = input_text.read()
    out_fig = Textarea(value=out_data)
    out_fig.layout = Layout(width="100%", height="400px")
    download = ui.Download(pdbFile[7:],  style = 'success',tooltip = 'DOWNLOAD PDB FILE',label = 'Download PDBfile', icon = 'arrow-circle-down')
    out_put = VBox([out_fig, download.w])
    return out_put

def fig_LAMMPSData(execution):
    dataFile = execution.read('LAMMPSDataFile', raw=True)
    out_data=''
    with open(dataFile[7:], 'r') as input_text:
        out_data = input_text.read()
    out_fig = Textarea(value=out_data)
    out_fig.layout = Layout(width="100%", height="400px")
    download = ui.Download(dataFile[7:],  style = 'success',tooltip = 'DOWNLOAD DATA FILE',label = 'Download Datafile', icon = 'arrow-circle-down')
    out_put = VBox([out_fig, download.w])
    return out_put

def fig_DataWarnings(execution):
    datawaringFile = execution.read('Datafile_warnings', raw=True)
    out_data=''
    with open(datawaringFile[7:], 'r') as input_text:
        out_data = input_text.read()
    out_fig = Textarea(value=out_data)
    out_fig.layout = Layout(width="100%", height="400px")
    return out_fig

def fig_LAMMPSInput(execution):
    inputFile = execution.read('LAMMPSinputfile', raw=True)
    out_data=''
    with open(inputFile[7:], 'r') as input_text:
        out_data = input_text.read()
    out_fig = Textarea(value=out_data)
    out_fig.layout = Layout(width="100%", height="400px")
    download = ui.Download(inputFile[7:],  style = 'success',tooltip = 'DOWNLOAD LAMMPS INPUT FILE',label = 'Download Inputfile', icon = 'arrow-circle-down')
    out_put = VBox([out_fig, download.w])
    return out_put

def fig_X6paircoeffs(execution):
    X6pairFile = execution.read('X6paircoeffs', raw=True)
    out_data=''
    with open(X6pairFile[7:], 'r') as input_text:
        out_data = input_text.read()
    out_fig = Textarea(value=out_data)
    out_fig.layout = Layout(width="100%", height="400px")
    download = ui.Download(X6pairFile[7:],  style = 'success',tooltip = 'DOWNLOAD X6paircoeffs.txt',label = 'Download X6paircoeffs.txt', icon = 'arrow-circle-down')
    out_put = VBox([out_fig, download.w])
    return out_put

def showFig(func, param, out, state=None):
    with out:
        clear_output()
        if state == None:
            display(func(param))
        else:
            display(func(param, state))


def on_button_clicked(b):
    run_folder = "RUNS"
    run_name = str(uuid.uuid4()) + "_run.ipynb"
    with QDOT_UI["sim_progress"]:
        clear_output()
        r = Run(nb, inputs, cache=True)

        buttons = []
        output = Output()
        but0 = Button(description="Parameters", layout=Layout(width="auto"))
        but0.on_click(
            lambda event, param=inputs, func=dict, out=output: showFig(func, param, out)
        )
        buttons.append(but0)

        but1 = Button(description="Structure", layout=Layout(width="auto"))
        but1.on_click(
            lambda event, param=r, func=fig_Structure, out=output: showFig(
                func, param, out
            )
        )
        buttons.append(but1)

        but2 = Button(description="PDB", layout=Layout(width="auto"))
        but2.on_click(
            lambda event, param=r, func=fig_PDB, out=output: showFig(
                func, param, out
            )
        )
        buttons.append(but2)

        if 'LAMMPSDataFile' in r.savedOutputs:
            but3 = Button(description="LAMMPS Data", layout=Layout(width="auto"))
            but3.on_click(lambda event, param=r, func=fig_LAMMPSData, out=output: showFig(func, param, out))
            buttons.append(but3)

        if 'Datafile_warnings' in r.savedOutputs:
            but4 = Button(description="Datafile Warnings", layout=Layout(width="auto"))
            but4.on_click(lambda event, param=r, func=fig_DataWarnings, out=output: showFig(func, param, out))
            buttons.append(but4)

        if 'LAMMPSinputfile' in r.savedOutputs:
            but5 = Button(description="LAMMPS Input", layout=Layout(width="auto"))
            but5.on_click(lambda event, param=r, func=fig_LAMMPSInput, out=output: showFig(func, param, out))
            buttons.append(but5)
            
            if 'X6paircoeffs' in r.savedOutputs:
                but6 = Button(description="X6paircoeffs.txt", layout=Layout(width="auto"))
                but6.on_click(lambda event, param=r, func=fig_X6paircoeffs, out=output: showFig(func, param, out))
                buttons.append(but6)
            
        new_box = HBox(
            [
                VBox(buttons, layout=Layout(width="150px")),
                VBox([output], layout=Layout(width="99%")),
            ],
            layout=Layout(width="99%"),
        )

        QDOT_UI["l2"].children = [*QDOT_UI["l2"].children, new_box]
        QDOT_UI["l2"].set_title(len(QDOT_UI["l2"].children) - 1, run_name)


QDOT_UI["button"].on_click(on_button_clicked)
