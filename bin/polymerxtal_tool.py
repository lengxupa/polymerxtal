# import os
# import sys
# import re
# import subprocess
# import time
# import hublib.tool as tool
# import math
import uuid
import numpy as np

# from plotly import __version__
# from plotly.offline import download_plotlyjs, init_notebook_mode, plot, iplot
from plotly.offline import init_notebook_mode

# from plotly.graph_objs import Figure, FigureWidget
# import plotly.figure_factory as FF
# import plotly.graph_objs as go
# from skimage import measure
# import matplotlib
# from simtool import findInstalledSimToolNotebooks, searchForSimTool, findSimToolNotebook
from simtool import searchForSimTool, getSimToolInputs, Run

# from simtool import getSimToolInputs, getSimToolOutputs, Run
# import papermill as pm
from ipywidgets import HBox, VBox, Button, Layout, Output, Accordion, Tab, Textarea

#    FloatProgress,
#    interactive,
from IPython.display import clear_output
import hublib.ui as ui

# import random, string
# from string import Template
import imolecule

init_notebook_mode(connected=True)

##nb = os.environ["TOOLDIR"] + "/" + 'qdot_basic.ipynb'
# nb = findSimToolNotebook(os.environ["TOOLDIR"] + "/" + 'Chain_Helix.ipynb')
simToolName = "polymerxtal"
nb = searchForSimTool("polymerxtal")
inputs = getSimToolInputs(nb)


def UI_SET_VALUE(widget, value):
    widget.value = value


PXTAL_UI = {}

PXTAL_UI["infinite"] = ui.Checkbox(
    name="Periodic Infinite Chain",
    value=True,
)
PXTAL_UI["infinite"].dd.observe(
    lambda obj: UI_SET_VALUE(inputs.infinite, obj.new), names="value"
)
PXTAL_UI["polymer_type"] = ui.Dropdown(
    name="Polymer Type",
    value=inputs.polymer_type.value,
    options=inputs.polymer_type.options,
)
PXTAL_UI["polymer_type"].dd.observe(
    lambda obj: UI_SET_VALUE(inputs.polymer_type, obj.new), names="value"
)
PXTAL_UI["helicity"] = ui.Dropdown(
    name="Helicity", value=inputs.helicity.value, options=inputs.helicity.options
)
PXTAL_UI["helicity"].dd.observe(
    lambda obj: UI_SET_VALUE(inputs.helicity, obj.new), names="value"
)
PXTAL_UI["monomers"] = ui.Integer(
    name="Monomers",
    value=inputs.monomers.value,
    min=inputs.monomers.min,
    max=inputs.monomers.max,
)
PXTAL_UI["monomers"].dd.observe(
    lambda obj: UI_SET_VALUE(inputs.monomers, obj.new), names="value"
)
PXTAL_UI["tacticity"] = ui.Dropdown(
    name="Tacticity", value=inputs.tacticity.value, options=inputs.tacticity.options
)
PXTAL_UI["tacticity"].dd.observe(
    lambda obj: UI_SET_VALUE(inputs.tacticity, obj.new), names="value"
)
PXTAL_UI["chiriality"] = ui.Dropdown(
    name="Chiriality", value=inputs.chiriality.value, options=inputs.chiriality.options
)
PXTAL_UI["chiriality"].dd.observe(
    lambda obj: UI_SET_VALUE(inputs.chiriality, obj.new), names="value"
)

PXTAL_UI["defect_checkbox"] = ui.Checkbox(
    name="Include Connectivity Defects",
    description="Include head-to-head & tail-to-tail connectivity defects in crystals",
    value=True,
)

PXTAL_UI["connectivity"] = ui.Number(
    name="Defect Ratio",
    description="Defect Ratio for head-to-head & tail-to-tail connections",
    value=inputs.head_tail_defect_ratio.value,
    min=inputs.head_tail_defect_ratio.min,
    max=inputs.head_tail_defect_ratio.max,
)
PXTAL_UI["connectivity"].dd.observe(
    lambda obj: UI_SET_VALUE(inputs.head_tail_defect_ratio, obj.new), names="value"
)
PXTAL_UI["configs"] = ui.Integer(
    name="Configs",
    description="Number of attempts to find a configuration that does not violate excluded region",
    value=inputs.configs.value,
    min=inputs.configs.min,
    max=inputs.configs.max,
)
PXTAL_UI["configs"].dd.observe(
    lambda obj: UI_SET_VALUE(inputs.configs, obj.new), names="value"
)

PXTAL_UI["defects"] = VBox([PXTAL_UI["connectivity"], PXTAL_UI["configs"]])

PXTAL_UI["defect"] = VBox([PXTAL_UI["defect_checkbox"], PXTAL_UI["defects"]])


def on_defect_checkbox(change):
    if change["new"]:
        PXTAL_UI["defects"].layout.display = ""
    else:
        PXTAL_UI["defects"].layout.display = "none"


PXTAL_UI["defect_checkbox"].dd.observe(on_defect_checkbox, names="value")

PXTAL_UI["create_lmpdata_file"] = ui.Checkbox(
    name="LAMMPS Data File",
    value=True,
)

ffieldB = ui.Dropdown(
    options=inputs.ffield.options,
    value=inputs.ffield.value,
    name="Force field",  # , 'PCFF'
)
ffieldB.dd.observe(lambda obj: UI_SET_VALUE(inputs.ffield, obj.new), names="value")

bondB = ui.Number(
    value=inputs.bondscale.value,
    name="Bondscale",
    description="Applied to equilibrium bond lengths",
    min=inputs.bondscale.min,
    max=inputs.bondscale.max,
)
bondB.dd.observe(lambda obj: UI_SET_VALUE(inputs.bondscale, obj.new), names="value")

chargeB = ui.Dropdown(
    options=inputs.charge.options, value=inputs.charge.value, name="Charge"
)
chargeB.dd.observe(lambda obj: UI_SET_VALUE(inputs.charge, obj.new), names="value")
PXTAL_UI["create_lmpinput_file"] = ui.Checkbox(
    name="LAMMPS Input File",
    value=True,
)
PXTAL_UI["create_lmpinput_file"].dd.observe(
    lambda obj: UI_SET_VALUE(inputs.create_lmpinput_file, obj.new), names="value"
)
box = VBox([ffieldB, bondB, chargeB, PXTAL_UI["create_lmpinput_file"]])


def on_lmpdata_checkbox(change):
    if change["new"]:
        box.layout.display = ""
    else:
        box.layout.display = "none"
    UI_SET_VALUE(inputs.create_lmpdata_file, change.new)


PXTAL_UI["create_lmpdata_file"].dd.observe(on_lmpdata_checkbox, names="value")

PXTAL_UI["lmpdata_file"] = VBox([PXTAL_UI["create_lmpdata_file"], box])


PXTAL_UI["button"] = Button(description="Build")
PXTAL_UI["button"].layout = Layout(width="99%")
PXTAL_UI["sim_progress"] = Output()

PXTAL_UI["a1"] = VBox(
    [
        PXTAL_UI["infinite"],
        PXTAL_UI["polymer_type"],
        PXTAL_UI["helicity"],
        PXTAL_UI["monomers"],
        PXTAL_UI["tacticity"],
        PXTAL_UI["chiriality"],
        PXTAL_UI["defect"],
    ]
)

PXTAL_UI["a2"] = VBox([PXTAL_UI["lmpdata_file"]])
PXTAL_UI["a"] = Accordion(children=[PXTAL_UI["a1"], PXTAL_UI["a2"]])
PXTAL_UI["a"].set_title(0, "Structure")
PXTAL_UI["a"].set_title(1, "Output")

PXTAL_UI["l1"] = VBox([PXTAL_UI["a"], PXTAL_UI["button"], PXTAL_UI["sim_progress"]])
PXTAL_UI["l2"] = Tab(children=[])
PXTAL_UI["bs"] = HBox([PXTAL_UI["l1"], PXTAL_UI["l2"]])
PXTAL_UI["l1"].layout = Layout(width="500px", border="1px")
PXTAL_UI["l2"].layout = Layout(width="100%", border="1px")
PXTAL_UI["bs"].layout = Layout(width="100%", border="1px")

PXTAL_UI["display"] = ui.Form(
    [PXTAL_UI["bs"]], name="PolymerXtal - Polymer Helix Chain Builder"
)


def fig_Structure(execution):
    pdbFile = execution.read("PDBview", raw=True)
    return imolecule.draw(pdbFile[7:])
    # helix_name = execution.db.read('Helix Name')
    # return imolecule.draw(execution.outdir + "/" + helix_name + ".pdb")  #, camera_type="orthographic")


def fig_PDB(execution):
    pdbFile = execution.read("PDB", raw=True)
    out_data = ""
    with open(pdbFile[7:], "r") as input_text:
        out_data = input_text.read()
    out_fig = Textarea(value=out_data)
    out_fig.layout = Layout(width="100%", height="400px")
    download = ui.Download(
        pdbFile[7:],
        style="success",
        tooltip="DOWNLOAD PDB FILE",
        label="Download PDBfile",
        icon="arrow-circle-down",
    )
    out_put = VBox([out_fig, download.w])
    return out_put


def fig_LAMMPSData(execution):
    dataFile = execution.read("LAMMPSDataFile", raw=True)
    out_data = ""
    with open(dataFile[7:], "r") as input_text:
        out_data = input_text.read()
    out_fig = Textarea(value=out_data)
    out_fig.layout = Layout(width="100%", height="400px")
    download = ui.Download(
        dataFile[7:],
        style="success",
        tooltip="DOWNLOAD DATA FILE",
        label="Download Datafile",
        icon="arrow-circle-down",
    )
    out_put = VBox([out_fig, download.w])
    return out_put


def fig_DataWarnings(execution):
    datawaringFile = execution.read("Datafile_warnings", raw=True)
    out_data = ""
    with open(datawaringFile[7:], "r") as input_text:
        out_data = input_text.read()
    out_fig = Textarea(value=out_data)
    out_fig.layout = Layout(width="100%", height="400px")
    return out_fig


def fig_LAMMPSInput(execution):
    inputFile = execution.read("LAMMPSinputfile", raw=True)
    out_data = ""
    with open(inputFile[7:], "r") as input_text:
        out_data = input_text.read()
    out_fig = Textarea(value=out_data)
    out_fig.layout = Layout(width="100%", height="400px")
    download = ui.Download(
        inputFile[7:],
        style="success",
        tooltip="DOWNLOAD LAMMPS INPUT FILE",
        label="Download Inputfile",
        icon="arrow-circle-down",
    )
    out_put = VBox([out_fig, download.w])
    return out_put


def fig_X6paircoeffs(execution):
    X6pairFile = execution.read("X6paircoeffs", raw=True)
    out_data = ""
    with open(X6pairFile[7:], "r") as input_text:
        out_data = input_text.read()
    out_fig = Textarea(value=out_data)
    out_fig.layout = Layout(width="100%", height="400px")
    download = ui.Download(
        X6pairFile[7:],
        style="success",
        tooltip="DOWNLOAD X6paircoeffs.txt",
        label="Download X6paircoeffs.txt",
        icon="arrow-circle-down",
    )
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
    with PXTAL_UI["sim_progress"]:
        clear_output()
        r = Run(nb, inputs, cache=True)

        if not hasattr(r, savedOutputs):
            r.savedOutputs = r.db.getSavedOutputs()

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
            lambda event, param=r, func=fig_PDB, out=output: showFig(func, param, out)
        )
        buttons.append(but2)

        if "LAMMPSDataFile" in r.savedOutputs:
            but3 = Button(description="LAMMPS Data", layout=Layout(width="auto"))
            but3.on_click(
                lambda event, param=r, func=fig_LAMMPSData, out=output: showFig(
                    func, param, out
                )
            )
            buttons.append(but3)

        if "Datafile_warnings" in r.savedOutputs:
            but4 = Button(description="Datafile Warnings", layout=Layout(width="auto"))
            but4.on_click(
                lambda event, param=r, func=fig_DataWarnings, out=output: showFig(
                    func, param, out
                )
            )
            buttons.append(but4)

        if "LAMMPSinputfile" in r.savedOutputs:
            but5 = Button(description="LAMMPS Input", layout=Layout(width="auto"))
            but5.on_click(
                lambda event, param=r, func=fig_LAMMPSInput, out=output: showFig(
                    func, param, out
                )
            )
            buttons.append(but5)

            if "X6paircoeffs" in r.savedOutputs:
                but6 = Button(
                    description="X6paircoeffs.txt", layout=Layout(width="auto")
                )
                but6.on_click(
                    lambda event, param=r, func=fig_X6paircoeffs, out=output: showFig(
                        func, param, out
                    )
                )
                buttons.append(but6)

        new_box = HBox(
            [
                VBox(buttons, layout=Layout(width="150px")),
                VBox([output], layout=Layout(width="99%")),
            ],
            layout=Layout(width="99%"),
        )

        PXTAL_UI["l2"].children = [*PXTAL_UI["l2"].children, new_box]
        PXTAL_UI["l2"].set_title(len(PXTAL_UI["l2"].children) - 1, run_name)


PXTAL_UI["button"].on_click(on_button_clicked)
