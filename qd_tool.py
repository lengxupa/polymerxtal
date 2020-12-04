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
import os
from simtool import findInstalledSimToolNotebooks, searchForSimTool, findSimToolNotebook
from simtool import getSimToolInputs, getSimToolOutputs, Run
import numpy as np
import papermill as pm
from ipywidgets import HBox, VBox, Button, Layout, FloatProgress, Output, Accordion, Tab, Textarea, interactive
from IPython.display import clear_output
import hublib.ui as ui
import random, string
from string import Template

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
    viridis_cmap = matplotlib.cm.get_cmap('viridis')
    norm = matplotlib.colors.Normalize(vmin=min_val, vmax=max_val)
    show_colorbar = True
    colormap = []
    #plotly hack to include custom colorbars
    trace = go.Scatter3d(
        x=[0 for i in range(num_iso + 2)],
        y=[0 for i in range(num_iso + 2)],
        z=[0 for i in range(num_iso + 2)],
        mode='markers',
        marker=dict(size=0.1,
                    color=np.linspace(min_val, max_val, num=num_iso + 2),
                    colorscale='Viridis',
                    opacity=1,
                    colorbar=dict(title="$\\Psi^2$")),
    )
    for value in np.linspace(min_val, max_val, num=num_iso + 2):
        color_rgb = matplotlib.colors.colorConverter.to_rgb(viridis_cmap(norm(value)))
        colormap.append(color_rgb)
        if value > min_val and value < max_val:
            vertices, simplices, faces, values = measure.marching_cubes_lewiner(
                vol, value, (float(X_DISTANCE), float(Y_DISTANCE), float(Z_DISTANCE)))
            x, y, z = zip(*vertices)
            fig = FF.create_trisurf(x=x,
                                    y=y,
                                    z=z,
                                    plot_edges=False,
                                    colormap=[color_rgb, color_rgb],
                                    simplices=simplices,
                                    show_colorbar=True,
                                    title="Isosurface")
            fig['data'][0]['opacity'] = 0.4
            show_colorbar = False
            data.append(fig['data'][0])
        data.append(trace)
    PLOT_ARRAY = data
    fig['layout']['showlegend'] = False
    PLOT_LAYOUT = fig['layout']
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
        "spaceframe": {
            "show": False,
            "fill": 0.0
        },
        "surface": {
            "show": True,
            "fill": 1.0,
            "count": NUM_CONTOURS
        },
        "caps": {
            "x": {
                "show": True,
                "fill": 1.0
            },
            "y": {
                "show": True,
                "fill": 1.0
            },
            "z": {
                "show": True,
                "fill": 1.0
            }
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
            }
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
    PLOT_LAYOUT = {"template": "plotly_white", "margin": {"t": 5, "b": 5, "l": 5, "r": 5}}
    #display(PLOT_ARRAY)
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
        "surface": {
            "show": False,
            "fill": 1.0
        },
        "spaceframe": {
            "show": False,
            "fill": 0.0
        },
        "slices": {
            "x": {
                "show": True,
                "fill": 0.1
            },
            "y": {
                "show": True,
                "fill": 0.1
            },
            "z": {
                "show": True,
                "fill": 0.1
            }
        },
        "caps": {
            "x": {
                "show": False,
                "fill": 1.0
            },
            "y": {
                "show": False,
                "fill": 1.0
            },
            "z": {
                "show": False,
                "fill": 1.0
            }
        },
        "contour": {
            "show": False,
            "width": 1
        },
        "isomin": n_value[1],
        "isomax": n_value[len(n_value) - 2],
        "value": V_COLUMN,
        "x": X_COLUMN,
        "y": Y_COLUMN,
        "z": Z_COLUMN,
        "opacity": 0.8
    }
    data.append(fig)
    PLOT_ARRAY = data
    PLOT_LAYOUT = {'template': 'plotly_white'}
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
    PLOT_LAYOUT['template'] = 'plotly_white'
    PLOT_FIG = FigureWidget(data=PLOT_ARRAY, layout=PLOT_LAYOUT)
    return PLOT_FIG


##nb = os.environ["TOOLDIR"] + "/" + 'qdot_basic.ipynb'
nb = findSimToolNotebook(os.environ["TOOLDIR"] + "/" + 'qdot_basic.ipynb')
inputs = getSimToolInputs(nb)


def UI_SET_VALUE(widget, value):
    widget.value = value


QDOT_UI = {}

QDOT_UI['nstates'] = ui.Integer(name='Number of States',
                                value=inputs.nstates.value,
                                min=inputs.nstates.min,
                                max=inputs.nstates.max,
                                units=inputs.nstates.units)
QDOT_UI['nstates'].dd.observe(lambda obj: UI_SET_VALUE(inputs.nstates, obj.new), names='value')
QDOT_UI['shape'] = ui.Dropdown(name='Shape', value=inputs.shape.value, options=inputs.shape.options)
QDOT_UI['shape'].dd.observe(lambda obj: UI_SET_VALUE(inputs.shape, obj.new), names='value')
QDOT_UI['dotx'] = ui.Number(name='X Dimensions',
                            value=inputs.dotx.value,
                            min=inputs.dotx.min,
                            max=inputs.dotx.max,
                            units=inputs.dotx.units)
QDOT_UI['dotx'].dd.observe(lambda obj: UI_SET_VALUE(inputs.dotx, obj.new), names='value')
QDOT_UI['doty'] = ui.Number(name='Y Dimensions',
                            value=inputs.doty.value,
                            min=inputs.doty.min,
                            max=inputs.doty.max,
                            units=inputs.doty.units)
QDOT_UI['doty'].dd.observe(lambda obj: UI_SET_VALUE(inputs.doty, obj.new), names='value')
QDOT_UI['dotz'] = ui.Number(name='Z Dimensions',
                            value=inputs.dotz.value,
                            min=inputs.dotz.min,
                            max=inputs.dotz.max,
                            units=inputs.dotz.units)
QDOT_UI['dotz'].dd.observe(lambda obj: UI_SET_VALUE(inputs.dotz, obj.new), names='value')
QDOT_UI['lattice'] = ui.Number(name='Lattice Constant',
                               value=inputs.lattice.value,
                               min=inputs.lattice.min,
                               max=inputs.lattice.max,
                               units=inputs.lattice.units)
QDOT_UI['lattice'].dd.observe(lambda obj: UI_SET_VALUE(inputs.lattice, obj.new), names='value')
QDOT_UI['emass'] = ui.Number(name='Effective Mass',
                             value=inputs.emass.value,
                             min=inputs.emass.min,
                             max=inputs.emass.max,
                             units=inputs.emass.units)
QDOT_UI['emass'].dd.observe(lambda obj: UI_SET_VALUE(inputs.emass, obj.new), names='value')
QDOT_UI['Eg'] = ui.Number(name='Energy Gap',
                          value=inputs.Eg.value,
                          min=inputs.Eg.min,
                          max=inputs.Eg.max,
                          units=inputs.Eg.units)
QDOT_UI['Eg'].dd.observe(lambda obj: UI_SET_VALUE(inputs.Eg, obj.new), names='value')

QDOT_UI['theta'] = ui.Number(name='Angle Tetha',
                             value=inputs.theta.value,
                             min=inputs.theta.min,
                             max=inputs.theta.max,
                             units=inputs.theta.units)
QDOT_UI['theta'].dd.observe(lambda obj: UI_SET_VALUE(inputs.theta, obj.new), names='value')
QDOT_UI['phi'] = ui.Number(name='Angle Phi',
                           value=inputs.phi.value,
                           min=inputs.phi.min,
                           max=inputs.phi.max,
                           units=inputs.phi.units)
QDOT_UI['phi'].dd.observe(lambda obj: UI_SET_VALUE(inputs.phi, obj.new), names='value')
#QDOT_UI['absolute_fermi'] = ui.Number(name='Absolute Fermi Level',value=inputs.Eg.value,min=inputs.Eg.min,max=inputs.Eg.max,units=inputs.Eg.units )
QDOT_UI['Ef'] = ui.Number(name='Electron Fermi Level',
                          value=inputs.Ef.value,
                          min=inputs.Ef.min,
                          max=inputs.Ef.max,
                          units=inputs.Ef.units)
QDOT_UI['Ef'].dd.observe(lambda obj: UI_SET_VALUE(inputs.Ef, obj.new), names='value')
QDOT_UI['temperature'] = ui.Number(name='Temperature',
                                   value=inputs.temperature.value,
                                   min=inputs.temperature.min,
                                   max=inputs.temperature.max,
                                   units=inputs.temperature.units)
QDOT_UI['temperature'].dd.observe(lambda obj: UI_SET_VALUE(inputs.temperature, obj.new), names='value')
#QDOT_UI['state_broadening'] = ui.Number(name='State Broadening',value=inputs.Eg.value,min=inputs.Eg.min,max=inputs.Eg.max,units=inputs.Eg.units )
QDOT_UI['sweep'] = ui.Dropdown(name='Sweep Parameter', value=inputs.sweep.value, options=inputs.sweep.options)
QDOT_UI['sweep'].dd.observe(lambda obj: UI_SET_VALUE(inputs.sweep, obj.new), names='value')
QDOT_UI['smin'] = ui.Number(name='Minimum',
                            value=inputs.smin.value,
                            min=inputs.smin.min,
                            max=inputs.smin.max,
                            units=inputs.smin.units)
QDOT_UI['smin'].dd.observe(lambda obj: UI_SET_VALUE(inputs.smin, obj.new), names='value')
QDOT_UI['smax'] = ui.Number(name='Maximum',
                            value=inputs.smax.value,
                            min=inputs.smax.min,
                            max=inputs.smax.max,
                            units=inputs.smax.units)
QDOT_UI['smax'].dd.observe(lambda obj: UI_SET_VALUE(inputs.smax, obj.new), names='value')
QDOT_UI['spts'] = ui.Integer(name='Number of Points',
                             value=inputs.spts.value,
                             min=inputs.spts.min,
                             max=inputs.spts.max)
QDOT_UI['spts'].dd.observe(lambda obj: UI_SET_VALUE(inputs.spts, obj.new), names='value')

QDOT_UI['button'] = Button(description='Simulate')
QDOT_UI['button'].layout = Layout(width='99%')
QDOT_UI['sim_progress'] = Output()

QDOT_UI['a1'] = VBox([
    QDOT_UI['nstates'],
    QDOT_UI['shape'],
    QDOT_UI['dotx'],
    QDOT_UI['doty'],
    QDOT_UI['dotz'],
    QDOT_UI['lattice'],
    QDOT_UI['emass'],
    QDOT_UI['Eg'],
])
QDOT_UI['a2'] = VBox([
    QDOT_UI['theta'],
    QDOT_UI['phi'],
    QDOT_UI['Ef'],
    QDOT_UI['temperature'],
    QDOT_UI['sweep'],
    QDOT_UI['smin'],
    QDOT_UI['smax'],
    QDOT_UI['spts'],
])
QDOT_UI['a'] = Accordion(children=[QDOT_UI['a1'], QDOT_UI['a2']])
QDOT_UI['a'].set_title(0, 'Structure')
QDOT_UI['a'].set_title(1, 'Optical')

QDOT_UI['l1'] = VBox([QDOT_UI['a'], QDOT_UI['button'], QDOT_UI['sim_progress']])
QDOT_UI['l2'] = Tab(children=[])
QDOT_UI['bs'] = HBox([QDOT_UI['l1'], QDOT_UI['l2']])
QDOT_UI['l1'].layout = Layout(width='500px', border='1px')
QDOT_UI['l2'].layout = Layout(width='100%', border='1px')
QDOT_UI['bs'].layout = Layout(width='100%', border='1px')

QDOT_UI['display'] = ui.Form([QDOT_UI['bs']], name='Quantum Dot Lab')


def fig_Energy_values(execution, states):
    out_data = np.array(execution.read('Eigenfunctions'))
    list_child = [Output() for i in range(states)]
    new_accord = Tab(children=list_child)
    new_accord.layout = Layout(width='100%')
    out_data_e = np.array(execution.read('Energy states'))
    for i, o in enumerate(out_data_e):
        new_accord.set_title(i, "l" + str(o))
    for i in range(states):
        out_data2 = out_data[i]
        out_fig = PLUGIN_PLOT_ISOVOLUMES(out_data2[:, 0], out_data2[:, 1], out_data2[:, 2], out_data2[:, 3], 5)
        with new_accord.children[i]:
            display(out_fig)
    return new_accord


def fig_Energy_value(execution, state):
    out_data = np.array(execution.read('Eigenfunctions'))
    out_data2 = out_data[state]
    out_fig = PLUGIN_PLOT_ISOVOLUMES(out_data2[:, 0], out_data2[:, 1], out_data2[:, 2], out_data2[:, 3], 5)
    #out_fig = PLUGIN_PLOT_ISOSURFACES(out_data2[:,0], out_data2[:,1], out_data2[:,2], out_data2[:,3], 5)
    return out_fig


def fig_Energy_states(execution):
    layout = {
        'xaxis': {
            'title': ''
        },
        'yaxis': {
            'title': 'Energy (eV)',
        },
        'title': 'Energy State',
    }
    out_data = np.array(execution.read('Energy states'))
    out_fig = PLUGIN_PLOT_2D([[0, 1] for o in out_data], [[o, o] for o in out_data], layout)
    return out_fig


def fig_X_polarized(execution):
    layout = {
        'xaxis': {
            'title': 'Energy (eV)',
        },
        'yaxis': {
            'title': 'Light/Dark Transition Strength',
            'type': 'log',
            'exponentformat': 'e',
            'showexponent': 'all'
        },
        'title': 'Light and Dark Transitions (X-polarized))',
    }
    out_data = np.array(execution.read('X-polarized'))
    out_fig = PLUGIN_PLOT_2D([out_data[:, 0]], [out_data[:, 1]], layout)
    return out_fig


def fig_X_polarized(execution):
    layout = {
        'xaxis': {
            'title': 'Energy (eV)',
        },
        'yaxis': {
            'title': 'Light/Dark Transition Strength',
            'type': 'log',
            'exponentformat': 'e',
            'showexponent': 'all'
        },
        'title': 'Light and Dark Transitions (X-polarized))',
    }
    out_data = np.array(execution.read('X-polarized'))
    out_fig = PLUGIN_PLOT_2D([out_data[:, 0]], [out_data[:, 1]], layout)
    return out_fig


def fig_Y_polarized(execution):
    layout = {
        'xaxis': {
            'title': 'Energy (eV)',
        },
        'yaxis': {
            'title': 'Light/Dark Transition Strength',
            'type': 'log',
            'exponentformat': 'e',
            'showexponent': 'all'
        },
        'title': 'Light and Dark Transitions (Y-polarized))',
    }
    out_data = np.array(execution.read('Y-polarized'))
    out_fig = PLUGIN_PLOT_2D([out_data[:, 0]], [out_data[:, 1]], layout)
    return out_fig


def fig_Z_polarized(execution):
    layout = {
        'xaxis': {
            'title': 'Energy (eV)',
        },
        'yaxis': {
            'title': 'Light/Dark Transition Strength',
            'type': 'log',
            'exponentformat': 'e',
            'showexponent': 'all'
        },
        'title': 'Light and Dark Transitions (Z-polarized))',
    }
    out_data = np.array(execution.read('Z-polarized'))
    out_fig = PLUGIN_PLOT_2D([out_data[:, 0]], [out_data[:, 1]], layout)
    return out_fig


def fig_Angle_polarized(execution):
    layout = {
        'xaxis': {
            'title': 'Energy (eV)',
        },
        'yaxis': {
            'title': 'Light/Dark Transition Strength',
            'type': 'log',
            'exponentformat': 'e',
            'showexponent': 'all'
        },
        'title':
        'Light and Dark Transitions (phi=' + str(inputs.phi.value) + ' theta=' + str(inputs.theta.value) + '))',
    }
    out_data = np.array(execution.read('Angle-polarized'))
    out_fig = PLUGIN_PLOT_2D([out_data[:, 0]], [out_data[:, 1]], layout)
    return out_fig


def fig_Absorption(execution):
    layout = {
        'xaxis': {
            'title': 'Energy (eV)',
        },
        'yaxis': {
            'title': 'Absorption',
            'type': 'log',
            'exponentformat': 'e',
            'showexponent': 'all'
        },
        'title': 'Absorption (phi=' + str(inputs.phi.value) + ' theta=' + str(inputs.theta.value) + '))',
    }
    out_data = np.array(execution.read('Absorption'))
    out_fig = PLUGIN_PLOT_2D([out_data[:, 0]], [out_data[:, 1]], layout)
    return out_fig


def fig_Absorption_sweep(execution):
    layout = {
        'xaxis': {
            'title': 'Energy (eV)',
        },
        'yaxis': {
            'title': 'Absorption',
            'type': 'log',
            'exponentformat': 'e',
            'showexponent': 'all'
        },
        'title': 'Absorption Sweep of ' + str(inputs.sweep.value),
    }
    out_data = np.array(execution.read('Absorption Sweep'))
    out_fig = PLUGIN_PLOT_2D([o[:, 0] for o in out_data], [o[:, 1] for o in out_data], layout)
    return out_fig


def fig_Integrated_absorption(execution):
    layout = {
        'xaxis': {
            'title': str(inputs.sweep.value),
        },
        'yaxis': {
            'title': 'Absorption',
            'type': 'log',
            'exponentformat': 'e',
            'showexponent': 'all'
        },
        'title': 'Integrated Absortion',
    }
    out_data = np.array(execution.read('Integrated Absortion'))
    out_fig = PLUGIN_PLOT_2D([out_data[:, 0]], [out_data[:, 1]], layout)
    return out_fig


def fig_Inputdeck(execution):
    out_data = execution.read('input_deck')
    out_fig = Textarea(value=out_data)
    out_fig.layout = Layout(width='100%', height='400px')
    return out_fig


def fig_Log(execution):
    out_data = execution.read('Output Log')
    out_fig = Textarea(value=out_data)
    out_fig.layout = Layout(width='100%', height='400px')
    return out_fig


def showFig(func, param, out, state=None):
    with out:
        clear_output()
        if state == None:
            display(func(param))
        else:
            display(func(param, state))


def on_button_clicked(b):
    run_folder = 'RUNS'
    run_name = str(uuid.uuid4()) + '_run.ipynb'
    with QDOT_UI['sim_progress']:
        clear_output()
        r = Run(nb, inputs, cache=True)

        buttons = []
        output = Output()
        but0 = Button(description='Parameters', layout=Layout(width='auto'))
        but0.on_click(lambda event, param=inputs, func=dict, out=output: showFig(func, param, out))
        buttons.append(but0)

        but1 = Button(description='Energy Values', layout=Layout(width='auto'))
        but1.on_click(lambda event, param=r.db, func=fig_Energy_states, out=output: showFig(func, param, out))
        buttons.append(but1)

        for l in range(inputs.nstates.value):
            but2 = Button(description='Energy Level (' + str(l + 1) + ')', layout=Layout(width='auto'))
            but2.on_click(
                lambda event, state=l, param=r.db, func=fig_Energy_value, out=output: showFig(func, param, out, state))
            buttons.append(but2)

        but3 = Button(description='X_polarized', layout=Layout(width='auto'))
        but3.on_click(lambda event, param=r.db, func=fig_X_polarized, out=output: showFig(func, param, out))
        buttons.append(but3)

        but4 = Button(description='Y_polarized', layout=Layout(width='auto'))
        but4.on_click(lambda event, param=r.db, func=fig_Y_polarized, out=output: showFig(func, param, out))
        buttons.append(but4)

        but5 = Button(description='Z_polarized', layout=Layout(width='auto'))
        but5.on_click(lambda event, param=r.db, func=fig_Z_polarized, out=output: showFig(func, param, out))
        buttons.append(but5)

        but6 = Button(description='Angle_polarized', layout=Layout(width='auto'))
        but6.on_click(lambda event, param=r.db, func=fig_Angle_polarized, out=output: showFig(func, param, out))
        buttons.append(but6)

        but7 = Button(description='Absorption', layout=Layout(width='auto'))
        but7.on_click(lambda event, param=r.db, func=fig_Absorption, out=output: showFig(func, param, out))
        buttons.append(but7)

        but8 = Button(description='Absorption sweep', layout=Layout(width='auto'))
        but8.on_click(lambda event, param=r.db, func=fig_Absorption_sweep, out=output: showFig(func, param, out))
        buttons.append(but8)

        but9 = Button(description='Integrated Absorption', layout=Layout(width='auto'))
        but9.on_click(lambda event, param=r.db, func=fig_Integrated_absorption, out=output: showFig(func, param, out))
        buttons.append(but9)

        but10 = Button(description='input_deck', layout=Layout(width='auto'))
        but10.on_click(lambda event, param=r.db, func=fig_Inputdeck, out=output: showFig(func, param, out))
        buttons.append(but10)

        but11 = Button(description='Output Log', layout=Layout(width='auto'))
        but11.on_click(lambda event, param=r.db, func=fig_Log, out=output: showFig(func, param, out))
        buttons.append(but11)

        new_box = HBox([VBox(buttons, layout=Layout(width='150px')),
                        VBox([output], layout=Layout(width='99%'))],
                       layout=Layout(width='99%'))

        QDOT_UI['l2'].children = [*QDOT_UI['l2'].children, new_box]
        QDOT_UI['l2'].set_title(len(QDOT_UI['l2'].children) - 1, run_name)


QDOT_UI['button'].on_click(on_button_clicked)
