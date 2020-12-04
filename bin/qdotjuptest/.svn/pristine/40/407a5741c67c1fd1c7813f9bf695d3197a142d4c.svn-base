from plotly.graph_objs import FigureWidget
import plotly.figure_factory as FF
import plotly.graph_objs as go
from skimage import measure
import matplotlib
from hublib.tool import read
from ipywidgets import Layout, Tab, Textarea
import numpy as np

def PLUGIN_PLOT_ISOSURFACES(X_COLUMN, Y_COLUMN, Z_COLUMN, V_COLUMN, NUM_CONTOURS):
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
    colormap = []
    #plotly hack to include custom colorbars

    trace = go.Scatter3d(
        x=[0] * (num_iso + 2),
        y=[0] * (num_iso + 2),
        z=[0] * (num_iso + 2),
        mode='markers',
        marker=dict(
            size=0.1,
            color=np.linspace(min_val, max_val, num=num_iso+2),
            colorscale='Viridis',
            opacity=1,
            colorbar={"title": "$\\Psi^2$"},
        ),
    )

    for value in np.linspace(min_val, max_val, num=num_iso+2):
        color_rgb = matplotlib.colors.colorConverter.to_rgb(viridis_cmap(norm(value)))
        colormap.append(color_rgb)
        if value > min_val and value < max_val:
            vertices, simplices, faces, values = measure.marching_cubes_lewiner(vol, value, (float(X_DISTANCE), float(Y_DISTANCE), float(Z_DISTANCE)))
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
            data.append(fig['data'][0])
        data.append(trace)
    fig['layout']['showlegend'] = False
    fig['layout']['scene'] = {'aspectmode': 'data'}  # preserve
    return FigureWidget(data=data, layout=fig['layout'])


def PLUGIN_PLOT_2D(xcol, ycol, layout={}):
    data = []
    for indx in range(len(ycol)):
        trace = go.Scatter(
            x=xcol[indx],
            y=ycol[indx],
        )
        data.append(trace)
    return FigureWidget(data=data, layout=layout)

def fig_Energy_values(execution, states):
    out_data = np.array(read(execution, 'Eigenfunctions'))
    children = []
    for i in range(states):
        out_data2 = out_data[i]
        out_fig = PLUGIN_PLOT_ISOSURFACES(
            out_data2[:, 0], out_data2[:, 1], out_data2[:, 2], out_data2[:, 3], 5)
        children.append(out_fig)
    new_accord = Tab(children, layout=Layout(width='100%'))
    out_data_e = np.array(read(execution, 'Energy states'))
    for i, o in enumerate(out_data_e):
        new_accord.set_title(i, "l" + str(o))
    return new_accord

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
    out_data = np.array(read(execution, 'Energy states'))
    out_fig = PLUGIN_PLOT_2D([[0, 1] for o in out_data], [
                             [o, o] for o in out_data], layout)
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
    out_data = np.array(read(execution, 'X-polarized'))
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
    out_data = np.array(read(execution, 'Y-polarized'))
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
    out_data = np.array(read(execution, 'Z-polarized'))
    out_fig = PLUGIN_PLOT_2D([out_data[:, 0]], [out_data[:, 1]], layout)
    return out_fig


def fig_Angle_polarized(execution, phi, theta):
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
        'title': 'Light and Dark Transitions (phi='+str(phi)+' theta='+str(theta)+'))',
    }
    out_data = np.array(read(execution, 'Angle-polarized'))
    out_fig = PLUGIN_PLOT_2D([out_data[:, 0]], [out_data[:, 1]], layout)
    return out_fig


def fig_Absorption(execution, phi, theta):
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
        'title': 'Absorption (phi='+str(phi)+' theta='+str(theta)+'))',
    }
    out_data = np.array(read(execution, 'Absorption'))
    out_fig = PLUGIN_PLOT_2D([out_data[:, 0]], [out_data[:, 1]], layout)
    return out_fig


def fig_Absorption_sweep(execution, sweep):
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
        'title': 'Absorption Sweep of ' + str(sweep),
    }
    out_data = np.array(read(execution, 'Absorption Sweep'))
    out_fig = PLUGIN_PLOT_2D([o[:, 0] for o in out_data], [
                             o[:, 1] for o in out_data], layout)
    return out_fig


def fig_Integrated_absorption(execution, sweep):
    layout = {
        'xaxis': {
            'title': str(sweep),
        },
        'yaxis': {
            'title': 'Absorption',
            'type': 'log',
            'exponentformat': 'e',
            'showexponent': 'all'
        },
        'title': 'Integrated Absortion',
    }
    out_data = np.array(read(execution, 'Integrated Absortion'))
    out_fig = PLUGIN_PLOT_2D([out_data[:, 0]], [out_data[:, 1]], layout)
    return out_fig


def fig_Inputdeck(execution):
    out_data = read(execution, 'input_deck')
    out_fig = Textarea(value=out_data)
    out_fig.layout = Layout(width='100%', height='400px')
    return out_fig


def fig_Log(execution):
    out_data = read(execution, 'Output Log')
    out_fig = Textarea(value=out_data)
    out_fig.layout = Layout(width='100%', height='400px')
    return out_fig

