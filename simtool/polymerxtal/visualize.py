"""
Functions for visualization of molecules
"""

import numpy as np
import sys

from mpl_toolkits.mplot3d import Axes3D

try:
    from ovito.io import import_file
    from ovito.vis import Viewport

    use_ovito = True
except:
    use_ovito = False

try:
    import matplotlib.pyplot as plt
except:
    from polymerxtal.io import check_nanohub

    use_nanohub = check_nanohub()
    if not (use_ovito and use_nanohub):
        import matplotlib.pyplot as plt

from polymerxtal.data import atom_colors


def draw_molecule(coordinates, symbols, draw_bonds=None, save_location=None, dpi=300):

    # Draw a picture of a molecule using matplotlib.

    # Create figure
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")

    # Get colors - based on atom name
    colors = []
    for atom in symbols:
        colors.append(atom_colors[atom])

    size = np.array(plt.rcParams["lines.markersize"] ** 2) * 200 / (len(coordinates))

    ax.scatter(
        coordinates[:, 0],
        coordinates[:, 1],
        coordinates[:, 2],
        marker="o",
        edgecolors="k",
        facecolors=colors,
        alpha=1,
        s=size,
    )

    # Draw bonds
    if draw_bonds:
        for atoms, bond_length in draw_bonds.items():
            atom1 = atoms[0]
            atom2 = atoms[1]

            ax.plot(
                coordinates[[atom1, atom2], 0],
                coordinates[[atom1, atom2], 1],
                coordinates[[atom1, atom2], 2],
                color="k",
            )

    # Save figure
    if save_location:
        plt.savefig(save_location, dpi=dpi, graph_min=0, graph_max=2)

    return ax


def bond_histogram(bond_list, save_location=None, dpi=300, graph_min=0, graph_max=2):
    # Draw a histogram of bond lengths based on a bond_list (output from build_bond_list function)

    lengths = []
    for atoms, bond_length in bond_list.items():
        lengths.append(bond_length)

    bins = np.linspace(graph_min, graph_max)

    fig = plt.figure()
    ax = fig.add_subplot(111)

    plt.xlabel("Bond Length (angstrom)")
    plt.ylabel("Number of Bonds")

    ax.hist(lengths, bins=bins)

    # Save figure
    if save_location:
        plt.savefig(save_location, dpi=dpi)

    return ax


def ovito_view(sample_path, filename, view="Perspective"):
    """
    Use the package ovito to make visualizaitons of molecules.

    Parameters
    ----------
    sample_path : str
        The path of the file to visualize
    filename : str
        The name of the output file image
    view : str (optional)
        The view to use
    """

    if use_ovito:
        # Import the sample file.
        pipeline = import_file(sample_path)
        pipeline.source.data.cell.vis.enabled = False
        pipeline.source.data.particles.vis.radius = 0.5
        pipeline.add_to_scene()

        vp = Viewport()

        if view == "Perspective":
            vp.type = Viewport.Type.Perspective
            vp.camera_dir = (-1, -1, -1)
        elif view == "Ortho":
            vp.type = Viewport.Type.Ortho
            vp.camera_dir = (-1, -1, -1)
        elif view == "Top":
            vp.type = Viewport.Type.Top
        elif view == "Bottom":
            vp.type = Viewport.Type.Bottom
        elif view == "Front":
            vp.type = Viewport.Type.Front
        elif view == "Back":
            vp.type = Viewport.Type.Back
        elif view == "Left":
            vp.type = Viewport.Type.Left
        elif view == "Right":
            vp.type = Viewport.Type.Right

        vp.zoom_all()
        vp.render_image(
            size=(800, 600), filename=filename, background=(0, 0, 0), frame=0
        )

        pipeline.remove_from_scene()

    else:
        print(
            "Cannot use function ovito_view - the package ovito is not installed or cannot be found."
        )


def main(sample_path, filename, view):
    ovito_view(sample_path, filename, view=view)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3])
