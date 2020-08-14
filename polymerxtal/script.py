from polymerxtal import Translator, Rotator, Sphere, Gay_Berne, Cluster

def main():

    translator = Translator()
    rotator = Rotator()

    origin = [0.0, 0.0]
    vector = [1.0, 0.0]

    argon = Sphere(origin)
    krypton = Sphere(origin)
    liquid_crystal = Gay_Berne(origin, vector)

    cluster = Cluster()
    cluster.add_particle(argon)
    cluster.add_particle(krypton)
    cluster.add_particle(liquid_crystal)

    argon.move(translator)
    argon.move(rotator)

    krypton.move(translator)
    krypton.move(rotator)

    liquid_crystal.move(translator)
    liquid_crystal.move(rotator)

    cluster.move(translator)
    cluster.move(rotator)

if __name__ == "__main__":
    main()
