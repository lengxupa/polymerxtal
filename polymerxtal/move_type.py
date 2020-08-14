class Translator:
    def move_sphere(self, sphere):
        print('Translating sphere should be very straightforward')
        sphere.center[0] += 1.0

    def move_gay_berne(self, gay_berne):
        print('Translating solid should involve COM computation and subsequent translation')

    def move_cluster(self, cluster):
        print('Translating cluster should involve COM computation and subsequent translation')
        for particle in cluster.particles:
            particle.move(self)

class Rotator:
    def move_sphere(self, sphere):
        print('Rotating sphere is nonsense due to rotation invariance')

    def move_gay_berne(self, gay_berne):
        print('Rotating Gay Berne using Eulerian angles')
        gay_berne.vector[0] = 1

    def move_cluster(self, cluster):
        print('Rotating cluster using Eulerian angles')
        for particle in cluster.particles:
            particle.move(self)