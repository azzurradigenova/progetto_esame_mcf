import numpy as np
import matplotlib.pyplot as plt

#costanti fisiche utili
X0 = 7e4 #lunghezza di radiazione (cm)
Ec = 85.0 #energia critica in aria (MeV)
me = 0.511 #massa elettrone (MeV)
dE_ion = 2.0e-3 #perdita per ionizzazione (MeV/cm)

H_START = 20_000 #quota di inizio dello sciame (m)
H_DET = 4_000 #quota del rivelatore (m)


#creo la classe particella
class Particle:
    def __init__(self, ptype, energy):
        self.type = ptype
        self.energy = energy


#definisco la funzione per la simulazione di uno sciame
def simulate_shower(E0_TeV, step, theta_deg):
    """
    E0_TeV : energia primaria in TeV
    step   : frazione di X0 (0 < step ≤ 1)
    theta_deg : angolo rispetto alla verticale
    """
    E0 = E0_TeV * 1e6 #da TeV a MeV
    theta = np.radians(theta_deg)

    #tengo conto della geometria dello sciame
    path_total = (H_START - H_DET) * 100 / np.cos(theta)
    step_len = step * X0
    n_steps = int(path_total / step_len)

    #inizializzazione dello sciame
    particles = [Particle("gamma", E0)]

    #evoluzione dello sciame
    for _ in range(n_steps):

        new_particles = []

        for p in particles:

            if p.type == "e":

                p.energy -= dE_ion * step_len

                if p.energy <= dE_ion * step_len:
                    continue

                if p.energy > Ec:
                    if np.random.rand() < (1 - np.exp(-step)):
                        e2 = p.energy / 2
                        new_particles.append(Particle("e", e2))
                        new_particles.append(Particle("gamma", e2))
                    else:
                        new_particles.append(p)
                else:
                    new_particles.append(p)

            elif p.type == "gamma":

                if p.energy <= 2 * me:
                    continue

                if np.random.rand() < (1 - np.exp(-(7/9) * step)):
                    e = p.energy / 2
                    new_particles.append(Particle("e", e))
                    new_particles.append(Particle("e", e))
                else:
                    new_particles.append(p)

        particles = new_particles

        if not particles:
            break

    n_e = sum(1 for p in particles if p.type == "e")
    n_g = sum(1 for p in particles if p.type == "gamma")

    return len(particles), n_e, n_g


#risposta statistica del rivelatore
def detector_response(E_TeV, theta_deg, step=0.2, N=300):

    hits = []

    for _ in range(N):
        h, _, _ = simulate_shower(E_TeV, step, theta_deg)
        hits.append(h)

    hits = np.array(hits)

    mean = np.mean(hits) #media, risposta attesa
    std = np.std(hits) #fluttuazioni dello sciame
    rel_fluct = std / mean if mean > 0 else 0 #rumore intrinseco del rivelatore

    return mean, std, rel_fluct, hits


#trovo la profondità teorica massima dello sciame
def rossi_tmax(E0_TeV):
    E0 = E0_TeV * 1e6
    return np.log(E0 / Ec) / np.log(2)


#istogramma delle fluttuazioni
mean, std, rel, hits = detector_response(10, 30)

plt.figure()
plt.hist(hits, bins=30, density=True)
plt.xlabel("Numero di hit")
plt.ylabel("Probabilità")
plt.title("Fluttuazioni statistiche dello sciame")
plt.show()

#mappa energia-angolo
#energies = np.logspace(0, 2, 8) #1-100 TeV
#angles = np.linspace(0, 45, 8) #0-45 gradi

#response = np.zeros((len(energies), len(angles)))

#for i, E in enumerate(energies):
    #for j, theta in enumerate(angles):
        #mean, _, _, _ = detector_response(E, theta)
        #response[i, j] = mean

#plt.figure()
#plt.imshow(response, origin="lower", aspect="auto")
#plt.colorbar(label="Hit medi")
#plt.xlabel("Indice angolo")
#plt.ylabel("Indice energia")
#plt.title("Risposta del rivelatore")
#plt.show()

from scipy.stats import lognorm

# fit lognormale della distribuzione degli hit
shape, loc, scale = lognorm.fit(hits, floc=0)

x = np.linspace(min(hits), max(hits), 200)
pdf = lognorm.pdf(x, shape, loc, scale)

plt.figure()
plt.hist(hits, bins=30, density=True, alpha=0.6, label="Monte Carlo")
plt.plot(x, pdf, 'r-', label="Fit lognormale")
plt.xlabel("Numero di hit")
plt.ylabel("Probabilità")
plt.title("Distribuzione degli hit e fit lognormale")
plt.legend()
plt.show()

energies = np.logspace(0, 2, 6)  # 1–100 TeV
means = []

for E in energies:
    mean, _, _, _ = detector_response(E, 30, N=100)
    means.append(mean)

energies = np.array(energies)
means = np.array(means)

# fit log-log
coeff = np.polyfit(np.log(energies), np.log(means), 1)
alpha = coeff[0]

plt.figure()
plt.loglog(energies, means, 'o', label="Monte Carlo")
plt.loglog(energies, np.exp(coeff[1]) * energies**alpha,
           label=f"Fit: N ∝ E^{alpha:.2f}")
plt.xlabel("Energia primaria [TeV]")
plt.ylabel("Hit medi")
plt.title("Scaling del numero di hit con l’energia")
plt.legend()
plt.show()

rel_flucts = []

for E in energies:
    _, _, rel, _ = detector_response(E, 30, N=100)
    rel_flucts.append(rel)

plt.figure()
plt.semilogx(energies, rel_flucts, 'o-')
plt.xlabel("Energia primaria [TeV]")
plt.ylabel("Fluttuazioni relative σ/μ")
plt.title("Stabilità statistica dello sciame")
plt.show()

