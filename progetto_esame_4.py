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
    step   : frazione di X0 (0 < step <= 1)
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

    return len(particles)


#risposta statistica del rivelatore
def detector_response(E_TeV, theta_deg, step=0.2, N=300):

    hits = []

    for _ in range(N):
        h = simulate_shower(E_TeV, step, theta_deg)
        hits.append(h)

    hits = np.array(hits)

    mean = np.mean(hits) #media, risposta attesa
    std = np.std(hits) #fluttuazioni dello sciame
    rel_fluct = std / mean if mean > 0 else 0 #rumore intrinseco del rivelatore

    return mean, std, rel_fluct, hits


#istogramma delle fluttuazioni
mean, std, rel, hits = detector_response(10, 30)

print(f"Il numero medio di hit è: {mean:.2f}")
print(f"La deviazione standard del numero di hit è: {std:.2f}")

plt.figure()
plt.hist(hits, bins=30, density=True, color='pink')
plt.xlabel("Numero di hit")
plt.ylabel("Densità di probabilità")
plt.title("Fluttuazioni statistiche dello sciame")
plt.show()

from scipy.optimize import curve_fit

#definisco una funzione opportuna per il fit
def lognormal(x, mu, sigma):
    """
    mu=media del logaritmo della variabile
    sigma=deviazione standard del logaritmo della variabile
    """
    return 1.0 / (x * sigma * np.sqrt(2*np.pi)) * np.exp(
        - (np.log(x) - mu)**2 / (2 * sigma**2)
    )
#cerco i parametri che meglio descrivono i dati
hits_pos = hits[hits > 0]
hist, bin_edges = np.histogram(hits_pos, bins=30, density=True) 
bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])

popt, pcov = curve_fit(
    lognormal,
    bin_centers,
    hist,
    p0=[np.log(np.mean(hits)), 0.5]
)
mu_fit, sigma_fit=popt

x = np.linspace(min(hits_pos), max(hits_pos), 200)

plt.figure()
plt.hist(hits_pos, bins=30, density=True, color='pink', alpha=0.6, label="Fluttuazioni statistiche dello sciame")
plt.plot(x, lognormal(x, mu_fit, sigma_fit), '-', color='purple', label="Fit lognormale")
plt.xlabel("Numero di hit")
plt.ylabel("Densità di probabilità")
plt.legend()
plt.show()

#dimostro che il numero di hit è proporzionale a una potenza dell'energia iniziale
energies = np.logspace(0, 2, 6)  # 1–100 TeV
means = []

for E in energies:
    mean, _, _, _ = detector_response(E, 30, N=100)
    means.append(mean)

energies = np.array(energies)
means = np.array(means)

def power_law(E, A, alpha):
    """
    E=energia della particella madre
    A=costante moltiplicativa
    alpha=potenza dell'energia E
    """
    return A * E**alpha

popt2, pcov2 = curve_fit(power_law, energies, means)

A_fit, alpha_fit = popt2
sigma_alpha = np.sqrt(pcov2[1,1])

E_plot = np.logspace(0, 2, 200)

plt.figure()
plt.loglog(energies, means, 'o')
plt.loglog(E_plot, power_law(E_plot, A_fit, alpha_fit), label=rf"Fit: $N = A E^{{\alpha}}$, $\alpha={alpha_fit:.2f}\pm{sigma_alpha:.2f}$")

plt.xlabel("Energia primaria [TeV]")
plt.ylabel("Hit medi")
plt.legend()
plt.show()


#dimostro che le fluttuazioni relative diminuiscono con l'energia iniziale
rel_flucts = []

for E in energies:
    _, _, rel, _ = detector_response(E, 30, N=100)
    rel_flucts.append(rel)

plt.figure()
plt.semilogx(energies, rel_flucts, 'o-')
plt.xlabel("Energia primaria [TeV]")
plt.ylabel(r"Fluttuazioni relative $\frac{\sigma}{\mu}$")
plt.title("Stabilità statistica dello sciame")
plt.show()

#fit lognormale degli hit con iminuit
from iminuit import Minuit
from iminuit.cost import UnbinnedNLL

cost = UnbinnedNLL(hits_pos, lognormal)

#inizializzazione e minimizzazione
m = Minuit(
    cost,
    mu=np.log(np.mean(hits_pos)),
    sigma=0.5
)

m.limits["sigma"] = (1e-3, None)  # vincolo fisico
m.migrad()   # minimizzazione
m.hesse()    # errori sui parametri


#confronto grafico Monte Carlo vs Fit
x = np.linspace(min(hits_pos), max(hits_pos), 300)
pdf_fit = lognormal(x, m.values["mu"], m.values["sigma"])

plt.figure()
plt.hist(hits_pos, bins=30, density=True, alpha=0.6, label="Fluttuazioni statistiche dello sciame")
plt.plot(x, pdf_fit, 'r-', lw=2, label="Fit lognormale (MINUIT)")
plt.xlabel("Numero di hit")
plt.ylabel("Densità di probabilità")
plt.title("Distribuzione degli hit – Fit con iminuit")
plt.legend()
plt.show()

# confronto grafico dei due fit
x = np.linspace(min(hits_pos), max(hits_pos), 300)
pdf_curve_fit = lognormal(x, mu_fit, sigma_fit)   # fit "a mano" con curve_fit
pdf_minuit    = lognormal(x, m.values["mu"], m.values["sigma"])  # fit Minuit

plt.figure()
plt.hist(hits_pos, bins=30, density=True, alpha=0.6, label="Fluttuazioni statistiche dello sciame", color='pink')
plt.plot(x, pdf_curve_fit, 'b-', lw=2, label="Fit lognormale (curve_fit)")
plt.plot(x, pdf_minuit, 'r--', lw=2, label="Fit lognormale (Minuit)")
plt.xlabel("Numero di hit")
plt.ylabel("Densità di probabilità")
plt.title("Distribuzione degli hit – Confronto dei fit")
plt.legend()
plt.show()

# errore sui parametri calcolato da curve_fit
sigma_mu_fit = np.sqrt(pcov[0,0])
sigma_sigma_fit = np.sqrt(pcov[1,1])

# confronto numerico dei parametri
print("Confronto parametri dei due fit:")
print(f"Curve_fit: mu = {mu_fit:.3f} ± {sigma_mu_fit:.3f}, "
      f"sigma = {sigma_fit:.3f} ± {sigma_sigma_fit:.3f}")
print(f"Minuit   : mu = {m.values['mu']:.3f} ± {m.errors['mu']:.3f}, "
      f"sigma = {m.values['sigma']:.3f} ± {m.errors['sigma']:.3f}")    