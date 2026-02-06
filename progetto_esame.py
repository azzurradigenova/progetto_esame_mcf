import numpy as np
import matplotlib.pyplot as plt

#costanti fisiche utili
X0 = 7e4          # lunghezza di radiazione (cm)
Ec = 85.0         # energia critica in aria (MeV)
me = 0.511        # massa elettrone (MeV)
dE_ion = 2.0e-3   # perdita per ionizzazione (MeV/cm)

H_START = 20_000  # quota inizio sciame (m)
H_DET = 4_000     # quota rivelatore (m)

#creo la classe particella
class Particle:
    def __init__(self, ptype, energy):
        self.type = ptype    # "e" oppure "gamma"
        self.energy = energy # MeV

#definisco la funzione per la simulazione di uno sciame
def simulate_shower(E0_TeV, step, theta_deg):
	"""
	E0_TeV : energia primaria in TeV
	step   : frazione di X0 (0 < step ≤ 1)
	theta_deg : angolo rispetto alla verticale
    	"""
	E0 = E0_TeV * 1e6      # da TeV a MeV
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
		if p.type == "e":
    			# perdita per ionizzazione (continua)
    			p.energy -= dE_ion * step_len
			
			#se l'energia è troppo bassa la particella muore
			if p.energy <= dE_ion * step_len:
    				continue
			#se l'energia è maggiore di quella critica la particella fa bremsstrahlung
			if p.energy > Ec:
    				if np.random.rand() < (1 - np.exp(-step)):
        			e2 = p.energy / 2
        			new_particles.append(Particle("e", e2))
        			new_particles.append(Particle("gamma", e2))
		elif p.type == "gamma":
    			if p.energy <= 2 * me:
        			continue #se l'energia è sotto soglia niente coppie
			if np.random.rand() < (1 - np.exp(-(7/9)*step)):
    				e = p.energy / 2
    				new_particles.append(Particle("e", e))
    				new_particles.append(Particle("e", e))


	particles = new_particles
	if not particles:
    		break

n_e = sum(1 for p in particles if p.type == "e")
n_g = sum(1 for p in particles if p.type == "gamma")

return len(particles), n_e, n_g

#simulo N sciami indipendenti per discutere statisticamente i risultati
def detector_response(E_TeV, theta_deg, step=0.2, N=300):
hits, electrons, gammas = [], [], []
for _ in range(N):
    h, ne, ng = simulate_shower(E_TeV, step, theta_deg)
    hits.append(h)
    electrons.append(ne)
    gammas.append(ng)

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

plt.hist(hits, bins=30, density=True)
plt.xlabel("Numero di hit")
plt.ylabel("Probabilità")
plt.title("Fluttuazioni statistiche dello sciame")
plt.show()

#mappa energia-angolo
energies = np.logspace(0, 2, 8)  # 1–100 TeV
angles = np.linspace(0, 45, 8)

response = np.zeros((len(energies), len(angles)))

for i, E in enumerate(energies):
    for j, theta in enumerate(angles):
        mean, _, _, _ = detector_response(E, theta)
        response[i, j] = mean

#firma del rivelatore
plt.imshow(response, origin="lower", aspect="auto")
plt.colorbar(label="Hit medi")
plt.xlabel("Angolo")
plt.ylabel("Energia")
plt.title("Risposta del rivelatore")
plt.show()
	

