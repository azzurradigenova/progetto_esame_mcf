# progetto_esame_mcf
Analisi statistica di sciami elettromagnetici ad alta quota.

Questo modulo implementa una simulazione semplificata di uno sciame elettromagnetico prodotto da una particella primaria ad alta energia che attraversa l’atmosfera. Il modello segue l’evoluzione discreta delle particelle secondarie (fotoni ed elettroni) tenendo conto dei principali processi fisici. 

Inizialmente sono state definite delle costanti fisiche utili durante tutto lo sviluppo del programma. 

La classe Particle rappresenta una particella dello sciame ed è caratterizzata da:

-type: il tipo di particella ("gamma" o "e")

-energy: l’energia della particella (in MeV)

La funzione simulate_shower esegue la simulazione dello sciame a partire da un fotone primario.

Parametri:

-E0_TeV: energia iniziale della particella primaria (in TeV)

-step: numero compreso fra 0 e1 per definire la frazione di X0 

-theta_deg: angolo di incidenza dello sciame rispetto alla verticale (in gradi)

Descrizione del metodo:

1)Conversione delle unità:

L’energia primaria viene convertita da TeV a MeV e l’angolo in radianti.

2)Geometria dello sciame:

La lunghezza totale del cammino atmosferico è calcolata tenendo conto dell’inclinazione dello sciame. Il percorso è suddiviso in passi di ampiezza step · X0.

4)Inizializzazione:

Lo sciame è inizialmente composto da un singolo fotone (gamma) con energia E0.

5)Evoluzione dello sciame:

A ogni passo:

-Gli elettroni perdono energia per ionizzazione.

-Se l’energia dell’elettrone è maggiore dell’energia critica Ec, può emettere un fotone tramite bremsstrahlung con una probabilità dipendente dal passo.

-I fotoni con energia sufficiente possono produrre una coppia elettrone–positrone tramite produzione di coppie.

-Le particelle con energia troppo bassa vengono rimosse dalla simulazione.

6)Condizione di arresto:

La simulazione termina quando non restano più particelle nello sciame oppure quando viene raggiunta la fine del percorso atmosferico.

7)Output:

La funzione restituisce il numero totale di particelle presenti nello sciame al termine della simulazione.

La funzione detector_response simula più volte uno sciame con le stesse condizioni iniziali e analizza statisticamente il numero di particelle rivelate (“hit”).

Parametri:

-E_TeV: energia della particella primaria (in TeV)

-theta_deg: angolo di incidenza dello sciame rispetto alla verticale (in gradi)

-step: numero compreso fra zero e uno per definire la frazione di X0 (default: 0.2)

-N: numero di simulazioni indipendenti dello sciame (default: 300)

Descrizione del metodo:

Simulazioni multiple dello sciame.

La funzione richiama ripetutamente simulate_shower per N volte, mantenendo fissi energia e angolo di incidenza.

Ogni simulazione restituisce il numero totale di particelle generate, interpretato come numero di hit nel rivelatore.

Distribuzione statistica degli hit:

I valori ottenuti vengono raccolti in un array NumPy per l’analisi statistica.

Quantità calcolate:

-Valore medio (mean): rappresenta la risposta attesa del rivelatore

-Deviazione standard (std): misura le fluttuazioni statistiche dello sciame

-Fluttuazione relativa (std / mean): stima del rumore intrinseco del rivelatore

-Distribuzione completa degli hit (hits)

Output:

La funzione restituisce:

mean, std, rel_fluct, hits

dove hits è l’insieme dei risultati delle singole simulazioni.

Dopo aver calcolato la risposta del rivelatore per un’energia primaria di 10 TeV e un angolo di 30°, il codice:

-stampa il numero medio di hit e la deviazione standard

-costruisce un istogramma normalizzato della distribuzione degli hit

L’istogramma rappresenta la distribuzione di densità probabilità delle fluttuazioni statistiche dello sciame, evidenziando la natura stocastica del processo di sviluppo e il limite di risoluzione del rivelatore.

Segue una sezione del codice in cui viene analizzata la forma della distribuzione statistica del numero di hit prodotti dallo sciame, verificando se le fluttuazioni possono essere descritte da una distribuzione lognormale.

Funzione di fit: distribuzione lognormale.

Viene definita una funzione di densità di probabilità lognormale:

dove:

-mu è la media del logaritmo della variabile casuale

-sigma è la deviazione standard del logaritmo della variabile

La scelta della lognormale è motivata dal fatto che lo sviluppo dello sciame è un processo moltiplicativo stocastico, per il quale questa distribuzione rappresenta una buona approssimazione.

Preparazione dei dati:

Vengono considerati solo i valori positivi del numero di hit (hits > 0), poiché la distribuzione lognormale è definita per x>0.

Si costruisce un istogramma normalizzato della distribuzione degli hit.

I centri dei bin dell’istogramma vengono utilizzati come punti sperimentali per il fit.

Fit dei dati:

Il fit viene eseguito utilizzando la funzione curve_fit di scipy.optimize, che determina i parametri mu e sigma che meglio descrivono i dati sperimentali: popt, pcov.

I valori iniziali (p0) sono scelti in modo fisicamente ragionevole per favorire la convergenza.

Il risultato del fit fornisce:

-mu_fit: media del logaritmo del numero di hit

-sigma_fit: ampiezza delle fluttuazioni logaritmiche

Visualizzazione dei risultati:

Infine, il codice:

-visualizza l’istogramma delle fluttuazioni statistiche dello sciame

-sovrappone la curva lognormale ottenuta dal fit

Il confronto grafico permette di valutare la bontà della descrizione lognormale delle fluttuazioni del rivelatore.

Interpretazione fisica:

Il buon accordo con una distribuzione lognormale indica che:

-le fluttuazioni del numero di hit non sono puramente poissoniane.

-il processo di sviluppo dello sciame è dominato da meccanismi moltiplicativi casuali.

-la risposta del rivelatore riflette la natura stocastica dell’evoluzione dello sciame elettromagnetico.

Successivamente il codice studia come la risposta media del rivelatore dipenda dall’energia iniziale della particella primaria, verificando una legge di tipo power-law.

Metodo:

Selezione delle energie:

L’energia primaria viene variata tra 1 e 100 TeV utilizzando una scala logaritmica.

Risposta media del rivelatore:

Per ciascun valore di energia, la funzione detector_response viene eseguita più volte (N = 100) per stimare il numero medio di hit prodotti dallo sciame.

Modello di fit:

I dati vengono descritti tramite una legge di potenza: N(E)=A * E ** alpha

dove:

-A è una costante di normalizzazione

-alpha è l’esponente che descrive la crescita dello sciame con l’energia

Fit dei dati:

Il fit viene effettuato con curve_fit, ottenendo una stima di A e dell’esponente alpha, insieme alla sua incertezza statistica derivata dalla matrice di covarianza.

Visualizzazione e interpretazione:

I dati e il fit vengono rappresentati in un grafico log–log, che permette di evidenziare chiaramente l’andamento a legge di potenza.

Il valore dell’esponente alpha quantifica la relazione tra energia primaria e moltiplicazione dello sciame, confermando che il numero medio di hit cresce come una potenza dell’energia iniziale.

Stabilità statistica dello sciame:

Oltre al valore medio, viene analizzata la stabilità statistica della risposta del rivelatore al variare dell’energia.

Metodo:

Per ciascun valore di energia primaria viene calcolata la fluttuazione relativa: sigma/mu, dove mu è il numero medio di hit e sigma la deviazione standard.

I risultati vengono visualizzati in funzione dell’energia su asse logaritmico.

Risultato:

Il grafico mostra che le fluttuazioni relative diminuiscono all’aumentare dell’energia primaria.

Questo comportamento indica che:

-gli sciami ad alta energia sono statisticamente più stabili

-l’effetto delle fluttuazioni casuali diventa meno rilevante per grandi moltiplicazioni dello sciame

-la risoluzione del rivelatore migliora con l’energia

Conclusione fisica:

L’analisi conferma che:

-il numero di hit cresce come una legge di potenza dell’energia primaria

-la risposta del rivelatore diventa più precisa e meno rumorosa all’aumentare dell’energia

-il comportamento osservato è coerente con i modelli teorici degli sciami elettromagnetici

Nell'ultima parte viene eseguito un fit unbinned della distribuzione del numero di hit utilizzando la libreria iminuit, ampiamente impiegata in analisi dati in fisica delle alte energie.

L’obiettivo è confrontare questo approccio con il fit basato su istogramma effettuato in precedenza tramite curve_fit.

Metodo di massima verosimiglianza:

Funzione di costo:

Il fit utilizza una negative log-likelihood non binning costruita direttamente sui dati Monte Carlo: cost = UnbinnedNLL(hits_pos, lognormal).

Questo approccio:

-non dipende dalla scelta dei bin

-sfrutta tutta l’informazione statistica contenuta nei dati

-fornisce stime più robuste dei parametri

Inizializzazione e minimizzazione:

Il minimizzatore Minuit viene inizializzato con valori di partenza fisicamente ragionevoli:

-mu ≈ ln(⟨N⟩)

-sigma ≈ 0.5

Viene inoltre imposto un vincolo fisico su sigma per garantire positività: m.limits["sigma"] = (1e-3, None).

La procedura di minimizzazione segue due passaggi standard:

m.migrad(): ricerca del minimo della funzione di costo

m.hesse(): stima delle incertezze sui parametri tramite matrice di Hess

Confronto Monte Carlo vs fit Minuit:

I parametri ottimali ottenuti vengono utilizzati per costruire la funzione di densità lognormale, che viene sovrapposta all’istogramma dei dati Monte Carlo.

Questo confronto grafico permette di verificare visivamente la bontà del fit e la capacità del modello lognormale di descrivere le fluttuazioni statistiche dello sciame.

Confronto tra i due metodi di fit:

Per una validazione incrociata, vengono confrontati:

-il fit basato su istogramma (curve_fit)

-il fit non binning basato su massima verosimiglianza (Minuit)

Entrambe le curve vengono sovrapposte alla distribuzione degli hit, evidenziando il buon accordo tra i due metodi.

Confronto numerico dei parametri:

Infine, vengono confrontati anche i valori numerici dei parametri e delle rispettive incertezze:

-curve_fit: errori derivati dalla matrice di covarianza del fit sui bin

-Minuit: errori stimati tramite la matrice di Hess della likelihood

Questo confronto mostra che:

-i valori centrali di mu e sigma sono compatibili entro gli errori

-il metodo basato su massima verosimiglianza non binning fornisce una stima più affidabile e meno dipendente dalle scelte di binning

Considerazioni finali:

L’uso di iminuit conferma che:

-la distribuzione del numero di hit è ben descritta da una lognormale

-le fluttuazioni osservate sono coerenti con un processo moltiplicativo stocastico

-l’analisi è robusta rispetto al metodo di fit utilizzato

Questa sezione conclude l’analisi statistica del progetto, fornendo una validazione incrociata dei risultati ottenuti.
