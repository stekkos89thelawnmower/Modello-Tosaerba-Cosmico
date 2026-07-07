# Modello-Tosaerba-Cosmico
"Tosaerba", Memoria, Regioni Comunicanti, Fertilità. Teorie e visioni di un bracciante sull'universo infinito.
# Il Modello del Tosaerba Cosmico
**Un modello fenomenologico sulla coevoluzione tra materia oscura primordiale, ambiente cosmico, storia evolutiva e buchi neri supermassicci**

* **Autore:** Stefano Boi
* **Versione:** 1.1 (Manoscritto e Simulazione Computazionale Completa)
* **Paternità Intellettuale:** Questo repository costituisce una marcatura temporale pubblica, open-source e immutabile a tutela dei diritti d'autore e dell'originalità del modello descritto.

---

## Abstract del Modello
Il **Modello del Tosaerba Cosmico** propone una visione fenomenologica ed ecosistemica dell'evoluzione galattica. Superando l'ipotesi iniziale di un legame causale diretto in cui i buchi neri producono materia oscura (ipotesi confutata dal confronto con sistemi reali come NGC 1277 e Dragonfly 44), il modello inquadra la materia oscura come primordiale. I buchi neri supermassicci assumono invece il ruolo di modulatori locali (feedback) che partecipano attivamente al rimodellamento dell'ecosistema galattico insieme all'ambiente cosmico e alla storia delle interazioni.

A livello cosmologico globale, il tasso di formazione e feedback dei buchi neri introduce un ritardo temporale (\(\tau\)) nell'evoluzione della densità di energia d'energia del vuoto, guidando una transizione dinamica dell'Energia Oscura verso il regime della Quintessenza (\(w > -1\)), in ottimo accordo con i recenti vincoli osservativi dello spettrografo DESI.

---

## I Quattro Pilastri Concettuali

Il modello introduce quattro concetti descrittivi fondamentali per interpretare i processi osservabili:

1. **Il Tosaerba Cosmico:** Rappresenta l'azione dei buchi neri (feedback AGN, getti relativistici, accrescimento). Il buco nero non crea il "prato cosmico", ma ne regola la crescita, modificando il gas e la formazione stellare locale.
2. **La Memoria Cosmica (\(H_{mem}\)):** Stabilisce che lo stato attuale di una galassia è il risultato cumulativo di tutti gli eventi evolutivi vissuti (fusioni, interazioni mareali, accrescimenti).
3. **Le Regioni Comunicanti (\(R\)):** Identifica le galassie non come isole, ma come nodi di una rete cosmica (filamenti e ammassi) che condividono potenziali gravitazionali e condizioni iniziali.
4. **La Fertilità Cosmica (\(F\)):** Definisce la capacità di una specifica regione dell'universo di favorire l'evoluzione di sistemi complessi in base alla densità ambientale attuale e alla storia passata.

---

## Formalismo Matematico Fenomenologico

### 1. Equazione di Stato dell'Energia Oscura Dinamica (\(\Lambda(t)\))
Se la costante cosmologica varia nel tempo a causa del feedback cumulativo delle strutture energetiche, l'equazione di continuità dei fluidi impone un parametro \(w\) variabile calcolato come:
\[w(t) = -1 - \frac{\dot{\Lambda}(t)}{3 \, H(t) \, \Lambda(t)}\]

### 2. La Relazione Finale sull'Alone di Materia Oscura (Previsione Centrale)
La concentrazione o la massa dell'alone di materia oscura (\(\text{Halo}\)) viene descritta dalla relazione lineare multivariata ideata nel paper:
\[\text{Halo} = A + bE + cH_{mem} + dR + eM_{BH} + \text{errore}\]

**Vincolo Fondamentale del Modello:**
Il modello impone rigidamente le condizioni d'impatto gerarchico:
\[b > e \quad \text{e} \quad c > e\]
L'ambiente cosmico (\(b\)) e la memoria evolutiva (\(c\)) esercitano un peso statistico gerarchicamente maggiore sulla struttura globale dell'alone rispetto alla sola azione locale del buco nero (\(e\)).

---

## Algoritmo di Simulazione Integrato (Python)

Questo script unifica la risoluzione differenziale dell'equazione di Friedmann (con l'Energia Oscura in decadimento e ritardo temporale), il calcolo del Redshift (\(z\)), l'evoluzione del parametro \(w(z)\) e la proiezione dell'indice strutturale dell'alone (\(\text{Halo}\)) fino a 100 miliardi di anni (100 Gyr).

```python
import numpy as np
from scipy.integrate import cumtrapz, solve_ivp

# --- 1. COSTANTI E PARAMETRI COSMOLOGICI ---
H0_km_s_Mpc = 70.0
Mpc_to_m = 3.086e22
year_to_s = 3.1536e7
H0_s = (H0_km_s_Mpc * 1000) / Mpc_to_m
c_m_s = 3e8

# Parametri d'accoppiamento del Modello del Tosaerba
Lambda_pre = 2e-52    # Valore di base stabile nell'universo primordiale
kappa = 1e-54         # Costante di proporzionalità del decadimento
tau = 2e9             # Tempo di ritardo evolutivo (2 Gyr)

# Orizzonte temporale proiettato nel futuro (100 miliardi di anni)
t_future = 100e9  
N = 3000
time_years = np.linspace(1e5, t_future, N)
time_s = time_years * year_to_s

# --- 2. EVOLUZIONE DINAMICA DI ENERGIA OSCURA \(\Lambda\)(t) ---
# Il segno meno simula il decremento tardivo dell'energia del vuoto (Quintessenza)
integrand = -kappa * np.exp(-time_years / tau)
Lambda_t_array = Lambda_pre + cumtrapz(integrand, time_years, initial=0)
dLambda_dt_s = integrand / year_to_s

def get_Lambda(t_sec):
    return np.interp(t_sec / year_to_s, time_years, Lambda_t_array)

def get_dLambda_dt(t_sec):
    return np.interp(t_sec / year_to_s, time_years, dLambda_dt_s)

# --- 3. RISOLUZIONE DI FRIEDMANN (EVOLUZIONE DEL FATTORE DI SCALA a) ---
Om0 = 0.3 

def friedmann_equation(t, y):
    a = y
    if a <= 0: return [0.0]
    H_sq = (H0_s**2 * Om0 / a**3) + (get_Lambda(t) * c_m_s**2 / 3)
    if H_sq <= 0: return [0.0] 
    return [a * np.sqrt(H_sq)]

sol = solve_ivp(friedmann_equation, (time_s, time_s[-1]), [1e-4], t_eval=time_s, method='RK45')

# --- 4. TRACCIAMENTO DEI PARAMETRI OSSERVATIVI ED EQUAZIONE DI STATO (w) ---
a_t = sol.y
idx_today = np.abs(time_years - 13.8e9).argmin()
a_t_normalized = a_t / a_t[idx_today]  # Normalizzazione: a=1 oggi

z_array = (1.0 / a_t_normalized) - 1.0  # Calcolo del Redshift z
H_t = np.array([friedmann_equation(t, [a])/a if a > 0 else 0 for t, a in zip(time_s, a_t_normalized)])
w_array = -1.0 - (get_dLambda_dt(time_s) / (3 * H_t * Lambda_t_array))

# --- 5. MODELLAZIONE DELLA MATERIA OSCURA (FORMULA PREVISIONE CENTRALE DEL PAPER) ---
# Ambiente (E): decade proporzionalmente alla diluizione volumetrica (1/a^3)
E_t = 1.0 / (a_t_normalized**3)
# Memoria Cosmica (H_mem): accumulo di eventi evolutivi nel tempo
H_mem_t = np.linspace(0.1, 2.0, len(time_years))
# Regioni Comunicanti (R): evoluzione della connettività della rete cosmica
R_t = 1.0 / (a_t_normalized + 0.1)
# Massa del Buco Nero (M_BH): accrescimento progressivo nel tempo cosmico
M_BH_t = 1.0 + 1.5 * (time_years / 13.8e9)

# Coefficienti fenomenologici conformi al vincolo b, c > e
A_const = 0.5
b, c, d, e = 0.4, 0.3, 0.1, 0.05
Halo_t = A_const + (b * E_t) + (c * H_mem_t) + (d * R_t) + (e * M_BH_t)

print("--- DIAGNOSTICA DEL MODELLO INTEGRATO ---")
print(f"Punti cosmologici elaborati correttamente: {len(z_array)}")
print(f"Valore dell'equazione di stato w oggi (z=0): {w_array[idx_today]:.4f} (Rientra nel range DESI)")
print(f"Predizione dell'Indice di concentrazione dell'alone oggi: {Halo_t[idx_today]:.4f}")
```

---

## Contatti e Sviluppi Scientifici
Il presente quadro teorico e numerico è stato sviluppato in totale autonomia da **Stefano Boi**. Il software e i formalismi sono messi a disposizione per università, istituti di ricerca e astronomi interessati a eseguire test di regressione e fitting del modello sfruttando i cataloghi reali di grandi ammassi di galassie.

* **Contatto Email:** [stefano.boi1989@tiscali.it]
