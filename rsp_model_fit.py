"""
Fit statisticamente onesto per un modello esteso di R_sp/R200m.

Modello base (consolidato in letteratura, Diemer & Kravtsov 2014 / More+2015):
    R_sp/R200m = A * (1 + B * exp(-Gamma / C))

Modello esteso proposto (da verificare, NON ancora validato in letteratura):
    R_sp/R200m = A * (1 + B * exp(-Gamma / C))
                 + D * log10(1 + N_vicini)
                 + E * t_form_orth
                 + residuo

dove t_form_orth e' il residuo di t_form dopo aver rimosso la componente
gia' spiegata da Gamma (per evitare collinearita' tra i due predittori:
un alone che accresce rapidamente di recente ha quasi per definizione
un tempo di formazione piu' giovane, quindi Gamma e t_form NON sono
indipendenti e vanno ortogonalizzati prima di sommarli in un modello
additivo).

IMPORTANTE: questo script assume che R_sp sia stato MISURATO direttamente
dalla simulazione (es. via SPARTA o profilo di densita'), non calcolato
da Gamma tramite Colossus -- altrimenti il fit torna perfetto per
costruzione (circolarita'), come discusso nel repository "Cosmic Lawnmower".
"""

import numpy as np

try:
    import statsmodels.api as sm
    from statsmodels.stats.outliers_influence import variance_inflation_factor
    HAS_STATSMODELS = True
except ImportError:
    HAS_STATSMODELS = False
    print("statsmodels non disponibile: pip install statsmodels --break-system-packages")

from scipy.optimize import curve_fit
from scipy import stats


# ============================================================
# 1. Modello base consolidato (baseline di riferimento)
# ============================================================

def baseline_model(Gamma, A, B, C):
    """Diemer & Kravtsov 2014 / More, Diemer & Kravtsov 2015, forma
    esponenziale (NON tanh -- vedi discussione)."""
    return A * (1 + B * np.exp(-Gamma / C))


# valori tipici pubblicati (More+2015, eq. 2, z~0):
BASELINE_PARAMS_LITERATURE = dict(A=0.58, B=1.08, C=2.26)


# ============================================================
# 2. Ortogonalizzazione di t_form rispetto a Gamma
#    (rimuove la collinearita' prima di sommare i due termini)
# ============================================================

def orthogonalize(t_form, Gamma):
    """Regressione lineare semplice t_form ~ Gamma, restituisce il
    residuo (la parte di t_form NON spiegata da Gamma)."""
    slope, intercept, r_value, p_value, std_err = stats.linregress(Gamma, t_form)
    t_form_predicted = slope * Gamma + intercept
    t_form_orth = t_form - t_form_predicted
    print(f"Correlazione Gamma vs t_form: r={r_value:.3f}, p={p_value:.2e}")
    if abs(r_value) > 0.3:
        print("  ATTENZIONE: correlazione non trascurabile -- l'ortogonalizzazione "
              "e' necessaria, non opzionale, per questo campione.")
    return t_form_orth, (slope, intercept, r_value)


# ============================================================
# 3. Modello esteso completo
# ============================================================

def extended_model(X, A, B, C, D, E):
    Gamma, log_n_vicini, t_form_orth = X
    return baseline_model(Gamma, A, B, C) + D * log_n_vicini + E * t_form_orth


def fit_extended_model(Gamma, N_vicini, t_form, R_sp_over_R200m):
    """Fit completo con diagnostica. R_sp_over_R200m deve essere MISURATO
    direttamente dalla simulazione (SPARTA o profilo di densita'),
    non calcolato da Gamma -- altrimenti il fit e' circolare."""

    log_n_vicini = np.log10(1 + N_vicini)
    t_form_orth, ortho_info = orthogonalize(t_form, Gamma)

    # --- controllo di collinearita' (VIF) tra i tre predittori finali ---
    if HAS_STATSMODELS:
        X_design = np.column_stack([Gamma, log_n_vicini, t_form_orth])
        X_design_const = sm.add_constant(X_design)
        vifs = [variance_inflation_factor(X_design_const, i)
                for i in range(1, X_design_const.shape[1])]
        print("\nVariance Inflation Factor per predittore (VIF>5 e' un problema, >10 e' grave):")
        for name, vif in zip(["Gamma", "log(1+N_vicini)", "t_form_orth"], vifs):
            flag = "  <-- ATTENZIONE" if vif > 5 else ""
            print(f"  {name:20s}: VIF={vif:.2f}{flag}")

    # --- fit non lineare ---
    p0 = [0.58, 1.08, 2.26, 0.0, 0.0]  # baseline letteratura + estensioni a zero
    X = (Gamma, log_n_vicini, t_form_orth)
    popt, pcov = curve_fit(extended_model, X, R_sp_over_R200m, p0=p0, maxfev=10000)
    perr = np.sqrt(np.diag(pcov))

    pred = extended_model(X, *popt)
    resid = R_sp_over_R200m - pred
    ss_res = np.sum(resid**2)
    ss_tot = np.sum((R_sp_over_R200m - R_sp_over_R200m.mean())**2)
    r2 = 1 - ss_res / ss_tot

    # --- confronto con la sola baseline (per capire se D,E aggiungono davvero informazione) ---
    popt_base, _ = curve_fit(baseline_model, Gamma, R_sp_over_R200m,
                              p0=[0.58, 1.08, 2.26])
    pred_base = baseline_model(Gamma, *popt_base)
    ss_res_base = np.sum((R_sp_over_R200m - pred_base)**2)
    r2_base = 1 - ss_res_base / ss_tot

    n = len(R_sp_over_R200m)
    k_full, k_base = 5, 3
    f_stat = ((ss_res_base - ss_res) / (k_full - k_base)) / (ss_res / (n - k_full))
    p_value_f = 1 - stats.f.cdf(f_stat, k_full - k_base, n - k_full)

    print(f"\n=== Risultati fit ===")
    print(f"Baseline (solo Gamma): R^2={r2_base:.4f}")
    print(f"Modello esteso (Gamma + ambiente + eta'): R^2={r2:.4f}")
    print(f"Test F per l'aggiunta di D,E: F={f_stat:.2f}, p={p_value_f:.4f}")
    if p_value_f < 0.05:
        print("  -> i termini aggiuntivi migliorano il fit in modo statisticamente significativo")
    else:
        print("  -> NON c'e' evidenza che D,E aggiungano potere predittivo reale "
              "(l'estensione non e' giustificata da questi dati)")

    param_names = ["A", "B", "C", "D", "E"]
    print("\nParametri stimati:")
    for name, val, err in zip(param_names, popt, perr):
        print(f"  {name} = {val:.4f} +/- {err:.4f}")

    return dict(popt=popt, perr=perr, r2=r2, r2_baseline=r2_base,
                f_stat=f_stat, p_value=p_value_f, residuals=resid,
                orthogonalization_info=ortho_info,
                predictors=dict(Gamma=Gamma, log_n_vicini=log_n_vicini,
                                 t_form_orth=t_form_orth))


# ============================================================
# 4. Termine AGN (M_BH) -- OPZIONALE, SPENTO DI DEFAULT
# ============================================================
# ATTENZIONE: dopo quanto trovato in letteratura (FLAMINGO, arXiv 2312.05126),
# il raggio di splashback della materia oscura risulta sostanzialmente
# indifferente al modello di feedback barionico (~5-6%, dentro il rumore).
# Questo termine va quindi trattato come SPERIMENTALE: va incluso solo se
# il proprio dataset mostra evidenza statistica reale della sua utilita',
# non aggiunto per abitudine o perche' "potrebbe servire". Il default e'
# delta=0, cioe' spento -- il modello standard (extended_model, sopra)
# resta quello raccomandato finche' non c'e' una giustificazione precisa.

def extended_model_with_bh(X, A, B, C, D, E, delta):
    Gamma, log_n_vicini, t_form_orth, log_m_bh = X
    return extended_model((Gamma, log_n_vicini, t_form_orth), A, B, C, D, E) \
        + delta * log_m_bh


def test_bh_term(Gamma, N_vicini, t_form, M_BH, R_sp_over_R200m, base_result=None):
    """Testa SE aggiungere il termine M_BH e' giustificato, con lo stesso
    approccio di confronto/test-F usato per gli altri termini -- non lo
    include mai per default, lo aggiunge solo se il test lo giustifica.

    base_result: opzionalmente, il risultato gia' calcolato da
    fit_extended_model sugli stessi dati (evita di rifittare due volte).
    """
    log_n_vicini = np.log10(1 + N_vicini)
    t_form_orth, _ = orthogonalize(t_form, Gamma)
    log_m_bh = np.log10(1 + M_BH)

    if base_result is None:
        base_result = fit_extended_model(Gamma, N_vicini, t_form, R_sp_over_R200m)

    ss_res_ext = np.sum(base_result["residuals"] ** 2)

    p0 = list(base_result["popt"]) + [0.0]  # parte da delta=0
    X_bh = (Gamma, log_n_vicini, t_form_orth, log_m_bh)
    popt_bh, pcov_bh = curve_fit(extended_model_with_bh, X_bh, R_sp_over_R200m,
                                  p0=p0, maxfev=10000)
    perr_bh = np.sqrt(np.diag(pcov_bh))

    pred_bh = extended_model_with_bh(X_bh, *popt_bh)
    resid_bh = R_sp_over_R200m - pred_bh
    ss_res_bh = np.sum(resid_bh ** 2)

    n = len(R_sp_over_R200m)
    k_ext, k_bh = 5, 6
    f_stat = ((ss_res_ext - ss_res_bh) / (k_bh - k_ext)) / (ss_res_bh / (n - k_bh))
    p_value_f = 1 - stats.f.cdf(f_stat, k_bh - k_ext, n - k_bh)

    delta_val, delta_err = popt_bh[-1], perr_bh[-1]
    delta_sigma = abs(delta_val) / delta_err if delta_err > 0 else np.nan

    print(f"\n=== Test del termine sperimentale delta*log(1+M_BH) ===")
    print(f"(atteso: FLAMINGO suggerisce un effetto piccolo o nullo -- "
          f"questo NON e' un test a favore del termine, e' un test contro di esso)")
    print(f"delta = {delta_val:.4f} +/- {delta_err:.4f}  ({delta_sigma:.1f} sigma da zero)")
    print(f"Test F per l'aggiunta di delta: F={f_stat:.2f}, p={p_value_f:.4f}")
    if p_value_f < 0.05 and delta_sigma > 2:
        print("  -> il termine M_BH migliora il fit in modo statisticamente significativo.")
        print("     Prima di accettarlo: verificare che non sia un confondimento con Gamma/massa")
        print("     (vedi il problema del terzo test nel README 'Cosmic Lawnmower' -- un segnale")
        print("     fortissimo li' si e' rivelato un artefatto strutturale, non fisica reale).")
    else:
        print("  -> NESSUNA evidenza che il termine M_BH aggiunga potere predittivo reale.")
        print("     Coerente con FLAMINGO. Il modello standard (senza questo termine) resta preferibile.")

    return dict(popt=popt_bh, perr=perr_bh, delta=delta_val, delta_err=delta_err,
                delta_sigma=delta_sigma, f_stat=f_stat, p_value=p_value_f,
                residuals=resid_bh)


if __name__ == "__main__":
    # ESEMPIO SINTETICO -- sostituire con dati reali (R_sp MISURATO, non calcolato)
    rng = np.random.default_rng(0)
    n = 200
    Gamma = rng.uniform(0.5, 6, n)
    # t_form correlato con Gamma di proposito, per mostrare l'effetto
    # dell'ortogonalizzazione su dati realistici
    t_form = 0.6 - 0.05 * Gamma + rng.normal(0, 0.08, n)
    N_vicini = rng.poisson(3, n)
    # M_BH correlato con la massa dell'alone (quindi indirettamente con Gamma) --
    # per costruzione qui NON c'e' nessun effetto diretto di M_BH su R_sp,
    # solo la correlazione spuria che ci si aspetta da un confondimento di massa
    M_BH = 1e6 * 10 ** (rng.normal(2, 0.5, n))
    R_sp_true = baseline_model(Gamma, 0.58, 1.08, 2.26) + rng.normal(0, 0.02, n)

    print("=== 1. MODELLO STANDARD (raccomandato di default) ===\n")
    result = fit_extended_model(Gamma, N_vicini, t_form, R_sp_true)

    print("\n\n=== 2. Termine AGN sperimentale (solo su richiesta esplicita, non di default) ===")
    result_bh = test_bh_term(Gamma, N_vicini, t_form, M_BH, R_sp_true, base_result=result)
