# =============== READ-ONLY: LIFT + auto-sök utan baspercentiler (mark_klass × torv_klass) ===============
#  1) Delar in i L/M/H via terciler (Elev_mean).
#  2) Stratifierar inom varje band på mark_klass x torv_klass.
#  3) Bygger SÖKRYMDEN AUTOMATISKT per stratum direkt från datan (inga bas-interval).
#     - För H, C (högre = bättre) börjar vi söket vid den percentil som motsvarar ett "golv" i ENHETER (H_FLOOR, C_FLOOR),
#       så att trösklarna inte kan bli orimligt låga, och sträcker oss upp mot p≈99.
#     - För G (och ev. F om aktiv; lägre = bättre) slutar vi söket vid den percentil som motsvarar ett "tak" i ENHETER (G_CEILING),
#       så att trösklarna inte kan bli orimligt generösa.
#  4) Kör en grid-sök över alla kombinationer av H/C/G/(F)-percentiler (konverterade till ABS-trösklar i dm/m/dm/år).
#  5) Printar trösklar (ABS & PCTL) under gång. INGEN skrivning till is_old.
# ========================================================================================================

import arcpy, math, numpy as np
arcpy.env.overwriteOutput = True

# ------------------ INDATA ------------------
# Feature class / lager (namn i Contents eller full sökväg)
FC = r"C:\Your\File_Path\To\Tradkronor_PreProcessed"

# Fält MÅSTE finnas i FC:
ELEV_FIELD = "Elev_mean" # m (Bandindelning via terciler)
H_FIELD    = "O2_PCT95"   # dm  (Högre=Bättre)
G_FIELD    = "gr_p95yr"   # dm/år (Lägre=Bättre, >=0)
C_FIELD    = "cr_dia"     # m   (Högre=Bättre)
F_FIELD    = "cr_flat"    # dm  (Lägre=Bättre, valfri)
USE_F = True              # F har stökat med thresholds tidigare, sätt F False om struliga resultat.

# Kategorier för stratifiering
MARK_FIELD = "mark_klass"   # t.ex. "torr-frisk", "frisk-fuktig"...
TORV_FIELD = "torv_klass"   # t.ex. "mineraljord", ">= 40cm torvdjup"...

# ------------------ STRIKTHETS-SWITCH ------------------
# "normal" | "strikt" | "superstrikt"
STRICTNESS = "normal"

# Profiler: hur tätt vi provar percentiler samt hur hårda räckena är.
# nH/nC/nG/nF = antal jämnt fördelade percentilpnkter per dimension (fler = långsammare, mer noggrant)
# Lift_Min = minsta tillåtna enrichment (observed / expected under oboeroendeantagandet)
# Delta_z_min = "signalstyrka" mellan valda objekt och referensen i standardavvikelser. Lågt värde -> tillåter svaga trender, högt värde -> kräver tydlig biologisk signal (Höga H och låga G samtidigt)
# Min_Selected  minst antal som måste passera törsklarna (support)
PROFILES = {
    "normal":      {"nH":7, "nC":7, "nG":7, "nF":7, "LIFT_MIN":2.0, "DELTA_Z_MIN":0.30, "MIN_SELECTED":5},
    "strikt":      {"nH":8, "nC":8, "nG":8, "nF":8, "LIFT_MIN":2.3, "DELTA_Z_MIN":0.40, "MIN_SELECTED":15},
    "superstrikt": {"nH":10,"nC":10,"nG":10,"nF":10,"LIFT_MIN":2.6, "DELTA_Z_MIN":0.50, "MIN_SELECTED":25},
}

"""
För Idre: "nH":7, "nC":7, "nG":7, "nF":7, "LIFT_MIN":2.0, "DELTA_Z_MIN":0.30, "MIN_SELECTED":50
För Idre: "nH":8, "nC":8, "nG":8, "nF":8, "LIFT_MIN":2.3, "DELTA_Z_MIN":0.40, "MIN_SELECTED":60
För Idre: "nH":10,"nC":10,"nG":10,"nF":10,"LIFT_MIN":2.6, "DELTA_Z_MIN":0.50, "MIN_SELECTED":70

För Lunsen: "nH":7, "nC":7, "nG":7, "nF":7, "LIFT_MIN":2.0, "DELTA_Z_MIN":0.30, "MIN_SELECTED":5
För Lunsen: "nH":8, "nC":8, "nG":8, "nF":8, "LIFT_MIN":2.3, "DELTA_Z_MIN":0.40, "MIN_SELECTED":15
För Lunsen: "nH":10,"nC":10,"nG":10,"nF":10,"LIFT_MIN":2.6, "DELTA_Z_MIN":0.50, "MIN_SELECTED":25
"""
if STRICTNESS not in PROFILES:
    raise SystemExit("STRICTNESS måste vara 'normal', 'strikt' eller 'superstrikt'.")
P = PROFILES[STRICTNESS]


# Dessa används både för att bygga sökytans percentilspann och när vi slår om percentil -> ABS.
# H/C: ABS-tröskeln får inte understiga golvet.  G: ABS-tröskeln får inte överstiga taket.
H_FLOOR   = {"L": 155.0, "M": 150.0, "H": 135.0}  # dm (P95 höjd) - lägsta rimliga nivå
C_FLOOR   = {"L": 7.00,  "M": 7.00,  "H": 6.90 }  # m (Krondiameter) - lägsta rimliga nivå
G_CEILING = {"L": 1.50,  "M": 1.50,  "H": 1.50 }  # dm/år (growth p95 year) - högsta rimliga nivå
# (F inget golv/tak i enheter här utan begränsas F via percentilspannet i auto-sökandet

# Guardrails (kandidaten MÅSTE klara dessa)
LIFT_MIN      = P["LIFT_MIN"]       # "Enrichment": hur mycket "oftare än väntat" kombinationen uppträder
DELTA_Z_MIN   = P["DELTA_Z_MIN"]    # Effektstyrka: H ska vara högt (positiv z), G lågt (negativ z -> minustecken i delta_z)
MIN_SELECTED  = P["MIN_SELECTED"]   # Support: Minst så här många måste passera tröskeln.
MIN_N_PER_STRATUM = 0   # sätt >0 om du vill hoppa små strata

# ------------------ Hjälpfunktioner ------------------
# Konvertera värde till float på ett säkert sätt
def to_float(x, default=None):
    try:
        v = float(x)
        return v if math.isfinite(v) else default  # Returnera tal om det är ändligt, annars d (t.ex. None)
    except:
        return default  # Misslyckad konvertering → d (default None)
    

# Beräkna percentil (linjär interpolation, type-7)
def percentile_sorted(vlist, p):
    if not vlist: 
        return None                 # Tom lista → saknar percentil
    v = sorted(vlist)               # Sortera stigande
    if p <= 0:   return v[0]        # p ≤ 0 → minvärde
    if p >= 100: return v[-1]       # p ≥ 100 → maxvärde
    k = (len(v)-1)*(p/100.0)        # Fraktionellt index k
    f = int(k); c = min(f+1, len(v)-1)  # Golvindex f, takindex c (klippt mot slutet)
    return v[f] if f==c else v[f] + (k-f)*(v[c]-v[f])  # Interpolera mellan v[f] och v[c]


def percentile_rank(vals, thr):

    """
    Omvänd percentil för tröskeln 'thr' givet rådata 'vals'.
    - Returnerar percentil i intervallet 0..100 (float)
    - Linjär interpolation mellan närliggande värden (matchar 'type-7' tänket)
    - Returnerar None för tom lista
    """
   
    arr = np.asarray(vals, dtype=float)     # Gör om till numpy-array av floats (möjliggör användandet av numpy-operations)
    arr = arr[np.isfinite(arr)]             # Tar bort icke-siffror
    n = arr.size                            # Antal gilitiga värden efter resning.
    if n == 0:
        return None                         # Tom lista -> Går inte att beräkna percentilrank = none
    if n == 1:
        return 100.0                        # Ett enda värde: definiera som rank 100. Går att välja 0 percentilen...
    
    arr.sort()                              # Sortera stigande så index blir meningsfulla.

    # Snabba gränsfall. 
    if thr <= arr[0]:
        return 0.0                          # Om tröskeln ligger under minvärdet -> 0:e percentilen
    if thr >= arr[-1]:
        return 100.0                        # Om tröskeln ligger över maxvärdet -> 100:e percentilen

    # Hitta insättningsposition (vänster) för thr i den sorterade arrayen (arr.sort())
    j = np.searchsorted(arr, thr, side='left')      # "left" betyder första platsen där thr kan sättas UTAN att bryta sorteringen via (searchsorted)

    # Välj det vänstra indexet i paret (i, i+1) som omger thr
    i = max(0, min(j - 1, n - 2))   # clamp till [0, n-2] för att undvika indexfel vid kanter

    # Värdena som omger thr (vänster och höger om thr)
    x0 = arr[i]
    x1 = arr[i + 1]

    # Linjär interpolation till "kontinuerligt index" k
    if x1 == x0:
        k = float(i)                    # Om x0 == x1 (dubbelvärde) kan vi inte interpolera, k = i
    else:
        k = i + (thr - x0) / (x1 - x0)  # Fraktion mellan x0 ochy x1 som thr representerar.

    # Översätt till percentil i 0..100.
    # Matchar vår skala som percentil-funktionen använder.
    rank = 100.0 * (k / (n - 1))
    return float(max(0.0, min(100.0, rank)))    # Klipp resultatet till [0, 100] för säkerhets skull. se [0, 100] som %.


def read_rows(fc, fields):
    # Läser ut fälten + koordinat till en lista av tuples.
    # Tuplen blir som denna: (elev, H, G, C, F, (x, y))
    out=[]
    with arcpy.da.SearchCursor(fc, fields+["SHAPE@XY"]) as cur:
        for *vals, xy in cur:
            out.append(tuple(vals) + (xy,))
    return out

def tercile_cuts(elevs):
    # Beräknar percentil gränser.
    return percentile_sorted(elevs, 33.33), percentile_sorted(elevs, 66.67)

def norm_cat(v):
    """Normalisera kategorivärde till icke-tomt (skapa - om None/blankt)"""
    if v is None: return "-"
    s = str(v).strip()
    return s if s != "" else "-"

def make_stratum_key(elev_band, mark, torv):
    """Nykcel: elev_band|mark_klass=...|torv_klass=... (stabil och läsbar)
    Överväg att ändra till något mer användarvänligt"""
    return f"{elev_band}|mark_klass={norm_cat(mark)}|torv_klass={norm_cat(torv)}"

def split_strata(rows, cutL, cutH):
    """Delar upp rader i strata (grupper) enligt elev_band x mark x torv"""
    strata = {}
    for row in rows:
        elevation, h, g, c, f, mark, torv, _xy = row   # Packa upp
        evf = to_float(elevation, None)                # Skapa float
        if evf is None:                             
            continue                                   # Hoppa okända höjder
        elev_band = 'L' if evf < cutL else ('M' if evf <= cutH else 'H')    # Band
        key = make_stratum_key(elev_band, mark, torv)                       # Strata-nyckel
        strata.setdefault(key, []).append({
            "H": to_float(h, None),
            "G": to_float(g, None),
            "C": to_float(c, None),
            "F": to_float(f, None) if USE_F else None
        })
    return strata

# --------- AUTO-sökrymd: bygg kandidat-percentiler från datan ----------
def _linspace_int(a, b, n):
    """
    Skapar n jämnt fördelade heltalspercentiler mellan a och b, avrundar varje punkt till närmaste heltal och klipper till [0, 100].
    Om n <= 1 returnerar den mittpunkten som heltal.
    Används för att bygga ett rutnät av percentiler för grid-söket
    """
    if n <= 1:
        mid = int(round((a + b) / 2.0))     # Mittpunkt mellan a och b
        return [mid]

    # Generera n jämnt ördelade flyttal mellan a och b
    xs = np.linspace(a, b, n)


    ints = []
    for x in xs:
        val = int(round(x))                 # Runda till närmsta heltalspercentil
        if 0 <= val <= 100:                 # Håll oss inom giltigt percentilintervall
            ints.append(val)
    return ints

def build_auto_percentiles_for_band(band, Hs, Cs, Gs, Fs):
    # Beräkna p_lo/p_hi utifrån golv/tak i ENHETER (per elevband).
    def p_of_floor(vals, floor):
        if not vals:
            return None
        return percentile_rank(vals, floor)     # Percentil som motsvarar golvet

    def p_of_ceiling(vals, ceiling):
        if not vals:
            return None
        return percentile_rank(vals, ceiling)   # Percentil som motsvarar taket

    # H (högre bättre) - börja över golvet i enheter
    pH_lo = 0
    if Hs:
        p_floor_H = p_of_floor(Hs, H_FLOOR[band])                
        if p_floor_H is None:
            pH_lo = 0       # Ingen data -> börja från percentil 0
        else:
            # Ceil() -> välj första heltals precentil som inte ligger under golvet.
            pH_lo = int(min(100, max(0, math.ceil(p_floor_H))))   
    # Undvik p = 100 (max-värdet) för att inte låsa söket på enstaka extremvärden
    pH_hi = 99
    if Hs:
        pH_vals = _linspace_int(pH_lo, pH_hi, P["nH"])      # Jämt spridda heltalspunkter i [pH_lo, 99]
    else:
        pH_vals = []

    # C (högre bättre) - börja över golvet i enheter
    pC_lo = 0
    if Cs:
        p_floor_C = p_of_floor(Cs, C_FLOOR[band])
        if p_floor_C is None:
            pC_lo = 0
        else:
            pC_lo = int(min(100, max(0, math.ceil(p_floor_C))))
    pC_hi = 99
    if Cs:
        pC_vals = _linspace_int(pC_lo, pC_hi, P["nC"])
    else:
        pC_vals = []

    # G (lägre bättre) - sluta vid percentilen som motsvarar taket i enheter
    pG_hi = 100
    if Gs:
        p_ceil_G = p_of_ceiling(Gs, G_CEILING[band])        # Percentilen där G-taket landar i datan
        if p_ceil_G is None:
            pG_hi = 100
        else:
            # Floor () -> välj sista heltals percentil som inte ligger över taket
            pG_hi = int(min(100, max(0, math.floor(p_ceil_G))))
    pG_lo = 0
    if Gs:
        pG_vals = _linspace_int(pG_lo, pG_hi, P["nG"])      # Kandidater i [0, pG_hi]
    else:
        pG_vals = []

    # F (lägre bättre) - standard 0..95
    if USE_F and Fs:
        pF_lo = 0
        pF_hi = 95
        pF_vals = _linspace_int(pF_lo, pF_hi, P["nF"])
    else:
        pF_vals = [None]  # inaktiv dimension

    # Säkerställ att varje aktiv dimension har åtminstone två punkter att prova
    # Om _linspace_int skulle ge distnikta p-värden p.g.a. avrundning / smalt intervall
    if len(pH_vals) < 2:
        pH_vals = sorted(set((pH_lo, pH_hi)))
    if len(pC_vals) < 2:
        pC_vals = sorted(set((pC_lo, pC_hi)))
    if Gs and len(pG_vals) < 2:
        pG_vals = sorted(set((pG_lo, pG_hi)))
    return pH_vals, pC_vals, pG_vals, pF_vals

def to_abs(vals, p, key, band):
    # Gör om percentil p till en absolut tröskel och applciera skyddsräcken.
    thr = percentile_sorted(vals, p) if vals else None
    if thr is None: return None
    # räcken i enheter
    if key == "H":
        thr = max(thr, H_FLOOR[band])   # H får inte gå under golv
    if key == "C":
        thr = max(thr, C_FLOOR[band])   # C får inte gå under golv
    if key == "G":
        thr = min(thr, G_CEILING[band]) # G får inte gå över tak
    return thr

def prevalence(items, thr):
    # Beräkna andel (prevalens) av items som klarar samtliga villkor:
    # H >= h_thr
    # 0 <= G <= G_thr
    # C >= C_thr
    # om USE_F = True: F <= F_thr

    sel = 0; N = 0
    for r in items:
        H = r["H"]; G = r["G"]; C = r["C"]; F = r["F"]
        if H is None or G is None or C is None:
            continue        # Kräver H/G/C för att kunna bedöma
        ok = (H >= thr["H"] and C >= thr["C"] and 0.0 <= G <= thr["G"])
        if USE_F and thr["F"] is not None:
            ok = ok and (F is not None and F <= thr["F"])
        if ok: sel += 1
        N += 1
    return (sel / float(N)) if N else 0.0

def andel_ge(array, thr):
    """
    Andel giltiga värden i array som är greater-or-equal (ge) >= träskeln 'thr'
    """
    a = np.asarray([x for x in array if x is not None], float)  # Filtrera bort None och gör om till float-array
    a = a[np.isfinite(a)]                                       # Ta bort Null eller NaN eller infinite
    return 0.0 if a.size==0 else float(np.mean(a >= thr))       # Andel som uppfiller villkoret (>= thr)

def andel_le(array, thr):
    """
    Andel giltiga värden u array som är less-or-equal (le) <= tröskeln 'thr'
    Används för marginalsannolikheten p(x <= thr).
    """
    a = np.asarray([x for x in array if x is not None], float)  # Filtrera bort None och gör om till float-array
    a = a[np.isfinite(a)]                                       # Ta bort null eller NaN eller infinite
    return 0.0 if a.size==0 else float(np.mean(a <= thr))       # Andel som uppfyuller villkoret (<= thr)

def zmean(values, ref_values):
    """
    Medel-z för "values" relativt fördelningen i "ref_values"
    Används för att uppskatta "effektstorlek": hur långt över/under referensfördelningen ligger de valda
    observationerna i genomsnitt, mätt i standardavvikelser (z-enheter)
    """
    ref = np.asarray([v for v in ref_values if v is not None], float)   # Referens som float-array
    ref = ref[np.isfinite(ref)]                                         # Ta bort Null eller NaN eller infinite
    if ref.size < 2:
        return None                                                     # Kräver minst 2 observationer för standard deviation
    
    mu, sd = ref.mean(), ref.std()
    if sd <= 0:
        return None                                                     # Ingen variation -> ingen z-skalning
    
    arr = np.asarray([v for v in values if v is not None], float)       # Målvärden
    arr = arr[np.isfinite(arr)]                                         # Ta bort Null....
    if arr.size == 0:
        return None                                                     # Om inget finns, finns inget att bedöma
   
    return float(((arr - mu)/sd).mean())                                # Medel-z relativt referensen

# Grid-sök på auto-byggda percentiler
def pick_best_auto(band, items):
    """
    Loopar över auto-genererade percentilrutor och väljer bästa kandidat:
        - Beräknas observed/expected, lift, support och delta_z
        - Filtrerar på LIFT_MIN, MIN_SELECTED och DELTA_Z_MIN
        - Prioriteringar: högst lift -> flest valda -> största delta_z
    Returnerar (best, grids) där grids är provade percentiler per dimension
    """
    Hs=[r["H"] for r in items if r["H"] is not None]    # rålistor per dimension
    Gs=[r["G"] for r in items if r["G"] is not None]
    Cs=[r["C"] for r in items if r["C"] is not None]
    Fs=[r["F"] for r in items if (USE_F and r["F"] is not None)]
    if not Hs or not Gs or not Cs:
        return None, None  # behöver minst H/G/C för att gå vidare

    pH_vals, pC_vals, pG_vals, pF_vals = build_auto_percentiles_for_band(band, Hs, Cs, Gs, Fs)  # Sökrymd

    best=None
    for pH in pH_vals:
        H_thr = to_abs(Hs, pH, "H", band)   # slå om percentil till ABS + räcken
        if H_thr is None:
            continue
        for pC in pC_vals:
            C_thr = to_abs(Cs, pC, "C", band)
            if C_thr is None:
                continue
            for pG in pG_vals:
                G_thr = to_abs(Gs, pG, "G", band)
                if G_thr is None:
                    continue
                for pF in pF_vals:
                    if USE_F and pF is not None:
                        F_thr = to_abs(Fs, pF, "F", band)
                    else:
                        F_thr = None

                    thr = {             # Samlad tröskel
                        "H": H_thr,
                        "C": C_thr,
                        "G": G_thr,
                        "F": F_thr,
                    }

                    # Observed vs. Expected (under oberoendeantagande)
                    obs = prevalence(items, thr)        # Verklig samsannolikhet
                    pH_m = andel_ge(Hs, H_thr)          # Marginal för exp
                    pC_m = andel_ge(Cs, C_thr)
                    pG_m = andel_le(Gs, G_thr)
                    pF_m = andel_le(Fs, F_thr) if (USE_F and F_thr is not None) else 1.0
                    exp_ind = pH_m * pC_m * pG_m * pF_m
                    lift = (obs/exp_ind) if exp_ind>0 else 0.0  # Enrichment (robust mot 0)

                    # Räkna antal valda + spara H/G för delta_z
                    sel_n = 0; N = 0
                    H_sel, G_sel = [], []
                    for r in items:
                        H=r["H"]; G=r["G"]; C=r["C"]; F=r["F"]
                        if H is None or G is None or C is None:
                            continue
                        ok = (H>=H_thr and C>=C_thr and 0.0<=G<=G_thr)      # Kärnfilter
                        if USE_F and F_thr is not None:
                            ok = ok and (F is not None and F <= F_thr)      # Lägg på F-filter om aktivt
                        if ok:
                            sel_n += 1
                            H_sel.append(H); G_sel.append(G)
                        N += 1

                    # Effektstorlek: högre H och lägre G ska “peka uppåt” tillsammans
                    zH = zmean(H_sel, Hs)
                    zG = zmean(G_sel, Gs)
                    delta_z = (zH - (-zG)) if (zH is not None and zG is not None) else -1e9     # Straffa om ej beräkningsbart

                    # Guardrails: kandidat måste klara alla tre
                    if not (lift >= LIFT_MIN and sel_n >= MIN_SELECTED and delta_z >= DELTA_Z_MIN):
                        continue

                    cand = {"thr":thr, "lift":lift, "sel_n":sel_n, "delta_z":delta_z, "N":N,
                            "pH":pH, "pC":pC, "pG":pG, "pF":pF}     # Packa ihop kandidat
                    # Prioritering: 1) LIFT, 2) sel_n, 3) Δz
                    if best is None:
                        best = cand
                    else:
                        higher_lift = cand["lift"] > best["lift"]
                        equal_lift = abs(cand["lift"] - best["lift"]) < 1e-12
                        more_selected = cand["sel_n"] > best["sel_n"]
                        equal_selected = cand["sel_n"] == best["sel_n"]
                        better_delta = cand["delta_z"] > best["delta_z"]

                        if higher_lift:
                            best = cand
                        elif equal_lift and more_selected:
                            best = cand
                        elif equal_lift and equal_selected and better_delta:
                            best = cand

    # returnera även de auto-genererade kandidatlistorna för utskrift
    return best, {"H":pH_vals, "C":pC_vals, "G":pG_vals, "F":pF_vals}

# ------------------ Körning ------------------
print("\n- READ-ONLY: LIFT med auto-sök (utan baspercentiler) -")
print(f"Läge: {STRICTNESS} | LIFT_MIN={LIFT_MIN}, MIN_SELECTED={MIN_SELECTED}, delta_z_min={DELTA_Z_MIN}")

# Kontroll fält
need = [ELEV_FIELD,H_FIELD,G_FIELD,C_FIELD,F_FIELD,MARK_FIELD,TORV_FIELD]   # obligatoriska fält
have = [fld.name for fld in arcpy.ListFields(FC)]                   
for f in need:
    if f not in have:
        raise SystemExit("Saknar fält: " + f)   # Fail om schema saknas
    
# 2) Läs alla rader till minnet: (elev, H, G, C, F, mark, torv)
rows = read_rows(
    FC,
    [ELEV_FIELD, H_FIELD, G_FIELD, C_FIELD, F_FIELD, MARK_FIELD, TORV_FIELD]
)

N = len(rows)
if N == 0:
    raise SystemExit("Inga rader i FC.")    # Inget att göra

print(f"Läste {N} kronor.")     # Visar att koden tuggar

# 3) Skapa höjdband via terciler (L/M/H)
elevs = []
for r in rows:
    val = to_float(r[0], None)      # r[0] = elev
    if val is not None:
        elevs.append(val)

cutL, cutH = tercile_cuts(elevs)    # Bandgränser

cutL_txt = f"{cutL:.2f}"
cutH_txt = f"{cutH:.2f}"
print(f"Bandgränser (terciler): L<{cutL_txt}, M:{cutL_txt}-{cutH_txt}, H>{cutH_txt}")   # Visar att koden tuggar

# 4) Stratifiera på band × mark_klass × torv_klass
strata = split_strata(rows, cutL, cutH)
keys = sorted(strata.keys())    # Stabil utskriftsordning

print(f"Antal strata: {len(keys)}")     # Visar att koden tuggar

for k in keys:
    n_k = len(strata[k])
    print(f"  {k}: N={n_k}")    # Visar storlekar, visar att koden tuggar

# 5) Kör auto-grid per stratum och välj bästa kandidat (om någon klarar räckena)
best_by_stratum = {}

for i, key in enumerate(keys, start=1):
    items = strata[key]


    if MIN_N_PER_STRATUM and len(items) < MIN_N_PER_STRATUM:        # Hoppar över strata med för låg population
        print(
            f"\n({i}/{len(keys)}) {key}: hoppar (för få, "
            f"N={len(items)} < {MIN_N_PER_STRATUM})"
        )
        continue

    # 'L' / 'M' / 'H' - behövs för bandvisa räcken
    elev_band = key.split("|", 1)[0]    # Extrahera första segmentet

    print(f"\n({i}/{len(keys)}) Söker trösklar för stratum: {key}, n = {len(items)}")     # Visar att koden tuggar

    best, grids = pick_best_auto(elev_band, items)      # Kör gridden

    # Logga auto-sökrutnätet (vilka percentiler som faktiskt provades)
    if grids:
        # Hjälpfunktion för att summera en lista av percentiler:
        # Tom lista -> "-" (visar att dimensionen saknade kandidater)
        # Annars: skriv ut lägsta-högsta samt hur många punkter som provades, t.ex. "72-99 (10p)"
        def span(vs):
            if not vs:
                return "-"                  # inga provade percentiler i denna dimension
            lo = min(vs)                    # Lägst provad percentil
            hi = max(vs)                    # Högst provad percentil
            n = len(vs)                     # Antal provade punkter
            return f"{lo}-{hi} ({n}p)"      # Formatera kompakt "lo-hi (n-percentil)"

        # Bygg en läsbar statusrad för varje dimension som var aktiv i gridden
        # Exempel: "H 72-99 (10p), C 65-99 (10p), G 0-85 (10p), F 0-95 (10p)"
        parts = [
            f"H {span(grids['H'])}",    # H (högre är bättre) - vilka p testades?
            f"C {span(grids['C'])}",    # C *(högre är bättre)
            f"G {span(grids['G'])}",    # G (lägre är bättre)
        ]
        if USE_F:
            f_vals = [x for x in grids['F'] if x is not None]
            parts.append(f"F {span(f_vals)}")       # Lägg till F om relevant

        # Samlad utskrift: ger snabb överblick över auto-genererade sökintervall per dimension
        # Visar att koden tuggar
        print("  Auto-sökpercentiler: " + ", ".join(parts))

    # om ingen kandidat passerade guardails (LIFT_MIN, MIN_SELECTED, DELTA_Z_MIN)
    # så printar "ingen evidens" för stratumet och går vidare
    if not best:
        print("  Ingen evidens (ingen kandidat klarade LIFT/Support/Effekt) — "
              "inga trösklar för detta stratum.")
        continue

    # Bästa ABS-trösklarna i enheter
    thr = best["thr"]

    # 6) Översätt tillbaka till ungefärliga percentiler (för läsbar utskrift)
    Hs = [r["H"] for r in items if r["H"] is not None]
    Cs = [r["C"] for r in items if r["C"] is not None]
    Gs = [r["G"] for r in items if r["G"] is not None]
    Fs = [r["F"] for r in items if (USE_F and r["F"] is not None)]

    # För varje dimension: beräkna vilken percentil (p) den valda ABS-tröskeln hamnar på i just detta stratum
    # Detta gör varje utskrift mer intuitiv (p85 etc.) än rena enheter och hjälper vid jämförelser mellan strata.
    if Hs:
        pH_eff = percentile_rank(Hs, thr["H"])      # p för H-tröskeln
    else:
        pH_eff = None

    if Cs:
        pC_eff = percentile_rank(Cs, thr["C"])      # p för C-tröskeln
    else:
        pC_eff = None

    if Gs:
        pG_eff = percentile_rank(Gs, thr["G"])      # p för G-tröskeln
    else:
        pG_eff = None

    if USE_F and Fs and (thr["F"] is not None):
        pF_eff = percentile_rank(Fs, thr["F"])      # p för F-tröskeln
    else:
        pF_eff = None

    # Utskrift: stöd (sel/total), enrichment (LIFT), effektstorek (delta_z)
    sel_n   = best['sel_n']                                            # Antal som passerade trösklarna
    total_n = best['N']                                                # Antal bedömda i stratumet
    lift_val = best['lift']                                            # Observed/Expected (oberoendeapprox)
    dz_val   = best['delta_z']                                         # Signalstyrka i z-enheter (H upp, G ner)

    print(f"  sel={sel_n}/{total_n}  LIFT={lift_val:.2f}  delta_z≈{dz_val:.2f}")

    # Utskrift av valda trösklar i absoluta enheter (dm/m/dm/år)
    abs_parts = [
        f"H={thr['H']:.2f} dm",                                        # H-tröskel i dm
        f"C={thr['C']:.2f} m",                                         # C-tröskel i m
        f"G={thr['G']:.2f} dm/år",                                     # G-tröskel i dm/år
    ]
    if USE_F and thr["F"] is not None:
        abs_parts.append(f"F={thr['F']:.2f} dm")                       # F-tröskel (om aktiv)

    print("     ABS:  " + ",  ".join(abs_parts))

    # för en kompakt percentile utskrift: används 0.0 som fallback när p är None (för formateringens skull)
    p_h = 0.0 if pH_eff is None else pH_eff
    p_c = 0.0 if pC_eff is None else pC_eff
    p_g = 0.0 if pG_eff is None else pG_eff

    # utskrift av trösklarna uttryckte i percentiler (lättare att "läsa av" relativ position i fördelningen)
    pctl_txt = f"     PCTL: H≈p{p_h:.1f}, C≈p{p_c:.1f}, G≈p{p_g:.1f}"
    if USE_F and pF_eff is not None:
        pctl_txt += f", F≈p{pF_eff:.1f}"

    print(pctl_txt)

    # Spara sammanfattning för detta stratum och använd det senare i "kopieringsblocket" av koden
    best_by_stratum[key] = {
        "thr": thr,                                                         # ABS-trösklar (används av andra skript)
        "pctl_eff": {"H": pH_eff, "C": pC_eff, "G": pG_eff, "F": pF_eff},   # Ungefärliga percentiler
        "lift": best["lift"],                                               # För ev. diagnostik/rapport
        "sel_n": best["sel_n"],                                             # Support
    }


# 7) Kopieringsblock (direkt körbara i andra skript)
if best_by_stratum:
    print("\nABS_THRESHOLDS_BY_STRATUM = {")
    for key in sorted(best_by_stratum.keys()):
        thr = best_by_stratum[key]["thr"]
        lift = best_by_stratum[key]["lift"]

        if USE_F and (thr["F"] is not None):
            F_out = thr["F"]
        else:
            F_out = None

        # ABS-trösklar i enheter - detta block kan copy-pastas in i is_old = 0/1 skriptet
        print(
            f"    '{key}': {{'H': {thr['H']:.2f}, 'G': {thr['G']:.2f}, "
            f"'C': {thr['C']:.2f}, 'F': {repr(F_out)}, 'LIFT': {lift:.3f}}},"
        )
    print("}")

    print("\nPCTS_BY_STRATUM = {  # ungefärliga percentiler som motsvarar ABS inom stratum")
    for key in sorted(best_by_stratum.keys()):
        pp = best_by_stratum[key]["pctl_eff"]

        # rundad percentil för snygg utskrift. Saknas p -> använd rimlig "neutral" fallback (H/G/C = 0.0, F = 100.0)
        # H/C tolkas som "högre är bättre", så p = 0.0 visar "ingen/låg position", F inaktiv -> visa som 100
        if pp['H'] is None:
            Hp = 0.0
        else:
            Hp = round(pp['H'], 1)

        if pp['G'] is None:
            Gp = 0.0
        else:
            Gp = round(pp['G'], 1)

        if pp['C'] is None:
            Cp = 0.0
        else:
            Cp = round(pp['C'], 1)

        if USE_F and pp['F'] is not None:
            Fp = round(pp['F'], 1)
        else:
            Fp = 100.0

        print(f"    '{key}': {{'H': {Hp}, 'G': {Gp}, 'C': {Cp}, 'F': {Fp}, 'LIFT': {round(lift,3)}}},")
    print("}")


print("\nKLART (READ-ONLY, AUTO-SÖK). Justera STRICTNESS, LIFT_MIN/MIN_SELECTED/DELTA_Z_MIN vid behov.")
