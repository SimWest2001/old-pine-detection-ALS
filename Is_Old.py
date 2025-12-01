# -*- coding: utf-8 -*-
# =============================================================================
# Sätter/uppdaterar is_old baserat på per-stratum PCTS_BY_STRATUM (percentiler)
#  - Beräknar L/M/H (elev_band) från Elev_mean (terciler)
#  - Slår om percentiler till absoluta trösklar per stratum (H,G,C,F)
#  - Använder G-cap per band (kan stängas av)
#  - Uppdaterar is_old = 1/0
#  - STRATUM: elev_band × mark_klass × torv_klass
# =============================================================================

import arcpy, math

# ------------------------------ INSTÄLLNINGAR ------------------------------
# Feature class / lager (namn i Contents eller full sökväg)
FC = r"C:\Your\File_Path\To\Tradkronor_PreProcessed"   # <- ändra vid behov

# Fält i FC
ELEV_FIELD   = "Elev_mean"       # m, används för att skapa L/M/H
BAND_FIELD   = "elev_band"       # skapas/uppdateras (TEXT, 1 tecken: L/M/H)
IS_OLD_FIELD = "is_old"          # skapas/uppdateras (SHORT, 0/1)
H_FIELD      = "O2_PCT95"        # dm — högre bättre (höjd P95)
G_FIELD      = "gr_p95yr"        # dm/år — lägre bättre (0..)
C_FIELD      = "cr_dia"          # m — högre bättre (kron-diameter)
F_FIELD      = "cr_flat"         # dm — lägre bättre (frivillig)

# Stratifieringsfält (ANVÄNDS NU)
MARK_FIELD   = "mark_klass"
TORV_FIELD   = "torv_klass"

# Använd F i regeln om percentil finns?
USE_F = True

# G-cap (tak för G per elevband). Sätt USE_G_CAP=False för att stänga av.
USE_G_CAP = True
G_CAP = {"L": 1.35, "M": 1.55, "H": 1.55}

# Minsta antal kompletta rader i ett stratum för att beräkna percentil-trösklarna
MIN_VALID_PER_STRATUM = 25

# ---------------------- DINA PERCENTILER PER STRATUM -----------------------
# Dessa percentiler omvandlas automatiskt till ABS-trösklar per stratum.
# Nyckelns format: 'Band|mark_klass=<text>|torv_klass=<text>'
PCTS_BY_STRATUM = {  # ungefärliga percentiler som motsvarar ABS inom stratum
    'H|mark_klass=frisk-fuktig|torv_klass=>= 30 cm torvdjup': {'H': 86.9, 'G': 42.8, 'C': 99.0, 'F': 95.0, 'LIFT': 7.877},
    'H|mark_klass=frisk-fuktig|torv_klass=>= 40 cm torvdjup': {'H': 81.5, 'G': 47.7, 'C': 93.5, 'F': 95.0, 'LIFT': 7.877},
    'H|mark_klass=frisk-fuktig|torv_klass=>= 50 cm torvdjup': {'H': 90.7, 'G': 68.2, 'C': 98.9, 'F': 95.0, 'LIFT': 7.877},
    'H|mark_klass=frisk-fuktig|torv_klass=mineraljord': {'H': 91.9, 'G': 35.7, 'C': 98.9, 'F': 95.0, 'LIFT': 7.877},
    'H|mark_klass=fuktig-blöt|torv_klass=>= 50 cm torvdjup': {'H': 95.8, 'G': 87.9, 'C': 95.7, 'F': 95.0, 'LIFT': 7.877},
    'H|mark_klass=torr-frisk|torv_klass=mineraljord': {'H': 94.7, 'G': 21.4, 'C': 99.0, 'F': 95.0, 'LIFT': 7.877},
    'L|mark_klass=frisk-fuktig|torv_klass=>= 30 cm torvdjup': {'H': 81.7, 'G': 36.7, 'C': 92.9, 'F': 95.0, 'LIFT': 7.877},
    'L|mark_klass=frisk-fuktig|torv_klass=mineraljord': {'H': 73.9, 'G': 32.0, 'C': 99.0, 'F': 95.0, 'LIFT': 7.877},
    'L|mark_klass=torr-frisk|torv_klass=mineraljord': {'H': 93.8, 'G': 25.2, 'C': 98.9, 'F': 95.0, 'LIFT': 7.877},
    'M|mark_klass=frisk-fuktig|torv_klass=>= 30 cm torvdjup': {'H': 85.8, 'G': 30.9, 'C': 94.0, 'F': 95.0, 'LIFT': 7.877},
    'M|mark_klass=frisk-fuktig|torv_klass=>= 40 cm torvdjup': {'H': 89.8, 'G': 38.6, 'C': 90.7, 'F': 95.0, 'LIFT': 7.877},
    'M|mark_klass=frisk-fuktig|torv_klass=>= 50 cm torvdjup': {'H': 86.7, 'G': 54.2, 'C': 93.9, 'F': 95.0, 'LIFT': 7.877},
    'M|mark_klass=frisk-fuktig|torv_klass=mineraljord': {'H': 90.8, 'G': 33.1, 'C': 93.6, 'F': 95.0, 'LIFT': 7.877},
    'M|mark_klass=torr-frisk|torv_klass=mineraljord': {'H': 93.7, 'G': 15.5, 'C': 98.9, 'F': 95.0, 'LIFT': 7.877},
}


# ----------------------------- HJÄLPFUNKTIONER -----------------------------
def add_msg(m): arcpy.AddMessage(m)
def add_warn(m): arcpy.AddWarning(m)
def add_err(m): arcpy.AddError(m)

def ensure_field(fc, name, ftype="TEXT", length=4):
    names = [f.name.lower() for f in arcpy.ListFields(fc)]
    if name.lower() not in names:
        if ftype.upper() == "TEXT":
            arcpy.management.AddField(fc, name, ftype, field_length=length)
        else:
            arcpy.management.AddField(fc, name, ftype)

def to_float(x, default=None):
    try:
        v = float(x)
        return v if math.isfinite(v) else default
    except:
        return default

def percentile(values, p):
    vals = [to_float(v, None) for v in values if to_float(v, None) is not None]
    if not vals:
        return None
    v = sorted(vals)
    if p <= 0:   return v[0]
    if p >= 100: return v[-1]
    k = (len(v)-1)*(p/100.0)
    f = int(k); c = min(f+1, len(v)-1)
    return v[f] if f == c else v[f] + (k-f)*(v[c]-v[f])

def update_elev_band(fc):
    """Skapar/uppdaterar L/M/H i BAND_FIELD från Elev_mean (terciler)."""
    ensure_field(fc, BAND_FIELD, "TEXT", length=1)
    elevs=[]
    with arcpy.da.SearchCursor(fc, [ELEV_FIELD]) as cur:
        for (ev,) in cur:
            evf = to_float(ev, None)
            if evf is not None:
                elevs.append(evf)
    if not elevs:
        raise SystemExit(f"{ELEV_FIELD} saknar giltiga värden.")
    cutL = percentile(elevs, 33.33)
    cutH = percentile(elevs, 66.67)
    add_msg(f"Bandgränser (terciler): L<{cutL:.2f}, M:{cutL:.2f}-{cutH:.2f}, H>{cutH:.2f}")

    cL=cM=cH=0
    with arcpy.da.UpdateCursor(fc, [ELEV_FIELD, BAND_FIELD]) as cur:
        for ev, _ in cur:
            b=None
            evf = to_float(ev, None)
            if evf is not None:
                if evf < cutL: b='L'; cL+=1
                elif evf <= cutH: b='M'; cM+=1
                else: b='H'; cH+=1
            cur.updateRow((ev, b))
    add_msg(f"Bandfördelning: L={cL}, M={cM}, H={cH}")

def _norm_cat(v):
    if v is None:
        return "—"
    s = str(v).strip()
    return s if s else "—"

def parse_stratum_key(key):
    """
    'H|mark_klass=frisk-fuktig|torv_klass=mineraljord'
    -> ('H', 'frisk-fuktig', 'mineraljord')
    """
    parts = key.split("|")
    band = parts[0]
    def val_of(p):
        if "=" in p:
            return p.split("=",1)[1]
        return ""
    mark = val_of(parts[1]) if len(parts)>1 else ""
    torv = val_of(parts[2]) if len(parts)>2 else ""
    return band, mark, torv

def make_stratum_key(band, mark, torv):
    return f"{band}|mark_klass={_norm_cat(mark)}|torv_klass={_norm_cat(torv)}"

def fmt_num(x, nd=2):
    if x is None: return "—"
    if isinstance(x, float) and math.isinf(x): return "inf"
    try:
        return f"{float(x):.{nd}f}"
    except:
        return str(x)

# ------------------------------- HUVUDLOGIK -------------------------------
def main():
    arcpy.env.overwriteOutput = True

    # 1) Validera fält (om F saknas -> stäng av USE_F och fortsätt)
    need = [ELEV_FIELD, H_FIELD, G_FIELD, C_FIELD, MARK_FIELD, TORV_FIELD]
    have = [f.name for f in arcpy.ListFields(FC)]
    missing = [f for f in need if f not in have]
    if F_FIELD not in have:
        add_warn("Fältet för F saknas -> USE_F=False (kör utan F).")
        use_f_now = False
    else:
        use_f_now = USE_F

    if missing:
        raise SystemExit("Saknar fält: " + ", ".join(missing))

    if not PCTS_BY_STRATUM:
        raise SystemExit("PCTS_BY_STRATUM är tomt. Inget att tillämpa.")

    # 2) Uppdatera/skriv elev_band (L/M/H)
    add_msg("\n--- Sätter/uppdaterar elev_band från Elev_mean ---")
    update_elev_band(FC)

    # 3) Läs in alla rader till minne för att räkna percentiler per stratum
    add_msg("\nLäser data...")
    fields_for_read = [BAND_FIELD, MARK_FIELD, TORV_FIELD, H_FIELD, G_FIELD, C_FIELD]
    if use_f_now:
        fields_for_read.append(F_FIELD)
    rows = []
    with arcpy.da.SearchCursor(FC, fields_for_read) as cur:
        for r in cur:
            # packa upp beroende på om F finns
            if use_f_now:
                b, mark, torv, H, G, C, F = r
            else:
                b, mark, torv, H, G, C = r
                F = None
            rows.append((b, mark, torv, H, G, C, F))

    # 4) Beräkna ABS-trösklar per stratum från PCTS_BY_STRATUM
    add_msg("\n--- Beräknar trösklar (ABS) per stratum från dina percentiler ---")
    abs_by_stratum = {}
    skipped = 0
    size_by_key = {}

    for key, pcts in PCTS_BY_STRATUM.items():
        band, mark_txt, torv_txt = parse_stratum_key(key)

        # Plocka ut värden i just detta stratum
        Hs, Gs, Cs, Fs = [], [], [], []
        n_all = 0
        for b, mark, torv, H, G, C, F in rows:
            if b != band:
                continue
            if _norm_cat(mark) != _norm_cat(mark_txt): continue
            if _norm_cat(torv) != _norm_cat(torv_txt): continue
            n_all += 1
            Hs.append(to_float(H, None))
            Gs.append(to_float(G, None))
            Cs.append(to_float(C, None))
            if use_f_now:
                Fs.append(to_float(F, None))
            else:
                Fs.append(None)

        size_by_key[key] = n_all  # spara storlek för logg

        # Rensa None inför percentiler
        Hs = [x for x in Hs if x is not None]
        Gs = [x for x in Gs if x is not None]
        Cs = [x for x in Cs if x is not None]
        Fs_valid = [x for x in Fs if (x is not None)] if use_f_now else []

        if min(len(Hs), len(Gs), len(Cs)) < MIN_VALID_PER_STRATUM:
            add_warn(f"Hoppar {key}: för få kompletta värden (H/G/C) i stratum (N={n_all}).")
            skipped += 1
            continue

        # Slå om percentil -> absolut tröskel
        H_thr = percentile(Hs, pcts.get("H", 80))
        G_thr = percentile(Gs, pcts.get("G", 80))
        C_thr = percentile(Cs, pcts.get("C", 80))
        F_thr = percentile(Fs_valid, pcts.get("F", 100)) if (use_f_now and Fs_valid) else None

        # G-cap per band (valfritt)
        if USE_G_CAP and G_thr is not None:
            cap = G_CAP.get(band)
            if cap is not None:
                try:
                    G_thr = min(float(G_thr), float(cap))
                except:
                    pass
            if G_thr < 0:
                G_thr = 0.0

        abs_by_stratum[key] = {"H": H_thr, "G": G_thr, "C": C_thr, "F": F_thr}
        add_msg(f"  {key}: N={n_all}  ->  H={fmt_num(H_thr)} dm, G={fmt_num(G_thr)} dm/år, C={fmt_num(C_thr)} m" +
                (f", F={fmt_num(F_thr)} dm" if (use_f_now and F_thr is not None) else ""))

    if not abs_by_stratum:
        raise SystemExit("Inga stratum fick beräknade trösklar – avbryter.")

    # 5) Uppdatera is_old med ABS-trösklar per stratum
    add_msg("\n--- Uppdaterar is_old ---")
    ensure_field(FC, IS_OLD_FIELD, "SHORT")

    total = pos = miss = 0
    fields_for_update = [BAND_FIELD, MARK_FIELD, TORV_FIELD, H_FIELD, G_FIELD, C_FIELD]
    if use_f_now:
        fields_for_update.append(F_FIELD)
    fields_for_update.append(IS_OLD_FIELD)

    with arcpy.da.UpdateCursor(FC, fields_for_update) as cur:
        for row in cur:
            if use_f_now:
                band, mark, torv, Hv, Gv, Cv, Fv, _old = row
            else:
                band, mark, torv, Hv, Gv, Cv, _old = row
                Fv = None

            key = make_stratum_key(band, mark, torv)
            thr = abs_by_stratum.get(key)
            val = 0
            if thr is None:
                miss += 1
            else:
                if None not in (Hv, Gv, Cv):
                    try:
                        h_ok = float(Hv) >= float(thr["H"])
                        g_ok = (0.0 <= float(Gv) <= float(thr["G"]))
                        c_ok = float(Cv) >= float(thr["C"])
                        if use_f_now and (thr.get("F") is not None):
                            f_ok = (Fv is not None) and (float(Fv) <= float(thr["F"]))
                        else:
                            f_ok = True
                        val = 1 if (h_ok and g_ok and c_ok and f_ok) else 0
                    except:
                        val = 0

            # skriv tillbaka
            if use_f_now:
                cur.updateRow((band, mark, torv, Hv, Gv, Cv, Fv, val))
            else:
                cur.updateRow((band, mark, torv, Hv, Gv, Cv, val))

            total += 1
            if val == 1: pos += 1

    rate = (100.0 * pos / total) if total else 0.0
    add_msg(f"\nResultat: is_old=1 för {pos}/{total} rader ({rate:.3f}%)")
    add_msg(f"Rader utan matchande stratum i PCTS_BY_STRATUM: {miss}")

    # 6) Summering per band
    add_msg("\n--- Summering per band ---")
    for b in ("L","M","H"):
        where = f"{arcpy.AddFieldDelimiters(FC, BAND_FIELD)} = '{b}'"
        tot = posb = 0
        with arcpy.da.SearchCursor(FC, [IS_OLD_FIELD], where_clause=where) as cur:
            for (v,) in cur:
                tot += 1
                if v == 1: posb += 1
        if tot:
            add_msg(f"{b}: {posb}/{tot} = {100.0*posb/tot:.3f}%")

    # Tydlig konsolrad som visar att skriptet är klart
    print(f"KLART - is_old uppdaterat i {FC}.  Flagged: {pos}/{total} ({rate:.3f}%).  Strata utan match: {miss}")

# ------------------------------- KÖRNING -------------------------------
if __name__ == "__main__":
    try:
        main()
        add_msg("\nKLART.")
        # Extra print som backup om AddMessage inte syns:
        print("Skriptet är klart.")
    except Exception as e:
        add_err(str(e))
        print(f"Fel: {e}")
        raise
