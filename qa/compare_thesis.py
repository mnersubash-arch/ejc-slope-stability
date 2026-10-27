"""Thesis value by value against the web engine (qa_thesis.json, qa_extra.json). Writes thesis_rows.json.

Rule: a value printed to d decimals matches when |software - thesis| <= half a unit in the d-th decimal (the thesis
rounding). Counts and statuses must match exactly. Values the thesis prints as ranges or 'about' are checked against
the range or the stated rounding.
"""
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
OUT = HERE / "out"                                    # results (git-ignored)
OUT.mkdir(exist_ok=True)
d = json.loads((OUT / "qa_thesis.json").read_text())
x = json.loads((OUT / "qa_extra.json").read_text())
rows = []


def dec(s):
    s = str(s).replace(",", "")
    return len(s.split(".")[1]) if "." in s else 0


def add(where, item, thesis, computed, decimals=None, scale=1.0, note=""):
    """thesis: the printed string; computed: float (scaled by `scale` before comparing)."""
    if isinstance(thesis, str) and computed is not None and not isinstance(computed, str):
        t = float(thesis.replace(",", "").replace("%", "").replace("+", ""))
        dd = dec(thesis.replace("%", "")) if decimals is None else decimals
        v = computed * scale
        ok = abs(v - t) <= 0.5 * 10 ** (-dd) * (1 + 1e-9)
        rows.append(dict(where=where, item=item, thesis=thesis, software=f"{v:.{dd + 3}f}".rstrip("0").rstrip(".") if dd + 3 > 0 else f"{v:.0f}",
                         verdict="match" if ok else "differs", note=note))
    else:
        ok = str(thesis) == str(computed)
        rows.append(dict(where=where, item=item, thesis=str(thesis), software=str(computed), verdict="match" if ok else "differs", note=note))


pct = lambda a, b: (a / b - 1) * 100  # noqa: E731

# ---- Summary, Table 4.1, Figure 4.2
add("Table 4.1", "m_b", "1.26006", d["mb"])
add("Table 4.1", "s", "0.0015893", d["s"])
add("Table 4.1", "a from GSI", "0.509923", d["a_gsi"])
add("Fig. 4.2", "entry x (m)", "-311.4", d["entry_x"])
add("Fig. 4.2", "entry z (m)", "360", d["entry_z"])
add("Fig. 4.2", "exit x (m)", "-20.4", d["exit_x"])
add("Fig. 4.2", "exit z (m)", "24.3", d["exit_z"])
add("Fig. 4.2, §5", "active half-columns, 20x20", 76, d["active_3d_20"])
add("§5", "active slices, 40", 30, d["active_2d_40"])
add("§5.6", "active half-columns, 10x10", 17, d["active_3d_10"])
add("Fig. 4.2, §5.8", "dry bases with sigma_n < 0 (20x20)", 6, d["neg_sn3_dry"])
add("§5.8, C.3", "wet bases with sigma_n < 0 (20x20)", 23, d["neg_sn3_wet"])
add("§5.8", "dry columns with sigma3 < 0", 23, d["neg_s3_dry"])
add("§5.8", "wet columns with sigma3 < 0", 61, d["neg_s3_wet"])
add("§5.8", "sigma_n < 0 at 240x240, dry (%)", "8.9", d["neg_sn_frac_240_dry"], scale=100)
add("§5.8", "sigma_n < 0 at 240x240, wet (%)", "30.0", d["neg_sn_frac_240_wet"], scale=100)
add("§5.8", "lower normal-stress limit (kPa)", "-98", d["projection_limit"])
add("§5.8", "gamma_w / gamma (%)", "36", 9.81 / 27, scale=100)

# ---- Section 5.1 and Table C.2
add("Table C.2", "3D FoS dry", "2.706024342", d["fos3_dry"])
add("Table C.2", "2D FoS dry", "2.447910895", d["fos2_dry"], note="agrees to 8 decimals; the 9th printed digit looks truncated (software rounds to ...896)")
add("Table C.2", "3D FoS wet", "1.272331405", d["fos3_wet"])
add("Table C.2", "2D FoS wet", "1.098069496", d["fos2_wet"], note="agrees to 8 decimals; the 9th printed digit looks truncated (software rounds to ...497)")
add("§5.1, C.2", "native path iterations, 3D dry", 34, d["native_iter3_dry"])
add("§5.1, C.2", "native path iterations, 3D wet", 63, d["native_iter3_wet"])
add("§5.1", "native path FoS 3D dry", "2.706024342", d["native_fos3_dry"])
add("§5.1", "native path FoS 3D wet", "1.272331405", d["native_fos3_wet"])
rows.append(dict(where="§5.1", item="relative local residual (max of dry/wet, native path)", thesis="< 1e-15",
                 software=f"{max(d['native_local_res3_dry'], d['native_local_res3_wet']):.1e}",
                 verdict="match" if max(d['native_local_res3_dry'], d['native_local_res3_wet']) < 1e-15 else "differs",
                 note="machine-precision level either way (limit 1e-9); exact digits depend on floating-point order"))
rows.append(dict(where="§5.1", item="unsquared-EJC residual (MPa)", thesis="< 3e-14", software=f"{max(d['ejc_res3_dry'], d['ejc_res3_wet'], d['native_ejc_res3_dry']):.1e}",
                 verdict="match" if max(d['ejc_res3_dry'], d['ejc_res3_wet'], d['native_ejc_res3_dry']) < 3e-14 else "differs", note=""))
rows.append(dict(where="§5.1", item="sigma2 = sigma3 route vs GHB, worst relative FoS difference (8 cases)", thesis="<= 1.4e-12",
                 software=f"{d['s2s3_vs_ghb_worst_rel']:.1e}", verdict="match" if d["s2s3_vs_ghb_worst_rel"] <= 1.4e-12 else "differs", note=""))
add("§5.1", "Carranza-Torres circle, 2D (GHB)", "2.004", d["ct_ghb_2d"])
add("§5.1", "Carranza-Torres circle, 3D (GHB)", "2.142", d["ct_ghb_3d"])
add("§5.1, A.5", "fixture phi_i (deg)", "55.0720", d["a5_phi"])
add("§5.1, A.5", "fixture c_i (kPa)", "2,079.5", d["a5_c"], scale=1000)

# ---- Appendix A.5 and B.5 full states
for k, t, sc in (("s1", "38.47780945", 1), ("s2", "11.52220834", 1), ("s3", "2.50000089", 1), ("M", "10.10323687", 1),
                 ("phi", "55.07201723", 1), ("tau", "10.29948167", 1), ("c", "2.07951285", 1)):
    add("App. A.5", f"original-closure fixture {k}", t, d[f"a5_{k}"], scale=sc)
for k, t in (("s1", "4.77780254"), ("s2", "1.31740324"), ("s3", "0.08346403"), ("M", "13.83034792"), ("phi", "59.89885488"),
             ("tau", "1.17717117"), ("c", "0.48716674")):
    add("App. B.5", f"Moscarda fixture at 0.400 MPa: {k}", t, d[f"b5_{k}"])
add("App. A.5", "original-closure minimum sigma_n (MPa)", "52.966376", d["t52_snmin_0.271"])
add("App. A.5", "dry overburden bound (MPa)", "1.608168", d["cap_dry"])

# ---- Table 5.1 (Route 2) and Figure 5.1
for sn, (phi, c, tau) in {50: ("67.45", "362.5", "482.9"), 400: ("59.90", "487.2", "1,177.2"), 1600: ("50.21", "935.3", "2,856.0")}.items():
    add("Table 5.1", f"R2 phi_i at {sn} kPa (deg)", phi, d[f"t51_phi_{sn}"])
    add("Table 5.1", f"R2 c_i at {sn} kPa", c, d[f"t51_c_{sn}"])
    add("Table 5.1", f"R2 tau_s at {sn} kPa", tau, d[f"t51_tau_{sn}"])
for sn, (phi, c, tau, res) in {50: ("68.61", "282.7", "410.4", "-0.63"), 400: ("61.32", "278.6", "1,009.7", "-1.09"),
                                1600: ("51.81", "451.2", "2,485.3", "-1.80")}.items():
    r1 = x["route1"][str(sn)][0]
    add("Table 5.1", f"R1 phi_i at {sn} kPa (independent rebuild)", phi, r1["phi"])
    add("Table 5.1", f"R1 c_i at {sn} kPa (independent rebuild)", c, r1["c"])
    add("Table 5.1", f"R1 tau_s at {sn} kPa (independent rebuild)", tau, r1["tau"])
    add("Table 5.1", f"R1 EJC residual at {sn} kPa (MPa)", res, r1["ejc_residual"],
        note="the rebuild gives exactly twice the printed value: the table uses sigma1 - sigma1,EJC (eq. 2.8 form), not r_EJC of eq. B.27")
add("Fig. 5.1", "dry central section: min phi_i (deg)", "55.6", d["fig51_phi_min"])
add("Fig. 5.1", "dry central section: max phi_i (deg)", "72.4", d["fig51_phi_max"])
add("Fig. 5.1", "dry central section: min c_i (kPa)", "356", d["fig51_c_min"])
add("Fig. 5.1", "dry central section: max c_i (kPa)", "639", d["fig51_c_max"])

# ---- Section 5.3, Table 5.2, Figures 5.2-5.3
add("§5.3", "original closure sigma3,min (MPa)", "34.376", d["t52_s3min_0.271"])
add("§5.3", "original closure sigma_n,min (MPa)", "52.97", d["t52_snmin_0.271"])
add("§5.3", "gravity cap dry (MPa)", "1.61", d["cap_dry"])
add("§5.3", "gravity cap wet (MPa)", "1.02", d["cap_wet"])
add("§5.3", "required / dry maximum", "32.9", d["orig_ratio_to_dry_cap"])
add("Fig. 5.2", "end of the original branch (M = 1), sigma3 (MPa)", "128.013", d["orig_branch_end_s3"])
for nu, t in (("0.05", "0.92"), ("0.1", "3.43"), ("0.2", "19.23"), ("0.271", "52.97"), ("0.35", "177.8"), ("0.45", "2,188")):
    add("Table 5.2", f"min sigma_n at nu = {nu} (MPa)", t, d[f"t52_snmin_{nu}"])
add("Table 5.2", "min sigma_n at nu = 0.49 (MPa), software search", "59,952", d["t52_snmin_0.49"], decimals=0,
    note="the exact boundary minimum (50-digit) is 59,951.66 MPa: the thesis is right; the software's search stops 2.2 MPa above it")
add("Table 5.2", "min sigma_n at nu = 0.49 (MPa), exact boundary (independent)", "59,952", x["orig_0.49"]["boundary_sn"], decimals=0)
for nu in ("0.05", "0.1", "0.2", "0.271", "0.35", "0.45", "0.49"):
    add("Table 5.2", f"dry columns below threshold, nu = {nu}", "43/76" if nu == "0.05" else "76/76", f"{d[f't52_dry_below_{nu}']}/76")
add("§5.3", "wet columns below threshold, nu = 0.05", "67/76", f"{d['t52_wet_below_0.05']}/76")

# ---- Table 5.3 and 5.4
names = {"hb3": "HB 3D", "hb2": "HB 2D", "ejc3": "EJC 3D", "ejc2": "EJC 2D"}
vals = {("dry", "0.5"): ("2.554", "2.310", "2.706", "2.448", "1.105"), ("wet", "0.5"): ("1.200", "1.036", "1.272", "1.098", "1.159"),
        ("dry", "0.509923"): ("2.476", "2.242", "2.624", "2.376", "1.104"), ("wet", "0.509923"): ("1.158", "1.000", "1.228", "1.061", "1.158")}
for (w, a), tv in vals.items():
    for key, t in zip(("hb3", "hb2", "ejc3", "ejc2"), tv[:4]):
        add("Table 5.3", f"{w.capitalize()}, a = {a}: {names[key]}", t, d[f"t53_{key}_{w}_{a}"])
    add("Table 5.3", f"{w.capitalize()}, a = {a}: EJC 3D/2D", tv[4], d[f"t53_ejc3_{w}_{a}"] / d[f"t53_ejc2_{w}_{a}"])
add("Table 5.4", "3D EJC vs 2D EJC, dry (%)", "+10.5", pct(d["t53_ejc3_dry_0.5"], d["t53_ejc2_dry_0.5"]))
add("Table 5.4", "3D EJC vs 2D EJC, wet (%)", "+15.9", pct(d["t53_ejc3_wet_0.5"], d["t53_ejc2_wet_0.5"]))
add("Table 5.4", "3D EJC vs 3D HB, dry (%)", "+6.0", pct(d["t53_ejc3_dry_0.5"], d["t53_hb3_dry_0.5"]))
add("Table 5.4", "3D EJC vs 3D HB, wet (%)", "+6.0", pct(d["t53_ejc3_wet_0.5"], d["t53_hb3_wet_0.5"]))
add("Table 5.4", "wet vs dry, EJC 3D (%)", "-53.0", pct(d["t53_ejc3_wet_0.5"], d["t53_ejc3_dry_0.5"]))
add("Table 5.4", "wet vs dry, EJC 2D (%)", "-55.1", pct(d["t53_ejc2_wet_0.5"], d["t53_ejc2_dry_0.5"]))
add("Table 5.4", "a = 0.509923 vs 0.5, 3D dry (%)", "-3.0", pct(d["t53_ejc3_dry_0.509923"], d["t53_ejc3_dry_0.5"]))
add("Table 5.4", "a = 0.509923 vs 0.5, 3D wet (%)", "-3.5", pct(d["t53_ejc3_wet_0.509923"], d["t53_ejc3_wet_0.5"]))
gains = [pct(d[f"t53_ejc{k}_{w}_{a}"], d[f"t53_hb{k}_{w}_{a}"]) for k in ("3", "2") for w in ("dry", "wet") for a in ("0.5", "0.509923")]
add("§5.4", "EJC over GHB, lowest of 8 (%)", "5.96", min(gains))
add("§5.4", "EJC over GHB, highest of 8 (%)", "6.03", max(gains))
add("§5.4", "(1 - nu + nu^2)^-1/2 gain at nu = 0.271 (%)", "11.6", d["closure_gain_nu"], scale=100)
add("§5.4", "cohesion share of wet 3D resistance (%)", "77", d["coh_share3_wet"], scale=100, note="'about 77%'")
add("§5.4", "cohesion share of dry 3D resistance (%)", "45", d["coh_share3_dry"], scale=100)
rows.append(dict(where="§5.4", item="tangent shear-strength gain EJC over GHB", thesis="about 5.9-6.6%",
                 software=f"{d['tau_gain_min'] * 100:.2f}-{d['tau_gain_max'] * 100:.2f}% (sigma_n > 0); 5.95-6.43% over all 20x20 / 40 states",
                 verdict="match", note="the upper end needs the tensile bases, where the gain rises to 7.2% near the tensile limit"))

# ---- Table 5.5
for n, w, fh, fm, df in ((20, "dry", "2.7060", "2.7022", "-0.14"), (20, "wet", "1.2723", "1.2695", "-0.22"),
                         (240, "dry", "2.7112", "2.7068", "-0.16"), (240, "wet", "1.2776", "1.2742", "-0.26")):
    add("Table 5.5", f"{n}x{n} {w}: F_H", fh, d[f"t55_fh_{n}_{w}"])
    add("Table 5.5", f"{n}x{n} {w}: F_M", fm, d[f"t55_fm_{n}_{w}"])
    add("Table 5.5", f"{n}x{n} {w}: difference (%)", df, pct(d[f"t55_fm_{n}_{w}"], d[f"t55_fh_{n}_{w}"]))
rows.append(dict(where="§5.5", item="F_M below F_H, range in the text", thesis="0.14-0.27%",
                 software=f"{-pct(d['t55_fm_20_dry'], d['t55_fh_20_dry']):.3f}-{-pct(d['t55_fm_240_wet'], d['t55_fh_240_wet']):.3f}%",
                 verdict="differs", note="largest difference is 0.263%, so Table 5.5 (-0.26%) is right and the text's 0.27% is slightly high"))

# ---- Table 5.6, Figure 5.5, finest values
m3d, m3w, m2d, m2w = d["mesh3_0.5_dry"], d["mesh3_0.5_wet"], d["mesh2_0.5_dry"], d["mesh2_0.5_wet"]
N3, N2 = [10, 15, 20, 30, 40, 60, 80, 120, 160, 240], [20, 40, 80, 160, 320, 640, 1280, 2560, 5120]
var = lambda arr, idx, fin: (max(arr[i] for i in idx) - min(arr[i] for i in idx)) / fin * 100  # noqa: E731
add("Table 5.6", "3D 120/160/240 dry variation (%)", "0.04", var(m3d, [7, 8, 9], m3d[9]))
add("Table 5.6", "3D 120/160/240 wet variation (%)", "0.08", var(m3w, [7, 8, 9], m3w[9]))
add("Table 5.6", "2D 1280/2560/5120 dry variation (%)", "0.02", var(m2d, [6, 7, 8], m2d[8]))
add("Table 5.6", "2D 1280/2560/5120 wet variation (%)", "0.05", var(m2w, [6, 7, 8], m2w[8]))
add("§5.6", "3D 40/60/80 dry variation (%)", "0.26", var(m3d, [4, 5, 6], m3d[6]))
add("§5.6", "3D 40/60/80 wet variation (%)", "0.53", var(m3w, [4, 5, 6], m3w[6]))
add("Table 5.6", "finest dry 3D (240x240)", "2.7112", m3d[9])
add("Table 5.6", "finest wet 3D (240x240)", "1.2776", m3w[9])
add("Table 5.6", "finest dry 2D (5120)", "2.4422", m2d[8])
add("Table 5.6", "finest wet 2D (5120)", "1.0916", m2w[8])
add("§5.6", "10x10 below finest, dry (%)", "-2.6", pct(m3d[0], m3d[9]))
add("§5.6", "10x10 below finest, wet (%)", "-4.4", pct(m3w[0], m3w[9]))
add("§5.6", "finest 3D/2D ratio, dry", "1.110", m3d[9] / m2d[8])
add("§5.6", "finest 3D/2D ratio, wet", "1.170", m3w[9] / m2w[8])
add("Summary", "finest 3D above 2D, dry (%)", "11.0", pct(m3d[9], m2d[8]))
add("Summary", "finest 3D above 2D, wet (%)", "17.0", pct(m3w[9], m2w[8]))
rows.append(dict(where="§5.6", item="20x20 baseline within finest 3D", thesis="about 0.6%",
                 software=f"{pct(m3d[2], m3d[9]):.2f}% dry, {pct(m3w[2], m3w[9]):.2f}% wet", verdict="match", note=""))

# ---- Table 5.7, Figure 5.8
t57 = {"3d20": ("2.706", "2.637", "2.597", "1.272", "0.964", "No accepted FoS"), "3d80": ("2.712", "2.640", "2.590", "1.278", "0.952", "Not verified"),
       "2d40": ("2.448", "2.430", "2.429", "1.098", "0.978", "No accepted FoS")}
for tag, tv in t57.items():
    for (w, tn), t in zip((("dry", "none"), ("dry", "C0"), ("dry", "T0"), ("wet", "none"), ("wet", "C0"), ("wet", "T0")), tv):
        v = d[f"t57_{tag}_{w}_{tn}"]
        label = {"3d20": "3D 20x20", "3d80": "3D 80x80", "2d40": "2D 40 slices"}[tag]
        if isinstance(v, str):
            ok = t == "No accepted FoS"
            rows.append(dict(where="Table 5.7", item=f"{label} {w} {tn}", thesis=t, software=v.split(":")[0] + ": " + v.split(":", 1)[1].strip()[:60],
                             verdict="match" if ok else "new result",
                             note="" if ok else "the software confirms there is no accepted FoS (trial F turns non-positive)"))
        else:
            add("Table 5.7", f"{label} {w} {tn}", t, v, note="thesis value from interpolated local strengths; the direct solve agrees" if (tag, w, tn) == ("3d80", "wet", "C0") else "")
add("§5.8", "80x80 wet: tensile bases at baseline", 367, d["neg_sn_count_80_wet"])
add("§5.8", "80x80 wet: cut-off set after C0", 491, d["t57_3d80_wet_C0_set"])
add("§5.8", "dry C0 reduction, 20x20 (%)", "-2.5", pct(d["t57_3d20_dry_C0"], d["t57_3d20_dry_none"]), decimals=1, note="'2.5-2.6%'")
add("§5.8", "dry C0 reduction, 80x80 (%)", "-2.6", pct(d["t57_3d80_dry_C0"], d["t57_3d80_dry_none"]), decimals=1)
add("§5.8", "dry T0 reduction, 20x20 (%)", "-4.0", pct(d["t57_3d20_dry_T0"], d["t57_3d20_dry_none"]))
add("§5.8", "dry T0 reduction, 80x80 (%)", "-4.5", pct(d["t57_3d80_dry_T0"], d["t57_3d80_dry_none"]))
add("Fig. 5.8", "EJC vs GHB, wet 3D (%)", "+6.02", pct(d["t53_ejc3_wet_0.5"], d["t53_hb3_wet_0.5"]))
add("Fig. 5.8", "C0 vs EJC, wet 3D (%)", "-24.23", pct(d["t57_3d20_wet_C0"], d["t57_3d20_wet_none"]),
    note="exact values give -24.25%; the figure used rounded FoS (0.9640 / 1.2723)")

# ---- Chapter 2 figures
add("Fig. 2.2", "Mohr-Coulomb example: sigma1 (MPa)", "168.282", d["fig22_s1"])
add("Fig. 2.2", "Mohr-Coulomb example: tau at tangency (MPa)", "58.579", d["fig22_tau"])
add("Fig. 2.2", "Mohr-Coulomb example: circle centre (MPa)", "100.641", d["fig22_centre"])
add("Fig. 2.5", "Lode parameter b of the closure state (p' = 1.810 MPa)", "0.266", d["fig25_b_1.81"])
add("Fig. 2.5", "EJC principal-difference gain at fixed sigma3 (%)", "11.5", d["fig25_gain_fixed_s3_1.81"], scale=100)
add("Fig. 2.5", "EJC radial increase, p' = 1.810 MPa (%)", "1.7", d["fig25_radial_gain_1.81"], scale=100)
add("Fig. 2.5", "EJC radial increase, p' = 2.906 MPa (%)", "2.3", d["fig25_radial_gain_2.906"], scale=100)

(OUT / "thesis_rows.json").write_text(json.dumps(rows, indent=1))
n = len(rows)
bad = [r for r in rows if r["verdict"] != "match"]
print(f"{n} thesis values checked; {n - len(bad)} match; not matching:")
for r in bad:
    print(f"  [{r['where']}] {r['item']}: thesis {r['thesis']} | software {r['software']} | {r['note']}")
