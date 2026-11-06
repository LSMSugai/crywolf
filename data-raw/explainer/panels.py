"""The panels of the article “How crywolf works”, as SVG strings."""
import math, random
from svg_helpers import T, R, C, L, P, svg, arrow, wave, note, NOTE, MONO

SM = 'style="font-size:10px"'

def legend(b, x, y, items):
    for cls, lab in items:
        b.append(R(x, y - 9, 12, 12, cls, 2)); b.append(T(x + 17, y + 1, lab, "ts", "start", SM))
        x += 40 + len(lab) * 6.2

# ---------------- 1 · the problem
def problem():
    b = []; Wd = 900
    random.seed(4)
    sites = [(1, [1, 0, 1, 1, 0, 1, 0, 1, 1, 0]), (1, [0, 1, 0, 0, 1, 0, 0, 1, 0, 0]),
             (1, [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]), (0, [0, 0, 2, 0, 0, 0, 0, 2, 0, 0]),
             (0, [0, 0, 0, 0, 0, 0, 2, 0, 0, 0]), (0, [0, 0, 0, 0, 0, 0, 0, 0, 0, 0])]
    b.append(T(0, 16, "site", "tc tm")); b.append(T(64, 16, "truth", "tc tm", "middle"))
    b.append(T(110, 16, "what the classifier says, night by night (occasions)", "tc tm"))
    b.append(T(560, 16, "naive reading", "tc tm")); b.append(T(720, 9, "crywolf: P(occupied),", "tc f4t")); b.append(T(720, 21, "from the occupancy model", "tc f4t"))
    for i, (z, ys) in enumerate(sites):
        y = 30 + i * 30
        b.append(T(0, y + 14, f"S{i+1}", "tc tink"))
        if z: b.append(C(64, y + 10, 8, "f3"))
        else: b.append(C(64, y + 10, 8, "fsurf st3 w15"))
        for k, v in enumerate(ys):
            cls = {0: "fline", 1: "f3", 2: "fbad"}[v]
            b.append(R(110 + k * 40, y, 34, 20, cls, 3))
        det = any(ys)
        b.append(T(560, y + 14, "occupied" if det else "empty", "ts " + ("tink" if det else "tm")))
        wrong = det != bool(z)
        if wrong: b.append(T(640, y + 14, "✗ false alarm" if det else "✗ missed", "ts fbadt"))
        n1 = sum(1 for v in ys if v); n0 = len(ys) - n1
        lo = 0.5 * 0.6 ** n1 * 0.4 ** n0; le = 0.5 * 0.05 ** n1 * 0.95 ** n0
        b.append(T(720, y + 14, f"{lo / (lo + le):.2f}", "tc f4t"))
    legend(b, 110, 222, [("f3", "true detection"), ("fbad", "false positive: the classifier cried wolf"), ("fline", "nothing detected")])
    b.append(C(116, 240, 6, "f3")); b.append(T(126, 244, "present", "ts", "start", SM))
    b.append(C(186, 240, 6, "fsurf st3 w15")); b.append(T(196, 244, "absent (the truth, unknown to you)", "ts", "start", SM))
    b.append(L(0, 258, Wd, 258, "sgrid"))
    note(b, 0, 276, "Counting a site as occupied after one detection turns S4 and S5 into occupied sites (ψ = 4/6 instead of 3/6), and misses S3.")
    note(b, 0, 290, "Ordinary occupancy models fix misses by repeating nights, but assume every detection is real. crywolf also lets detections be false,")
    note(b, 0, 304, "and gives each site a probability of being occupied (right column: computed with ψ = 0.5, p11 = 0.6, p10 = 0.05; see panel 2b).")
    return svg(Wd, 312, "".join(b), "Classifier detections over nights at six sites, with true and false detections")

# ---------------- 2 · the model
def model():
    b = []; W = 440
    b.append(R(170, 6, 100, 30, "fs4 st4", 8)); b.append(T(220, 26, "a site", "ts tink", "middle"))
    b.append(P("M220 36 V50 M100 50 H340 M100 50 V64 M340 50 V64", "nofill sline"))
    b.append(T(150, 46, "ψ", "tb f4t", "middle")); b.append(T(296, 46, "1 − ψ", "tb f4t", "middle"))
    b.append(R(30, 66, 140, 30, "fs3 st3", 8)); b.append(T(100, 86, "occupied (z = 1)", "ts tink", "middle"))
    b.append(R(270, 66, 140, 30, "fsurf st3 w15", 8)); b.append(T(340, 86, "empty (z = 0)", "ts tink", "middle"))
    b.append(arrow(100, 98, 100, 118)); b.append(arrow(340, 98, 340, 118))
    b.append(T(100, 134, "each night: detection", "ts", "middle")); b.append(T(100, 150, "with p11 (e.g. 0.6)", "tb f3t", "middle"))
    b.append(T(340, 134, "each night: detection", "ts", "middle")); b.append(T(340, 150, "with p10 (e.g. 0.05)", "tb fbadt", "middle"))
    for k in range(8):
        b.append(R(30 + k * 18, 162, 14, 12, "f3" if k in (0, 2, 3, 5, 7) else "fline", 2))
        b.append(R(270 + k * 18, 162, 14, 12, "fbad" if k == 4 else "fline", 2))
    b.append(L(0, 190, W, 190, "sgrid"))
    note(b, 0, 208, "z is hidden: the model weighs both stories for every site")
    note(b, 0, 222, "y (night) ~ Bernoulli(z · p11 + (1 − z) · p10)")
    note(b, 0, 236, "p11 = true-detection rate, p10 = false-positive rate (crying wolf)")
    note(b, 0, 250, "constraint: p11 > p10 at an average site (the gap Δ), so the")
    note(b, 0, 264, "two states cannot swap labels")
    return svg(W, 272, "".join(b), "The false-positive occupancy model as two branches")

# ---------------- 3 · how it tells them apart
def mixture():
    b = []; W = 440
    x0, x1, y0, y1 = 36, 430, 190, 40
    n = 20
    def binom(k, p): return math.comb(n, k) * p ** k * (1 - p) ** (n - k)
    occ = [60 * binom(k, .6) for k in range(n + 1)]
    emp = [60 * binom(k, .05) for k in range(n + 1)]
    mx = max(max(occ), max(emp))
    bw = (x1 - x0) / (n + 1)
    for k in range(n + 1):
        x = x0 + k * bw
        he = emp[k] / mx * (y0 - y1); ho = occ[k] / mx * (y0 - y1)
        b.append(R(x + 1, y0 - he, bw - 2, he, "fbad", 1, 'opacity=".75"'))
        b.append(R(x + 1, y0 - ho, bw - 2, ho, "f3", 1, 'opacity=".75"'))
    b.append(L(x0, y0, x1, y0, "sgrid"))
    for k in (0, 5, 10, 15, 20): b.append(T(x0 + (k + .5) * bw, y0 + 13, str(k), "tc", "middle"))
    b.append(T((x0 + x1) / 2, y0 + 28, "nights with a detection, out of 20", "ts", "middle"))
    b.append(T(x0, y1 - 10, "number of sites", "ts"))
    b.append(T(x0 + 3 * bw, y1 + 6, "empty sites:", "ts fbadt")); b.append(T(x0 + 3 * bw, y1 + 19, "a few false alarms", "ts fbadt"))
    b.append(T(x0 + 14.5 * bw, y1 + 6, "occupied sites:", "ts f3t")); b.append(T(x0 + 14.5 * bw, y1 + 19, "detected most nights", "ts f3t"))
    b.append(R(x0 + 0.5 * bw, y0 - 128, 4.5 * bw, 132, "nofill st4 w15 dash", 4))
    b.append(T(x0 + 5.4 * bw, y0 - 112, "1–4 detections: a rare caller", "tny f4t")); b.append(T(x0 + 5.4 * bw, y0 - 100, "or a false alarm?", "tny f4t"))
    b.append(L(0, 230, W, 230, "sgrid"))
    note(b, 0, 248, "With repeated nights the two groups form two humps; the model")
    note(b, 0, 262, "estimates ψ, p11 and p10 as the mixture that fits them best. Where")
    note(b, 0, 276, "the humps overlap, the data alone cannot say how many detections")
    note(b, 0, 290, "are false: that is what the known truth (next panels) settles")
    return svg(W, 298, "".join(b), "Distribution of detections per site: empty and occupied sites")

# ---------------- 4 · the known truth: two designs
def designs():
    b = []; Wd = 900
    def site_row(x, y, cells, extra=None):
        for k, v in enumerate(cells):
            b.append(R(x + k * 22, y, 18, 14, {0: "fline", 1: "f3", 2: "fbad"}[v], 2))
    # left: multiple methods
    b.append(R(0, 0, 436, 290, "fs2", 10, 'opacity=".45"'))
    b.append(T(14, 22, "Multiple methods", "tb f2t")); b.append(T(140, 22, "Miller et al. 2011 · site confirmation", "ts"))
    note(b, 14, 42, "an expert surveys some sites by a second method that cannot")
    note(b, 14, 56, "cry wolf: whole recordings annotated by ear, at random")
    b.append(T(14, 82, "automated", "tc tm")); site_row(90, 72, [0, 1, 0, 0, 2, 0, 1, 0])
    b.append(T(14, 106, "verified", "tc f2t"))
    for k, v in enumerate([0, 1, 0]):
        b.append(R(90 + k * 44, 96, 40, 14, "f2" if v else "fsurf st2 w15", 2))
    b.append(T(232, 107, "a heard call ⇒ occupied, for sure", "ts f2t"))
    note(b, 14, 140, "verified detections at occupied and empty sites anchor p10:")
    note(b, 14, 154, "the model learns how often the classifier fires at sites known")
    note(b, 14, 168, "to be occupied and at sites where the expert heard nothing")
    b.append(L(14, 184, 422, 184, "sgrid"))
    b.append(T(14, 204, "in crywolf", "tc f2t"))
    note(b, 14, 222, "fp_data(verified = ): one row per annotated")
    note(b, 14, 236, "recording, with its site and y (1 if the species")
    note(b, 14, 250, "was heard)")
    note(b, 14, 274, "independent of the classifier: no score, no threshold")
    # right: observation confirmation
    x0 = 464
    b.append(R(x0, 0, 436, 290, "fs4", 10, 'opacity=".55"'))
    b.append(T(x0 + 14, 22, "Observation confirmation", "tb f4t")); b.append(T(x0 + 196, 22, "Chambert et al. 2015", "ts"))
    note(b, x0 + 14, 42, "an expert checks the detections of some occasions,")
    note(b, x0 + 14, 56, "drawn at random, and labels each one true or false")
    states = [("0", "nothing flagged", "fline", None), ("1", "flagged, all false", "fbad", None),
              ("2", "flagged, all true", "f3", None), ("3", "true and false", "f3", "fbad")]
    for k, (st, lab, c1, c2) in enumerate(states):
        x = x0 + 14 + k * 104
        b.append(R(x, 70, 30, 22, c1, 3))
        if c2: b.append(R(x + 15, 70, 15, 22, c2, 0))
        b.append(T(x + 15, 85, st, "ts finv" if c1 != "fline" else "ts tink", "middle"))
        b.append(T(x, 108, lab, "ts", "start", 'style="font-size:9.6px"'))
    note(b, x0 + 14, 140, "false detections are seen directly, at occupied sites too,")
    note(b, x0 + 14, 154, "so p10 is estimated even without a separate survey;")
    note(b, x0 + 14, 168, "fp_data(confirmed = c(\"n_true\", \"n_false\"))")
    b.append(L(x0 + 14, 184, x0 + 422, 184, "sgrid"))
    b.append(T(x0 + 14, 204, "in crywolf", "tc f4t"))
    note(b, x0 + 14, 222, "reviewed clips (fp_read_validation()), drawn at the")
    note(b, x0 + 14, 236, "threshold of the analysis → fp_aggregate(validation = )")
    note(b, x0 + 14, 250, "→ n_true, n_false → fp_data(confirmed = )")
    b.append(T(x0 + 14, 274, "clips picked for high scores understate p10", "ts fbadt", "start", 'style="font-size:10px"'))
    b.append(T(450, 312, "the two can be combined; without either, p10 rests only on the priors and the gap (fp_occu() warns)", "ts", "middle"))
    return svg(Wd, 320, "".join(b), "Two designs that give the model known truth")

# ---------------- 5 · where the threshold enters
def threshold():
    b = []; W = 440
    note(b, 0, 14, "fp_aggregate(): an occasion is a detection when any of its")
    note(b, 0, 28, "predictions scores at or above the threshold")
    occs = [[.35, .52, .81], [.22, .4], [.67], [.12, .3, .45, .58], []]
    t = 0.5
    x0, x1, y0, y1 = 40, 420, 160, 66
    Y = lambda s: y0 - s * (y0 - y1)
    b.append(L(x0, Y(t), x1, Y(t), "st4 w15 dash")); b.append(T(x1, y1 - 12, "- - threshold 0.5", "tc f4t", "end"))
    for v in (0, 1): b.append(T(x0 - 6, Y(v) + 4, str(v), "tc", "end"))
    b.append(T(x0 - 30, y1 - 12, "score", "ts"))
    for i, ss in enumerate(occs):
        x = x0 + 20 + i * 76
        b.append(R(x - 18, y1 - 4, 36, y0 - y1 + 8, "fsurf sline", 4))
        for j, s in enumerate(ss): b.append(C(x - 8 + j * 6, Y(s), 3.2, "f4" if s >= t else "fmuted"))
        det = any(s >= t for s in ss)
        b.append(R(x - 10, y0 + 10, 20, 14, "f4" if det else "fline", 2))
        b.append(T(x, y0 + 38, "y = 1" if det else "y = 0", "tc " + ("f4t" if det else "tm"), "middle"))
    b.append(T(x0 - 30, y0 + 21, "occasion", "tny tm"))
    b.append(L(0, 214, W, 214, "sgrid"))
    note(b, 0, 232, "a lower threshold lets in more true and more false detections;")
    note(b, 0, 246, "p10 is estimated at that threshold and absorbs the extra false")
    note(b, 0, 260, "ones. So the confirmed clips must be drawn at the same threshold.")
    note(b, 0, 274, "fp_threshold_sweep() refits at several thresholds to check")
    return svg(W, 282, "".join(b), "Predictions become occasions at a threshold")

# ---------------- 6 · covariates
def covs():
    b = []; W = 440
    note(b, 0, 14, "three parts, each with its own covariates:")
    rows = [("ψ", "f4t", "occ = ~ x", "occupancy: habitat, gradient (the ecology)"),
            ("p11", "f3t", "det = ~ noise", "true detection: noise, recorder, night"),
            ("p10", "fbadt", "fp = ~ w", "false positives: what makes it cry wolf")]
    for k, (sym, c, f, lab) in enumerate(rows):
        y = 36 + k * 22
        b.append(T(0, y, sym, "tb " + c)); b.append(T(40, y, f, "tc tink")); b.append(T(150, y, lab, "ts", "start", SM))
    x0, x1, y0, y1 = 40, 420, 230, 112
    b.append(L(x0, y0, x1, y0, "sgrid")); b.append(L(x0, y0, x0, y1, "sgrid"))
    X = lambda v: x0 + (v + 2) / 4 * (x1 - x0); Yp = lambda p: y0 - p * (y0 - y1)
    xs = [-2 + i * .1 for i in range(41)]
    lg = lambda z: 1 / (1 + math.exp(-z))
    b.append(P("M" + " L".join(f"{X(v):.1f} {Yp(lg(-1.6 * v)):.1f}" for v in xs), "nofill st4 w2"))
    b.append(P("M" + " L".join(f"{X(v):.1f} {Yp(lg(.4 - .5 * v)):.1f}" for v in xs), "nofill stbad w15 dash"))
    b.append(T(X(-1.9), Yp(.97), "true ψ along x", "tny f4t")); b.append(T(X(-1.9), Yp(.3), "model with constant p10, while false alarms", "tny fbadt"))
    b.append(T(X(-1.9), Yp(.2), "rise where the species is rare: slope flattened", "tny fbadt"))
    b.append(T(x1, y0 + 13, "x (e.g. forest cover) →", "tc", "end")); b.append(T(x0 - 6, y1 + 4, "ψ", "tc", "end"))
    b.append(L(0, 250, W, 250, "sgrid"))
    note(b, 0, 268, "if false positives follow a gradient tied to occupancy (a soundalike")
    note(b, 0, 282, "species where the target is rare), let p10 depend on it (fp = ~ x)")
    return svg(W, 290, "".join(b), "Covariates of occupancy, true detection and false positives")

# ---------------- 7 · from detections to results (wide)
def flow():
    b = []; Wd = 900
    def box(x, y, w, h, cls, lines, tc="tink"):
        b.append(R(x, y, w, h, cls, 8))
        for k, (t, c) in enumerate(lines): b.append(T(x + 10, y + 18 + k * 14, t, c))
    box(0, 30, 172, 64, "fsurf sline", [("BirdNET output", "ts tink"), ("fp_read_birdnet()", "tc tink"), ("+ fp_read_validation()", "tc tink")])
    box(0, 110, 172, 64, "fs1 st1", [("a sheepdog export", "ts f1t"), ("fp_read_sheepdog()", "tc tink"), ("occasions, scores, reviews", "ts")])
    b.append(arrow(174, 62, 196, 92)); b.append(arrow(174, 142, 196, 112))
    box(198, 70, 140, 64, "fs4 st4", [("occasions", "ts f4t"), ("with scores and", "ts tink"), ("known truth", "ts tink")])
    b.append(P("M340 102 H356 V40 H372 M356 102 V166 H372", "nofill sline"))
    b.append(arrow(360, 40, 372, 40)); b.append(arrow(360, 166, 372, 166))
    box(374, 10, 226, 60, "fsurf sline", [("one score threshold", "ts tink"), ("fp_aggregate(threshold = )", "tc f4t"), ("→ fp_data() → fp_occu()", "tc f4t")])
    box(374, 136, 226, 60, "fsurf sline", [("one version per draw", "ts tink"), ("fp_draws_data(rule = )", "tc f4t"), ("→ fp_occu_mc(): one fit per draw", "tc f4t")])
    b.append(arrow(602, 40, 630, 92)); b.append(arrow(602, 166, 630, 112))
    box(632, 40, 268, 126, "fs4 st4", [("results", "ts f4t"), ("ψ and its slopes (coef, fp_effects)", "ts tink"), ("p11, p10 and their covariates", "ts tink"),
                                         ("P(occupied) of each site (fp_sites)", "ts tink"), ("checks: fp_ppc(), loo(), fp_calibrate()", "ts tink"),
                                         ("with draws: sd_within + sd_between", "ts tink")])
    b.append(T(487, 92, "any score threshold", "tny tm", "middle")); b.append(T(487, 214, "calibration draws from sheepdog (brms)", "tny tm", "middle"))
    b.append(L(0, 232, Wd, 232, "sgrid"))
    note(b, 0, 250, "Either way crywolf estimates p10 from the known truth. Calibration draws only add the uncertainty about which")
    note(b, 0, 264, "detections are true; the occupancy fit does not feed back into the calibration (a “cut” model).")
    return svg(Wd, 272, "".join(b), "From detections to the crywolf results")


# ---------------- 2b · how P(occupied) is computed for one site (wide)
def psite():
    b = []; Wd = 900
    psi, p11, p10 = .5, .6, .05
    note(b, 0, 14, "two stories for the same nights; P(occupied) is the share of the evidence that the “occupied” story gets (Bayes’ rule)")
    # general formulas, with what each term means
    boxes = [(0, "fs3 st3", "f3t", "Occupied  =  ψ  ×  p11^n₁  ×  (1 − p11)^n₀",
              ["ψ: chance that the site is occupied at all", "p11^n₁: each detected night is a true detection (p11 each)", "(1 − p11)^n₀: each silent night is a miss (1 − p11 each)"]),
             (456, "fsbad stbad", "fbadt", "Empty  =  (1 − ψ)  ×  p10^n₁  ×  (1 − p10)^n₀",
              ["1 − ψ: chance that the site is empty", "p10^n₁: each detected night is a false alarm (p10 each)", "(1 − p10)^n₀: each silent night has no false alarm"])]
    for x, cls, tc, f, ls in boxes:
        b.append(R(x, 28, 444, 92, cls, 8))
        b.append(T(x + 14, 50, f, "tc " + tc, "start", 'style="font-size:12px"'))
        for k, l in enumerate(ls): note(b, x + 14, 70 + k * 15, l)
    b.append(R(170, 132, 560, 30, "fs4 st4", 8))
    b.append(T(450, 152, "P(occupied)  =  Occupied  ÷  (Occupied + Empty)", "tc f4t", "middle", 'style="font-size:12.5px"'))
    note(b, 0, 182, "n₁ = nights with a detection, n₀ = nights without; the nights are multiplied because each is an independent yes/no outcome (Bernoulli).")
    note(b, 0, 196, "The two stories are the only possibilities, so their shares add up to 1. Each number is the chance of one exact sequence of nights: only their ratio matters.")
    cases = [("S2", [0, 1, 0, 0, 1, 0, 0, 1, 0, 0]), ("S4", [0, 0, 1, 0, 0, 0, 0, 1, 0, 0])]
    for i, (nm, ys) in enumerate(cases):
        y0 = 220 + i * 108
        n1 = sum(ys); n0 = len(ys) - n1
        lo = psi * p11 ** n1 * (1 - p11) ** n0; le = (1 - psi) * p10 ** n1 * (1 - p10) ** n0
        P_ = lo / (lo + le)
        b.append(T(0, y0 + 13, nm, "tc tink"))
        for k, v in enumerate(ys): b.append(R(30 + k * 18, y0 + 2, 14, 14, "f4" if v else "fline", 2))
        b.append(T(30, y0 + 32, f"n₁ = {n1}, n₀ = {n0}", "ts", "start", SM))
        b.append(T(230, y0 + 13, "occupied:", "ts f3t")); b.append(T(300, y0 + 13, f"0.5 × 0.6^{n1} × 0.4^{n0} = {lo:.6f}", "tc tink", "start", MONO))
        b.append(T(230, y0 + 32, "empty:", "ts fbadt")); b.append(T(300, y0 + 32, f"0.5 × 0.05^{n1} × 0.95^{n0} = {le:.6f}", "tc tink", "start", MONO))
        b.append(T(560, y0 + 22, f"→ {lo:.6f} ÷ ({lo:.6f} + {le:.6f})", "tc f4t", "start", MONO))
        bw = 520; x0 = 230
        b.append(R(x0, y0 + 44, bw * P_, 16, "f3", 2)); b.append(R(x0 + bw * P_, y0 + 44, bw * (1 - P_), 16, "fbad", 2))
        b.append(T(x0 + bw + 10, y0 + 57, f"P(occupied) = {P_:.2f}", "ts f4t"))
        note(b, x0, y0 + 78, ["three false alarms in ten nights would be rare at p10 = 0.05: “occupied” wins",
                               "an occupied site would rarely be this quiet (8 misses), two false alarms are plausible: “empty” wins"][i])
    b.append(L(0, 442, Wd, 442, "sgrid"))
    note(b, 0, 460, "In practice ψ, p11 and p10 are not known: crywolf estimates them from all the sites and the known truth, repeats this for every")
    note(b, 0, 474, "posterior draw and averages (fp_sites()). Covariates give each site its own ψ and p10. A verified or confirmed true detection sets P = 1.")
    note(b, 0, 488, "Values here: ψ = 0.5, p11 = 0.6, p10 = 0.05.")
    return svg(Wd, 496, "".join(b), "How P(occupied) is computed for one site")
