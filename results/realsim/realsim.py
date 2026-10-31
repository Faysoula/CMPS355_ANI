"""Real-genome simulations: genome pairs built from REAL chromosomes with a known true ANI.

Genome A = a real chromosome. Genome B = A mutated with a known regional rate pattern (5 kb regions),
plus small insertions/deletions at the same local rate, plus (in one setting) whole regions replaced by
real sequence from another genus ("genes only B has"). True ANI = identity over shared, non-deleted sites.

Every pair is scored by: Mash, sourmash, skani, dnadiff (published tools) and by the project's
estimators, copied unchanged from the notebooks: uniform control, original multi-k gamma,
C (window s) and G (window s + aligned-only).

Usage:  python realsim.py PROJECT_DIR TOOLS_DIR WORKER N_WORKERS [--quick]
  - resumable: finished jobs are skipped; each worker appends to its own CSV in PROJECT_DIR/results/realsim/
  - job j always uses random seed 10000 + j, so results do not depend on how the work is split
"""
import os, sys, json, time, subprocess
import numpy as np, pandas as pd
from scipy.optimize import curve_fit
from scipy.special import gammainc
import pyskani, sourmash

PROJECT, TOOLS, WID, NW = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
QUICK = "--quick" in sys.argv
OUT = f"{PROJECT}/results/realsim"; os.makedirs(OUT, exist_ok=True)
import shutil
MASH_BIN = f"{TOOLS}/mash-Linux64-v2.3/mash"
if not os.path.exists(MASH_BIN): MASH_BIN = shutil.which("mash")      # e.g. on a cluster, Mash installed with conda
assert MASH_BIN and shutil.which("dnadiff"), "Mash or dnadiff not found"

# ---------------- design ----------------
GENOMES = {"E. coli": "GCF_000597845.1",            # 4.76 Mb, GC 51%
           "K. pneumoniae": "GCF_000364385.3",      # 5.44 Mb chromosome, GC 57%
           "S. pneumoniae": "GCF_001095405.2"}      # 2.15 Mb, GC 40%
FOREIGN = "GCF_000635895.2"                         # Bacillus cereus: source of foreign genes
REGION_LEN, WIN = 5_000, 1000
KS = np.arange(15, 32, 2); K_MASH = 21; QSTAR = 0.78
KS_SAVE = np.arange(11, 32, 2)                      # raw counts saved for k = 11..31 (estimators use 15..31)
WIN_KS = (15, 21)                                   # per-1kb-window hit counts saved at these k (for later estimators)
TARGETS = [0.80, 0.84, 0.88, 0.91, 0.93, 0.95, 0.97]
REPS = 10
DNADIFF_REPS = 3                                    # dnadiff (slowest step, ~70% of the time) only on replicates 0-2:
                                                    # it is a reference tool, not an estimator; 378 pairs are plenty
INDEL_RATIO = 0.1                                   # indel events per substitution, at the same local rate

# ---------------- reading / writing ----------------
LOOKUP = np.full(256, 255, dtype=np.uint8)
for i, b in enumerate(b"ACGT"): LOOKUP[b] = i
for i, b in enumerate(b"acgt"): LOOKUP[b] = i
DEC = np.frombuffer(b"ACGT", dtype=np.uint8)

def chromosome(acc):
    contigs, cur = [], []
    for line in open(f"{PROJECT}/data/genomes/{acc}.fna", "rb"):
        if line.startswith(b">"):
            if cur: contigs.append(b"".join(cur)); cur = []
        else: cur.append(line.strip())
    contigs.append(b"".join(cur))
    s = LOOKUP[np.frombuffer(max(contigs, key=len), dtype=np.uint8)]
    return s[s < 4]

def write_fa(path, seq):
    with open(path, "wb") as f: f.write(b">seq\n" + DEC[seq].tobytes() + b"\n")

# ---------------- rate patterns (as in the notebooks, plus the generator passed in) ----------------
def rates_uniform(rng, n): return np.ones(n)
def rates_gamma(alpha): return lambda rng, n: rng.gamma(alpha, 1 / alpha, n)
def rates_two_class(rng, n, frac_fast=0.3, fold=5):
    r = np.where(rng.random(n) < frac_fast, fold, 1.0); return r / r.mean()
def rates_wave(rng, n, fold=4, waves=2):
    x = np.arange(n) / n; amp = (fold - 1) / (fold + 1)
    r = 1 + amp * np.sin(2 * np.pi * waves * x + rng.uniform(0, 2 * np.pi)); return r / r.mean()

SETTINGS = {"Uniform": (rates_uniform, 0.0),
            "Gamma a=2": (rates_gamma(2.0), 0.0),
            "Gamma a=0.5": (rates_gamma(0.5), 0.0),
            "Two-class": (rates_two_class, 0.0),
            "Wave": (rates_wave, 0.0),
            "Gamma a=2 + 20% foreign genes": (rates_gamma(2.0), 0.20)}

def make_pair(rng, A, mean_p, rate_fn, gene_diff, foreign):
    L = len(A); n_reg = -(-L // REGION_LEN)
    reg = np.arange(L) // REGION_LEN
    p_site = np.clip(mean_p * rate_fn(rng, n_reg)[reg], 0, 0.75)
    sub = rng.random(L) < p_site
    B = A.copy(); B[sub] = (A[sub] + rng.integers(1, 4, sub.sum(), dtype=np.uint8)) % 4
    rep = rng.choice(n_reg - 1, int(round(gene_diff * n_reg)), replace=False)       # never the short last region
    is_rep = np.isin(reg, rep)
    for r in rep:
        o = rng.integers(0, len(foreign) - REGION_LEN)
        B[r * REGION_LEN:(r + 1) * REGION_LEN] = foreign[o:o + REGION_LEN]
    ev = np.flatnonzero((rng.random(L) < p_site * INDEL_RATIO) & ~is_rep)
    lens = np.minimum(rng.geometric(1 / 3, ev.size), 20)
    is_ins = rng.random(ev.size) < 0.5
    deleted = np.zeros(L, bool)
    for p, l in zip(ev[~is_ins], lens[~is_ins]): deleted[p:p + l] = True
    pieces, last = [], 0
    for p, l in zip(ev[is_ins], lens[is_ins]):
        pieces.append(B[np.flatnonzero(~deleted[last:p]) + last]); pieces.append(rng.integers(0, 4, l, dtype=np.uint8)); last = p
    pieces.append(B[np.flatnonzero(~deleted[last:]) + last])
    keep = ~is_rep & ~deleted
    rid = 1 - np.bincount(reg, weights=sub & keep, minlength=n_reg) / np.maximum(np.bincount(reg, weights=keep, minlength=n_reg), 1)
    shared_regions = ~np.isin(np.arange(n_reg), rep)
    return np.concatenate(pieces), 1 - sub[keep].mean(), 1 - is_rep.mean(), rid[shared_regions]

# ---------------- k-mers ----------------
def pack(seq, k):
    n = len(seq) - k + 1; v = np.zeros(n, dtype=np.uint64)
    for j in range(k): v = (v << np.uint64(2)) | seq[j:j + n].astype(np.uint64)
    return v
def uniq(v):
    v = np.sort(v); return v[np.concatenate(([True], v[1:] != v[:-1]))]
def member(x, sorted_set):
    i = np.searchsorted(sorted_set, x); i[i == sorted_set.size] = 0
    return sorted_set[i] == x

# ---------------- estimators (copied unchanged from the notebooks) ----------------
def est_mash(C):
    C21 = C[KS == K_MASH][0]; J = C21 / (2 - C21); return 1 + np.log(2 * J / (1 + J)) / K_MASH
def est_uniform(C):
    slope, intercept = np.polyfit(KS, np.log(C), 1); return np.exp(slope), np.exp(intercept), -slope
def gamma_model(k, s, d, alpha): return np.log(s) - alpha * np.log1p(k * d / alpha)
def est_gamma(C):
    y = np.log(C); _, s0, d0 = est_uniform(C); best = None
    for a0 in [0.3, 1, 5, 50]:
        try:
            p, _ = curve_fit(gamma_model, KS, y, p0=[min(s0, 0.99), max(d0, 1e-4), a0],
                             bounds=([1e-6, 1e-6, 0.01], [1, 3, 1000]), maxfev=20000)
            err = np.sum((gamma_model(KS, *p) - y) ** 2)
            if best is None or err < best[1]: best = (p, err)
        except RuntimeError: pass
    s, d, alpha = best[0]; return (1 + d / alpha) ** (-alpha), s, alpha
def fit_gamma_s_fixed(C, s_fixed):
    y = np.log(np.maximum(C, 1e-12)); slope, _ = np.polyfit(KS, y, 1); d0 = max(-slope, 1e-4); best = None
    f = lambda k, d, a: gamma_model(k, s_fixed, d, a)
    for a0 in [0.3, 1, 5, 50]:
        try:
            q, _ = curve_fit(f, KS, y, p0=[d0, a0], bounds=([1e-6, 0.01], [3, 1000]), maxfev=20000)
            sse = np.sum((f(KS, *q) - y) ** 2)
            if best is None or sse < best[1]: best = (q, sse)
        except (RuntimeError, ValueError): pass
    return best[0]
def ani_gamma(d, a): return (1 + d / a) ** (-a)
def ani_gamma_aligned(d, a, qstar=QSTAR):
    lam, scale = -np.log(qstar), d / a
    p_al = gammainc(a, lam / scale)
    return (1 + scale) ** (-a) * gammainc(a, lam * (1 + scale) / scale) / max(p_al, 1e-300)

# ---------------- published tools ----------------
def parse_dnadiff(path):
    ani = al = None; one = False
    for line in open(path):
        if line.startswith("AlignedBases") and al is None: al = float(line.split()[1].split("(")[1].rstrip("%)")) / 100
        if line.startswith("1-to-1"): one = True
        if one and line.startswith("AvgIdentity"): ani = float(line.split()[1]) / 100; break
    return ani, al

def frac_minhash(seq_bytes):
    mh = sourmash.MinHash(n=0, ksize=21, scaled=1000); mh.add_sequence(seq_bytes.decode(), force=True); return mh

# ---------------- main loop ----------------
def main():
    tmp = f"{os.environ.get('TMPDIR', '/tmp')}/realsim_w{WID}"; os.makedirs(tmp, exist_ok=True)   # node-local temp files
    tag = "_quick" if QUICK else ""
    CURVES = f"{OUT}/curves{tag}"; os.makedirs(CURVES, exist_ok=True)   # raw k-mer counts per pair (job_<j>.npz)
    out_csv = f"{OUT}/raw_w{WID}{tag}.csv"
    done = set()                                     # finished jobs from ANY worker file, so the split can change safely
    for f in os.listdir(OUT):
        if f.startswith("raw_w") and f.endswith(f"{tag}.csv") and (QUICK or "_quick" not in f):
            done |= set(pd.read_csv(f"{OUT}/{f}").job)
    foreign = chromosome(FOREIGN)
    # replicate is the OUTER loop: if a run stops early, the finished replicates are complete
    jobs = [(g, s, t, r) for r in range(REPS) for g in GENOMES for s in SETTINGS for t in TARGETS]
    sel = [j for j in range(len(jobs)) if not QUICK or (jobs[j][2] in (0.88, 0.95) and jobs[j][3] == 0)]
    mine = [j for j in sel[WID::NW] if j not in done]
    print(f"worker {WID}/{NW}: {len(mine)} jobs to do ({len(done)} finished overall)", flush=True)
    cur = None
    for n, j in enumerate(mine, 1):
        gname, sname, target, rep = jobs[j]; t0 = time.time()
        if cur is None or cur[0] != gname:
            A = chromosome(GENOMES[gname]); faA = f"{tmp}/A.fa"; write_fa(faA, A)
            subprocess.run([MASH_BIN, "sketch", "-k", "21", "-s", "10000", "-o", f"{tmp}/A", faA], capture_output=True, check=True)
            skdb = pyskani.Database(); skdb.sketch("A", DEC[A].tobytes())
            cur = (gname, A, faA, {k: uniq(pack(A, k)) for k in KS_SAVE}, {k: pack(A, k) for k in WIN_KS}, skdb,
                   frac_minhash(DEC[A].tobytes()))
        _, A, faA, kA, Apos, skdb, smA = cur
        A21 = Apos[21]
        rng = np.random.default_rng(10_000 + j)
        fn, gd = SETTINGS[sname]
        B, true_ani, true_s, rid = make_pair(rng, A, 1 - target, fn, gd, foreign)
        Bbytes = DEC[B].tobytes()
        # shared-k-mer curve (same chance correction as the main notebook)
        C, raw, kBs = [], {"nA": [], "nB": [], "shared": []}, {}
        for k in KS_SAVE:
            kB = uniq(pack(B, k)); sh = int(member(kB, kA[k]).sum())
            raw["nA"].append(kA[k].size); raw["nB"].append(kB.size); raw["shared"].append(sh)
            if k in WIN_KS: kBs[k] = kB
            if k in KS: C.append(max(sh - kA[k].size * kB.size / 4.0 ** k, 1) / kA[k].size)
        C = np.array(C)
        win = {k: np.bincount(np.arange(Apos[k].size) // WIN, weights=member(Apos[k], kBs[k])).astype(np.uint16) for k in WIN_KS}
        sw = (win[21] >= 1).mean()
        np.savez_compressed(f"{CURVES}/job_{j}.npz", ks=KS_SAVE, nA=np.array(raw["nA"]), nB=np.array(raw["nB"]),
                            shared=np.array(raw["shared"]), win15=win[15], win21=win[21],
                            region_identity=rid.astype(np.float32), win_len=WIN, region_len=REGION_LEN)
        g_ani, g_s, g_a = est_gamma(C)
        dC, aC = fit_gamma_s_fixed(C, min(sw, 1))
        # published tools
        faB = f"{tmp}/B.fa"; write_fa(faB, B)
        mash = 1 - float(subprocess.run([MASH_BIN, "dist", f"{tmp}/A.msh", faB], capture_output=True, text=True,
                                        check=True).stdout.split("\t")[2])
        h = skdb.query("B", Bbytes); skani = h[0].identity if h else np.nan
        r = smA.containment_ani(frac_minhash(Bbytes)); smash = r.ani if r.ani is not None else np.nan
        dd = dd_al = np.nan
        if rep < DNADIFF_REPS:
            subprocess.run(["dnadiff", "-p", f"{tmp}/dd", faA, faB], capture_output=True, check=True, cwd=tmp)
            dd, dd_al = parse_dnadiff(f"{tmp}/dd.report")
        row = {"job": j, "genome": gname, "setting": sname, "target": target, "rep": rep, "true": true_ani, "true_s": true_s,
               "Mash (tool)": mash, "sourmash": smash, "skani": skani, "dnadiff": dd, "dnadiff_aligned": dd_al,
               "Mash formula": est_mash(C), "uniform (control)": est_uniform(C)[0], "original": g_ani,
               "C: window s": ani_gamma(dC, aC), "G: window s + aligned-only": ani_gamma_aligned(dC, aC),
               "s_orig": g_s, "alpha_orig": g_a, "s_window": sw, "alpha_C": aC,
               "region_sd": float(rid.std()), "region_min": float(rid.min()), "seconds": round(time.time() - t0, 1)}
        pd.DataFrame([row]).to_csv(out_csv, mode="a", header=not os.path.exists(out_csv), index=False)
        print(f"[w{WID} {n}/{len(mine)}] {gname} | {sname} | {target}: true {true_ani:.4f}  Mash {mash:.4f}  "
              f"skani {skani:.4f}  dnadiff {dd:.4f}  ours {g_ani:.4f}  ({row['seconds']:.0f}s)", flush=True)
    print(f"worker {WID} finished", flush=True)

if __name__ == "__main__":
    main()
