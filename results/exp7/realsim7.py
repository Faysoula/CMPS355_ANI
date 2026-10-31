"""Experiment 7: genome pairs simulated from real chromosomes using REAL mutation-rate profiles.

Each rate profile is the identity of consecutive 5 kb regions along a real genome pair, measured from its
whole-genome alignment (dnadiff), turned into relative rates (mean 1) and laid onto chromosome A in genomic
order, so the real regional structure (neighbouring regions with similar rates) is kept. Everything else is
as in Experiment 6: small indels at the local rate, true ANI over non-deleted sites, same estimators.
Also saved per pair: FracMinHash sketches (scaled 100 and 1000) at every k, to test a fast sketched version.

Usage:  python realsim7.py PROJECT_DIR TOOLS_DIR WORKER N_WORKERS [--quick]     (resumable; seed 20000 + job)
"""
import os, sys, json, time, subprocess, shutil
import numpy as np, pandas as pd
from scipy.optimize import curve_fit
import pyskani, sourmash

PROJECT, TOOLS, WID, NW = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
QUICK = "--quick" in sys.argv
OUT = f"{PROJECT}/results/exp7"; os.makedirs(OUT, exist_ok=True)
MASH_BIN = f"{TOOLS}/mash-Linux64-v2.3/mash"
if not os.path.exists(MASH_BIN): MASH_BIN = shutil.which("mash")

GENOMES = {"E. coli": "GCF_000597845.1", "K. pneumoniae": "GCF_000364385.3", "S. pneumoniae": "GCF_001095405.2"}
REGION_LEN, WIN = 5_000, 1000
KS = np.arange(15, 32, 2); KS_SAVE = np.arange(11, 32, 2); WIN_KS = (15, 21)
TARGETS = [0.85, 0.88, 0.91, 0.94, 0.97]
INDEL_RATIO = 0.1
SCALED = (100, 1000)                                 # FracMinHash: keep k-mers whose hash < 2^64 / scaled

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

def profile_rates(rel, n):
    """real relative rates in genomic order, repeated (tiled) or cut to n regions, mean 1"""
    r = np.resize(np.asarray(rel, float), n); return r / r.mean()

def make_pair(rng, A, mean_p, rel):
    L = len(A); n_reg = -(-L // REGION_LEN); reg = np.arange(L) // REGION_LEN
    off = rng.integers(0, len(rel))                                     # random start point along the profile
    p_site = np.clip(mean_p * profile_rates(np.roll(rel, -off), n_reg)[reg], 0, 0.75)
    sub = rng.random(L) < p_site
    B = A.copy(); B[sub] = (A[sub] + rng.integers(1, 4, sub.sum(), dtype=np.uint8)) % 4
    ev = np.flatnonzero(rng.random(L) < p_site * INDEL_RATIO)
    lens = np.minimum(rng.geometric(1 / 3, ev.size), 20); is_ins = rng.random(ev.size) < 0.5
    deleted = np.zeros(L, bool)
    for p, l in zip(ev[~is_ins], lens[~is_ins]): deleted[p:p + l] = True
    pieces, last = [], 0
    for p, l in zip(ev[is_ins], lens[is_ins]):
        pieces.append(B[np.flatnonzero(~deleted[last:p]) + last]); pieces.append(rng.integers(0, 4, l, dtype=np.uint8)); last = p
    pieces.append(B[np.flatnonzero(~deleted[last:]) + last])
    keep = ~deleted
    return np.concatenate(pieces), 1 - sub[keep].mean()

def pack(seq, k):
    n = len(seq) - k + 1; v = np.zeros(n, dtype=np.uint64)
    for j in range(k): v = (v << np.uint64(2)) | seq[j:j + n].astype(np.uint64)
    return v
def uniq(v):
    v = np.sort(v); return v[np.concatenate(([True], v[1:] != v[:-1]))]
def member(x, s):
    i = np.searchsorted(s, x); i[i == s.size] = 0; return s[i] == x
def mix64(x):                                                           # splitmix64 finaliser: a good 64-bit hash
    with np.errstate(over="ignore"):
        z = x + np.uint64(0x9E3779B97F4A7C15)
        z = (z ^ (z >> np.uint64(30))) * np.uint64(0xBF58476D1CE4E5B9)
        z = (z ^ (z >> np.uint64(27))) * np.uint64(0x94D049BB133111EB)
        return z ^ (z >> np.uint64(31))
def sketch(kset, scaled):
    return kset[mix64(kset) < np.uint64((2 ** 64 - 1) // scaled)]

def gamma_model(k, s, d, alpha): return np.log(s) - alpha * np.log1p(k * d / alpha)
def est_gamma(C):
    y = np.log(C); slope, inter = np.polyfit(KS, y, 1); s0, d0 = np.exp(inter), max(-slope, 1e-4); best = None
    for a0 in [0.3, 1, 5, 50]:
        try:
            p, _ = curve_fit(gamma_model, KS, y, p0=[min(s0, 0.99), d0, a0], bounds=([1e-6, 1e-6, 0.01], [1, 3, 1000]), maxfev=20000)
            err = np.sum((gamma_model(KS, *p) - y) ** 2)
            if best is None or err < best[1]: best = (p, err)
        except RuntimeError: pass
    s, d, a = best[0]; return (1 + d / a) ** (-a), s, a

def frac_minhash(b):
    mh = sourmash.MinHash(n=0, ksize=21, scaled=1000); mh.add_sequence(b.decode(), force=True); return mh

def main():
    prof = pd.read_csv(f"{OUT}/real_profiles.csv")
    rels = {r.profile_id: np.array(json.loads(r.relative_rates)) for r in prof.itertuples()}
    tmp = f"{os.environ.get('TMPDIR', '/tmp')}/exp7_w{WID}"; os.makedirs(tmp, exist_ok=True)
    tag = "_quick" if QUICK else ""
    CURVES = f"{OUT}/curves{tag}"; os.makedirs(CURVES, exist_ok=True)
    out_csv = f"{OUT}/raw_w{WID}{tag}.csv"
    done = set()
    for f in os.listdir(OUT):
        if f.startswith("raw_w") and f.endswith(f"{tag}.csv") and (QUICK or "_quick" not in f):
            done |= set(pd.read_csv(f"{OUT}/{f}").job)
    jobs = [(g, pid, t) for g in GENOMES for pid in sorted(rels) for t in TARGETS]
    sel = [j for j in range(len(jobs)) if not QUICK or (jobs[j][1] in sorted(rels)[:2] and jobs[j][2] in (0.88, 0.94))]
    mine = [j for j in sel[WID::NW] if j not in done]
    print(f"worker {WID}/{NW}: {len(mine)} jobs to do ({len(done)} finished overall)", flush=True)
    cur = None
    for n, j in enumerate(mine, 1):
        gname, pid, target = jobs[j]; t0 = time.time()
        if cur is None or cur[0] != gname:
            A = chromosome(GENOMES[gname]); faA = f"{tmp}/A.fa"; write_fa(faA, A)
            subprocess.run([MASH_BIN, "sketch", "-k", "21", "-s", "10000", "-o", f"{tmp}/A", faA], capture_output=True, check=True)
            skdb = pyskani.Database(); skdb.sketch("A", DEC[A].tobytes())
            kA = {k: uniq(pack(A, k)) for k in KS_SAVE}
            cur = (gname, A, faA, kA, {k: pack(A, k) for k in WIN_KS}, skdb, frac_minhash(DEC[A].tobytes()),
                   {(k, sc): sketch(kA[k], sc) for k in KS for sc in SCALED})
        _, A, faA, kA, Apos, skdb, smA, skA = cur
        rng = np.random.default_rng(20_000 + j)
        B, true_ani = make_pair(rng, A, 1 - target, rels[pid]); Bb = DEC[B].tobytes()
        raw = {"nA": [], "nB": [], "shared": []}; kBs = {}; C = []
        sk = {sc: {"nA": [], "nB": [], "shared": []} for sc in SCALED}
        for k in KS_SAVE:
            kB = uniq(pack(B, k)); sh = int(member(kB, kA[k]).sum())
            raw["nA"].append(kA[k].size); raw["nB"].append(kB.size); raw["shared"].append(sh)
            if k in WIN_KS: kBs[k] = kB
            if k in KS:
                C.append(max(sh - kA[k].size * kB.size / 4.0 ** k, 1) / kA[k].size)
                for sc in SCALED:
                    sa, sb = skA[(k, sc)], sketch(kB, sc)
                    sk[sc]["nA"].append(sa.size); sk[sc]["nB"].append(sb.size); sk[sc]["shared"].append(int(member(sb, sa).sum()))
        win = {k: np.bincount(np.arange(Apos[k].size) // WIN, weights=member(Apos[k], kBs[k])).astype(np.uint16) for k in WIN_KS}
        np.savez_compressed(f"{CURVES}/job_{j}.npz", ks=KS_SAVE, nA=np.array(raw["nA"]), nB=np.array(raw["nB"]),
                            shared=np.array(raw["shared"]), win15=win[15], win21=win[21], sketch_ks=KS,
                            **{f"s{sc}_{x}": np.array(sk[sc][x]) for sc in SCALED for x in ("nA", "nB", "shared")})
        g_ani, g_s, g_a = est_gamma(np.array(C))
        faB = f"{tmp}/B.fa"; write_fa(faB, B)
        mash = 1 - float(subprocess.run([MASH_BIN, "dist", f"{tmp}/A.msh", faB], capture_output=True, text=True, check=True).stdout.split("\t")[2])
        h = skdb.query("B", Bb); skani = h[0].identity if h else np.nan
        r = smA.containment_ani(frac_minhash(Bb)); smash = r.ani if r.ani is not None else np.nan
        row = {"job": j, "genome": gname, "profile_id": pid, "target": target, "true": true_ani,
               "Mash (tool)": mash, "sourmash": smash, "skani": skani, "original": g_ani, "s_orig": g_s, "alpha_orig": g_a,
               "seconds": round(time.time() - t0, 1)}
        pd.DataFrame([row]).to_csv(out_csv, mode="a", header=not os.path.exists(out_csv), index=False)
        print(f"[w{WID} {n}/{len(mine)}] {gname} | {pid} | {target}: true {true_ani:.4f}  Mash {mash:.4f}  skani {skani:.4f}  "
              f"ours {g_ani:.4f}  ({row['seconds']:.0f}s)", flush=True)
    print(f"worker {WID} finished", flush=True)

if __name__ == "__main__":
    main()
