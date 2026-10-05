# synthetic species-structured collection: gen.py <outdir> <species> <genomes_per_species> <genome_len> <seed>
import numpy as np, sys, os
out, S, G, L, seed = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
rng = np.random.default_rng(seed)
os.makedirs(out, exist_ok=True)
alpha = np.frombuffer(b'ACGT', dtype=np.uint8)
lines = []
for s in range(S):
    base = alpha[rng.integers(0, 4, L)]
    # within-species phylogeny: each genome derives from a random earlier genome of the species
    pool = [base]
    for g in range(G):
        par = pool[rng.integers(0, len(pool))].copy()
        nsnp = rng.poisson(len(par) * 0.0005)
        pos = rng.integers(0, len(par), nsnp)
        par[pos] = alpha[rng.integers(0, 4, nsnp)]
        for _ in range(rng.poisson(3)):   # indels
            p = rng.integers(0, len(par) - 100)
            if rng.random() < 0.5: par = np.delete(par, slice(p, p + rng.integers(1, 50)))
            else: par = np.insert(par, p, alpha[rng.integers(0, 4, rng.integers(1, 50))])
        pool.append(par)
        fn = f'{out}/sp{s}_g{g}.fa'
        with open(fn, 'wb') as f:
            f.write(b'>g\n'); f.write(par.tobytes()); f.write(b'\n')
        lines.append(f'species{s}\t{fn}')
open(f'{out}/list.tsv', 'w').write('\n'.join(lines) + '\n')
