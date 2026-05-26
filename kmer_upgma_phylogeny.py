#!/usr/bin/env python3

import sys


# ==========================================================
# STEP 1: Read FASTA
# ==========================================================

def read_fasta(file):
    name = None
    seq = []

    with open(file) as f:
        for line in f:
            line = line.strip()

            if line.startswith(">"):
                if name:
                    yield name, ''.join(seq)

                name = line[1:].strip().replace(" ", "_")
                seq = []

            else:
                seq.append(line.upper())

        if name:
            yield name, ''.join(seq)


# ==========================================================
# STEP 2: Generate k-mers
# ==========================================================

def kmers(seq, k=13):
    return set(
        seq[i:i+k]
        for i in range(len(seq)-k+1)
        if 'N' not in seq[i:i+k]
    )


# ==========================================================
# STEP 3: Build distance matrix
# ==========================================================

def build_distance_matrix(fasta_file, k=13):

    print("Reading genomes...")

    seqs = list(read_fasta(fasta_file))

    print(f"Loaded {len(seqs)} genomes")

    kmer_sets = [
        (name, kmers(seq, k))
        for name, seq in seqs
    ]

    names = [name for name, _ in kmer_sets]

    dist = {}

    print("Calculating pairwise distances...")

    for i, (n1, k1) in enumerate(kmer_sets):

        dist[n1] = {}

        for j, (n2, k2) in enumerate(kmer_sets):

            if i == j:
                d = 0.0

            else:
                intersection = len(k1 & k2)
                union = len(k1 | k2)

                d = 1 - (intersection / union)

            dist[n1][n2] = d

    return names, dist


# ==========================================================
# STEP 4: Save distance matrix
# ==========================================================

def save_distance_matrix(names, dist, outfile):

    with open(outfile, "w") as out:

        out.write("\t" + "\t".join(names) + "\n")

        for n1 in names:

            row = [n1]

            for n2 in names:
                row.append(f"{dist[n1][n2]:.6f}")

            out.write("\t".join(row) + "\n")

    print(f"Distance matrix saved: {outfile}")


# ==========================================================
# STEP 5: Build UPGMA tree
# ==========================================================

def build_upgma_tree(names, dist):

    print("Constructing UPGMA tree...")

    clusters = {name: name for name in names}
    sizes = {name: 1 for name in names}

    D = {
        a: {b: dist[a][b] for b in names}
        for a in names
    }

    while len(clusters) > 1:

        keys = list(clusters.keys())

        min_d = float("inf")
        pair = None

        # Find closest pair
        for i in range(len(keys)):
            for j in range(i + 1, len(keys)):

                a = keys[i]
                b = keys[j]

                if D[a][b] < min_d:
                    min_d = D[a][b]
                    pair = (a, b)

        a, b = pair

        # Build Newick node
        newick = (
            f"({clusters[a]}:{min_d/2:.6f},"
            f"{clusters[b]}:{min_d/2:.6f})"
        )

        new_name = a + "|" + b

        clusters[new_name] = newick
        sizes[new_name] = sizes[a] + sizes[b]

        D[new_name] = {}

        # Update average distances
        for k in list(D.keys()):

            if k not in (a, b, new_name):

                d = (
                    D[a][k] * sizes[a]
                    + D[b][k] * sizes[b]
                ) / (sizes[a] + sizes[b])

                D[new_name][k] = d
                D[k][new_name] = d

        D[new_name][new_name] = 0.0

        # Remove merged clusters
        del clusters[a]
        del clusters[b]

        del sizes[a]
        del sizes[b]

        del D[a]
        del D[b]

        for k in D:
            D[k].pop(a, None)
            D[k].pop(b, None)

    final_cluster = list(clusters.keys())[0]
    tree = clusters[final_cluster] + ";"

    return tree


# ==========================================================
# STEP 6: Save tree
# ==========================================================

def save_tree(tree, outfile):

    with open(outfile, "w") as out:
        out.write(tree)

    print(f"Newick tree saved: {outfile}")


# ==========================================================
# MAIN
# ==========================================================

if len(sys.argv) != 2:
    print("Usage:")
    print("python3 build_upgma_tree.py input.fasta")
    sys.exit(1)

fasta = sys.argv[1]

names, dist = build_distance_matrix(fasta, k=13)

save_distance_matrix(
    names,
    dist,
    "distances.tsv"
)

tree = build_upgma_tree(
    names,
    dist
)

save_tree(
    tree,
    "phage_tree_final.nwk"
)

print("Done.")