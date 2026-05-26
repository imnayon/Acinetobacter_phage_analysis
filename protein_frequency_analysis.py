import csv
from collections import Counter

input_file = "phage_database_master.faa"
output_file = "all_protein_frequency.csv"

# We will count the protein descriptions
protein_counts = Counter()

print("Analyzing 34,525 proteins...")

with open(input_file, 'r') as f:
    for line in f:
        if line.startswith('>'):
            # Prokka format: >ID Description
            # We want to skip the ID and just count the Description
            parts = line.strip().split(maxsplit=1)
            if len(parts) > 1:
                description = parts[1]
                # Remove the unique Prokka ID from the description if it's there
                # (e.g. "PCCACAJL_00003 hypothetical protein" becomes "hypothetical protein")
                clean_desc = " ".join([word for word in description.split() if "_" not in word])
                if not clean_desc: clean_desc = "hypothetical protein"
                protein_counts[clean_desc] += 1

# Write the master spreadsheet
with open(output_file, 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(["Protein Annotation", "Frequency across 341 Phages"])
    for prot, count in protein_counts.most_common():
        writer.writerow([prot, count])

print(f"Master spreadsheet created: {output_file}")