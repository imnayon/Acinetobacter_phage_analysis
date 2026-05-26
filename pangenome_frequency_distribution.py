import csv
from collections import Counter

freq_counter = Counter()

with open("true_protein_frequency.csv") as f:
    reader = csv.DictReader(f)
    for row in reader:
        freq_counter[int(row["Frequency"])] += 1

with open("pangenome_curve_data.csv", "w", newline="") as out:
    writer = csv.writer(out)
    writer.writerow(["Number_of_Phages", "Protein_Count"])
    for freq in sorted(freq_counter):
        writer.writerow([freq, freq_counter[freq]])

print("DONE: pangenome_frequency_data.csv created")