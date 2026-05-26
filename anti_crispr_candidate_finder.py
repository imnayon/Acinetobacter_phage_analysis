import csv
import os

input_file = "phage_database_master.faa"
output_file = "acr_candidate_report.csv"

def get_length(seq):
    return len(seq.replace("\n", "").strip())

candidates = []
current_phage = ""
protein_buffer = []

with open(input_file, 'r') as f:
    content = f.read().split('##')

for section in content:
    if not section.strip(): continue
    lines = section.strip().split('\n')
    phage_name = lines[0].strip()
    
    # Parse proteins in this phage
    proteins = []
    current_prot = None
    
    for line in lines[1:]:
        if line.startswith('>'):
            if current_prot: proteins.append(current_prot)
            header = line[1:].strip()
            current_prot = {'id': header.split()[0], 'desc': header, 'seq': ''}
        elif current_prot:
            current_prot['seq'] += line.strip()
    if current_prot: proteins.append(current_prot)

    # Analyze Neighborhood
    for i, p in enumerate(proteins):
        p_len = len(p['seq'])
        desc_lower = p['desc'].lower()
        
        # Criteria 1: Prokka already found it
        if "acr" in desc_lower or "anti-crispr" in desc_lower:
            candidates.append([phage_name, p['id'], p_len, p['desc'], "DIRECT HIT"])
            continue
            
        # Criteria 2: Small protein near an HTH/Aca regulator
        if 30 <= p_len <= 150:
            # Check neighbors (one up, one down)
            neighbors = []
            if i > 0: neighbors.append(proteins[i-1]['desc'].lower())
            if i < len(proteins)-1: neighbors.append(proteins[i+1]['desc'].lower())
            
            for n in neighbors:
                if any(k in n for k in ["hth", "helix-turn-helix", "aca", "regulator", "binding"]):
                    candidates.append([phage_name, p['id'], p_len, p['desc'], "POTENTIAL (Near Regulator)"])
                    break

# Write to Excel CSV
with open(output_file, 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(["Phage Name", "Protein ID", "Length (AA)", "Prokka Annotation", "Discovery Method"])
    writer.writerows(candidates)

print(f"Success! Analyzed 341 phages. Candidates saved to {output_file}")