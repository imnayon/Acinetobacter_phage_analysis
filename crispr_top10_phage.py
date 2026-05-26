import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Load CSV with manual column names
df = pd.read_csv(
    "top10_phages.csv",
    header=None,
    names=["Spacer_Match_Count", "Phage_ID", "Full_Phage_Name"]
)

# Sort for better visualization
df = df.sort_values("Spacer_Match_Count", ascending=True)

plt.figure(figsize=(8,6))
plt.barh(df["Full_Phage_Name"], df["Spacer_Match_Count"])

plt.xlabel("Number of Matching CRISPR Spacers")
plt.ylabel("Phage Genome")
plt.title("Top 10 Acinetobacter Phages with CRISPR Spacer Matches")

plt.tight_layout()
plt.savefig("Top10_CRISPR_Phage_Matches.png", dpi=300)
