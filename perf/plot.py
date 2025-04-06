import matplotlib.pyplot as plt
import json
import sys
import os
from math import log2, ceil, sqrt

# Check for command-line argument
if len(sys.argv) < 2:
    print("Usage: python plot_query.py <file.json>")
    sys.exit(1)

# Read the JSON file path from command-line
json_file = sys.argv[1]

# Load JSON data from each line
datasets = []
with open(json_file, 'r') as f:
    for line in f:
        if line.strip():  # Skip empty lines
            datasets.append(json.loads(line))

num_plots = len(datasets)
cols = 2 # ceil(sqrt(num_plots))
rows = ceil(num_plots / cols)

# Set up subplots
fig, axs = plt.subplots(rows, cols, figsize=(cols * 6, rows * 4))
axs = axs.flatten()

# Helper function for powers of 2 labels
def format_as_power_of_two(x):
    return f"$2^{{{int(log2(x))}}}$"

# Generate each subplot
for i, data in enumerate(datasets):
    ax = axs[i]
    x = data["sequence_lengths"]
    y = data["avg_ns_per_query"]

    ax.plot(x, y, marker='o', linestyle='-', color='blue')
    ax.set_xscale('log', base=2)
    ax.set_xlabel("sequence length")
    ax.set_ylabel("avg. ns/query")
    ax.set_title(data.get("query", f"Query {i+1}"))
    ax.grid(True, which="both", linestyle='--', linewidth=0.5)
    ax.set_xticks(x)
    ax.set_xticklabels([format_as_power_of_two(val) for val in x], rotation=0)

# Hide any unused subplots
for j in range(i + 1, len(axs)):
    fig.delaxes(axs[j])

plt.tight_layout()

# Save as PDF
pdf_file = os.path.splitext(json_file)[0] + ".pdf"
plt.savefig(pdf_file)
print(f"Plot saved to {pdf_file}")