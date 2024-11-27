import matplotlib.pyplot as plt

data = [
    {"Gene Name": "shuttlecraft", "Start": 3983998, "Stop": 3999923, "Strand": "+", "Color": "gray"},
    {"Gene Name": "HEAT2", "Start": 4268556, "Stop": 4285716, "Strand": "+", "Color": "gray"},
    {"Gene Name": "putative ivory", "Start": 4338269, "Stop": 4565668, "Strand": "-", "Color": "red"},
    {"Gene Name": "cortex", "Start": 4574204, "Stop": 4589927, "Strand": "+", "Color": "green"},
    {"Gene Name": "GPC", "Start": 4697935, "Stop": 4720018, "Strand": "-", "Color": "gray"},
    {"Gene Name": "Rab-ggt", "Start": 4965639, "Stop": 4969451, "Strand": "-", "Color": "gray"}
]

# Normalize the start and stop positions
min_pos = min(gene["Start"] for gene in data)
max_pos = max(gene["Stop"] for gene in data)
for gene in data:
    gene["Normalized Start"] = (gene["Start"] - min_pos) / (max_pos - min_pos)
    gene["Normalized Stop"] = (gene["Stop"] - min_pos) / (max_pos - min_pos)

# Create a plot
fig, ax = plt.subplots(figsize=(12, 6))
ax.set_xlim(-0.05, 1.05)
ax.set_ylim(-1, 1)
ax.set_xlabel("Relative Position")
ax.set_yticks([])

# Draw arrows for each gene
for gene in data:
    start = gene["Normalized Start"]
    stop = gene["Normalized Stop"]
    color = gene["Color"]
    gene_name = gene["Gene Name"]
    if gene_name == "Rab-ggt":
        gene_name = r"Rab-ggt$\beta$"
    if gene["Strand"] == "+":
        ax.annotate(
            "",
            xy=(stop, 0),
            xytext=(start, 0),
            arrowprops=dict(facecolor=color, edgecolor=color, arrowstyle="->", lw=1.5)
        )
    else:
        ax.annotate(
            "",
            xy=(stop, 0),
            xytext=(start, 0),
            arrowprops=dict(facecolor=color, edgecolor=color, arrowstyle="<-", lw=1.5)
        )
    ax.text((start + stop) / 2, 0.05, gene_name, horizontalalignment='center', verticalalignment='bottom')

# Add scale-bar
kb_in_range = 100000
scale_length_normalized = kb_in_range / (max_pos - min_pos)
ax.plot([0.05, 0.05 + scale_length_normalized], [-0.2, -0.2], color="black", lw=2)
ax.text(0.05 + scale_length_normalized / 2, -0.17, "100 kb", horizontalalignment='center', verticalalignment='bottom', fontsize=10)
ax.plot([0.05, 0.05], [-0.22, -0.18], color="black", lw=2)  # Move vertical lines to bottom-left at y = -0.2
ax.plot([0.05 + scale_length_normalized, 0.05 + scale_length_normalized], [-0.22, -0.18], color="black", lw=2)  # Move vertical lines to bottom-left at y = -0.2


# Show the plot
plt.show()
