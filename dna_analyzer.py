
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
from Bio.SeqUtils import molecular_weight
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import matplotlib.pyplot as plt
import os

# ─────────────────────────────────────
# SETTINGS
# ─────────────────────────────────────
FASTA_FILE = "ncbi1.fasta"
MOTIF = "ATGC"
OUTPUT_FILE = "dna_report.txt"

# ─────────────────────────────────────
# FUNCTION: Species specification
# ─────────────────────────────────────
def identify_species(description):
    known_species = {
        "homo sapiens": "Human",
        "mus musculus": "Mouse",
        "rattus norvegicus": "Rat",
        "drosophila melanogaster": "Fruit Fly",
        "danio rerio": "Zebrafish",
        "sars-cov-2": "SARS-CoV-2 Virus (COVID-19)",
        "escherichia coli": "E. Coli Bacteria",
        "saccharomyces cerevisiae": "Yeast",
    }
    desc_lower = description.lower()
    for key, value in known_species.items():
        if key in desc_lower:
            return value
    return "Unknown Species"

# ─────────────────────────────────────
# FUNCTION: STR Analysis
# ─────────────────────────────────────
def find_str(seq, min_repeat=3):
    lines = []
    lines.append(f"\n--- STR Analysis (Short Tandem Repeats) ---")
    found_any = False
    for unit_len in range(2, 7):
        i = 0
        while i < len(seq) - unit_len:
            unit = seq[i:i+unit_len]
            count = 1
            j = i + unit_len
            while j + unit_len <= len(seq) and seq[j:j+unit_len] == unit:
                count += 1
                j += unit_len
            if count >= min_repeat:
                lines.append(f"  Repeat unit: {unit} | Count: {count} | Position: {i}")
                found_any = True
            i += 1
    if not found_any:
        lines.append("  No significant STR patterns found.")
    return lines

# ─────────────────────────────────────
# FUNCTION: Mutation Detection
# ─────────────────────────────────────
def detect_mutations(seq1, seq2, id1, id2):
    lines = []
    lines.append(f"\n--- Mutation Detection: {id1} vs {id2} ---")
    min_len = min(len(seq1), len(seq2))
    mutations = []
    for i in range(min_len):
        if seq1[i] != seq2[i]:
            mutations.append(f"  Position {i+1}: {seq1[i]} → {seq2[i]}")
    if mutations:
        lines.append(f"  Total mutations found: {len(mutations)}")
        for m in mutations:
            lines.append(m)
    else:
        lines.append("  No mutations detected between these sequences.")
    if len(seq1) != len(seq2):
        lines.append(f"  Length difference: {abs(len(seq1) - len(seq2))} bp")
    return lines

# ─────────────────────────────────────
# FUNCTION: Nucleotide Chart
# ─────────────────────────────────────
def plot_nucleotide_chart(records):
    names = []
    a_counts = []
    t_counts = []
    g_counts = []
    c_counts = []

    for record in records:
        names.append(record.id)
        a_counts.append(record.seq.count('A'))
        t_counts.append(record.seq.count('T'))
        g_counts.append(record.seq.count('G'))
        c_counts.append(record.seq.count('C'))

    x = range(len(names))
    width = 0.2

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.bar([i - 1.5*width for i in x], a_counts, width, label='A', color='steelblue')
    ax.bar([i - 0.5*width for i in x], t_counts, width, label='T', color='coral')
    ax.bar([i + 0.5*width for i in x], g_counts, width, label='G', color='mediumseagreen')
    ax.bar([i + 1.5*width for i in x], c_counts, width, label='C', color='mediumpurple')

    ax.set_xlabel('Sequence ID')
    ax.set_ylabel('Nucleotide Count')
    ax.set_title('Nucleotide Composition per Sequence')
    ax.set_xticks(list(x))
    ax.set_xticklabels(names)
    ax.legend()
    plt.tight_layout()
    plt.savefig("nucleotide_chart.png")
    plt.show()
    print("\nChart saved as: nucleotide_chart.png")

# ─────────────────────────────────────
# FUNCTION: Analyze single sequence
# ─────────────────────────────────────
def analyze_sequence(record):
    lines = []
    lines.append(f"{'='*50}")
    lines.append(f"SEQUENCE ID   : {record.id}")
    lines.append(f"DESCRIPTION   : {record.description}")

    # Species
    species = identify_species(record.description)
    lines.append(f"SPECIES       : {species}")
    lines.append(f"{'='*50}")

    seq = record.seq

    # Basic info
    lines.append(f"Sequence      : {seq}")
    lines.append(f"Length        : {len(seq)} bp")

    # Nucleotide counts
    lines.append(f"\n--- Nucleotide Composition ---")
    lines.append(f"A : {seq.count('A')}")
    lines.append(f"T : {seq.count('T')}")
    lines.append(f"G : {seq.count('G')}")
    lines.append(f"C : {seq.count('C')}")

    # GC Content
    gc = gc_fraction(seq) * 100
    lines.append(f"\n--- GC Content ---")
    lines.append(f"GC Content    : {gc:.2f}%")

    # Molecular weight
    try:
        mw = molecular_weight(seq, "DNA")
        lines.append(f"\n--- Molecular Weight ---")
        lines.append(f"Mol. Weight   : {mw:.2f} Da")
    except Exception:
        lines.append("Molecular weight could not be calculated.")

    # Complement and Reverse Complement
    lines.append(f"\n--- Forensic Strand Analysis ---")
    lines.append(f"Complement    : {seq.complement()}")
    lines.append(f"Rev Complement: {seq.reverse_complement()}")

    # Motif search
    lines.append(f"\n--- Motif Search: '{MOTIF}' ---")
    positions = []
    start = 0
    while True:
        pos = seq.find(MOTIF, start)
        if pos == -1:
            break
        positions.append(pos)
        start = pos + 1
    if positions:
        lines.append(f"Found at positions: {positions}")
    else:
        lines.append(f"Motif '{MOTIF}' not found in this sequence.")

    # STR Analysis
    str_results = find_str(seq)
    lines.extend(str_results)

    return lines

# ─────────────────────────────────────
# MAIN PROGRAM
# ─────────────────────────────────────
print(f"\n{'='*50}")
print(f"   DNA FORENSIC ANALYZER")
print(f"{'='*50}\n")

all_lines = []
records = list(SeqIO.parse(FASTA_FILE, "fasta"))

if not records:
    print("No sequences found in file!")
else:
    for record in records:
        result = analyze_sequence(record)
        for line in result:
            print(line)
        all_lines.extend(result)

    # Compare first two sequences if available
    if len(records) >= 2:
        # Sequence comparison
        comparison = []
        comparison.append(f"\n{'='*50}")
        comparison.append(f"FORENSIC SEQUENCE COMPARISON")
        comparison.append(f"{'='*50}")
        alignments = pairwise2.align.globalxx(records[0].seq, records[1].seq)
        if alignments:
            comparison.append(format_alignment(*alignments[0]))
            score = alignments[0].score
            max_len = max(len(records[0].seq), len(records[1].seq))
            similarity = (score / max_len) * 100
            comparison.append(f"Similarity    : {similarity:.2f}%")
            if similarity > 90:
                comparison.append("Verdict: HIGHLY SIMILAR — likely same individual or close relative")
            elif similarity > 70:
                comparison.append("Verdict: MODERATELY SIMILAR — possible distant relation")
            else:
                comparison.append("Verdict: LOW SIMILARITY — likely different individuals")
        for line in comparison:
            print(line)
        all_lines.extend(comparison)

        # Mutation detection
        mutations = detect_mutations(
            records[0].seq, records[1].seq,
            records[0].id, records[1].id
        )
        for line in mutations:
            print(line)
        all_lines.extend(mutations)

    # Plot chart
    plot_nucleotide_chart(records)

    # Export report
    with open(OUTPUT_FILE, "w", encoding="utf-8") as f:
        f.write("\n".join(all_lines))
    print(f"\nReport saved to: {OUTPUT_FILE}")