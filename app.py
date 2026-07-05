import streamlit as st
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
from Bio.SeqUtils import molecular_weight
import matplotlib.pyplot as plt
import io

st.set_page_config(page_title="DNA Forensic Analyzer", page_icon="🧬")

st.title("🧬 DNA Forensic Analyzer")
st.write("Upload a FASTA file to analyze your DNA sequence")

# ─────────────────────────────────────
# SPECIES IDENTIFICATION
# ─────────────────────────────────────
def identify_species(description):
    known_species = {
        "homo sapiens": "Human",
        "mus musculus": "Mouse",
        "rattus norvegicus": "Rat",
        "sars-cov-2": "SARS-CoV-2 Virus (COVID-19)",
        "escherichia coli": "E. Coli Bacteria",
    }
    desc_lower = description.lower()
    for key, value in known_species.items():
        if key in desc_lower:
            return value
    return "Unknown Species"

# ─────────────────────────────────────
# STR ANALYSIS
# ─────────────────────────────────────
def find_str(seq, min_repeat=3):
    results = []
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
                results.append(f"Repeat unit: `{unit}` | Count: {count} | Position: {i}")
            i += 1
    return results

# ─────────────────────────────────────
# FILE UPLOAD
# ─────────────────────────────────────
uploaded_file = st.file_uploader("Upload your FASTA file", type=["fasta", "fa", "txt"])

if uploaded_file:
    content = uploaded_file.read().decode("utf-8")
    records = list(SeqIO.parse(io.StringIO(content), "fasta"))

    if not records:
        st.error("No sequences found in the file!")
    else:
        st.success(f"Found {len(records)} sequence(s) in the file!")

        for record in records:
            st.markdown(f"---")
            st.subheader(f"🔬 {record.id}")

            species = identify_species(record.description)
            col1, col2, col3 = st.columns(3)
            col1.metric("Species", species)
            col2.metric("Length", f"{len(record.seq)} bp")

            gc = gc_fraction(record.seq) * 100
            col3.metric("GC Content", f"{gc:.2f}%")

            try:
                mw = molecular_weight(record.seq, "DNA")
                st.metric("Molecular Weight", f"{mw:.2f} Da")
            except:
                pass

            # Nucleotide counts
            st.markdown("### Nucleotide Composition")
            a = record.seq.count('A')
            t = record.seq.count('T')
            g = record.seq.count('G')
            c = record.seq.count('C')

            col1, col2, col3, col4 = st.columns(4)
            col1.metric("A", a)
            col2.metric("T", t)
            col3.metric("G", g)
            col4.metric("C", c)

            # Chart
            fig, ax = plt.subplots()
            ax.bar(['A', 'T', 'G', 'C'], [a, t, g, c],
                   color=['steelblue', 'coral', 'mediumseagreen', 'mediumpurple'])
            ax.set_title(f"Nucleotide Composition — {record.id}")
            ax.set_ylabel("Count")
            st.pyplot(fig)

            # Strands
            st.markdown("### 🔍 Forensic Strand Analysis")
            st.code(f"Complement    : {record.seq.complement()}")
            st.code(f"Rev Complement: {record.seq.reverse_complement()}")

            # STR
            st.markdown("### 🧪 STR Analysis")
            str_results = find_str(record.seq)
            if str_results:
                for r in str_results:
                    st.write(r)
            else:
                st.write("No significant STR patterns found.")