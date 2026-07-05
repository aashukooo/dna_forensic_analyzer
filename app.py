import streamlit as st
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
from Bio.SeqUtils import molecular_weight
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import io

st.set_page_config(
    page_title="DNA Forensic Analyzer",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# ─────────────────────────────────────
# CUSTOM CSS
# ─────────────────────────────────────
st.markdown("""
<style>
    .main { background-color: #0a0a0a; }
    .stApp { background: linear-gradient(135deg, #0a0a0a 0%, #0d1117 100%); }
    
    .hero-title {
        text-align: center;
        font-size: 3.5em;
        font-weight: 900;
        background: linear-gradient(90deg, #00ff88, #00bfff, #ff00ff);
        -webkit-background-clip: text;
        -webkit-text-fill-color: transparent;
        margin-bottom: 0;
        padding: 20px;
    }
    
    .hero-subtitle {
        text-align: center;
        font-size: 1.2em;
        color: #888;
        margin-bottom: 30px;
    }

    .feature-card {
        background: linear-gradient(135deg, #1a1a2e, #16213e);
        border: 1px solid #00ff8833;
        border-radius: 15px;
        padding: 20px;
        margin: 10px 0;
        box-shadow: 0 4px 15px rgba(0, 255, 136, 0.1);
    }

    .metric-card {
        background: linear-gradient(135deg, #1a1a2e, #16213e);
        border: 1px solid #00bfff33;
        border-radius: 12px;
        padding: 15px;
        text-align: center;
    }

    .result-header {
        font-size: 1.5em;
        font-weight: 700;
        color: #00ff88;
        border-bottom: 2px solid #00ff8833;
        padding-bottom: 10px;
        margin-bottom: 15px;
    }

    .forensic-badge {
        background: linear-gradient(90deg, #ff4444, #ff8800);
        color: white;
        padding: 4px 12px;
        border-radius: 20px;
        font-size: 0.8em;
        font-weight: 600;
    }

    .verdict-high {
        background: linear-gradient(135deg, #00ff88, #00bfff);
        color: black;
        padding: 15px;
        border-radius: 10px;
        font-weight: 700;
        text-align: center;
        font-size: 1.1em;
    }

    .verdict-low {
        background: linear-gradient(135deg, #ff4444, #ff8800);
        color: white;
        padding: 15px;
        border-radius: 10px;
        font-weight: 700;
        text-align: center;
        font-size: 1.1em;
    }

    .stButton > button {
        background: linear-gradient(90deg, #00ff88, #00bfff);
        color: black;
        font-weight: 700;
        border: none;
        border-radius: 25px;
        padding: 10px 30px;
        font-size: 1em;
    }

    div[data-testid="stFileUploader"] {
        border: 2px dashed #00ff8866;
        border-radius: 15px;
        padding: 20px;
        background: #0d1117;
    }
</style>
""", unsafe_allow_html=True)

# ─────────────────────────────────────
# SIDEBAR
# ─────────────────────────────────────
with st.sidebar:
    st.markdown("""
    <div style='text-align:center; padding: 20px 0;'>
        <span style='font-size: 4em;'>🧬</span>
        <h2 style='color: #00ff88; margin: 0;'>DNA Forensic<br>Analyzer</h2>
        <p style='color: #888; font-size: 0.85em;'>Bioinformatics Tool v1.0</p>
    </div>
    """, unsafe_allow_html=True)

    st.markdown("---")
    st.markdown("### 📋 How to Use")
    st.markdown("""
    1. 📁 Upload a `.fasta` file
    2. 🔬 View instant analysis
    3. 📊 Explore visual charts
    4. 💾 Download your report
    """)

    st.markdown("---")
    st.markdown("### 🌐 Get FASTA Files")
    st.markdown("[NCBI Nucleotide Database](https://www.ncbi.nlm.nih.gov/nuccore)")
    st.markdown("[Ensembl Genome Browser](https://www.ensembl.org)")
    st.markdown("[UCSC Genome Browser](https://genome.ucsc.edu)")

    st.markdown("---")
    st.markdown("### ✅ Supported Formats")
    st.markdown("`.fasta` `.fa` `.txt`")

    st.markdown("---")
    st.markdown("""
    <div style='text-align:center; color: #888; font-size: 0.85em;'>
        Built by <b style='color:#00ff88'>Aastha Sharma</b><br>
        Undergraduate Bioinformatics Project<br><br>
        <a href='https://github.com/aashukooo/dna_forensic_analyzer' 
           style='color:#00bfff'>GitHub Repository</a>
    </div>
    """, unsafe_allow_html=True)

# ─────────────────────────────────────
# HERO
# ─────────────────────────────────────
st.markdown("<div class='hero-title'>🧬 DNA Forensic Analyzer</div>", unsafe_allow_html=True)
st.markdown("<div class='hero-subtitle'>Professional bioinformatics tool for forensic DNA sequence analysis using real FASTA genomic data</div>", unsafe_allow_html=True)

col1, col2, col3, col4, col5 = st.columns(5)
col1.success("🔬 Nucleotide Analysis")
col2.info("🧪 STR Profiling")
col3.warning("🔍 Mutation Detection")
col4.error("🦠 Species ID")
col5.success("📊 Visual Charts")

st.markdown("---")

# ─────────────────────────────────────
# FUNCTIONS
# ─────────────────────────────────────
def identify_species(description):
    known_species = {
        "homo sapiens": "🧑 Human (Homo sapiens)",
        "mus musculus": "🐭 Mouse (Mus musculus)",
        "rattus norvegicus": "🐀 Rat (Rattus norvegicus)",
        "sars-cov-2": "🦠 SARS-CoV-2 (COVID-19)",
        "escherichia coli": "🧫 E. Coli Bacteria",
        "saccharomyces cerevisiae": "🍞 Yeast",
        "danio rerio": "🐟 Zebrafish",
    }
    desc_lower = description.lower()
    for key, value in known_species.items():
        if key in desc_lower:
            return value
    return "❓ Unknown Species"

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
                results.append({
                    "unit": str(unit),
                    "count": count,
                    "position": i
                })
            i += 1
    return results

def detect_mutations(seq1, seq2):
    mutations = []
    min_len = min(len(seq1), len(seq2))
    for i in range(min_len):
        if seq1[i] != seq2[i]:
            mutations.append({
                "position": i + 1,
                "from": seq1[i],
                "to": seq2[i]
            })
    return mutations

def plot_nucleotide_chart(record):
    seq = record.seq
    counts = {
        'A': seq.count('A'),
        'T': seq.count('T'),
        'G': seq.count('G'),
        'C': seq.count('C')
    }
    colors = ['#00ff88', '#00bfff', '#ff00ff', '#ff8800']
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    fig.patch.set_facecolor('#0d1117')

    # Bar chart
    ax1.set_facecolor('#0d1117')
    bars = ax1.bar(counts.keys(), counts.values(), color=colors, width=0.5, edgecolor='none')
    ax1.set_title('Nucleotide Composition', color='white', fontsize=14, fontweight='bold')
    ax1.set_ylabel('Count', color='white')
    ax1.tick_params(colors='white')
    for spine in ax1.spines.values():
        spine.set_edgecolor('#333')
    for bar, val in zip(bars, counts.values()):
        ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 5,
                str(val), ha='center', va='bottom', color='white', fontweight='bold')

    # Pie chart
    ax2.set_facecolor('#0d1117')
    wedges, texts, autotexts = ax2.pie(
        counts.values(),
        labels=counts.keys(),
        colors=colors,
        autopct='%1.1f%%',
        startangle=90,
        textprops={'color': 'white'}
    )
    for at in autotexts:
        at.set_color('black')
        at.set_fontweight('bold')
    ax2.set_title('Nucleotide Distribution', color='white', fontsize=14, fontweight='bold')

    plt.tight_layout()
    return fig

# ─────────────────────────────────────
# FILE UPLOAD
# ─────────────────────────────────────
st.markdown("### 📁 Upload Your FASTA File")
uploaded_file = st.file_uploader(
    "Drag and drop or click to upload",
    type=["fasta", "fa", "txt"],
    help="Upload a FASTA formatted genomic sequence file"
)

if uploaded_file:
    content = uploaded_file.read().decode("utf-8")
    records = list(SeqIO.parse(io.StringIO(content), "fasta"))

    if not records:
        st.error("❌ No sequences found in the file! Make sure it's a valid FASTA format.")
    else:
        st.success(f"✅ Successfully loaded **{len(records)}** sequence(s) from **{uploaded_file.name}**")

        for idx, record in enumerate(records):
            st.markdown("---")
            species = identify_species(record.description)
            gc = gc_fraction(record.seq) * 100

            st.markdown(f"<div class='result-header'>🔬 Sequence {idx+1}: {record.id}</div>", unsafe_allow_html=True)

            # Top metrics
            col1, col2, col3, col4 = st.columns(4)
            col1.metric("🦠 Species", species)
            col2.metric("📏 Length", f"{len(record.seq):,} bp")
            col3.metric("🧪 GC Content", f"{gc:.2f}%")
            try:
                mw = molecular_weight(record.seq, "DNA")
                col4.metric("⚖️ Molecular Weight", f"{mw:,.2f} Da")
            except:
                col4.metric("⚖️ Molecular Weight", "N/A")

            st.markdown("---")

            # Nucleotide counts
            st.markdown("#### 🔢 Nucleotide Composition")
            seq = record.seq
            c1, c2, c3, c4 = st.columns(4)
            c1.metric("Adenine (A)", f"{seq.count('A'):,}")
            c2.metric("Thymine (T)", f"{seq.count('T'):,}")
            c3.metric("Guanine (G)", f"{seq.count('G'):,}")
            c4.metric("Cytosine (C)", f"{seq.count('C'):,}")

            # Chart
            st.markdown("#### 📊 Visual Analysis")
            fig = plot_nucleotide_chart(record)
            st.pyplot(fig)

            # Forensic strand analysis
            st.markdown("#### 🔍 Forensic Strand Analysis")
            col1, col2 = st.columns(2)
            with col1:
                st.markdown("**Complement Strand:**")
                st.code(str(seq.complement()), language=None)
            with col2:
                st.markdown("**Reverse Complement:**")
                st.code(str(seq.reverse_complement()), language=None)

            # STR Analysis
            st.markdown("#### 🧪 STR Analysis (Short Tandem Repeats)")
            str_results = find_str(seq)
            if str_results:
                st.success(f"Found **{len(str_results)}** STR pattern(s)!")
                str_data = {
                    "Repeat Unit": [r["unit"] for r in str_results],
                    "Count": [r["count"] for r in str_results],
                    "Position": [r["position"] for r in str_results]
                }
                st.dataframe(str_data, use_container_width=True)
            else:
                st.info("No significant STR patterns found in this sequence.")

        # Comparison between first two sequences
        if len(records) >= 2:
            st.markdown("---")
            st.markdown("### ⚖️ Forensic Sequence Comparison")
            st.markdown(f"Comparing **{records[0].id}** vs **{records[1].id}**")

            seq1 = records[0].seq
            seq2 = records[1].seq
            min_len = min(len(seq1), len(seq2))
            max_len = max(len(seq1), len(seq2))

            matches = sum(1 for a, b in zip(seq1, seq2) if a == b)
            similarity = (matches / max_len) * 100

            st.metric("Sequence Similarity", f"{similarity:.2f}%")

            if similarity > 90:
                st.markdown("<div class='verdict-high'>✅ HIGHLY SIMILAR — Likely same individual or close relative</div>", unsafe_allow_html=True)
            elif similarity > 70:
                st.markdown("<div class='verdict-low'>⚠️ MODERATELY SIMILAR — Possible distant relation</div>", unsafe_allow_html=True)
            else:
                st.markdown("<div class='verdict-low'>❌ LOW SIMILARITY — Likely different individuals</div>", unsafe_allow_html=True)

            # Mutation detection
            st.markdown("#### 🔬 Mutation Detection")
            mutations = detect_mutations(seq1, seq2)
            if mutations:
                st.warning(f"Found **{len(mutations)}** mutation(s) between sequences!")
                mut_data = {
                    "Position": [m["position"] for m in mutations[:50]],
                    "From": [m["from"] for m in mutations[:50]],
                    "To": [m["to"] for m in mutations[:50]]
                }
                st.dataframe(mut_data, use_container_width=True)
                if len(mutations) > 50:
                    st.info(f"Showing first 50 of {len(mutations)} mutations.")
            else:
                st.success("No mutations detected between these sequences!")