import streamlit as st
from Bio.SeqUtils import MeltingTemp as mt


def clean_sequence(raw_sequence):
    """Remove FASTA headers, whitespace, and convert to uppercase."""
    lines = raw_sequence.strip().splitlines()

    sequence_lines = [
        line.strip()
        for line in lines
        if line.strip() and not line.strip().startswith(">")
    ]

    return "".join(sequence_lines).replace(" ", "").upper()


def gc_percent(sequence):
    """Calculate GC percentage."""
    gc_count = sequence.count("G") + sequence.count("C")
    return gc_count / len(sequence) * 100


def validate_sequence(sequence):
    """Check that the sequence contains only DNA bases."""
    if not sequence:
        return False, "Please enter a DNA sequence."

    invalid_bases = set(sequence) - set("ACGT")

    if invalid_bases:
        return False, f"Invalid base(s): {', '.join(invalid_bases)}"

    if len(sequence) < 40:
        return False, "Please enter a sequence at least 40 bases long."

    return True, ""


st.set_page_config(
    page_title="Primer Design Tool",
    page_icon="🧬",
)

st.title("🧬 Primer Design Tool")
st.write("Analyze a DNA target sequence.")

raw_sequence = st.text_area(
    "Paste your DNA sequence",
    height=180,
    placeholder="ATGCGTACGATCGATCGATCG...",
)

if st.button("Analyze sequence"):
    sequence = clean_sequence(raw_sequence)
    valid, message = validate_sequence(sequence)

    if not valid:
        st.error(message)
    else:
        gc = gc_percent(sequence)
        tm = mt.Tm_NN(sequence)

        st.success("Valid DNA sequence")

        column1, column2, column3 = st.columns(3)

        with column1:
            st.metric("Length", f"{len(sequence)} bp")

        with column2:
            st.metric("GC content", f"{gc:.2f}%")

        with column3:
            st.metric("Approximate Tm", f"{tm:.2f} °C")

        st.subheader("Cleaned sequence")
        st.code(sequence)
        