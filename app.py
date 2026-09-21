import streamlit as st
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction
from Bio.SeqUtils import MeltingTemp as mt


st.set_page_config(
    page_title="Primer Design Tool",
    page_icon="🧬",
    layout="wide"
)


# -----------------------------
# Sequence processing
# -----------------------------

def clean_sequence(raw_sequence):
    """Remove FASTA headers and whitespace."""
    lines = raw_sequence.strip().splitlines()
    sequence_parts = []

    for line in lines:
        line = line.strip()

        if line.startswith(">"):
            continue

        sequence_parts.append("".join(line.split()))

    return "".join(sequence_parts).upper()


def validate_sequence(sequence):
    """Validate that the sequence contains only DNA bases."""
    if len(sequence) < 40:
        return False, "Sequence must contain at least 40 bases."

    invalid_bases = set(sequence) - set("ACGT")

    if invalid_bases:
        return (
            False,
            f"Invalid DNA base(s) found: {', '.join(sorted(invalid_bases))}"
        )

    return True, ""


# -----------------------------
# Primer calculations
# -----------------------------

def calculate_gc(sequence):
    """Calculate GC percentage."""
    return gc_fraction(sequence) * 100


def calculate_tm(sequence):
    """Calculate approximate melting temperature."""
    return mt.Tm_Wallace(sequence)


def generate_forward_primers(sequence, min_length, max_length):
    """Generate forward primers from the beginning of the target sequence."""
    primers = []

    for length in range(min_length, max_length + 1):
        if length <= len(sequence):
            primer = sequence[:length]

            primers.append({
                "Type": "Forward",
                "Sequence": primer,
                "Length": len(primer),
                "GC%": calculate_gc(primer),
                "Tm (°C)": calculate_tm(primer),
            })

    return primers


def generate_reverse_primers(sequence, min_length, max_length):
    """Generate reverse primers as reverse complements of target slices."""
    primers = []

    for length in range(min_length, max_length + 1):
        if length <= len(sequence):
            target_slice = sequence[-length:]
            primer = str(Seq(target_slice).reverse_complement())

            primers.append({
                "Type": "Reverse",
                "Sequence": primer,
                "Length": len(primer),
                "GC%": calculate_gc(primer),
                "Tm (°C)": calculate_tm(primer),
            })

    return primers


# -----------------------------
# Sidebar controls
# -----------------------------

st.sidebar.header("Primer Parameters")

min_length = st.sidebar.number_input(
    "Minimum primer length",
    min_value=10,
    max_value=50,
    value=18,
    step=1
)

max_length = st.sidebar.number_input(
    "Maximum primer length",
    min_value=10,
    max_value=50,
    value=24,
    step=1
)

min_gc = st.sidebar.number_input(
    "Minimum GC%",
    min_value=0.0,
    max_value=100.0,
    value=40.0,
    step=1.0
)

max_gc = st.sidebar.number_input(
    "Maximum GC%",
    min_value=0.0,
    max_value=100.0,
    value=60.0,
    step=1.0
)

min_tm = st.sidebar.number_input(
    "Minimum Tm (°C)",
    min_value=0.0,
    max_value=100.0,
    value=50.0,
    step=1.0
)

max_tm = st.sidebar.number_input(
    "Maximum Tm (°C)",
    min_value=0.0,
    max_value=100.0,
    value=65.0,
    step=1.0
)


# -----------------------------
# Main application
# -----------------------------

st.title("🧬 Primer Design Tool")

st.write(
    "Paste a DNA sequence to generate forward and reverse primer candidates "
    "based on your selected primer parameters."
)

raw_sequence = st.text_area(
    "Paste your DNA sequence",
    height=250,
    placeholder=(
        ">Example_sequence\n"
        "ATGCGTACGATCGATCGATCGATCGATCGATCGATCGATCGATCG"
    )
)


if st.button("Generate Primers", type="primary"):

    if not raw_sequence.strip():
        st.warning("Please paste a DNA sequence first.")

    elif min_length > max_length:
        st.error("Minimum primer length cannot be greater than maximum length.")

    elif min_gc > max_gc:
        st.error("Minimum GC% cannot be greater than maximum GC%.")

    elif min_tm > max_tm:
        st.error("Minimum Tm cannot be greater than maximum Tm.")

    else:
        sequence = clean_sequence(raw_sequence)

        valid, error_message = validate_sequence(sequence)

        if not valid:
            st.error(error_message)

        else:
            st.success(f"Valid DNA sequence: {len(sequence)} bp")

            # Generate candidate primers
            forward_primers = generate_forward_primers(
                sequence,
                min_length,
                max_length
            )

            reverse_primers = generate_reverse_primers(
                sequence,
                min_length,
                max_length
            )

            # Filter candidates
            forward_filtered = [
                primer for primer in forward_primers
                if min_gc <= primer["GC%"] <= max_gc
                and min_tm <= primer["Tm (°C)"] <= max_tm
            ]

            reverse_filtered = [
                primer for primer in reverse_primers
                if min_gc <= primer["GC%"] <= max_gc
                and min_tm <= primer["Tm (°C)"] <= max_tm
            ]

            # Sequence information
            st.subheader("Target Sequence")

            col1, col2, col3 = st.columns(3)

            with col1:
                st.metric("Sequence Length", f"{len(sequence)} bp")

            with col2:
                st.metric("Target GC%", f"{calculate_gc(sequence):.2f}%")

            with col3:
                st.metric(
                    "Approx. Target Tm",
                    f"{calculate_tm(sequence):.2f} °C"
                )

            # Primer results
            st.subheader("Forward Primer Candidates")

            if forward_filtered:
                st.dataframe(
                    forward_filtered,
                    use_container_width=True
                )
            else:
                st.info(
                    "No forward primers matched the selected GC% and Tm ranges."
                )

            st.subheader("Reverse Primer Candidates")

            if reverse_filtered:
                st.dataframe(
                    reverse_filtered,
                    use_container_width=True
                )
            else:
                st.info(
                    "No reverse primers matched the selected GC% and Tm ranges."
                )

            # Clean sequence
            with st.expander("View Cleaned DNA Sequence"):
                st.code(sequence, language="text")
