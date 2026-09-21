import streamlit as st
import pandas as pd
import primer3

from functools import lru_cache

from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction
from Bio.SeqUtils import MeltingTemp as mt


st.set_page_config(
    page_title="Primer Design Tool",
    page_icon="🧬",
    layout="wide",
)


VALID_BASES = set("ACGT")


# -----------------------------
# Sequence processing
# -----------------------------

def clean_sequence(raw_sequence):
    """Remove FASTA headers and whitespace."""
    lines = raw_sequence.strip().splitlines()
    sequence_parts = []

    for line in lines:
        line = line.strip()

        if not line or line.startswith(">"):
            continue

        sequence_parts.append("".join(line.split()))

    return "".join(sequence_parts).upper()


def validate_sequence(sequence):
    """Validate that the sequence contains only DNA bases."""
    if not sequence:
        return False, "Please paste a DNA sequence."

    if len(sequence) < 40:
        return False, "Sequence must contain at least 40 bases."

    invalid_bases = set(sequence) - VALID_BASES

    if invalid_bases:
        return (
            False,
            f"Invalid DNA base(s) found: "
            f"{', '.join(sorted(invalid_bases))}",
        )

    return True, ""


# -----------------------------
# Primer calculations
# -----------------------------

def calculate_gc(sequence):
    """Calculate GC percentage."""
    return gc_fraction(sequence) * 100


def calculate_tm(sequence):
    """Calculate approximate melting temperature using Wallace rule."""
    return mt.Tm_Wallace(sequence)


def primer_record(
    primer_type,
    sequence,
    start,
    end,
):
    """Create a consistent primer record."""
    return {
        "type": primer_type,
        "sequence": sequence,
        "start": start,
        "end": end,
        "length": len(sequence),
        "gc_percent": calculate_gc(sequence),
        "tm": calculate_tm(sequence),
    }


def generate_forward_primers(sequence, min_length, max_length):
    """Generate forward primers from all target-sequence positions."""
    primers = []

    for start_index in range(len(sequence)):
        for length in range(min_length, max_length + 1):
            end_index = start_index + length

            if end_index > len(sequence):
                continue

            primer = sequence[start_index:end_index]

            primers.append(
                primer_record(
                    primer_type="Forward",
                    sequence=primer,
                    start=start_index + 1,
                    end=end_index,
                )
            )

    return primers


def generate_reverse_primers(sequence, min_length, max_length):
    """Generate reverse primers from reverse-complemented target regions."""
    primers = []

    for start_index in range(len(sequence)):
        for length in range(min_length, max_length + 1):
            end_index = start_index + length

            if end_index > len(sequence):
                continue

            target_region = sequence[start_index:end_index]
            primer = str(Seq(target_region).reverse_complement())

            primers.append(
                primer_record(
                    primer_type="Reverse",
                    sequence=primer,
                    start=start_index + 1,
                    end=end_index,
                )
            )

    return primers


def filter_primers(
    primers,
    min_gc,
    max_gc,
    min_tm,
    max_tm,
):
    """Filter individual primers by GC percentage and Tm."""
    return [
        primer
        for primer in primers
        if (
            min_gc <= primer["gc_percent"] <= max_gc
            and min_tm <= primer["tm"] <= max_tm
        )
    ]


# -----------------------------
# Primer-pair generation
# -----------------------------

def generate_primer_pairs(
    forward_candidates,
    reverse_candidates,
    min_product_size,
    max_product_size,
    max_tm_difference,
):
    """Generate compatible forward/reverse primer pairs."""
    pairs = []

    for forward in forward_candidates:
        for reverse in reverse_candidates:
            # The reverse-binding region must be downstream and
            # must not overlap the forward-binding region.
            if reverse["start"] <= forward["end"]:
                continue

            product_size = reverse["end"] - forward["start"] + 1
            tm_difference = abs(forward["tm"] - reverse["tm"])

            if not (
                min_product_size
                <= product_size
                <= max_product_size
            ):
                continue

            if tm_difference > max_tm_difference:
                continue

            pairs.append(
                {
                    "forward_sequence": forward["sequence"],
                    "reverse_sequence": reverse["sequence"],
                    "forward_start": forward["start"],
                    "forward_end": forward["end"],
                    "reverse_start": reverse["start"],
                    "reverse_end": reverse["end"],
                    "product_size": product_size,
                    "forward_length": forward["length"],
                    "reverse_length": reverse["length"],
                    "forward_gc_percent": forward["gc_percent"],
                    "reverse_gc_percent": reverse["gc_percent"],
                    "forward_tm": forward["tm"],
                    "reverse_tm": reverse["tm"],
                    "tm_difference": tm_difference,
                }
            )

    # Rank by Tm similarity first, then shorter products.
    pairs.sort(
        key=lambda pair: (
            pair["tm_difference"],
            pair["product_size"],
        )
    )

    return pairs


# -----------------------------
# Quality checks
# -----------------------------

@lru_cache(maxsize=10000)
def analyze_primer_quality(sequence):
    """Analyze hairpin and self-dimer structures for one primer."""
    hairpin = primer3.calc_hairpin(sequence)
    self_dimer = primer3.calc_homodimer(sequence)

    return {
        "hairpin_tm": float(hairpin.tm),
        "hairpin_found": bool(hairpin.structure_found),
        "self_dimer_tm": float(self_dimer.tm),
        "self_dimer_found": bool(self_dimer.structure_found),
    }


@lru_cache(maxsize=10000)
def analyze_cross_dimer(forward_sequence, reverse_sequence):
    """Analyze cross-dimer formation between two primers."""
    cross_dimer = primer3.calc_heterodimer(
        forward_sequence,
        reverse_sequence,
    )

    return {
        "cross_dimer_tm": float(cross_dimer.tm),
        "cross_dimer_found": bool(cross_dimer.structure_found),
    }


def add_quality_results(
    primer_pairs,
    max_hairpin_tm,
    max_self_dimer_tm,
    max_cross_dimer_tm,
):
    """Add hairpin, self-dimer, and cross-dimer results to pairs."""
    checked_pairs = []

    for pair in primer_pairs:
        forward_quality = analyze_primer_quality(
            pair["forward_sequence"]
        )

        reverse_quality = analyze_primer_quality(
            pair["reverse_sequence"]
        )

        cross_quality = analyze_cross_dimer(
            pair["forward_sequence"],
            pair["reverse_sequence"],
        )

        pair_with_quality = {
            **pair,
            "forward_hairpin_tm": forward_quality["hairpin_tm"],
            "forward_self_dimer_tm": (
                forward_quality["self_dimer_tm"]
            ),
            "reverse_hairpin_tm": reverse_quality["hairpin_tm"],
            "reverse_self_dimer_tm": (
                reverse_quality["self_dimer_tm"]
            ),
            "cross_dimer_tm": cross_quality["cross_dimer_tm"],
            "passes_quality_checks": (
                forward_quality["hairpin_tm"] <= max_hairpin_tm
                and reverse_quality["hairpin_tm"] <= max_hairpin_tm
                and forward_quality["self_dimer_tm"]
                <= max_self_dimer_tm
                and reverse_quality["self_dimer_tm"]
                <= max_self_dimer_tm
                and cross_quality["cross_dimer_tm"]
                <= max_cross_dimer_tm
            ),
        }

        checked_pairs.append(pair_with_quality)

    return checked_pairs


# -----------------------------
# Display helpers
# -----------------------------

def format_primer_dataframe(primers):
    """Format individual primer records for display."""
    dataframe = pd.DataFrame(primers)

    if dataframe.empty:
        return dataframe

    dataframe = dataframe.rename(
        columns={
            "type": "Type",
            "sequence": "Sequence",
            "start": "Start",
            "end": "End",
            "length": "Length",
            "gc_percent": "GC%",
            "tm": "Tm (°C)",
        }
    )

    dataframe["GC%"] = dataframe["GC%"].round(2)
    dataframe["Tm (°C)"] = dataframe["Tm (°C)"].round(2)

    return dataframe


# -----------------------------
# Sidebar controls
# -----------------------------

st.sidebar.header("Primer Parameters")

min_length = st.sidebar.number_input(
    "Minimum primer length",
    min_value=10,
    max_value=50,
    value=18,
    step=1,
)

max_length = st.sidebar.number_input(
    "Maximum primer length",
    min_value=10,
    max_value=50,
    value=24,
    step=1,
)

min_gc = st.sidebar.number_input(
    "Minimum GC%",
    min_value=0.0,
    max_value=100.0,
    value=40.0,
    step=1.0,
)

max_gc = st.sidebar.number_input(
    "Maximum GC%",
    min_value=0.0,
    max_value=100.0,
    value=60.0,
    step=1.0,
)

min_tm = st.sidebar.number_input(
    "Minimum Tm (°C)",
    min_value=0.0,
    max_value=100.0,
    value=50.0,
    step=1.0,
)

max_tm = st.sidebar.number_input(
    "Maximum Tm (°C)",
    min_value=0.0,
    max_value=100.0,
    value=65.0,
    step=1.0,
)

st.sidebar.header("Primer-Pair Settings")

min_product_size = st.sidebar.number_input(
    "Minimum amplicon size",
    min_value=40,
    max_value=5000,
    value=100,
    step=10,
)

max_product_size = st.sidebar.number_input(
    "Maximum amplicon size",
    min_value=40,
    max_value=5000,
    value=1000,
    step=10,
)

max_tm_difference = st.sidebar.number_input(
    "Maximum primer Tm difference",
    min_value=0.0,
    max_value=20.0,
    value=3.0,
    step=0.5,
)

maximum_pairs_to_display = st.sidebar.number_input(
    "Maximum pairs to display",
    min_value=1,
    max_value=500,
    value=50,
    step=1,
)

st.sidebar.header("Quality Checks")

max_hairpin_tm = st.sidebar.number_input(
    "Maximum hairpin Tm",
    min_value=0.0,
    max_value=100.0,
    value=45.0,
    step=0.5,
)

max_self_dimer_tm = st.sidebar.number_input(
    "Maximum self-dimer Tm",
    min_value=0.0,
    max_value=100.0,
    value=45.0,
    step=0.5,
)

max_cross_dimer_tm = st.sidebar.number_input(
    "Maximum cross-dimer Tm",
    min_value=0.0,
    max_value=100.0,
    value=45.0,
    step=0.5,
)

show_only_passing_pairs = st.sidebar.checkbox(
    "Show only pairs passing quality checks",
    value=True,
)


# -----------------------------
# Main application
# -----------------------------

st.title("🧬 Primer Design Tool")

st.write(
    "Paste a DNA sequence to generate and evaluate basic PCR primer pairs."
)

raw_sequence = st.text_area(
    "Paste your DNA sequence",
    height=250,
    placeholder=(
        ">Example_sequence\n"
        "ATGCGTACGATCGATCGATCGATCGATCGATCGATCGATCGATCG"
    ),
)


if st.button("Generate Primers", type="primary"):

    if not raw_sequence.strip():
        st.warning("Please paste a DNA sequence first.")
        st.stop()

    if min_length > max_length:
        st.error(
            "Minimum primer length cannot be greater than maximum length."
        )
        st.stop()

    if min_gc > max_gc:
        st.error("Minimum GC% cannot be greater than maximum GC%.")
        st.stop()

    if min_tm > max_tm:
        st.error("Minimum Tm cannot be greater than maximum Tm.")
        st.stop()

    if min_product_size > max_product_size:
        st.error(
            "Minimum amplicon size cannot exceed maximum amplicon size."
        )
        st.stop()

    sequence = clean_sequence(raw_sequence)
    valid, error_message = validate_sequence(sequence)

    if not valid:
        st.error(error_message)
        st.stop()

    # Generate and filter individual primers.
    forward_primers = generate_forward_primers(
        sequence,
        min_length,
        max_length,
    )

    reverse_primers = generate_reverse_primers(
        sequence,
        min_length,
        max_length,
    )

    forward_filtered = filter_primers(
        forward_primers,
        min_gc,
        max_gc,
        min_tm,
        max_tm,
    )

    reverse_filtered = filter_primers(
        reverse_primers,
        min_gc,
        max_gc,
        min_tm,
        max_tm,
    )

    # Generate compatible primer pairs.
    primer_pairs = generate_primer_pairs(
        forward_filtered,
        reverse_filtered,
        min_product_size,
        max_product_size,
        max_tm_difference,
    )

    # Run quality checks on every compatible pair.
    quality_checked_pairs = add_quality_results(
        primer_pairs,
        max_hairpin_tm,
        max_self_dimer_tm,
        max_cross_dimer_tm,
    )

    all_pairs_df = pd.DataFrame(quality_checked_pairs)

    if not all_pairs_df.empty:
        numeric_columns = [
            "forward_gc_percent",
            "reverse_gc_percent",
            "forward_tm",
            "reverse_tm",
            "tm_difference",
            "forward_hairpin_tm",
            "forward_self_dimer_tm",
            "reverse_hairpin_tm",
            "reverse_self_dimer_tm",
            "cross_dimer_tm",
        ]

        all_pairs_df[numeric_columns] = (
            all_pairs_df[numeric_columns].round(2)
        )

    passing_pairs_df = all_pairs_df

    if show_only_passing_pairs and not all_pairs_df.empty:
        passing_pairs_df = all_pairs_df[
            all_pairs_df["passes_quality_checks"]
        ].reset_index(drop=True)

    # Target information.
    st.success(f"Valid DNA sequence: {len(sequence)} bp")

    col1, col2, col3, col4 = st.columns(4)

    with col1:
        st.metric("Sequence Length", f"{len(sequence)} bp")

    with col2:
        st.metric("Target GC%", f"{calculate_gc(sequence):.2f}%")

    with col3:
        st.metric(
            "Approx. Target Tm",
            f"{calculate_tm(sequence):.2f} °C",
        )

    with col4:
        st.metric("Compatible Pairs", len(primer_pairs))

    # Individual primer results.
    st.subheader("Forward Primer Candidates")

    forward_df = format_primer_dataframe(forward_filtered)

    if forward_df.empty:
        st.info(
            "No forward primers matched the selected GC% and Tm ranges."
        )
    else:
        st.dataframe(
            forward_df,
            use_container_width=True,
            hide_index=True,
        )

    st.subheader("Reverse Primer Candidates")

    reverse_df = format_primer_dataframe(reverse_filtered)

    if reverse_df.empty:
        st.info(
            "No reverse primers matched the selected GC% and Tm ranges."
        )
    else:
        st.dataframe(
            reverse_df,
            use_container_width=True,
            hide_index=True,
        )

    # Pair results.
    st.subheader("Primer Pairs and Quality Checks")

    passing_count = 0

    if not all_pairs_df.empty:
        passing_count = int(
            all_pairs_df["passes_quality_checks"].sum()
        )

    st.write(
        f"{passing_count} of {len(primer_pairs)} compatible pairs "
        "passed the selected quality thresholds."
    )

    if passing_pairs_df.empty:
        if show_only_passing_pairs and not all_pairs_df.empty:
            st.warning(
                "No primer pairs passed the hairpin and dimer thresholds. "
                "Uncheck 'Show only pairs passing quality checks' to inspect "
                "all compatible pairs."
            )
        else:
            st.info(
                "No compatible primer pairs matched the selected settings."
            )
    else:
        pair_columns = [
            "forward_sequence",
            "reverse_sequence",
            "product_size",
            "forward_start",
            "forward_end",
            "reverse_start",
            "reverse_end",
            "forward_tm",
            "reverse_tm",
            "tm_difference",
            "forward_hairpin_tm",
            "forward_self_dimer_tm",
            "reverse_hairpin_tm",
            "reverse_self_dimer_tm",
            "cross_dimer_tm",
            "passes_quality_checks",
        ]

        st.dataframe(
            passing_pairs_df[pair_columns].head(
                int(maximum_pairs_to_display)
            ),
            use_container_width=True,
            hide_index=True,
        )

    with st.expander("View Cleaned DNA Sequence"):
        st.code(sequence, language="text")