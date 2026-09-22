import streamlit as st
import pandas as pd

from Bio.Data import CodonTable


st.set_page_config(
    page_title="RSCU Calculator",
    page_icon="🧬",
    layout="wide",
)


STANDARD_TABLE_ID = 1
VALID_BASES = set("ACGT")


def clean_sequence(raw_sequence):
    """Remove FASTA headers, whitespace, and convert to uppercase."""
    lines = raw_sequence.strip().splitlines()
    sequence_parts = []

    for line in lines:
        line = line.strip()

        if not line or line.startswith(">"):
            continue

        sequence_parts.append("".join(line.split()))

    return "".join(sequence_parts).upper()


def validate_cds(sequence, reading_frame):
    """Validate a coding DNA sequence for codon analysis."""
    if not sequence:
        return False, "Please paste a coding DNA sequence."

    invalid_bases = set(sequence) - VALID_BASES

    if invalid_bases:
        return (
            False,
            f"Invalid DNA base(s): "
            f"{', '.join(sorted(invalid_bases))}",
        )

    framed_sequence = sequence[reading_frame:]

    if not framed_sequence:
        return False, "The selected reading frame contains no sequence."

    if len(framed_sequence) % 3 != 0:
        return (
            False,
            "The sequence length in the selected reading frame "
            "must be divisible by 3.",
        )

    return True, ""


def get_codon_groups():
    """Return codons grouped by encoded amino acid."""
    table = CodonTable.unambiguous_dna_by_id[STANDARD_TABLE_ID]

    codon_to_amino_acid = {
        codon: amino_acid
        for codon, amino_acid in table.forward_table.items()
    }

    amino_acid_to_codons = {}

    for codon, amino_acid in codon_to_amino_acid.items():
        amino_acid_to_codons.setdefault(
            amino_acid,
            [],
        ).append(codon)

    for amino_acid in amino_acid_to_codons:
        amino_acid_to_codons[amino_acid].sort()

    return codon_to_amino_acid, amino_acid_to_codons, table.stop_codons


def calculate_rscu(sequence, reading_frame=0, include_stop_codons=False):
    """Calculate codon counts, frequencies, and RSCU values."""
    codon_to_amino_acid, amino_acid_to_codons, stop_codons = (
        get_codon_groups()
    )

    framed_sequence = sequence[reading_frame:]
    codons = [
        framed_sequence[index:index + 3]
        for index in range(0, len(framed_sequence), 3)
    ]

    codon_counts = {}

    for codon in codons:
        if codon in codon_to_amino_acid:
            codon_counts[codon] = codon_counts.get(codon, 0) + 1
        elif include_stop_codons and codon in stop_codons:
            codon_counts[codon] = codon_counts.get(codon, 0) + 1

    rows = []

    for amino_acid, synonymous_codons in amino_acid_to_codons.items():
        total_count = sum(
            codon_counts.get(codon, 0)
            for codon in synonymous_codons
        )

        number_of_synonymous_codons = len(synonymous_codons)
        expected_count = (
            total_count / number_of_synonymous_codons
            if total_count > 0
            else 0
        )

        for codon in synonymous_codons:
            observed_count = codon_counts.get(codon, 0)

            if expected_count > 0:
                rscu = observed_count / expected_count
            else:
                rscu = 0.0

            rows.append(
                {
                    "Amino acid": amino_acid,
                    "Codon": codon,
                    "Count": observed_count,
                    "Expected count": expected_count,
                    "RSCU": rscu,
                    "Synonymous codons": number_of_synonymous_codons,
                }
            )

    if include_stop_codons:
        stop_count = sum(
            codon_counts.get(codon, 0)
            for codon in stop_codons
        )

        expected_stop_count = stop_count / len(stop_codons)

        for codon in sorted(stop_codons):
            observed_count = codon_counts.get(codon, 0)

            if expected_stop_count > 0:
                rscu = observed_count / expected_stop_count
            else:
                rscu = 0.0

            rows.append(
                {
                    "Amino acid": "Stop",
                    "Codon": codon,
                    "Count": observed_count,
                    "Expected count": expected_stop_count,
                    "RSCU": rscu,
                    "Synonymous codons": len(stop_codons),
                }
            )

    dataframe = pd.DataFrame(rows)

    if not dataframe.empty:
        dataframe["Expected count"] = dataframe[
            "Expected count"
        ].round(3)

        dataframe["RSCU"] = dataframe["RSCU"].round(3)

    return dataframe, codons


st.title("🧬 Codon Usage Bias Calculator")
st.write(
    "Calculate codon counts, expected counts, and Relative "
    "Synonymous Codon Usage (RSCU) from a coding DNA sequence."
)

reading_frame = st.selectbox(
    "Reading frame",
    options=[0, 1, 2],
    format_func=lambda frame: f"Frame {frame + 1}",
)

include_stop_codons = st.checkbox(
    "Include stop codons",
    value=False,
    help=(
        "Stop codons are normally excluded from RSCU calculations "
        "because they do not encode an amino acid."
    ),
)

raw_sequence = st.text_area(
    "Paste a coding DNA sequence",
    height=220,
    placeholder=(
        ">example_CDS\n"
        "ATGGCTGCTGCTAAAGGCGGTTAA"
    ),
)

if st.button("Calculate RSCU", type="primary"):
    sequence = clean_sequence(raw_sequence)

    valid, message = validate_cds(
        sequence,
        reading_frame,
    )

    if not valid:
        st.error(message)
        st.stop()

    dataframe, codons = calculate_rscu(
        sequence,
        reading_frame,
        include_stop_codons,
    )

    if dataframe.empty:
        st.warning(
            "No complete coding codons were found in the sequence."
        )
        st.stop()

    total_coding_codons = int(dataframe["Count"].sum())
    used_codons = int((dataframe["Count"] > 0).sum())

    st.success(
        f"Sequence accepted: {len(sequence[reading_frame:])} bases "
        f"({len(codons)} codons)"
    )

    column1, column2, column3 = st.columns(3)

    with column1:
        st.metric("Codons analyzed", len(codons))

    with column2:
        st.metric("Observed codon count", total_coding_codons)

    with column3:
        st.metric("Codons observed", used_codons)

    st.subheader("RSCU results")

    st.dataframe(
        dataframe,
        use_container_width=True,
        hide_index=True,
    )

    st.subheader("RSCU chart")

    chart_data = dataframe.set_index("Codon")[["RSCU"]]
    st.bar_chart(chart_data)

    csv_data = dataframe.to_csv(index=False).encode("utf-8")

    st.download_button(
        label="Download RSCU results as CSV",
        data=csv_data,
        file_name="rscu_results.csv",
        mime="text/csv",
    )