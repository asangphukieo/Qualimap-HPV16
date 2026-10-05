#!/usr/bin/env python3
"""
Generate mock test data for Qualimap HPV16 QC pipeline.

Creates:
  - Mock HPV16 reference genome BAM files (with gi|333031|lcl|HPV16REF.1| header)
  - Example output files (summary_table.txt format)
  - HPV16 gene annotation BED file
"""

import os
import random
import pysam

random.seed(42)

HPV16_LEN = 7906
HPV16_REF = "gi|333031|lcl|HPV16REF.1|"
READ_LEN = 150
NUM_PAIRS = 200
INSERT_SIZE = 300

# HPV16 gene coordinates
HPV16_GENES = {
    "E6":  (104, 559),
    "E7":  (562, 858),
    "E1":  (865, 2814),
    "E2":  (2756, 3853),
    "E5a": (3850, 4101),
    "L2":  (4237, 5658),
    "L1":  (5639, 7156),
}

OUTPUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "test_data")
BAM_DIR = os.path.join(OUTPUT_DIR, "input_bam")
EXAMPLE_DIR = os.path.join(OUTPUT_DIR, "example_output")

os.makedirs(BAM_DIR, exist_ok=True)
os.makedirs(EXAMPLE_DIR, exist_ok=True)


def random_seq(length):
    return "".join(random.choice("ACGT") for _ in range(length))


def reverse_complement(seq):
    comp = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}
    return "".join(comp.get(b, "N") for b in reversed(seq))


def simulate_reads(genome_len, num_pairs, read_len, insert_size):
    reads = []
    genome = random_seq(genome_len)
    for i in range(num_pairs):
        ins = max(read_len + 10, int(random.gauss(insert_size, 50)))
        if ins > genome_len:
            ins = genome_len
        start = random.randint(0, genome_len - ins)
        fragment = genome[start : start + ins]
        r1_seq = fragment[:read_len]
        r2_seq = reverse_complement(fragment[-read_len:])
        reads.append({
            "name": f"read_{i:04d}",
            "r1_seq": r1_seq,
            "r2_seq": r2_seq,
            "pos": start,
            "ins": ins,
        })
    return reads


def create_bam(bam_path, sample_name, ref_name, ref_len, reads, read_len):
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": ref_name, "LN": ref_len}],
        "RG": [{"ID": sample_name, "SM": sample_name, "PL": "ILLUMINA"}],
    }
    with pysam.AlignmentFile(bam_path, "wb", header=header) as outf:
        for r in reads:
            a1 = pysam.AlignedSegment()
            a1.query_name = r["name"]
            a1.query_sequence = r["r1_seq"]
            a1.flag = 99
            a1.reference_id = 0
            a1.reference_start = r["pos"]
            a1.mapping_quality = 60
            a1.cigar = [(0, read_len)]
            a1.next_reference_id = 0
            a1.next_reference_start = r["pos"] + r["ins"] - read_len
            a1.template_length = r["ins"]
            a1.query_qualities = pysam.qualitystring_to_array("I" * read_len)
            a1.set_tag("RG", sample_name)

            a2 = pysam.AlignedSegment()
            a2.query_name = r["name"]
            a2.query_sequence = r["r2_seq"]
            a2.flag = 147
            a2.reference_id = 0
            a2.reference_start = r["pos"] + r["ins"] - read_len
            a2.mapping_quality = 60
            a2.cigar = [(0, read_len)]
            a2.next_reference_id = 0
            a2.next_reference_start = r["pos"]
            a2.template_length = -r["ins"]
            a2.query_qualities = pysam.qualitystring_to_array("I" * read_len)
            a2.set_tag("RG", sample_name)

            outf.write(a1)
            outf.write(a2)

    pysam.sort("-o", bam_path + ".tmp", bam_path)
    os.rename(bam_path + ".tmp", bam_path)
    pysam.index(bam_path)


def write_bed(path):
    with open(path, "w") as f:
        for gene, (start, end) in HPV16_GENES.items():
            f.write(f"{HPV16_REF}\t{start}\t{end}\t{gene}\n")


def write_example_summary(path):
    header = (
        "Sample\ttotal_passed\ttotal_failed\tmapped_passed\tmapped_failed\t"
        "total_read_original\ttotal_read_qc_pri\tmost_read_map\t"
        "total_read_hpv16\ttotal_read_unmap\t"
        "E6\tE7\tE1\tE2\tE5_a\tL2\tL1\t"
        "E6_f\tE7_f\tE1_f\tE2_f\tE5_f\tL2_f\tL1_f\t"
        "Lineage\tSublineage\tLineage_counts\tSublineage_counts\t"
        "Average_Depth\tQ1_Depth\tQ3_Depth\tMedian_Depth\tStd_Depth"
    )
    row1 = (
        "sample1\t400\t0\t380\t0\t"
        "400\t350\tgi|333031|lcl|HPV16REF.1|\t"
        "320\t20\t"
        "45\t30\t120\t80\t25\t90\t95\t"
        "3164\t3243\t1967\t2330\t3187\t2025\t2004\t"
        "A\tA1\t320:A\t280:A1\t"
        "12.5\t8.0\t16.0\t11.0\t5.2"
    )
    row2 = (
        "sample2\t400\t0\t370\t0\t"
        "400\t340\tgi|333031|lcl|HPV16REF.1|\t"
        "310\t30\t"
        "40\t28\t115\t75\t22\t85\t88\t"
        "2812\t2946\t1816\t2104\t2699\t1841\t1786\t"
        "D\tD2\t310:D\t270:D2\t"
        "11.8\t7.5\t15.5\t10.5\t4.9"
    )
    with open(path, "w") as f:
        f.write(header + "\n")
        f.write(row1 + "\n")
        f.write(row2 + "\n")


def main():
    print("Generating mock test data for Qualimap HPV16 QC pipeline...")

    samples = ["sample1", "sample2"]
    for sample in samples:
        reads = simulate_reads(HPV16_LEN, NUM_PAIRS, READ_LEN, INSERT_SIZE)
        bam_path = os.path.join(BAM_DIR, f"{sample}.bam")
        create_bam(bam_path, sample, HPV16_REF, HPV16_LEN, reads, READ_LEN)
        print(f"  Created BAM: {bam_path} ({NUM_PAIRS} read pairs)")

    bed_path = os.path.join(OUTPUT_DIR, "hpv16_gene_annot.bed")
    write_bed(bed_path)
    print(f"  Created BED: {bed_path}")

    summary_path = os.path.join(EXAMPLE_DIR, "summary_table.txt")
    write_example_summary(summary_path)
    print(f"  Created example output: {summary_path}")

    print(f"\nTest data: {OUTPUT_DIR}")
    print(f"BAM files: {BAM_DIR}")
    print(f"Example output: {EXAMPLE_DIR}")


if __name__ == "__main__":
    main()
