import os
import subprocess
import argparse
import pandas as pd
import numpy as np

# command line argument gathering
def parse_args():
    parser = argparse.ArgumentParser(description="Create MSA of gene lines, find systematic differences and potential conserved motifs, generate logo plot")
    
    # Define mutually exclusive group for input options
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument("-i", "--input", help="Input FASTA file containing all sequences including the reference")
    input_group.add_argument("-d", "--input_dir", help="Input directory containing multiple FASTA files for each gene line")
    input_group.add_argument("--msa", type=str, help="Path to pre-computed MSA file in FASTA format")

    parser.add_argument("-r", "--reference", required=True, help="Reference sequence ID present in the input FASTA or MSA file")
    parser.add_argument("-o", "--output_dir", required=True, help="Output directory for the alignment results")
    parser.add_argument("--header_delim", type=str, default=" ", help="Delimiter to use for splitting the headers (default: space)")
    parser.add_argument("--element_pos", type=int, default=0, help="Position of the element indicating the gene line (1-based, default: last)")
    parser.add_argument("--generate_logos", action='store_true', help="Generate sequence logos and comparative logos")
    parser.add_argument("--generate_barplot", action='store_true', help="Generate stacked bar plots for each position using Plotly")
    parser.add_argument("--generate_motifs", action='store_true', help="Identify motifs and highlight them on sequence logos")
    parser.add_argument("--cons_cutoff", type=float, default=0.8, help="Conservation threshold for motif identification (default: 0.8)")
    parser.add_argument("--motif_len", type=int, default=4, help="minimum length of identified conservation-motifs (default: 4)")
    parser.add_argument("--num_aa", type=int, help="Number of amino acids to display in the logos (default: all)")
    parser.add_argument("--rank_differences", action='store_true', help="Rank sequence differences between gene lines. Results will also be indicate in logo plot")
    parser.add_argument("--diff_cutoff", type=float, default=1.0, help="Difference score cutoff for filtering positions (default: 1.0)")
    parser.add_argument("--mirror_yax", action='store_true', help="Mirror y-axis in comparative logo plot")
    parser.add_argument("--muscle", action='store_true', help="Use MUSCLE for multiple sequence alignment instead of MAFFT")
    
    # Add arguments for specifying binary paths
    parser.add_argument("--mafft_path", type=str, help="Path to the MAFFT binary")
    parser.add_argument("--muscle_path", type=str, help="Path to the MUSCLE binary")
    
    return parser.parse_args()


# if dir with different gene line fasta files is given, combine them for following steps, mark gene lines 
def combine_fastas(input_dir, combined_output_file):
    try:
        with open(combined_output_file, "w") as outfile:
            for filename in os.listdir(input_dir):
                if filename.endswith((".fa", ".fasta")):
                    file_path = os.path.join(input_dir, filename)
                    gene_line = os.path.splitext(filename)[0]
                    with open(file_path, "r") as infile:
                        for line in infile:
                            if line.startswith(">"):
                                outfile.write(f"{line.strip()} {gene_line}\n")
                            else:
                                outfile.write(line)
        print(f"Combined FASTA file created at {combined_output_file}")
        return combined_output_file
    except Exception as e:
        print(f"Error combining FASTA files: {e}")
        return None

# check whether reference file for indexing is in given input fasta
def check_reference(input_file, reference_id):
    sequences = {}
    reference_seq = None
    reference_header = None

    try:
        print(f"Reading sequences from {input_file}")

        with open(input_file, "r") as infile:
            header, seq = None, []
            for line in infile:
                if line.startswith(">"):
                    if header:
                        sequences[header] = ''.join(seq)
                        if reference_id in header:
                            reference_header = header
                            reference_seq = ''.join(seq)
                    header = line.strip()
                    seq = []
                else:
                    seq.append(line.strip())
            if header:
                sequences[header] = ''.join(seq)
                if reference_id in header:
                    reference_header = header
                    reference_seq = ''.join(seq)

        if not reference_seq:
            raise ValueError(f"Reference sequence '{reference_id}' not found in input file.")

        print(f"Reference sequence '{reference_header}' found.")
        return reference_header, reference_seq
    except Exception as e:
        print(f"Error reading sequences: {e}")
        return None, None

# check if alignment tool is installed
def check_tool_installed(tool):
    """Check if a tool is installed and available in PATH."""
    return subprocess.call(f"type {tool}", shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL) == 0

# perform multiple sequence alignment with specified tool
def perform_msa(input_file, output_file, use_muscle=False, mafft_path=None, muscle_path=None):
    try:
        if use_muscle:
            tool = muscle_path if muscle_path else "muscle"
            cmd = [tool, "-align", input_file, "-output", output_file]
        else:
            tool = mafft_path if mafft_path else "mafft"
            cmd = [tool, "--auto", input_file]
            with open(output_file, "w") as out_f:
                print(f"Running MAFFT on {input_file}...")
                subprocess.run(cmd, stdout=out_f, check=True)
                print(f"Multiple sequence alignment completed. Results saved to {output_file}")
                return

        if check_tool_installed(tool):
            print(f"Running {tool.upper()} on {input_file}...")
            subprocess.run(cmd, check=True)
            print(f"Multiple sequence alignment completed. Results saved to {output_file}")
        else:
            raise FileNotFoundError(f"{tool.upper()} is not installed. Please install {tool.upper()} or use the other tool.")
    except Exception as e:
        print(f"Error during MSA: {e}")

# read alignment file into dict
def read_alignment(alignment_file, delimiter, element_pos):
    sequences = {}  # Use a regular dict
    try:
        with open(alignment_file, "r") as infile:
            header, seq = None, []
            for line in infile:
                if line.startswith(">"):
                    if header:
                        gene_line = parse_gene_line(header, delimiter, element_pos)
                        if gene_line not in sequences:
                            sequences[gene_line] = {}  # Handle missing keys
                        sequences[gene_line][header] = ''.join(seq)
                    header = line.strip()
                    seq = []
                else:
                    seq.append(line.strip())
            if header:
                gene_line = parse_gene_line(header, delimiter, element_pos)
                if gene_line not in sequences:
                    sequences[gene_line] = {}  # Handle missing keys
                sequences[gene_line][header] = ''.join(seq)
        print(f"Read aligned sequences from {alignment_file}")
        return sequences
    except Exception as e:
        print(f"Error reading alignment: {e}")
        return None


# extract part of header that indicates gene line
def parse_gene_line(header, delimiter, element_pos):
    try:
        # 1-based index provided by user, convert to 0-based for Python
        element_index = element_pos - 1
        desc_parts = header[1:].split(delimiter)
        gene_line = desc_parts[element_index]
        return gene_line
    except IndexError:
        raise ValueError(f"Header does not have enough parts to extract the element at position {element_pos}.")
        
        
# trim msa with respect to reference seq
def trim_msa(sequences, reference_header, reference_seq):
    non_gap_positions = [i for i, aa in enumerate(reference_seq) if aa != '-']
    
    trimmed_sequences = {}
    for gene_line, seqs in sequences.items():
        trimmed_sequences[gene_line] = {}  # Initialize sub-dict
        for header, seq in seqs.items():
            trimmed_seq = ''.join([seq[i] for i in non_gap_positions])
            trimmed_sequences[gene_line][header] = trimmed_seq
    trimmed_reference_seq = ''.join([reference_seq[i] for i in non_gap_positions])
    
    print("Trimmed aligned sequences")
    return trimmed_sequences, trimmed_reference_seq

# generate consensus sequence for each gene line
def generate_consensus(sequences, threshold=0.8):
    consensus_sequences = {}
    for gene_line, seqs in sequences.items():
        if not seqs:
            continue
        seq_len = len(next(iter(seqs.values())))
        consensus = []
        for i in range(seq_len):
            column = [seq[i] for seq in seqs.values()]

            count_dict = {}
            for aa in column:
                if aa in count_dict:
                    count_dict[aa] += 1
                else:
                    count_dict[aa] = 1
            
            most_common_aa, count = max(count_dict.items(), key=lambda item: item[1])

            if len(count_dict) == 1:
                consensus.append(most_common_aa.upper())  
            elif count / len(column) >= threshold:
                consensus.append(most_common_aa.upper())  
            else:
                max_count = max(count_dict.values())
                alternatives = [aa for aa, cnt in count_dict.items() if cnt == max_count]
                if len(alternatives) > 1:
                    consensus.append('X')  
                else:
                    consensus.append(most_common_aa.lower())  
        consensus_sequences[gene_line] = ''.join(consensus)
    print("Generated consensus sequences")
    return consensus_sequences

# save consensus as preliminary output file
def save_consensus(output_file, reference_header, reference_seq, consensus_sequences):
    try:
        with open(output_file, "w") as out_handle:
            out_handle.write(f"{reference_header}\n")
            for i in range(0, len(reference_seq), 60):
                out_handle.write(f"{reference_seq[i:i+60]}\n")

            for gene_line, seq in consensus_sequences.items():
                out_handle.write(f">consensus_seq_{gene_line}\n")
                for i in range(0, len(seq), 60):
                    out_handle.write(f"{seq[i:i+60]}\n")
        print(f"Consensus sequences saved to {output_file}")
    except Exception as e:
        print(f"Error saving consensus sequences: {e}")

# filter top rank of amino acids for later plotting steps
def filter_top_aa(counts_df, num_aa):
    if num_aa is None or num_aa >= counts_df.shape[1]:
        return counts_df

    def top_n_filter(row):
        sorted_row = row.sort_values(ascending=False)
        if len(sorted_row) > num_aa:
            threshold = sorted_row.iloc[num_aa - 1]
            filtered_row = row.apply(lambda x: x if x >= threshold else 0)
        else:
            filtered_row = row
        return filtered_row
    
    return counts_df.apply(top_n_filter, axis=1)
    
# calculate sequence identity between two groups using consensus sequences
def calculate_consensus_identity(sequences, gene_line1, gene_line2):
    consensus1 = generate_consensus({gene_line1: sequences[gene_line1]})
    consensus2 = generate_consensus({gene_line2: sequences[gene_line2]})
    consensus_seq1 = consensus1[gene_line1]
    consensus_seq2 = consensus2[gene_line2]
    identity = calculate_pairwise_identity(consensus_seq1, consensus_seq2)
    return identity

# calculate sequence identity between two sequences
def calculate_pairwise_identity(seq1, seq2):
    identical_positions = 0
    valid_positions = 0
    for a, b in zip(seq1, seq2):
        if a != '-' and b != '-':
            valid_positions += 1
            if a == b:
                identical_positions += 1
    if valid_positions > 0:
        identity = identical_positions / valid_positions
    else:
        identity = 0
    return identity * 100  # Convert to percentage

# decide which blosum matrix to use based on sequence identity
def get_blosum_matrix(identity):
    if identity >= 80:
        print("Using BLOSUM80 matrix for high identity sequences.")
        return BLOSUM80
    elif identity >= 62:
        print("Using BLOSUM62 matrix for moderate identity sequences.")
        return BLOSUM62
    else:
        print("Using BLOSUM45 matrix for low identity sequences.")
        return BLOSUM45

# rank positions based on 1: difference score 2: Blosum score
def rank_sequence_differences(sequences, output_dir, diff_cutoff=1.0):
    try:
        gene_lines = list(sequences.keys())
        if len(gene_lines) < 2:
            raise ValueError("Need at least two gene lines to compare differences.")

        num_positions = len(next(iter(next(iter(sequences.values())).values())))
        difference_scores = {f"{gene_lines[i]}_vs_{gene_lines[j]}": np.zeros(num_positions) for i in range(len(gene_lines) - 1) for j in range(i + 1, len(gene_lines))}
        blosum_scores = {f"{gene_lines[i]}_vs_{gene_lines[j]}": [None] * num_positions for i in range(len(gene_lines) - 1) for j in range(i + 1, len(gene_lines))}
        filtered_positions_all = {}

        print("Calculating frequency matrices for each gene line...")

        # Compute frequency matrices for each gene line, including gaps
        frequency_matrices = {}
        for gene_line, seqs in sequences.items():
            count_matrix = {aa: np.zeros(num_positions) for aa in "ACDEFGHIKLMNPQRSTVWY-"}
            for seq in seqs.values():
                for pos, aa in enumerate(seq):
                    count_matrix[aa][pos] += 1
            total_counts = np.sum(list(count_matrix.values()), axis=0)
            frequency_matrices[gene_line] = {}
            for aa in count_matrix:
                with np.errstate(divide='ignore', invalid='ignore'):
                    frequencies = np.true_divide(count_matrix[aa], total_counts)
                    frequencies[~np.isfinite(frequencies)] = 0
                frequency_matrices[gene_line][aa] = frequencies
            if '-' in frequency_matrices[gene_line]:
                del frequency_matrices[gene_line]['-']

        print("Calculating difference scores and BLOSUM scores for tie-breaking...")
        for i in range(len(gene_lines) - 1):
            for j in range(i + 1, len(gene_lines)):
                gene_line1 = gene_lines[i]
                gene_line2 = gene_lines[j]

                rep_seq1 = next(iter(sequences[gene_line1].values()))
                rep_seq2 = next(iter(sequences[gene_line2].values()))
                identity = calculate_pairwise_identity(rep_seq1, rep_seq2)
                blosum_matrix = get_blosum_matrix(identity)

                print(f"Calculated sequence identity for {gene_line1} vs {gene_line2}: {identity:.2f}%")

                most_prevalent_aa1_list = []
                most_prevalent_aa2_list = []
                most_prevalent_aa1_single_list = []
                most_prevalent_aa2_single_list = []

                for pos in range(num_positions):
                    score = 0
                    # Calculate the difference score excluding gaps (since we deleted them)
                    for aa in "ACDEFGHIKLMNPQRSTVWY":
                        freq1 = frequency_matrices[gene_line1][aa][pos]
                        freq2 = frequency_matrices[gene_line2][aa][pos]
                        score += abs(freq1 - freq2) * max(freq1, freq2)
                    difference_scores[f"{gene_line1}_vs_{gene_line2}"][pos] = score

                    # Determine the most prevalent AAs with ties considered
                    aa_freq1 = {aa: frequency_matrices[gene_line1][aa][pos] for aa in "ACDEFGHIKLMNPQRSTVWY"}
                    aa_freq2 = {aa: frequency_matrices[gene_line2][aa][pos] for aa in "ACDEFGHIKLMNPQRSTVWY"}

                    max_freq1 = max(aa_freq1.values())
                    max_freq2 = max(aa_freq2.values())

                    if max_freq1 == 0:
                        most_prevalent_aa1_all = ["-"]
                        most_prevalent_aa1_single = "-"
                    else:
                        most_prevalent_aa1_all = [aa for aa, freq in aa_freq1.items() if freq == max_freq1]
                        most_prevalent_aa1_single = most_prevalent_aa1_all[0]  # Use the first one for BLOSUM

                    if max_freq2 == 0:
                        most_prevalent_aa2_all = ["-"]
                        most_prevalent_aa2_single = "-"
                    else:
                        most_prevalent_aa2_all = [aa for aa, freq in aa_freq2.items() if freq == max_freq2]
                        most_prevalent_aa2_single = most_prevalent_aa2_all[0]  # Use the first one for BLOSUM

                    most_prevalent_aa1_list.append(','.join(most_prevalent_aa1_all))
                    most_prevalent_aa2_list.append(','.join(most_prevalent_aa2_all))
                    most_prevalent_aa1_single_list.append(most_prevalent_aa1_single)
                    most_prevalent_aa2_single_list.append(most_prevalent_aa2_single)

                    if most_prevalent_aa1_single == "-" or most_prevalent_aa2_single == "-":
                        blosum_scores[f"{gene_line1}_vs_{gene_line2}"][pos] = None  # Indicate gaps
                    else:
                        blosum_scores[f"{gene_line1}_vs_{gene_line2}"][pos] = blosum_matrix[most_prevalent_aa1_single].get(most_prevalent_aa2_single, 0)

                # Prepare data for sorting and exporting
                positions = np.arange(1, num_positions + 1)
                scores = difference_scores[f"{gene_line1}_vs_{gene_line2}"]
                blosum_ties = blosum_scores[f"{gene_line1}_vs_{gene_line2}"]

                # Add most prevalent AAs to the DataFrame
                data = {
                    'Position': positions,
                    'Diff_Score': scores,
                    'BLOSUM_Score': blosum_ties,
                    f'Most_Prevalent_{gene_line1}': most_prevalent_aa1_list,
                    f'Most_Prevalent_{gene_line2}': most_prevalent_aa2_list
                }
                df = pd.DataFrame(data)

                df['BLOSUM_Score'] = df['BLOSUM_Score'].fillna(-float('inf'))  # Adjust for sorting

                sorted_df = df.sort_values(by=['Diff_Score', 'BLOSUM_Score'], ascending=[False, True])
                sorted_output_file = os.path.join(output_dir, f"{gene_line1}_vs_{gene_line2}_scores_sorted.csv")
                sorted_df.to_csv(sorted_output_file, index=False)
                print(f"Sorted scores for {gene_line1} vs {gene_line2} saved to {sorted_output_file}")

                filtered_df = sorted_df[sorted_df['Diff_Score'] >= diff_cutoff]
                filtered_positions = set(filtered_df['Position'])
                filtered_positions_all[f"{gene_line1}_vs_{gene_line2}"] = filtered_positions
                print(f"Filtered positions for {gene_line1} vs {gene_line2} (cutoff={diff_cutoff}) identified.")

        print("Ranking completed.")
        return filtered_positions_all

    except Exception as e:
        print(f"Error in ranking sequence differences: {e}")
        return None


# identify conserved parts of gene (pseudomotifs)
def identify_motifs(sequences, threshold, reference_len, output_dir, motif_min_len):
    for gene_line, seqs in sequences.items():
        motifs = []
        seq_count = len(seqs)
        position_freq = [{aa: 0 for aa in "ACDEFGHIKLMNPQRSTVWY"} for _ in range(reference_len)]

        for seq in seqs.values():
            for pos, aa in enumerate(seq):
                if aa != '-':
                    position_freq[pos][aa] += 1

        for pos in range(reference_len):
            for aa in position_freq[pos]:
                position_freq[pos][aa] /= seq_count

        current_pos = 0
        while current_pos < reference_len:
            dominant_aa, freq = max(position_freq[current_pos].items(), key=lambda x: x[1])
            if freq >= threshold:
                motif_start = current_pos
                motif = [dominant_aa]
                current_pos += 1
                motif_valid = True
                while current_pos < reference_len:
                    next_dominant_aa, next_freq = max(position_freq[current_pos].items(), key=lambda x: x[1])
                    if next_freq >= threshold:
                        motif.append(next_dominant_aa)
                        current_pos += 1
                    else:
                        break
                if current_pos - motif_start >= motif_min_len and motif_valid:
                    motif_seq = ''.join(motif)
                    motifs.append([motif_start + 1, current_pos, motif_seq])
            else:
                current_pos += 1

        motifs_df = pd.DataFrame(motifs, columns=['Start', 'End', 'Motif'])
        motifs_file = os.path.join(output_dir, f"identified_motifs_{gene_line}.csv")
        motifs_df.to_csv(motifs_file, index=False)
        print(f"Motifs identified for {gene_line} and saved to {motifs_file}")

# load positions above cutoff for logo plot creation
def load_significant_positions(output_dir, gene_line1, gene_line2):
    significant_positions_file = os.path.join(output_dir, f"{gene_line1}_vs_{gene_line2}_scores_filtered.csv")
    significant_positions = set()
    if os.path.isfile(significant_positions_file):
        significant_positions_df = pd.read_csv(significant_positions_file)
        significant_positions = set(significant_positions_df['Position'])
    return significant_positions

# create logo plots for each sequence comparison, if difference score was calculated indicate pairs with large diff score
def prepare_comparative_logos(sequences, output_dir, num_aa=None, mirror_yax=False, generate_motifs=False, filtered_positions_dict=None):
    import logomaker as lm
    import matplotlib.pyplot as plt
    try:
        gene_lines = list(sequences.keys())
        for i in range(len(gene_lines) - 1):
            for j in range(i + 1, len(gene_lines)):
                gene_line1 = gene_lines[i]
                gene_line2 = gene_lines[j]

                aligned_seqs1 = [seq.upper() for seq in sequences[gene_line1].values()]
                aligned_seqs2 = [seq.upper() for seq in sequences[gene_line2].values()]

                counts_df1 = lm.alignment_to_matrix(sequences=aligned_seqs1, to_type='counts')
                counts_df2 = lm.alignment_to_matrix(sequences=aligned_seqs2, to_type='counts')

                if num_aa:
                    counts_df1 = filter_top_aa(counts_df1, num_aa)
                    counts_df2 = filter_top_aa(counts_df2, num_aa)

                total_seqs1 = len(aligned_seqs1)
                total_seqs2 = len(aligned_seqs2)

                counts_df1 = (counts_df1 / total_seqs1) * 100
                counts_df2 = (counts_df2 / total_seqs2) * 100

                counts_df1 = counts_df1.replace([np.inf, -np.inf, np.nan], 0)
                counts_df2 = counts_df2.replace([np.inf, -np.inf, np.nan], 0)

                fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(200, 10), sharex=True, gridspec_kw={'height_ratios': [1, 1]})
                fig.subplots_adjust(hspace=0.05)

                logo1 = lm.Logo(counts_df1, ax=ax1, baseline_width=0.0, stack_order='small_on_top', color_scheme='chemistry')
                ax1.set_title(f"Comparative Logos: {gene_line1} vs {gene_line2}")
                ax1.set_ylabel(f"Frequency {gene_line1} (%)")

                if mirror_yax:
                    logo2 = lm.Logo(counts_df2, ax=ax2, baseline_width=0.0, stack_order='small_on_top', color_scheme='chemistry')
                    ax2.invert_yaxis()
                    ax2.set_ylabel(f"Frequency {gene_line2} (%)")
                else:
                    logo2 = lm.Logo(counts_df2, ax=ax2, baseline_width=0.0, stack_order='big_on_top', color_scheme='chemistry')
                    ax2.set_ylabel(f"Frequency {gene_line2} (%)")

                # Use filtered positions specific to each comparison
                filtered_positions = filtered_positions_dict.get(f"{gene_line1}_vs_{gene_line2}", set())

                total_positions = counts_df1.shape[0]
                for pos in range(total_positions):
                    if pos + 1 not in filtered_positions:
                        for aa in "ACDEFGHIKLMNPQRSTVWY":
                            logo1.style_single_glyph(p=pos, c=aa, color='grey')
                            logo2.style_single_glyph(p=pos, c=aa, color='grey')

                if generate_motifs:
                    motif_file1 = os.path.join(output_dir, f"identified_motifs_{gene_line1}.csv")
                    if os.path.isfile(motif_file1):
                        motifs_df1 = pd.read_csv(motif_file1)
                        for _, row in motifs_df1.iterrows():
                            start = row['Start'] - 1
                            end = row['End'] - 1
                            logo1.highlight_position_range(start, end, color='yellow', alpha=0.3)

                            rect1 = plt.Rectangle(
                                (start - 0.5, ax1.get_ylim()[0]),
                                end - start + 1,
                                ax1.get_ylim()[1] - ax1.get_ylim()[0],
                                linewidth=3, edgecolor='black', facecolor='none'
                            )
                            ax1.add_patch(rect1)

                    motif_file2 = os.path.join(output_dir, f"identified_motifs_{gene_line2}.csv")
                    if os.path.isfile(motif_file2):
                        motifs_df2 = pd.read_csv(motif_file2)
                        for _, row in motifs_df2.iterrows():
                            start = row['Start'] - 1
                            end = row['End'] - 1
                            if mirror_yax:
                                logo2.highlight_position_range(start, end, floor=ax2.get_ylim()[1], ceiling=ax2.get_ylim()[0], color='yellow', alpha=0.3)
                            else:
                                logo2.highlight_position_range(start, end, floor=ax2.get_ylim()[0], ceiling=ax2.get_ylim()[1], color='yellow', alpha=0.3)

                            rect2 = plt.Rectangle(
                                (start - 0.5, ax2.get_ylim()[0]),
                                end - start + 1,
                                ax2.get_ylim()[1] - ax2.get_ylim()[0] if mirror_yax else ax2.get_ylim()[0] - ax2.get_ylim()[1],
                                linewidth=3, edgecolor='black', facecolor='none'
                            )
                            ax2.add_patch(rect2)

                num_ticks = counts_df1.shape[0]
                ax2.set_xticks(range(num_ticks))
                ax2.set_xticklabels(range(1, num_ticks + 1))

                ax1.set_xticks(range(num_ticks))

                logo_output_file = os.path.join(output_dir, f"logo_comparison_{gene_line1}_vs_{gene_line2}.png")
                fig.tight_layout()
                fig.savefig(logo_output_file)
                plt.close(fig)
                print(f"Comparative sequence logo saved for {gene_line1} vs {gene_line2} to {logo_output_file}")

    except Exception as e:
        print(f"Error preparing comparative Logomaker data with highlights: {e}")

# generate stacked barplot with plotly
def generate_plotly_stack_barplot(sequences, output_dir, num_aa=None, mirror_yax=False, threshold=None, differing_positions=None):
    import plotly.graph_objects as go
    import plotly.express as px
    import numpy as np

    gene_lines = list(sequences.keys())
    palette = px.colors.qualitative.Dark24

    for i in range(len(gene_lines) - 1):
        for j in range(i + 1, len(gene_lines)):
            gene_line1 = gene_lines[i]
            gene_line2 = gene_lines[j]
            comparison_key = f"{gene_line1}_vs_{gene_line2}"

            aligned_seqs1 = [seq.upper() for seq in sequences[gene_line1].values()]
            aligned_seqs2 = [seq.upper() for seq in sequences[gene_line2].values()]

            aa_list = "ACDEFGHIKLMNPQRSTVWY"

            # Manually calculate counts for each position
            counts_df1 = {aa: np.zeros(len(aligned_seqs1[0])) for aa in aa_list}
            counts_df2 = {aa: np.zeros(len(aligned_seqs2[0])) for aa in aa_list}

            for seq in aligned_seqs1:
                for pos, aa in enumerate(seq):
                    if aa in aa_list:
                        counts_df1[aa][pos] += 1

            for seq in aligned_seqs2:
                for pos, aa in enumerate(seq):
                    if aa in aa_list:
                        counts_df2[aa][pos] += 1

            total_seqs1 = len(aligned_seqs1)
            total_seqs2 = len(aligned_seqs2)

            counts_df1 = {aa: (counts / total_seqs1) * 100 for aa, counts in counts_df1.items()}
            counts_df2 = {aa: (counts / total_seqs2) * 100 for aa, counts in counts_df2.items()}

            total_positions = len(aligned_seqs1[0])

            fig = go.Figure()

            for pos in range(total_positions):
                aa_freq1 = {aa: counts_df1[aa][pos] for aa in aa_list}
                aa_freq2 = {aa: counts_df2[aa][pos] for aa in aa_list}

                sorted_aa1 = sorted(aa_freq1.items(), key=lambda item: item[1], reverse=True)
                sorted_aa2 = sorted(aa_freq2.items(), key=lambda item: item[1], reverse=True)

                for aa, freq in sorted_aa1:
                    fig.add_trace(go.Bar(
                        x=[pos + 1],
                        y=[freq],
                        name=f'{gene_line1} {aa}',
                        marker_color=palette[aa_list.index(aa)],
                        hovertemplate=f'Pos.: {pos + 1}, Freq.: {freq:.2f}, Gene_line: {gene_line1}, AA: {aa}<extra></extra>'
                    ))

                for aa, freq in sorted_aa2:
                    fig.add_trace(go.Bar(
                        x=[pos + 1],
                        y=[-freq],
                        name=f'{gene_line2} {aa}',
                        marker_color=palette[aa_list.index(aa)],
                        hovertemplate=f'Pos.: {pos + 1}, Freq.: {abs(freq):.2f}, Gene_line: {gene_line2}, AA: {aa}<extra></extra>'
                    ))

            yaxis_tickvals = list(range(-100, 101, 20))
            yaxis_ticktext = [str(abs(v)) for v in yaxis_tickvals]

            fig.update_layout(
                barmode='relative',
                title=f"Stacked Barplot: {gene_line1} vs {gene_line2}",
                yaxis_title='Frequency (%)',
                xaxis_title='Position',
                yaxis=dict(
                    title='Frequency (%)',
                    tickmode='array',
                    tickvals=yaxis_tickvals,
                    ticktext=yaxis_ticktext
                )
            )

            # Add semi-transparent white overlay for positions below the threshold
            if differing_positions is not None and comparison_key in differing_positions:
                differing_pos_set = differing_positions[comparison_key]
                for pos in range(total_positions):
                    if pos + 1 not in differing_pos_set:
                        fig.add_shape(
                            type="rect",
                            x0=pos + 0.5, x1=pos + 1.5,
                            y0=-100, y1=100,
                            fillcolor="white", opacity=0.8,
                            line=dict(color="white", width=0)
                        )

            plot_output_file = os.path.join(output_dir, f"stacked_barplot_{gene_line1}_vs_{gene_line2}.html")
            fig.write_html(plot_output_file)
            print(f"Stacked bar plot saved for {gene_line1} vs {gene_line2} to {plot_output_file}")


# handle calls to all subfunctions
def main():
    args = parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    sequences = None
    alignment_file = None

    if args.msa:
        print("Using pre-computed MSA file...")
        alignment_file = args.msa

        print("Reading alignment file...")
        sequences = read_alignment(alignment_file, args.header_delim, args.element_pos)
        if not sequences:
            print("Failed to read alignment file. Exiting.")
            return

        # Check reference sequence in the provided MSA
        print("Checking for reference sequence in MSA...")
        reference_header, reference_seq = check_reference(alignment_file, args.reference)
        if not reference_header:
            print("Failed to find reference sequence in MSA. Exiting.")
            return

    else:
        input_file = args.input
        if args.input_dir:
            combined_fasta = os.path.join(args.output_dir, "combined_sequences.fasta")
            input_file = combine_fastas(args.input_dir, combined_fasta)
            if not input_file:
                print("Failed to combine FASTA files. Exiting.")
                return

        print("Checking for reference sequence in input file...")
        reference_header, reference_seq = check_reference(input_file, args.reference)
        if not reference_header:
            print("Failed to find reference sequence in input file. Exiting.")
            return

        alignment_file = os.path.join(args.output_dir, "aligned_sequences.fasta")

        print("Performing multiple sequence alignment...")
        perform_msa(input_file, alignment_file, use_muscle=args.muscle, mafft_path=args.mafft_path, muscle_path=args.muscle_path)

        print("Reading alignment file...")
        sequences = read_alignment(alignment_file, args.header_delim, args.element_pos)
        if not sequences:
            print("Failed to read alignment file. Exiting.")
            return

        reference_seq = None
        for gene_line, seqs in sequences.items():
            for header, seq in seqs.items():
                if args.reference in header:
                    reference_seq = seq
                    break
            if reference_seq:
                break

        if not reference_seq:
            print("Failed to find aligned reference sequence. Exiting.")
            return

    print("Trimming aligned sequences...")
    trimmed_sequences, trimmed_reference_seq = trim_msa(sequences, reference_header, reference_seq)

    print("Generating consensus sequences...")
    consensus_sequences = generate_consensus(trimmed_sequences)

    consensus_output_file = os.path.join(args.output_dir, "consensus_sequences.fasta")
    save_consensus(consensus_output_file, reference_header, trimmed_reference_seq, consensus_sequences)

    if args.generate_motifs:
        print("Identifying motifs...")
        reference_len = len(trimmed_reference_seq)
        identify_motifs(trimmed_sequences, args.cons_cutoff, reference_len, args.output_dir, args.motif_len)

    filtered_positions_all = None
    if args.rank_differences:
        print("Ranking sequence differences between gene lines...")
        filtered_positions_all = rank_sequence_differences(trimmed_sequences, args.output_dir, args.diff_cutoff)

    if args.generate_logos:
        print("Creating comparative logo plots...")
        prepare_comparative_logos(trimmed_sequences, args.output_dir, args.num_aa, args.mirror_yax, args.generate_motifs, filtered_positions_all)

    if args.generate_barplot:
        print("Creating stacked bar plots with Plotly...")
        differing_positions = filtered_positions_all if filtered_positions_all else set()
        print(differing_positions)
        generate_plotly_stack_barplot(trimmed_sequences, args.output_dir, args.num_aa, args.mirror_yax, args.diff_cutoff, differing_positions)


# blosum matrices for secondary diff ranking 
BLOSUM45 = {
    'A': {'A':  5, 'C': -1, 'D': -2, 'E': -1, 'F': -2, 'G':  0, 'H': -2, 'I': -1, 'K': -1, 'L': -2, 'M': -1, 'N': -1, 'P': -1, 'Q': -1, 'R': -2, 'S':  1, 'T':  0, 'V':  0, 'W': -3, 'Y': -2},
    'C': {'A': -1, 'C': 12, 'D': -3, 'E': -2, 'F': -2, 'G': -3, 'H': -3, 'I': -2, 'K': -3, 'L': -2, 'M': -2, 'N': -2, 'P': -4, 'Q': -3, 'R': -4, 'S':  0, 'T': -2, 'V': -2, 'W': -5, 'Y': -3},
    'D': {'A': -2, 'C': -3, 'D':  7, 'E':  2, 'F': -4, 'G': -1, 'H': -1, 'I': -3, 'K':  0, 'L': -4, 'M': -3, 'N':  2, 'P': -1, 'Q':  0, 'R': -1, 'S':  0, 'T': -1, 'V': -3, 'W': -4, 'Y': -2},
    'E': {'A': -1, 'C': -2, 'D':  2, 'E':  6, 'F': -3, 'G': -2, 'H':  0, 'I': -2, 'K':  1, 'L': -3, 'M': -2, 'N':  0, 'P': -1, 'Q':  2, 'R':  0, 'S': -1, 'T': -1, 'V': -3, 'W': -3, 'Y': -2},
    'F': {'A': -2, 'C': -2, 'D': -4, 'E': -3, 'F':  9, 'G': -4, 'H': -2, 'I':  0, 'K': -3, 'L':  0, 'M':  0, 'N': -3, 'P': -4, 'Q': -3, 'R': -3, 'S': -3, 'T': -2, 'V': -1, 'W':  1, 'Y':  5},
    'G': {'A':  0, 'C': -3, 'D': -1, 'E': -2, 'F': -4, 'G':  6, 'H': -2, 'I': -3, 'K': -2, 'L': -3, 'M': -3, 'N':  0, 'P': -2, 'Q': -2, 'R': -3, 'S':  0, 'T': -2, 'V': -3, 'W': -2, 'Y': -3},
    'H': {'A': -2, 'C': -3, 'D': -1, 'E':  0, 'F': -2, 'G': -2, 'H':  8, 'I': -3, 'K': -1, 'L': -2, 'M': -2, 'N':  1, 'P': -2, 'Q':  0, 'R':  0, 'S': -1, 'T': -2, 'V': -3, 'W': -2, 'Y':  2},
    'I': {'A': -1, 'C': -2, 'D': -3, 'E': -2, 'F':  0, 'G': -3, 'H': -3, 'I':  7, 'K': -3, 'L':  2, 'M':  1, 'N': -3, 'P': -3, 'Q': -2, 'R': -3, 'S': -2, 'T':  0, 'V':  4, 'W': -3, 'Y': -1},
    'K': {'A': -1, 'C': -3, 'D':  0, 'E':  1, 'F': -3, 'G': -2, 'H': -1, 'I': -3, 'K':  6, 'L': -2, 'M': -1, 'N':  0, 'P': -1, 'Q':  1, 'R':  3, 'S':  0, 'T': -1, 'V': -2, 'W': -3, 'Y': -2},
    'L': {'A': -2, 'C': -2, 'D': -4, 'E': -3, 'F':  0, 'G': -3, 'H': -2, 'I':  2, 'K': -2, 'L':  6, 'M':  2, 'N': -3, 'P': -3, 'Q': -2, 'R': -3, 'S': -3, 'T': -1, 'V':  2, 'W': -2, 'Y': -1},
    'M': {'A': -1, 'C': -2, 'D': -3, 'E': -2, 'F':  0, 'G': -3, 'H': -2, 'I':  1, 'K': -1, 'L':  2, 'M':  7, 'N': -2, 'P': -2, 'Q':  0, 'R': -2, 'S': -2, 'T': -1, 'V':  1, 'W': -1, 'Y': -1},
    'N': {'A': -1, 'C': -2, 'D':  2, 'E':  0, 'F': -3, 'G':  0, 'H':  1, 'I': -3, 'K':  0, 'L': -3, 'M': -2, 'N':  7, 'P': -2, 'Q':  0, 'R':  0, 'S':  1, 'T':  0, 'V': -3, 'W': -4, 'Y': -2},
    'P': {'A': -1, 'C': -4, 'D': -1, 'E': -1, 'F': -4, 'G': -2, 'H': -2, 'I': -3, 'K': -1, 'L': -3, 'M': -2, 'N': -2, 'P': 10, 'Q': -1, 'R': -3, 'S': -1, 'T': -1, 'V': -3, 'W': -4, 'Y': -3},
    'Q': {'A': -1, 'C': -3, 'D':  0, 'E':  2, 'F': -3, 'G': -2, 'H':  0, 'I': -2, 'K':  1, 'L': -2, 'M':  0, 'N':  0, 'P': -1, 'Q':  7, 'R':  1, 'S': -1, 'T': -1, 'V': -2, 'W': -2, 'Y': -1},
    'R': {'A': -2, 'C': -4, 'D': -1, 'E':  0, 'F': -3, 'G': -3, 'H':  0, 'I': -3, 'K':  3, 'L': -3, 'M': -2, 'N':  0, 'P': -3, 'Q':  1, 'R':  7, 'S': -1, 'T': -1, 'V': -3, 'W': -3, 'Y': -1},
    'S': {'A':  1, 'C':  0, 'D':  0, 'E': -1, 'F': -3, 'G':  0, 'H': -1, 'I': -2, 'K':  0, 'L': -3, 'M': -2, 'N':  1, 'P': -1, 'Q': -1, 'R': -1, 'S':  4, 'T':  2, 'V': -2, 'W': -4, 'Y': -2},
    'T': {'A':  0, 'C': -2, 'D': -1, 'E': -1, 'F': -2, 'G': -2, 'H': -2, 'I':  0, 'K': -1, 'L': -1, 'M': -1, 'N':  0, 'P': -1, 'Q': -1, 'R': -1, 'S':  2, 'T':  5, 'V':  0, 'W': -3, 'Y': -1},
    'V': {'A':  0, 'C': -2, 'D': -3, 'E': -3, 'F': -1, 'G': -3, 'H': -3, 'I':  4, 'K': -2, 'L':  2, 'M':  1, 'N': -3, 'P': -3, 'Q': -2, 'R': -3, 'S': -2, 'T':  0, 'V':  5, 'W': -3, 'Y': -1},
    'W': {'A': -3, 'C': -5, 'D': -4, 'E': -3, 'F':  1, 'G': -2, 'H': -2, 'I': -3, 'K': -3, 'L': -2, 'M': -1, 'N': -4, 'P': -4, 'Q': -2, 'R': -3, 'S': -4, 'T': -3, 'V': -3, 'W': 15, 'Y':  2},
    'Y': {'A': -2, 'C': -3, 'D': -2, 'E': -2, 'F':  5, 'G': -3, 'H':  2, 'I': -1, 'K': -2, 'L': -1, 'M': -1, 'N': -2, 'P': -3, 'Q': -1, 'R': -1, 'S': -2, 'T': -1, 'V': -1, 'W':  2, 'Y':  8},
}

BLOSUM62 = {
    'A': {'A':  4, 'C':  0, 'D': -2, 'E': -1, 'F': -2, 'G':  0, 'H': -2, 'I': -1, 'K': -1, 'L': -1, 'M': -1, 'N': -2, 'P': -1, 'Q': -1, 'R': -1, 'S':  1, 'T':  0, 'V':  0, 'W': -3, 'Y': -2},
    'C': {'A':  0, 'C':  9, 'D': -3, 'E': -4, 'F': -2, 'G': -3, 'H': -3, 'I': -1, 'K': -3, 'L': -1, 'M': -1, 'N': -3, 'P': -3, 'Q': -3, 'R': -3, 'S': -1, 'T': -1, 'V': -1, 'W': -2, 'Y': -2},
    'D': {'A': -2, 'C': -3, 'D':  6, 'E':  2, 'F': -3, 'G': -1, 'H': -1, 'I': -3, 'K': -1, 'L': -4, 'M': -3, 'N':  1, 'P': -1, 'Q':  0, 'R': -2, 'S':  0, 'T': -1, 'V': -3, 'W': -4, 'Y': -3},
    'E': {'A': -1, 'C': -4, 'D':  2, 'E':  5, 'F': -3, 'G': -2, 'H':  0, 'I': -3, 'K':  1, 'L': -2, 'M': -2, 'N':  0, 'P': -1, 'Q':  2, 'R':  0, 'S':  0, 'T': -1, 'V': -2, 'W': -3, 'Y': -2},
    'F': {'A': -2, 'C': -2, 'D': -3, 'E': -3, 'F':  6, 'G': -3, 'H': -1, 'I':  0, 'K': -3, 'L':  0, 'M':  0, 'N': -3, 'P': -4, 'Q': -3, 'R': -3, 'S': -2, 'T': -2, 'V': -1, 'W':  1, 'Y':  3},
    'G': {'A':  0, 'C': -3, 'D': -1, 'E': -2, 'F': -3, 'G':  6, 'H': -2, 'I': -4, 'K': -2, 'L': -4, 'M': -3, 'N':  0, 'P': -2, 'Q': -2, 'R': -2, 'S':  0, 'T': -2, 'V': -3, 'W': -2, 'Y': -3},
    'H': {'A': -2, 'C': -3, 'D': -1, 'E':  0, 'F': -1, 'G': -2, 'H':  8, 'I': -3, 'K': -1, 'L': -3, 'M': -2, 'N':  1, 'P': -2, 'Q':  0, 'R':  0, 'S': -1, 'T': -2, 'V': -3, 'W': -2, 'Y':  2},
    'I': {'A': -1, 'C': -1, 'D': -3, 'E': -3, 'F':  0, 'G': -4, 'H': -3, 'I':  4, 'K': -3, 'L':  2, 'M':  1, 'N': -3, 'P': -3, 'Q': -3, 'R': -3, 'S': -2, 'T': -1, 'V':  3, 'W': -3, 'Y': -1},
    'K': {'A': -1, 'C': -3, 'D': -1, 'E':  1, 'F': -3, 'G': -2, 'H': -1, 'I': -3, 'K':  5, 'L': -2, 'M': -1, 'N':  0, 'P': -1, 'Q':  1, 'R':  2, 'S':  0, 'T': -1, 'V': -2, 'W': -3, 'Y': -2},
    'L': {'A': -1, 'C': -1, 'D': -4, 'E': -2, 'F':  0, 'G': -4, 'H': -3, 'I':  2, 'K': -2, 'L':  4, 'M':  2, 'N': -3, 'P': -3, 'Q': -2, 'R': -2, 'S': -2, 'T': -1, 'V':  1, 'W': -2, 'Y': -1},
    'M': {'A': -1, 'C': -1, 'D': -3, 'E': -2, 'F':  0, 'G': -3, 'H': -2, 'I':  1, 'K': -1, 'L':  2, 'M':  5, 'N': -2, 'P': -2, 'Q':  0, 'R': -1, 'S': -1, 'T': -1, 'V':  1, 'W': -1, 'Y': -1},
    'N': {'A': -2, 'C': -3, 'D':  1, 'E':  0, 'F': -3, 'G':  0, 'H':  1, 'I': -3, 'K':  0, 'L': -3, 'M': -2, 'N':  6, 'P': -2, 'Q':  0, 'R':  0, 'S':  1, 'T':  0, 'V': -3, 'W': -4, 'Y': -2},
    'P': {'A': -1, 'C': -3, 'D': -1, 'E': -1, 'F': -4, 'G': -2, 'H': -2, 'I': -3, 'K': -1, 'L': -3, 'M': -2, 'N': -2, 'P':  7, 'Q': -1, 'R': -2, 'S': -1, 'T': -1, 'V': -2, 'W': -4, 'Y': -3},
    'Q': {'A': -1, 'C': -3, 'D':  0, 'E':  2, 'F': -3, 'G': -2, 'H':  0, 'I': -3, 'K':  1, 'L': -2, 'M':  0, 'N':  0, 'P': -1, 'Q':  5, 'R':  1, 'S':  0, 'T': -1, 'V': -2, 'W': -2, 'Y': -1},
    'R': {'A': -1, 'C': -3, 'D': -2, 'E':  0, 'F': -3, 'G': -2, 'H':  0, 'I': -3, 'K':  2, 'L': -2, 'M': -1, 'N':  0, 'P': -2, 'Q':  1, 'R':  5, 'S': -1, 'T': -1, 'V': -3, 'W': -3, 'Y': -2},
    'S': {'A':  1, 'C': -1, 'D':  0, 'E':  0, 'F': -2, 'G':  0, 'H': -1, 'I': -2, 'K':  0, 'L': -2, 'M': -1, 'N':  1, 'P': -1, 'Q':  0, 'R': -1, 'S':  4, 'T':  1, 'V': -2, 'W': -3, 'Y': -2},
    'T': {'A':  0, 'C': -1, 'D': -1, 'E': -1, 'F': -2, 'G': -2, 'H': -2, 'I': -1, 'K': -1, 'L': -1, 'M': -1, 'N':  0, 'P': -1, 'Q': -1, 'R': -1, 'S':  1, 'T':  5, 'V':  0, 'W': -2, 'Y': -2},
    'V': {'A':  0, 'C': -1, 'D': -3, 'E': -2, 'F': -1, 'G': -3, 'H': -3, 'I':  3, 'K': -2, 'L':  1, 'M':  1, 'N': -3, 'P': -2, 'Q': -2, 'R': -3, 'S': -2, 'T':  0, 'V':  4, 'W': -3, 'Y': -1},
    'W': {'A': -3, 'C': -2, 'D': -4, 'E': -3, 'F':  1, 'G': -2, 'H': -2, 'I': -3, 'K': -3, 'L': -2, 'M': -1, 'N': -4, 'P': -4, 'Q': -2, 'R': -3, 'S': -3, 'T': -2, 'V': -3, 'W': 11, 'Y':  2},
    'Y': {'A': -2, 'C': -2, 'D': -3, 'E': -2, 'F':  3, 'G': -3, 'H':  2, 'I': -1, 'K': -2, 'L': -1, 'M': -1, 'N': -2, 'P': -3, 'Q': -1, 'R': -2, 'S': -2, 'T': -2, 'V': -1, 'W':  2, 'Y':  7},
}

BLOSUM80 = {
    'A': {'A':  7, 'C': -1, 'D': -3, 'E': -1, 'F': -3, 'G':  0, 'H': -2, 'I': -2, 'K': -1, 'L': -3, 'M': -2, 'N': -2, 'P': -1, 'Q': -1, 'R': -2, 'S':  1, 'T':  0, 'V':  0, 'W': -5, 'Y': -3},
    'C': {'A': -1, 'C': 10, 'D': -5, 'E': -4, 'F': -3, 'G': -4, 'H': -4, 'I': -2, 'K': -4, 'L': -2, 'M': -2, 'N': -4, 'P': -4, 'Q': -4, 'R': -4, 'S':  0, 'T': -1, 'V': -2, 'W': -5, 'Y': -4},
    'D': {'A': -3, 'C': -5, 'D': 10, 'E':  1, 'F': -6, 'G': -2, 'H': -2, 'I': -4, 'K': -1, 'L': -5, 'M': -4, 'N':  2, 'P': -2, 'Q': -1, 'R': -2, 'S':  0, 'T': -1, 'V': -4, 'W': -7, 'Y': -4},
    'E': {'A': -1, 'C': -4, 'D':  1, 'E':  8, 'F': -5, 'G': -3, 'H':  0, 'I': -3, 'K':  1, 'L': -3, 'M': -2, 'N':  0, 'P': -1, 'Q':  2, 'R': -1, 'S': -1, 'T': -1, 'V': -3, 'W': -5, 'Y': -3},
    'F': {'A': -3, 'C': -3, 'D': -6, 'E': -5, 'F':  8, 'G': -5, 'H': -2, 'I':  0, 'K': -5, 'L': -1, 'M': -1, 'N': -4, 'P': -5, 'Q': -4, 'R': -4, 'S': -3, 'T': -2, 'V': -1, 'W':  1, 'Y':  4},
    'G': {'A':  0, 'C': -4, 'D': -2, 'E': -3, 'F': -5, 'G':  8, 'H': -3, 'I': -5, 'K': -2, 'L': -4, 'M': -3, 'N':  0, 'P': -3, 'Q': -3, 'R': -3, 'S':  0, 'T': -2, 'V': -4, 'W': -3, 'Y': -4},
    'H': {'A': -2, 'C': -4, 'D': -2, 'E':  0, 'F': -2, 'G': -3, 'H': 10, 'I': -4, 'K':  0, 'L': -3, 'M': -2, 'N':  2, 'P': -2, 'Q':  2, 'R':  0, 'S': -2, 'T': -2, 'V': -4, 'W': -3, 'Y':  3},
    'I': {'A': -2, 'C': -2, 'D': -4, 'E': -3, 'F':  0, 'G': -5, 'H': -4, 'I':  7, 'K': -3, 'L':  2, 'M':  2, 'N': -3, 'P': -4, 'Q': -2, 'R': -3, 'S': -3, 'T': -1, 'V':  4, 'W': -4, 'Y': -1},
    'K': {'A': -1, 'C': -4, 'D': -1, 'E':  1, 'F': -5, 'G': -2, 'H':  0, 'I': -3, 'K':  8, 'L': -3, 'M': -1, 'N':  0, 'P': -1, 'Q':  1, 'R':  3, 'S':  0, 'T': -1, 'V': -3, 'W': -4, 'Y': -3},
    'L': {'A': -3, 'C': -2, 'D': -5, 'E': -3, 'F': -1, 'G': -4, 'H': -3, 'I':  2, 'K': -3, 'L':  6, 'M':  2, 'N': -3, 'P': -4, 'Q': -2, 'R': -3, 'S': -3, 'T': -2, 'V':  1, 'W': -3, 'Y': -1},
    'M': {'A': -2, 'C': -2, 'D': -4, 'E': -2, 'F': -1, 'G': -3, 'H': -2, 'I':  2, 'K': -1, 'L':  2, 'M':  10, 'N': -3, 'P': -3, 'Q': -1, 'R': -2, 'S': -2, 'T': -1, 'V':  1, 'W': -1, 'Y': -1},
    'N': {'A': -2, 'C': -4, 'D':  2, 'E':  0, 'F': -4, 'G':  0, 'H':  2, 'I': -3, 'K':  0, 'L': -3, 'M': -3, 'N': 10, 'P': -3, 'Q':  0, 'R':  0, 'S':  1, 'T':  0, 'V': -3, 'W': -5, 'Y': -2},
    'P': {'A': -1, 'C': -4, 'D': -2, 'E': -1, 'F': -5, 'G': -3, 'H': -2, 'I': -4, 'K': -1, 'L': -4, 'M': -3, 'N': -3, 'P': 10, 'Q': -1, 'R': -3, 'S': -1, 'T': -1, 'V': -3, 'W': -5, 'Y': -4},
    'Q': {'A': -1, 'C': -4, 'D': -1, 'E':  2, 'F': -4, 'G': -3, 'H':  2, 'I': -2, 'K':  1, 'L': -2, 'M': -1, 'N':  0, 'P': -1, 'Q':  8, 'R':  1, 'S': -1, 'T': -1, 'V': -2, 'W': -3, 'Y': -2},
    'R': {'A': -2, 'C': -4, 'D': -2, 'E': -1, 'F': -4, 'G': -3, 'H':  0, 'I': -3, 'K':  3, 'L': -3, 'M': -2, 'N':  0, 'P': -3, 'Q':  1, 'R':  8, 'S': -1, 'T': -1, 'V': -3, 'W': -4, 'Y': -3},
    'S': {'A':  1, 'C':  0, 'D':  0, 'E': -1, 'F': -3, 'G':  0, 'H': -2, 'I': -3, 'K':  0, 'L': -3, 'M': -2, 'N':  1, 'P': -1, 'Q': -1, 'R': -1, 'S':  7, 'T':  2, 'V': -2, 'W': -4, 'Y': -3},
    'T': {'A':  0, 'C': -1, 'D': -1, 'E': -1, 'F': -2, 'G': -2, 'H': -2, 'I': -1, 'K': -1, 'L': -2, 'M': -1, 'N':  0, 'P': -1, 'Q': -1, 'R': -1, 'S':  2, 'T':  8, 'V':  0, 'W': -4, 'Y': -2},
    'V': {'A':  0, 'C': -2, 'D': -4, 'E': -3, 'F': -1, 'G': -4, 'H': -4, 'I':  4, 'K': -3, 'L':  1, 'M':  1, 'N': -3, 'P': -3, 'Q': -2, 'R': -3, 'S': -2, 'T':  0, 'V':  7, 'W': -4, 'Y': -1},
    'W': {'A': -5, 'C': -5, 'D': -7, 'E': -5, 'F':  1, 'G': -3, 'H': -3, 'I': -4, 'K': -4, 'L': -3, 'M': -1, 'N': -5, 'P': -5, 'Q': -3, 'R': -4, 'S': -4, 'T': -4, 'V': -4, 'W': 17, 'Y':  2},
    'Y': {'A': -3, 'C': -4, 'D': -4, 'E': -3, 'F':  4, 'G': -4, 'H':  3, 'I': -1, 'K': -3, 'L': -1, 'M': -1, 'N': -2, 'P': -4, 'Q': -2, 'R': -3, 'S': -3, 'T': -2, 'V': -1, 'W':  2, 'Y': 10},
}            

if __name__ == "__main__":
    main()

