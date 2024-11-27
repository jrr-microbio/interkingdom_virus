import os
import subprocess
from ete3 import Tree

def concatenate_files(fasta1, fasta2, output_fasta):
    """
    Concatenate two FASTA files into one.
    """
    with open(output_fasta, 'wb') as outfile:
        for f in [fasta1, fasta2]:
            if os.path.exists(f):
                with open(f, 'rb') as infile:
                    outfile.write(infile.read())
    print(f"Concatenated {fasta1} and {fasta2} into {output_fasta}")

def run_muscle(input_fasta, output_fasta):
    """
    Run MUSCLE on the input_fasta file and save the aligned sequences to output_fasta file.
    """
    cmd = f"muscle -super5 {input_fasta} -output {output_fasta}"
    subprocess.run(cmd, shell=True, check=True)
    print(f"MUSCLE has aligned the sequences and saved to: {output_fasta}")

def run_trimal(input_fasta, output_fasta, html_output):
    """
    Run trimAl on the input_fasta file and save the trimmed sequences to output_fasta file.
    Also generate an HTML report.
    """
    cmd = f"trimal -in {input_fasta} -out {output_fasta} -automated1 -htmlout {html_output}"
    subprocess.run(cmd, shell=True, check=True)
    print(f"trimAl has trimmed the sequences and saved to: {output_fasta} with HTML report: {html_output}")

def run_fasttree(input_fasta, output_tree):
    """
    Run FastTree on the input_fasta file and save the tree to output_tree file.
    """
    cmd = f"fasttree {input_fasta} > {output_tree}"
    subprocess.run(cmd, shell=True, check=True)
    print(f"FastTree has created the initial tree: {output_tree}")

def reroot_tree(input_tree, output_tree, outgroup_name):
    """
    Reroot the tree using the specified outgroup and save it to the output_tree file.
    """
    # Load the initial tree
    tree = Tree(input_tree)

    # Check if the outgroup is in the tree
    if outgroup_name in tree:
        tree.set_outgroup(tree & outgroup_name)
    else:
        raise ValueError(f"Outgroup {outgroup_name} not found in the tree")

    # Save the rooted tree
    tree.write(outfile=output_tree, format=1)
    print(f"Tree successfully rerooted with the outgroup and saved to: {output_tree}")

def extract_outgroup_id(outgroup_fasta):
    """
    Extract the outgroup ID from the first sequence in the outgroup FASTA file.
    """
    with open(outgroup_fasta, 'r') as file:
        for line in file:
            if line.startswith(">"):
                return line[1:].strip().split()[0]  # Get the sequence ID, excluding '>'
    raise ValueError(f"No valid sequence ID found in {outgroup_fasta}")

def find_nearest_neighbors(tree_file, sm10_prefix, output_file):
    """
    Find the nearest neighbor for sequences prefixed with `sm10_prefix` in the tree.

    Args:
    - tree_file (str): Path to the input tree file.
    - sm10_prefix (str): Prefix to identify target sequences (e.g., "SM10").
    - output_file (str): Path to the output file to save nearest neighbors.
    """
    # Load the tree
    tree = Tree(tree_file)

    # Separate sequences into those prefixed with `sm10_prefix` and others
    SM10_seqs = [leaf.name for leaf in tree if leaf.name.startswith(sm10_prefix)]
    ref_seqs = [leaf.name for leaf in tree if not leaf.name.startswith(sm10_prefix)]

    # Dictionary to keep the closest reference neighbors
    closest_ref_neighbors = {}

    # Get the distance matrix function
    dist_matrix = tree.get_distance

    # Find nearest neighbors
    for sm10_seq in SM10_seqs:
        sm10_node = tree&sm10_seq
        ref_distances = []

        for ref_seq in ref_seqs:
            ref_node = tree&ref_seq
            distance = dist_matrix(sm10_node, ref_node)
            ref_distances.append((ref_seq, distance))

        # Find the closest reference based on distance
        if ref_distances:
            closest_ref_neighbor = min(ref_distances, key=lambda x: x[1])[0]
            closest_ref_neighbors[sm10_seq] = closest_ref_neighbor

    # Write the results to the output file
    with open(output_file, "w") as file:
        for sequence, neighbor in closest_ref_neighbors.items():
            file.write(f"{sequence}: {neighbor}\n")

def process_directory(directory):
    # Location of the outgroup file directory
    outgroup_dir = os.path.join(directory, "_outgroup_faas")

    # Iterate through each file in the directory
    for filename in os.listdir(directory):
        if filename.endswith("_concatenated.faa"):
            base_name = filename.split("_concatenated.faa")[0]
            taxon_name, taxon_level = base_name.split("_")

            input_fasta = os.path.join(directory, filename)
            outgroup_file = os.path.join(outgroup_dir, f"{base_name}_root.faa")

            if not os.path.exists(outgroup_file):
                print(f"Outgroup file not found for {base_name}, skipping.")
                continue
            
            concatenated_fasta_with_outgroup = os.path.join(directory, f"{base_name}_with_outgroup.faa")
            aligned_fasta = os.path.join(directory, f"{base_name}_aligned.faa")
            trimmed_fasta = os.path.join(directory, f"{base_name}_trimmed.faa")
            initial_tree = os.path.join(directory, f"{base_name}_initial_tree.nwk")
            rooted_tree = os.path.join(directory, f"{base_name}_rooted_tree.nwk")
            html_output = os.path.join(directory, f"{base_name}_trim_report.html")
            nearest_neighbors_output = os.path.join(directory, f"{base_name}_nearest_neighbors.txt")

            # Step 1: Concatenate the input sequences with the outgroup sequences
            concatenate_files(input_fasta, outgroup_file, concatenated_fasta_with_outgroup)

            # Step 2: Run MUSCLE to align sequences
            run_muscle(concatenated_fasta_with_outgroup, aligned_fasta)

            # Step 3: Run trimAl to trim the alignment and generate HTML report
            run_trimal(aligned_fasta, trimmed_fasta, html_output)

            # Step 4: Run FastTree to generate the initial tree
            run_fasttree(trimmed_fasta, initial_tree)

            # Step 5: Extract the outgroup ID from the outgroup FASTA file
            try:
                outgroup_name = extract_outgroup_id(outgroup_file)
            except ValueError as e:
                print(e)
                continue

            # Step 6: Reroot the tree with the specified outgroup
            try:
                reroot_tree(initial_tree, rooted_tree, outgroup_name)
            except ValueError as e:
                print(e)

            # Step 7: Find and save nearest neighbors
            find_nearest_neighbors(rooted_tree, "SM10", nearest_neighbors_output)

if __name__ == "__main__":
    # Set your working directory path
    directory = "./"  # Current directory

    # Run the function
    process_directory(directory)
