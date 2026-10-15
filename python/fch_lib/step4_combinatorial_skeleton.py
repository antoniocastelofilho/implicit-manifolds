from fch_wrapper import FCH_CPU_Runner

runner = FCH_CPU_Runner("./lib/")

# Run the combinatorial skeleton to generate the .pol file (text file).
runner.skeleton(
    in_path="out/results.bin", 
    out_path="out/results.pol"
    )

# Run the combinatorial skeleton with simplicial decomposition to generate the .spx file (text file).
runner.skeleton_simplicial_decomp(
    in_path="out/results.bin", 
    out_path="out/results.spx"
    )

# Run the combinatorial skeleton with simplicial decomposition to generate the .spx_bin file (binary file).
runner.skeleton_simplicial_decomp(
    in_path="out/results.bin", 
    out_path="out/results.spx_bin",
    is_output_binary=True
    )
