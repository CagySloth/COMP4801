# tests/test_pipeline_integration.py

import os, json, subprocess, sys

def test_full_pipeline_end_to_end(tmp_path):
    """Run a mini end-to-end pipeline: simulate, phase with one algorithm, evaluate accuracy."""
    # Parameters for a tiny simulation
    outprefix = tmp_path / "mini"
    ploidy = 2
    num_variants = 50
    num_reads = 100
    # 1. Simulate data
    sim_cmd = [
        sys.executable, "dataset/simulate.py",
        "--ploidy", str(ploidy),
        "--num-variants", str(num_variants),
        "--num-reads", str(num_reads),
        "--read-length", "20",
        "--error-rate", "0.05",
        "--output", str(outprefix)
    ]
    subprocess.run(sim_cmd, check=True)
    # Ensure simulation outputs exist
    hap_file = f"{outprefix}.haplotypes.tsv"
    reads_tsv = f"{outprefix}.reads.sparse.tsv"
    assert os.path.exists(hap_file) and os.path.exists(reads_tsv), "Simulation failed to produce expected files."
    # 2. Convert TSV to NPZ
    conv_cmd = [
        sys.executable, "-m", "algorithms.cli.convert",
        str(reads_tsv), "--output", f"{outprefix}.reads.npz"
    ]
    subprocess.run(conv_cmd, check=True)
    npz_file = f"{outprefix}.reads.npz"
    assert os.path.exists(npz_file), "Conversion to NPZ failed."
    # 3. Phase with diploid EM (for example)
    phase_cmd = [
        sys.executable, "-m", "algorithms.cli.phase", "diploid-em",
        "-i", str(npz_file), "--output-prefix", str(outprefix)+".phased"
    ]
    subprocess.run(phase_cmd, check=True)
    phased_hap = f"{outprefix}.phased.haplotypes.tsv"
    phased_assign = f"{outprefix}.phased.assignments.tsv"
    assert os.path.exists(phased_hap) and os.path.exists(phased_assign), "Phasing algorithm did not produce outputs."
    # 4. Evaluate accuracy
    acc_json = f"{outprefix}.phased.accuracy.json"
    acc_cmd = [
        sys.executable, "benchmark/benchmark_accuracy.py",
        "--truth", str(hap_file),
        "--pred", str(phased_hap),
        "--output", str(acc_json)
    ]
    subprocess.run(acc_cmd, check=True)
    assert os.path.exists(acc_json), "Accuracy evaluation did not produce JSON output."
    # Load accuracy result and basic sanity check
    with open(acc_json) as f:
        acc_data = json.load(f)
    assert 0.5 <= acc_data["accuracy"] <= 1.0, "Accuracy out of expected range (should be between 0.5 and 1.0 for diploid)."
    # Print accuracy for informational purpose (not needed for test to pass)
    print(f"Accuracy for mini pipeline: {acc_data['accuracy']:.2f}")
