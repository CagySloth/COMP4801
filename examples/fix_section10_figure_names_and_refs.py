from pathlib import Path
import re
import shutil

ROOT = Path(".")
FIG_DIR = ROOT / "report" / "figures" / "10_experiments"
MD_PATH = ROOT / "report" / "10_experiments.md"
BACKUP_PATH = ROOT / "report" / "10_experiments.md.before_renumber_backup"

if not FIG_DIR.exists():
    raise FileNotFoundError(f"Missing figure directory: {FIG_DIR}")
if not MD_PATH.exists():
    raise FileNotFoundError(f"Missing markdown file: {MD_PATH}")

if not BACKUP_PATH.exists():
    shutil.copy2(MD_PATH, BACKUP_PATH)
    print(f"Created backup: {BACKUP_PATH}")

def new_name(name: str) -> str:
    # Keep baseline numbering as-is because Section 10.2 still maps cleanly.
    if name.startswith("fig_10_2_") or name == "fig_10_2_combined_oracle_vs_called.png":
        return name

    # Old 10.3 -> new 10.3.1
    m = re.match(r"fig_10_3_(\d+)_(.+)\.png$", name)
    if m:
        return f"fig_10_3_1_{m.group(1)}_{m.group(2)}.png"

    # Old 10.4 -> new 10.3.2
    m = re.match(r"fig_10_4_(\d+)_(.+)\.png$", name)
    if m:
        return f"fig_10_3_2_{m.group(1)}_{m.group(2)}.png"

    # Old 10.5 -> new 10.3.3
    m = re.match(r"fig_10_5_(\d+)_(.+)\.png$", name)
    if m:
        return f"fig_10_3_3_{m.group(1)}_{m.group(2)}.png"

    # Old 10.6 -> new 10.3.4
    m = re.match(r"fig_10_6_(\d+)_(.+)\.png$", name)
    if m:
        return f"fig_10_3_4_{m.group(1)}_{m.group(2)}.png"

    # Old 10.7 -> new 10.3.5
    m = re.match(r"fig_10_7_(\d+)_(.+)\.png$", name)
    if m:
        return f"fig_10_3_5_{m.group(1)}_{m.group(2)}.png"

    # Old 10.8 -> new 10.4
    m = re.match(r"fig_10_8_(\d+)_(.+)\.png$", name)
    if m:
        return f"fig_10_4_{m.group(1)}_{m.group(2)}.png"

    # Optimization screening figures -> new 10.5.2.*
    # Coarse max_coverage: old 10.9.1 -> new 10.5.2.1-4
    m = re.match(r"fig_10_9_1_(\d+)_(.+)\.png$", name)
    if m:
        return f"fig_10_5_2_{int(m.group(1))}_{m.group(2)}.png"

    # Phasing quality thresholds: old 10.9.2 -> new 10.5.2.5-8
    m = re.match(r"fig_10_9_2_(\d+)_(.+)\.png$", name)
    if m:
        return f"fig_10_5_2_{4 + int(m.group(1))}_{m.group(2)}.png"

    # Caller threshold grid: old 10.9.3 -> new 10.5.2.9-12
    m = re.match(r"fig_10_9_3_(\d+)_(.+)\.png$", name)
    if m:
        return f"fig_10_5_2_{8 + int(m.group(1))}_{m.group(2)}.png"

    # Runtime-focused lower max_coverage: old 10.9.4 -> new 10.5.2.13-16
    m = re.match(r"fig_10_9_4_(\d+)_(.+)\.png$", name)
    if m:
        return f"fig_10_5_2_{12 + int(m.group(1))}_{m.group(2)}.png"

    # Fine phasing min_baseq: old 10.9.5 -> new 10.5.2.17-20
    m = re.match(r"fig_10_9_5_(\d+)_(.+)\.png$", name)
    if m:
        return f"fig_10_5_2_{16 + int(m.group(1))}_{m.group(2)}.png"

    # Local joint search: old 10.9.7 -> new 10.5.2.21-24
    m = re.match(r"fig_10_9_7_(\d+)_(.+)\.png$", name)
    if m:
        return f"fig_10_5_2_{20 + int(m.group(1))}_{m.group(2)}.png"

    # Caller mapq rule-out: old 10.9.10 -> new 10.5.2.25-28
    m = re.match(r"fig_10_9_10_(\d+)_(.+)\.png$", name)
    if m:
        return f"fig_10_5_2_{24 + int(m.group(1))}_{m.group(2)}.png"

    # Representative configuration comparison -> new 10.5.3.*
    # Frontier comparison: old 10.9.8 -> new 10.5.3.1-4
    m = re.match(r"fig_10_9_8_(\d+)_(.+)\.png$", name)
    if m:
        return f"fig_10_5_3_{int(m.group(1))}_{m.group(2)}.png"

    # Default vs optimized confirmation: old 10.9.6 -> new 10.5.3.5-8
    m = re.match(r"fig_10_9_6_(\d+)_(.+)\.png$", name)
    if m:
        return f"fig_10_5_3_{4 + int(m.group(1))}_{m.group(2)}.png"

    # Robustness across scenarios -> new 10.5.4.*
    m = re.match(r"fig_10_9_9_(\d+)_(.+)\.png$", name)
    if m:
        return f"fig_10_5_4_{int(m.group(1))}_{m.group(2)}.png"

    # DP scaling -> new 10.5.5.*
    m = re.match(r"fig_10_9_11_(\d+)_(.+)\.png$", name)
    if m:
        return f"fig_10_5_5_{int(m.group(1))}_{m.group(2)}.png"

    return name

# ------------------------------------------------------------------
# Rename figure files
# ------------------------------------------------------------------
renames = {}
for path in sorted(FIG_DIR.glob("*.png")):
    nn = new_name(path.name)
    if nn != path.name:
        renames[path.name] = nn

for old, new in renames.items():
    src = FIG_DIR / old
    dst = FIG_DIR / new
    if dst.exists():
        print(f"Target exists, skipping rename: {dst.name}")
    else:
        src.rename(dst)
        print(f"Renamed: {old} -> {new}")

# ------------------------------------------------------------------
# Patch markdown paths and figure numbers
# ------------------------------------------------------------------
text = MD_PATH.read_text(encoding="utf-8")

for old, new in renames.items():
    text = text.replace(
        f"figures/10_experiments/{old}",
        f"figures/10_experiments/{new}"
    )

# Fix displayed caption numbering so it matches the new structure.
caption_replacements = {
    "**Figure 10.2.1 — Called effective phased recall vs sequencing depth (end-to-end performance).**":
        "**Figure 10.2.3 — Called effective phased recall vs sequencing depth (end-to-end performance).**",

    "**Figure 10.3.1 — Called effective phased recall vs dropout fraction.**":
        "**Figure 10.3.2.4 — Called effective phased recall vs dropout fraction.**",

    "**Figure 10.4.1 — Called effective phased recall across duplication × dropout conditions.**":
        "**Figure 10.4.2 — Called effective phased recall across duplication × dropout conditions.**",

    "**Figure 10.5.1 — Called effective phased recall across representative configurations.**":
        "**Figure 10.5.3.1 — Called effective phased recall across representative configurations.**",

    "**Figure 10.5.2 — Called effective phased recall across scenarios and representative configurations.**":
        "**Figure 10.5.4.1 — Called effective phased recall across scenarios and representative configurations.**",

    "**Figure 10.5.3 — Called solve time vs SNP count under default and optimized configurations.**":
        "**Figure 10.5.5.1 — Called solve time vs SNP count under default and balanced (optimized) configurations.**",
}
for old, new in caption_replacements.items():
    text = text.replace(old, new)

# Standardize wording in the runtime table / DP-scaling subsection.
text = text.replace("| Optimized |", "| Balanced (optimized) |")
text = text.replace("| Speed-oriented |", "| Runtime (speed-oriented) |")
text = text.replace(
    "default and optimized configurations",
    "default and balanced (optimized) configurations",
)

# "Local optimum" is too strong; use "stable local plateau".
text = text.replace("stable local optimum", "stable local plateau")

# ------------------------------------------------------------------
# Add figure-link sentences after results summaries / support blocks
# ------------------------------------------------------------------
summary_insertions = {
    "### 10.2 Baseline attribution study: depth scaling":
        "These summary numbers are reflected by Figures 10.2.1–10.2.6, especially Figure 10.2.3 for the oracle-vs-called gap and Figure 10.2.1 for the depth-limited rise in calling recall.",

    "#### 10.3.1 Duplicated regions":
        "These summary numbers are reflected by Figures 10.3.1.1–10.3.1.5, especially Figure 10.3.1.2 for the drop in called effective phased recall only at the strongest duplication setting and Figure 10.3.1.3 for the corresponding switch-error increase.",

    "#### 10.3.2 Coverage dropout":
        "These summary numbers are reflected by Figures 10.3.2.1–10.3.2.5, especially Figure 10.3.2.4 for the steep decline in called effective phased recall and Figure 10.3.2.1 for the strong rise in oracle phase-set count.",

    "#### 10.3.4 Truth indels and SNP-only policy":
        "These summary numbers are reflected by Figures 10.3.5.1–10.3.5.4, which show that the main effect is evaluation robustness rather than a systematic performance trend.",

    "### 10.4 Interaction study: duplicated regions × coverage dropout":
        "These summary numbers are reflected by Figures 10.4.1–10.4.4, especially Figure 10.4.2 for the compounded drop in called effective phased recall and Figure 10.4.3 for the dropout-driven fragmentation pattern.",

    "#### 10.5.3 Representative configuration comparison":
        "These summary numbers are reflected by Figures 10.5.3.1–10.5.3.4, with Figure 10.5.3.1 showing the small but consistent advantage of the balanced and runtime configurations.",

    "#### 10.5.4 Robustness across scenarios":
        "These summary numbers are reflected by Figures 10.5.4.1–10.5.4.4, where the tuned configurations remain above the default across all four scenarios.",

    "#### 10.5.5 Runtime breakdown and DP-scaling interpretation":
        "These summary numbers are reflected by Figures 10.5.5.1–10.5.5.4, especially Figure 10.5.5.1 for the steeper default solve-time growth and Figure 10.5.5.4 for the corresponding fragmentation difference.",
}

for heading, sentence in summary_insertions.items():
    pattern = re.compile(
        rf"({re.escape(heading)}.*?#### Results summary.*?)(\n#### Key observations)",
        re.DOTALL
    )
    m = pattern.search(text)
    if m and sentence not in m.group(1):
        replacement = m.group(1).rstrip() + "\n\n" + sentence + m.group(2)
        text = text[:m.start()] + replacement + text[m.end():]

# Add support lines for 10.3.3 burst and read-length mini-blocks
burst_marker = (
    "Interpretation: under the current parameterization, correlated error bursts are a weak and somewhat noisy called-regime correctness stressor rather than a strong practical failure mode."
)
burst_sentence = (
    "\n\nThese trends are reflected by Figures 10.3.3.1–10.3.3.5, especially Figures 10.3.3.1 and 10.3.3.2, where the intermediate-burst dip is visible but non-monotonic."
)
if burst_marker in text and burst_sentence.strip() not in text:
    text = text.replace(burst_marker, burst_marker + burst_sentence)

len_marker = (
    "Interpretation: the read length model is informative for mechanism analysis, but under the present parameterization it is a mixed continuity / correctness trade-off rather than a clear optimization lever."
)
len_sentence = (
    "\n\nThese trends are reflected by Figures 10.3.4.1–10.3.4.5, especially Figure 10.3.4.2 for the continuity gain and Figure 10.3.4.4 for the small net change in called effective phased recall."
)
if len_marker in text and len_sentence.strip() not in text:
    text = text.replace(len_marker, len_marker + len_sentence)

MD_PATH.write_text(text, encoding="utf-8")
print(f"Patched markdown: {MD_PATH}")
print("Done.")