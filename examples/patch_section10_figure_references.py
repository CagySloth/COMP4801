from pathlib import Path
import re
import shutil

MD_PATH = Path("report/10_experiments.md")
BACKUP_PATH = Path("report/10_experiments.md.before_reference_patch")

if not MD_PATH.exists():
    raise FileNotFoundError(MD_PATH)

text = MD_PATH.read_text(encoding="utf-8")

if not BACKUP_PATH.exists():
    shutil.copy2(MD_PATH, BACKUP_PATH)
    print(f"Backup created: {BACKUP_PATH}")

# ------------------------------------------------------------
# 1) Remove all previously auto-inserted support sentences
# ------------------------------------------------------------
text = re.sub(
    r"\nThese (summary numbers|trends) (?:are|is) (?:visually consistent with|reflected by).*?\n",
    "\n",
    text,
)

# ------------------------------------------------------------
# 2) Fix subsection numbering in 10.3 so it matches figure files
# ------------------------------------------------------------
text = text.replace(
    "#### 10.3.3 Secondary stressors: correlated bursts and read length model",
    "#### 10.3.3 Correlated error bursts",
)
text = text.replace("\n##### Correlated error bursts\n", "\n")
text = text.replace("\n##### Read length model\n", "\n#### 10.3.4 Read length model\n")
text = text.replace(
    "#### 10.3.4 Truth indels and SNP-only policy",
    "#### 10.3.5 Truth indels and SNP-only policy",
)

# ------------------------------------------------------------
# 3) Fix stale displayed figure number in DP-scaling caption
# ------------------------------------------------------------
text = text.replace(
    "**Figure 10.5.3 — Called solve time vs SNP count under default and balanced (optimized) configurations.**",
    "**Figure 10.5.5.1 — Called solve time vs SNP count under default and balanced (optimized) configurations.**",
)

# ------------------------------------------------------------
# 4) Helper to insert support sentence before Key observations
# ------------------------------------------------------------
def insert_before_key_obs(section_heading: str, sentence: str):
    global text
    pattern = re.compile(
        rf"({re.escape(section_heading)}.*?#### Results summary.*?)(\n#### Key observations)",
        re.DOTALL
    )
    m = pattern.search(text)
    if m and sentence not in m.group(1):
        replacement = m.group(1).rstrip() + "\n\n" + sentence + m.group(2)
        text = text[:m.start()] + replacement + text[m.end():]

# ------------------------------------------------------------
# 5) Add accurate support sentences to results-summary blocks
# ------------------------------------------------------------
insert_before_key_obs(
    "### 10.2 Baseline attribution study: depth scaling",
    "These trends are reflected by Figures 10.2.1–10.2.6, especially Figure 10.2.1 for the depth-limited rise in calling recall, Figure 10.2.3 for the large oracle-vs-called gap at low depth, and Figures 10.2.4–10.2.5 for the reduction in fragmentation with increasing coverage."
)

insert_before_key_obs(
    "#### 10.3.1 Duplicated regions",
    "These trends are reflected by Figures 10.3.1.1–10.3.1.5, especially Figure 10.3.1.2 for the drop in called effective phased recall only at the strongest duplication setting and Figure 10.3.1.3 for the corresponding switch-error increase."
)

insert_before_key_obs(
    "#### 10.3.2 Coverage dropout",
    "These trends are reflected by Figures 10.3.2.1–10.3.2.5, especially Figure 10.3.2.4 for the steep decline in called effective phased recall and Figure 10.3.2.1 for the strong rise in oracle phase-set count."
)

insert_before_key_obs(
    "#### 10.3.5 Truth indels and SNP-only policy",
    "These trends are reflected by Figures 10.3.5.1–10.3.5.4, which show that the main result is evaluation robustness rather than a monotonic degradation trend."
)

insert_before_key_obs(
    "### 10.4 Interaction study: duplicated regions × coverage dropout",
    "These trends are reflected by Figures 10.4.1–10.4.4, especially Figure 10.4.2 for the compounded drop in called effective phased recall and Figure 10.4.3 for the dropout-driven fragmentation pattern."
)

insert_before_key_obs(
    "#### 10.5.3 Representative configuration comparison",
    "These trends are reflected by Figures 10.5.3.1–10.5.3.4, and are reinforced by the final default-vs-balanced confirmation in Figures 10.5.3.5–10.5.3.8."
)

insert_before_key_obs(
    "#### 10.5.4 Robustness across scenarios",
    "These trends are reflected by Figures 10.5.4.1–10.5.4.4, where the tuned configurations remain above the default across all four scenarios."
)

insert_before_key_obs(
    "#### 10.5.5 Runtime breakdown and DP-scaling interpretation",
    "These trends are reflected by Figures 10.5.5.1–10.5.5.4, especially Figure 10.5.5.1 for the steeper default solve-time growth and Figure 10.5.5.4 for the corresponding fragmentation difference."
)

# ------------------------------------------------------------
# 6) Add support lines to burst / read-length mini-blocks
# ------------------------------------------------------------
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

# ------------------------------------------------------------
# 7) Add figure links inside 10.5.2 optimization screening
# ------------------------------------------------------------
replacements = {
    "A lower-coverage sweep shows that performance plateaus by `max_coverage = 8–10`; `max_coverage = 4` is slightly too aggressive, while `15` is unnecessarily expensive.":
        "A lower-coverage sweep shows that performance plateaus by `max_coverage = 8–10`; `max_coverage = 4` is slightly too aggressive, while `15` is unnecessarily expensive.\n\nThese trends are reflected by Figures 10.5.2.1–10.5.2.4 for the coarse sweep and Figures 10.5.2.13–10.5.2.16 for the lower-coverage refinement, which together show a clear performance plateau with increasing runtime at higher retained coverage.",

    "A fine sweep confirms that performance is effectively flat for `min_baseq = 0–15` and degrades only at `20`.":
        "A fine sweep confirms that performance is effectively flat for `min_baseq = 0–15` and degrades only at `20`.\n\nThese trends are reflected by Figures 10.5.2.5–10.5.2.8 for the coarse phasing-threshold grid and Figures 10.5.2.17–10.5.2.20 for the fine `min_baseq` sweep, which together show that performance is stable across moderate settings and degrades only when filtering becomes too strict.",

    "`call_min_mapq` has only weak effects and does not provide meaningful optimization headroom.":
        "`call_min_mapq` has only weak effects and does not provide meaningful optimization headroom.\n\nThese trends are reflected by Figures 10.5.2.9–10.5.2.12 for the caller-threshold grid and Figures 10.5.2.25–10.5.2.28 for the `call_min_mapq` rule-out sweep, showing that caller base-quality is the stronger optimization lever.",

    "The decisive local sensitivity remains caller base-quality; phasing base-quality is effectively flat once overly strict filtering has been avoided.":
        "The decisive local sensitivity remains caller base-quality; phasing base-quality is effectively flat once overly strict filtering has been avoided.\n\nThese trends are reflected by Figures 10.5.2.21–10.5.2.24, which show that the recommended threshold pair lies on a stable local performance plateau.",
}

for old, new in replacements.items():
    text = text.replace(old, new)

MD_PATH.write_text(text, encoding="utf-8")
print(f"Patched: {MD_PATH}")