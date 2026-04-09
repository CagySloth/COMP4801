## 11. Discussion

### 11.1 Discussion of main findings

This study set out to understand WhatsHap performance under realistic long-read conditions, identify where performance is lost in an end-to-end pipeline, and explore whether practical optimization opportunities exist. The experimental results show that these goals were achieved at three levels: attribution, stress characterization, and optimization.

First, the baseline depth sweep showed a consistent attribution pattern. When correct variant sites are supplied (oracle regime), WhatsHap phasing is already highly accurate, and the main limitation at lower depth is incompleteness due to missing connectivity. In contrast, in the called regime, a large proportion of the end-to-end performance loss is caused by the overlap of correctly called heterozygous variants rather than by the phasing algorithm itself. This establishes an important baseline conclusion: in practical settings, the phaser is not always the dominant bottleneck, and improvements in end-to-end phased recall often depend as much on variant recovery as on phasing quality.

Second, the realism-stressor studies revealed that different perturbations harm the pipeline in qualitatively different ways. Coverage dropout emerged as the clearest and strongest isolated stressor, primarily by fragmenting phase blocks and reducing phasing completeness, and at high severity by also reducing variant recovery. Duplicated regions, in contrast, were relatively mild under the current parameterization, with their main impact appearing only at the highest setting and mostly in the called regime. Correlated error bursts produced a weaker and more variable signal than expected, while the read-length model showed a mixed trade-off: better continuity and overlap, but no substantial end-to-end gain because correctness became noisier. The indel study served mainly as a methodology-validation result, showing that SNP-only safeguards preserve interpretable evaluation when truth indels are present.

Third, the interaction study demonstrated that the most informative failure modes may not come from isolated knobs alone. Duplication combined with dropout produced a substantially worse end-to-end condition than either stressor individually, especially in the called regime. This indicates that realistic phasing difficulty can arise from compounded weaknesses: ambiguous mapping worsens low-coverage conditions by simultaneously reducing callable overlap and weakening phasing completeness. This result was particularly useful because it linked the stress studies to the optimization studies: the key practical weaknesses are mixed, rather than purely phaser-limited.

Finally, the optimization experiments showed that useful but limited optimization headroom exists. The strongest positive signal came from caller-side base-quality tuning: reducing `call_min_baseq` from an overly strict value improved variant recovery, shared heterozygous overlap, and end-to-end phased recall without materially harming precision. On the phasing side, the clearest lesson was not to over-filter: `min_baseq = 20` was consistently harmful, while moderate values such as 10 lay on a stable performance plateau. Retained coverage (`max_coverage`) also showed a practical runtime trade-off, with values around 8–10 preserving performance while avoiding unnecessary compute. In contrast, both caller and phasing mapping-quality thresholds provided little useful headroom. Taken together, these results suggest that practical optimization is not a matter of one magic WhatsHap knob, but rather a combination of preserving useful evidence, recovering more overlapping heterozygous variants upstream, and avoiding overly conservative filtering.

An additional algorithm-facing result came from the SNP-density scaling experiment. As variant density increased, the downstream solve stage became more expensive, but the balanced tuned configuration scaled substantially better than the default configuration. This suggests that tuning does not only improve final phasing output; it also changes the difficulty of the residual wMEC-style problem presented to the solver. Although total runtime in the research pipeline remained dominated by the project-specific readset-construction adaptor, this scaling result still supports the conclusion that practical tuning can improve both output quality and solver-facing problem structure.

Overall, the findings support a balanced interpretation of WhatsHap in realistic long-read pipelines. WhatsHap itself is robust when correct variant evidence is available, but end-to-end performance depends strongly on how much useful evidence survives calling, filtering, and connectivity disruptions. Therefore, optimization opportunities exist, but they are fundamentally pipeline-level rather than purely phaser-internal.

### 11.2 Practical implications for WhatsHap-based phasing

From a practical perspective, the results suggest several concrete recommendations for WhatsHap-based long-read phasing pipelines.

The first implication is that end-to-end phasing quality should not be interpreted using phasing metrics alone. In many conditions tested here, especially at moderate depth and under compounded stress, the dominant loss arose before or alongside phasing, through reduced overlap of correctly called heterozygous variants. This means that practitioners should evaluate both variant recovery and phasing quality together. A strong oracle result does not guarantee equally strong end-to-end performance if the called VCF fails to preserve enough informative heterozygous sites.

The second implication is that conservative filtering can be counterproductive. On the calling side, overly strict base-quality filtering removed too much true variant evidence and substantially reduced downstream phased recall. On the phasing side, strict base-quality filtering similarly weakened connectivity and fragmentation metrics once the threshold became too aggressive. In the current pipeline, moderate settings performed better than highly conservative ones. This supports a practical strategy of preserving useful evidence unless there is a strong reason to discard it.

The third implication is that practical tuning should combine caller-side and phasing-side decisions. The best-performing general configuration was not achieved by optimizing the phaser in isolation, nor by optimizing the caller in isolation, but by combining both types of tuning. The frontier and transferability experiments further showed that this balanced tuned configuration generalizes across baseline, dropout, interaction, and hard conditions. Therefore, the study supports recommending a small tuned configuration set rather than a single-parameter adjustment.

A practical configuration supported by the experiments is:
- `call_min_baseq = 10`
- `call_min_mapq = 20` (or default, since the effect is negligible)
- `max_coverage = 8`
- `min_mapq = 20`
- `min_baseq = 10`

A slightly more runtime-biased alternative with `max_coverage = 6` performed almost identically in the tested scenarios and may be acceptable when efficiency is prioritized.

The fourth implication is that not every intuitive tuning direction is worth pursuing. Increasing `max_coverage` above the plateau did not help. Caller-side and phasing-side mapQ thresholds were also weak knobs. This is important because negative findings are still useful for practical tuning: they narrow the search space and show where optimization effort is unlikely to be rewarded.

Finally, the results suggest that some of the most important remaining opportunities lie upstream of the phasing solver itself. Whenever performance loss was primarily driven by missing shared heterozygous sites, improvements in variant recovery are likely to matter more than further refinement of phasing thresholds. This does not diminish the value of WhatsHap tuning, but it clarifies where future practical gains are most likely to come from.

### 11.3 Limitations of the present study

Several limitations should be acknowledged when interpreting the results.

First, the study uses a synthetic benchmarking pipeline rather than a real biological truth set. This was a deliberate design choice because it enabled precise control over stressors such as duplication, dropout, and correlated error bursts, and made oracle-vs-called attribution possible. However, synthetic truth still simplifies some aspects of real sequencing data and genome structure. As a result, the absolute metric values reported here should be interpreted as controlled experimental measurements rather than direct predictions of production performance on all real datasets.

Second, some realism knobs produced weaker or noisier signals than expected. In particular, correlated error bursts did not show a strong monotonic degradation trend under the tested parameterization, and duplication alone was milder than expected except at the strongest setting. This does not invalidate the experiments, but it means that some conclusions are better interpreted as ruling out large effects under the current settings rather than as definitive evidence that the stressor is unimportant in general.

Third, the phasing pipeline used a custom adaptor layer to construct the `ReadSet` consumed by the vendored WhatsHap core. This adaptor was necessary for the project’s experimental interface, including support for oracle/called benchmarking, unified JSON provenance, and controlled SNP-only evaluation. Runtime profiling showed that the dominant runtime cost in the research pipeline lies in this adaptor rather than in the solver itself. For this reason, the runtime results should be interpreted carefully: they are informative about the project’s benchmarking infrastructure, but they do not necessarily represent the standard production-facing WhatsHap front-end directly.

Fourth, although the optimization sweeps were broad enough to identify practical recommendations, they do not exhaust the full parameter space. The recommended settings should therefore be interpreted as strong experimentally supported candidates rather than mathematically proven optima. The local search reduced this uncertainty substantially by showing that the selected threshold pair is locally robust, but broader or more adaptive searches could still be explored in future.

Fifth, the transferability and scaling studies were intentionally kept compact to remain tractable within the project scope. They are sufficient to support the main conclusions about robustness and solver-facing scaling, but they do not replace a full industrial-scale benchmarking campaign.

### 11.4 Future work

Several natural extensions follow from this study.

A first direction is to test the balanced tuned configuration on broader and more realistic datasets, including real long-read data or more complex synthetic references. This would help determine how far the present recommendations transfer beyond the controlled pipeline used here.

A second direction is to study upstream variant recovery more directly. Since many end-to-end losses were shown to be caller-limited, a natural extension would be to compare different variant callers, different calling models, or more adaptive caller-side filtering strategies. The present work suggests that this may be a more promising route to further gains than continued refinement of phasing thresholds alone.

A third direction is to deepen the algorithm-facing analysis of the residual phasing problem. The SNP-density scaling study already suggests that tuning can change solver difficulty, but a richer study could investigate how variant density, read connectivity, and overlap structure jointly affect the residual wMEC problem instance presented to the DP solver.

A fourth direction is to compare the project-specific adaptor against the standard WhatsHap front-end more directly. This would help separate properties of the benchmarking interface from properties of the production-facing workflow. Such a comparison is especially relevant for runtime interpretation.

A fifth direction is to investigate region-aware optimization strategies under compounded stress. The duplication × dropout interaction suggests that ambiguous regions and low-coverage gaps may require different handling from simpler regions. Future work could therefore explore more adaptive filtering, region-aware read selection, or other strategies that respond to local evidence quality rather than using globally fixed thresholds.

### 11.5 Conclusion

This study showed that WhatsHap is highly effective when correct variant evidence is available, but that end-to-end phasing performance in realistic long-read pipelines depends strongly on variant recovery, connectivity, and the interaction of stress conditions. The most important isolated stressor in this study was coverage dropout, while the most informative compounded weakness arose from duplication combined with dropout. Optimization experiments showed that practical gains are possible, especially through moderate caller-side and phasing-side base-quality thresholds and a lower retained-coverage cap. These gains generalize across multiple representative scenarios and also improve the scaling behaviour of the downstream solve stage. Overall, the work supports a pipeline-level view of optimization: the strongest practical improvements arise not from tuning the phasing algorithm in isolation, but from preserving useful evidence and improving the quality of the phasing problem that is ultimately passed to WhatsHap.