# blocks.py (new, under whatshap)
import logging
from typing import Sequence, Mapping, Optional, Set

from .core import ReadSet, NumericSampleIds
from .graph import ComponentFinder
from .utils import plural_s

def find_components(
    phased_positions: Sequence[int],
    reads: ReadSet,
    master_block: Optional[Sequence[int]] = None,
    heterozygous_positions: Optional[Mapping[int, Set[int]]] = None,
) -> Mapping[int, int]:
    """
    Return a dict that maps each variant position to the component it is in.
    Variants are considered to be in the same component if a read exists that
    covers both. A component is identified by the position of its leftmost
    variant.
    master_block -- List of positions in a "master block", i.e. all blocks containing
                    any of these positions are merged into one block.
    heterozygous_positions -- A dictionary mapping numeric sample ids to sets of
                              positions. Component building is then restricted to variants
                              at these positions. If None, all variants are used.
    """
    logger = logging.getLogger(__name__)
    logger.debug("Finding connected components ...")
    assert phased_positions == sorted(phased_positions)

    # Find connected components.
    # A component is identified by the position of its leftmost variant.
    component_finder = ComponentFinder(phased_positions)
    phased_positions_set = set(phased_positions)
    for read in reads:
        if heterozygous_positions is None:
            positions = [
                variant.position for variant in read if variant.position in phased_positions_set
            ]
        else:
            positions = [
                variant.position
                for variant in read
                if (variant.position in phased_positions_set)
                and (variant.position in heterozygous_positions[read.sample_id])
            ]
        for position in positions[1:]:
            component_finder.merge(positions[0], position)
    if master_block is not None:
        for position in master_block[1:]:
            component_finder.merge(master_block[0], position)
    components = {position: component_finder.find(position) for position in phased_positions_set}
    return components

def compute_overall_components(
    accessible_positions: Sequence[int],
    all_reads: ReadSet,
    distrust_genotypes: bool,
    family: Sequence[str],
    genetic_haplotyping: bool,
    homozygous_positions: Sequence[int],
    numeric_sample_ids: NumericSampleIds,
    superreads_list: Sequence[ReadSet],
) -> Mapping[int, int]:
    master_block = None
    heterozygous_positions_by_sample: Optional[Dict[int, Set[int]]] = None
    accessible_positions_set = set(accessible_positions)
    # If we distrusted genotypes, we need to re-determine which sites are homo-/heterozygous after phasing
    if distrust_genotypes:
        hom_in_any_sample = set()
        heterozygous_positions_by_sample = {}
        heterozygous_gts = frozenset({(0, 1), (1, 0)})
        homozygous_gts = frozenset({(0, 0), (1, 1)})
        for sample, sample_superreads in zip(family, superreads_list):
            hets = set()
            for v1, v2 in zip(*sample_superreads):
                assert v1.position == v2.position
                if v1.position not in accessible_positions_set:
                    continue
                gt = (v1.allele, v2.allele)
                if gt in heterozygous_gts:
                    hets.add(v1.position)
                elif gt in homozygous_gts:
                    hom_in_any_sample.add(v1.position)
            heterozygous_positions_by_sample[numeric_sample_ids[sample]] = hets
        if len(family) > 1 and genetic_haplotyping:
            master_block = sorted(hom_in_any_sample)
    else:
        if len(family) > 1 and genetic_haplotyping:
            master_block = sorted(set(homozygous_positions).intersection(accessible_positions_set))
    return find_components(
        accessible_positions, all_reads, master_block, heterozygous_positions_by_sample
    )
