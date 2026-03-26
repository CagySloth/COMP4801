import numpy as np

from algorithms.io.reads_data import ReadsData
from algorithms.vendor.whatshap_vendor import import_whatshap_vendor

_, core, _ = import_whatshap_vendor()


def build_readset_from_readsdata(
    data: ReadsData,
    sort_readset: bool = False,
    *,
    use_positions: bool = False,
    mapq: int = 60,
    source_id: int = 0,
    sample_id: int = 0,
    default_qual: int = 40,
) -> core.ReadSet:
    """
    Convert ReadsData (dense matrix of alleles) into a WhatsHap ReadSet.

    Default behavior (backwards-compatible with your current tests/pipeline):
      - Variant positions are column indices 0..N-1.

    If use_positions=True:
      - Variant positions are taken from data.positions[row, col].
        This is what you want in VCF-mode after subsetting to heterozygous sites,
        because you can keep the *original* VCF positions instead of renumbering.

    Alleles: 0/1; ignore -1 (missing).
    """
    A = data.alleles  # shape (R, N), values in {0,1,-1}
    R, N = A.shape
    P = getattr(data, "positions", None)

    readset = core.ReadSet()
    for i in range(R):
        row = A[i]
        if np.all(row < 0):
            continue  # skip entirely missing reads

        read = core.Read(
            name=f"r{i}",
            mapq=int(mapq),
            source_id=int(source_id),
            sample_id=int(sample_id),
        )

        for col in range(N):
            allele = int(row[col])
            if allele < 0:
                continue

            if use_positions and P is not None:
                pos = int(P[i, col])
                if pos < 0:
                    continue
            else:
                pos = col

            read.add_variant(pos, allele, int(default_qual))

        # WhatsHap expects variants within a read to be sorted
        read.sort()
        readset.add(read)

    if sort_readset:
        readset.sort()

    return readset
