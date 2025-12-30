import numpy as np

from algorithms.io.reads_data import ReadsData
from whatshap import core

def build_readset_from_readsdata(data: ReadsData) -> core.ReadSet:
    """
    Convert ReadsData (dense matrix of alleles) into a WhatsHap ReadSet.

    - Treat each row as a 'read' with name = f"r{i}"
    - Positions are column indices (0..N-1)
    - Alleles: 0/1; ignore -1 (missing)
    - Use uniform mapq/sample_id/source_id for now (you can refine later).
    """
    A = data.alleles      # shape (R, N), values in {0,1,-1}
    R, N = A.shape

    readset = core.ReadSet()
    for i in range(R):
        row = A[i]
        if np.all(row < 0):
            continue  # skip entirely missing reads

        # Build a Read
        read = core.Read(
            name=f"r{i}",
            mapq=60,          # arbitrary 'good' MAPQ
            source_id=0,      # all from same source; extend later if needed
            sample_id=0,      # single diploid sample
        )

        # Add variants (positions 0..N-1)
        for pos in range(N):
            allele = int(row[pos])
            if allele < 0:
                continue
            read.add_variant(pos, allele, 40)  # default: Q=40

        # make sure variants in this read are sorted by position (belt-and-suspenders)
        read.sort()

        readset.add(read)

    # VERY IMPORTANT: sort reads in ReadSet by first variant position
    readset.sort()

    return readset


# import numpy as np

# from algorithms.io.reads_data import ReadsData
# from whatshap import core

# def build_readset_from_readsdata(data: ReadsData) -> core.ReadSet:
#     """
#     Convert ReadsData (dense matrix of alleles) into a WhatsHap ReadSet.

#     - Treat each row as a 'read' with name = f"r{i}"
#     - Positions are column indices (0..N-1)
#     - Alleles: 0/1; ignore -1 (missing)
#     - Use uniform mapq/sample_id/source_id for now (you can refine later).
#     """
#     A = data.alleles      # shape (R, N), values in {0,1,-1}
#     R, N = A.shape

#     readset = core.ReadSet()
#     for i in range(R):
#         row = A[i]
#         if np.all(row < 0):
#             continue  # skip entirely missing reads

#         # Build a Read
#         read = core.Read(
#             name=f"r{i}",
#             mapq=60,          # arbitrary 'good' MAPQ
#             source_id=0,      # all from same source; extend later if needed
#             sample_id=0,      # single diploid sample
#         )

#         # Add variants
#         for pos in range(N):
#             allele = int(row[pos])
#             if allele < 0:
#                 continue
#             # In WhatsHap, a read stores variant positions + alleles; in Cython it's usually:
#             # read.add_variant(position, allele, quality)
#             read.add_variant(pos, allele, 40)  # default: Q=40

#         readset.add(read)

#     return readset