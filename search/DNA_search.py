from __future__ import annotations
from enum import IntEnum
from typing import Tuple, List, TypeVar, Iterable, Sequence, Generic, List, Callable, Set, Deque, Dict, Any, Optional
from typing_extensions import Protocol
from heapq import heappush, heappop

Nucleotide: IntEnum = IntEnum('Nucleotide', ('A', 'C', 'G', 'T'))
Codon = Tuple[Nucleotide, Nucleotide, Nucleotide]
Gene = List[Codon]

def string_to_gene(input: str) -> Gene:
    gene: Gene = []

    for i in range(0, len(input), 3):
        if (i + 2) >= len(input):
            return gene

        codon: Codon = (Nucleotide[input[i]], Nucleotide[input[i+1]], Nucleotide[input[i+2]])
        gene.append(codon)

    return gene

class SpecificContains():
    def linear_contains(self, gene: Gene, key_codon: Codon) -> bool:
        for codon in gene:
            if codon == key_codon:
                return True

                return False

    def binary_contains(self, gene: Gene, key_codon: Codon) -> bool:
        low: int = 0
        high: int = len(gene) - 1

        while low <= high:
            mid: int = (low + high)//2

            if gene[mid] < key_codon:
                low = mid + 1
            elif gene[mid] > key_codon:
                high = mid - 1
            else:
                return True

        return False

C = TypeVar("C", bound="Comparable")

class Comparable(Protocol):
    def __eq__(self, other: Any) -> bool:
        pass

    def __lt__(self: C, other: C) -> bool:
        pass

    def __gt__(self: C, other:C) -> bool:
        pass

    def __le__(self: C, other: C) -> bool:
        return self < other or self == other

    def __ge__(self: C, other: C) -> bool:
        return not self < other

T = TypeVar('T')

class GenericContains():
    def linear_contains(self, iterable: Iterable[T], key: T) -> bool:
        for item in iterable:
            if item == key:
                return True

        return False

    def binary_contains(self, sequence: Sequence[C], key: C) -> bool:
        low: int = 0
        high: int = len(sequence) - 1

        while low<= high:
            mid: int = (low + high)//2

            if sequence[mid] < key:
                low = mid + 1
            elif sequence[mid] > key:
                high = mid - 1
            else:
                return True

        return False


if __name__ == "__main__":
    gene_str: str = "ACGTGGCTCTCTAACGTACGTACGTACGGGGTTTATATATACCCTAGGACTCCCTTT"
    my_gene: Gene = string_to_gene(gene_str)

    acg: Codon = (Nucleotide.A, Nucleotide.C, Nucleotide.G)
    gat: Codon = (Nucleotide.G, Nucleotide.A, Nucleotide.T)
    print(f"specific linear_contains: {SpecificContains().linear_contains(my_gene, acg)}")
    print(f"specific linear_contains: {SpecificContains().linear_contains(my_gene, gat)}")

    sorted_gene: Gene = sorted(my_gene)
    print(f"specific binary_contains: {SpecificContains().binary_contains(sorted_gene, acg)}")
    print(f"specific binary_contains: {SpecificContains().binary_contains(sorted_gene, gat)}")

    print(f"Generic linear_contains: {GenericContains().linear_contains([1,2,3,4,5], 2)}")
    print(f"Generic linear_contains: {GenericContains().linear_contains(['a','b','c'], 'd')}")
    print(f"Generic binary_contains: {GenericContains().binary_contains([1,2,3,4,5], 2)}")
    print(f"Generic binary_contains: {GenericContains().binary_contains(['a','b','c'], 'd')}")
