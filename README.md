# grpnat
Exploration of abelian groups and naturality

This is a code to try to find examples of permutations that violate the naturality of automorphisms of abelian groups of the type (Z/2Z)^n.

The original purpose was to find counterexamples showing that some abelian groups (Z/2Z)^n were not naturally isomorphic to their character groups Hom((Z/2Z)^n, S1).

However, a mistake led to only test permutations as *h* functions, instead of all of the linear combinations.

For the original problem with all linear combinations as *h*, a counterexample was manually found later, for n = 3.

The mistake led to an interesting conjecture, that I intend to prove:
that if we restrict the homomorphisms, in Hom((Z/2Z)^n, (Z/2Z)^n), to only the permutations,
then there exists a natural isomorphism between (Z/2Z)^n and its character group.
