# Two finite-order pairs with equal invariants and different embeddings

I have two examples where a fixed list of measurements gives the same answers for two different finite orders and only one member embeds in the specified Minkowski space. This is a mathematical question about those measurements and those orders. Novelty not established; this repository exists to ask.

Two fixed lists of invariants, the link-and-interval checklist and the published spectral toolkit with directed link degrees, do not decide whether a finite order has a strict order embedding in 2+1 Minkowski space; an explicit connected pair has 30 elements.

A second connected pair has 70 elements and the same agreement on the two lists and exactly one member has a strict order embedding in 3+1 Minkowski space. The files contain the cover relations and exact coordinate certificates and the written reasons the other members cannot embed.

Both pairs here are two-layer orders (related constructions may have three layers). Nothing is claimed about faithful embeddings, sprinkled or manifoldlike sets. The perturbation method of arXiv:2008.02291 was not tested. These lists are not all possible invariants and no smallest-size claim is made.

## Run the check

Python 3 with its standard library is sufficient. From this folder run:

```sh
python3 check/check_pairs.py
```

The check rebuilds both lists for each member and verifies every coordinate relation and the supplied obstruction data. It prints PASS lines only after its known answers and deliberately changed inputs are checked. It uses exact integers and fractions and polynomial block identities and no floating eigenvalue comparison. The written non-embedding proofs also use mathematical theorems; the program is not a proof assistant.

## Read the examples

- [Definitions and the exact polynomial formulas](check/METHOD.md).
- [30-element pair and the non-embedding reason](pairs/pair30.txt). Its JSON contains both cover lists and coordinates for P.
- [70-element pair and the non-embedding reason](pairs/pair70.txt). Its JSON contains both cover lists and coordinates for B.
- [Two-page evidence sheet](EVIDENCE_SHEET.txt).

The strict convention is i<j exactly when t_j>t_i and (t_j-t_i)^2 exceeds the squared spatial distance. The certificates have no null-separated pairs.

## Credit

The question builds on the causal-set census and spectral tests of Yazdi and Kempf (arXiv:1611.09947) and Yazdi, Letizia and Kempf (arXiv:2008.02291). Their fixed-toolkit census is the reference for the small-order reproductions. The 30-element example uses the cospectral planar and nonplanar graphs in van Dam and Haemers (2003), Figure 2. Its non-embedding reason uses Scheinerman's circle-order theorem and the circle-order literature including Sidney, Sidney and Urrutia (1988). These prior results are not claimed as new. Full citations and access limits are in [SOURCES.md](SOURCES.md).

Code uses the MIT license. Original text and data use CC BY 4.0. Existing published material retains its own terms; see LICENSE.
