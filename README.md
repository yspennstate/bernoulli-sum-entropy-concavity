# Entropy concavity for at most six positively weighted Bernoulli variables

Jörg Neunhäuserer conjectures that a family of recursively overlapped Bernoulli distributions has concave entropy. In Amir's formulation the question is whether the Shannon entropy of a sum of independent Bernoulli(p) variables with positive real weights is concave in p; Amir proves this for at most three summands. This paper proves strict concavity for four, five and six summands, which settles the conjecture through six summands.

The entropy depends on the weights only through which subset sums coincide. An integer representative in the style of Muroga, Toda and Takasu keeps every coincidence; there are 23, 216 and 6,767 entropy profiles for four, five and six summands, and for six summands they are enumerated through the spaces of subset-sum relations (23,273 orbits). Exact rational Bernstein enclosures prove negative curvature for each pattern away from the endpoints, and an analytic estimate covers the endpoints. For every number of summands the paper also bounds the curvature near the endpoints and finds its leading term. Arbitrary positive weights from seven summands onward remain open.

## Contents

- `main.tex`, `main.pdf`: the article.
- `code/`: the programs and the complete proof data; `code/README.md` gives the commands and the expected outputs. Python 3.9 or later (standard library only) and a C++17 compiler suffice. The six-summand certificate is about 11 MB.
- `CREDITS.md`: what is known, what is new, and where each piece comes from.
- `SHA256SUMS`: checksums of every file in the repository.

## License

The text is released under CC BY 4.0 and the code under the MIT license; see `LICENSE`.
