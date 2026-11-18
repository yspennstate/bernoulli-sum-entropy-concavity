# Reproduction

Use Python 3.9 or later, with assertions enabled (no `-O` or `PYTHONOPTIMIZE`), and a C++17 compiler. All Python programs use only the standard library. Run sequentially on one core from this directory. Commands use `nice -n 10`; compilers and interpreters may be replaced by compatible installations.

## Verify the supplied proof data

```sh
nice -n 10 c++ -O2 -std=c++17 enumerate_flats.cpp -o enumerate_flats
nice -n 10 ./enumerate_flats 4 flats4
nice -n 10 ./enumerate_flats 5 flats5
nice -n 10 ./enumerate_flats 6 flats6
nice -n 10 c++ -O2 -std=c++17 verify_geometry.cpp -o verify_geometry
nice -n 10 ./verify_geometry
nice -n 10 python3 -B validate_geometry.py
nice -n 10 python3 -B replay.py
nice -n 10 python3 -B replay_six.py
nice -n 10 python3 -B verify_six.py
nice -n 10 python3 -B summarize.py
nice -n 10 python3 -B check_latex.py
```

The orbit program writes the corresponding `flats*_flats.txt` and `flats*_profiles.json` files. Output includes elapsed seconds, which vary. Its final profile counts are 23, 216 and 6,767. Six-summand rank counts, with cumulative profile counts, are:

```text
rank=0 orbits=1 profiles=1
rank=1 orbits=9 profiles=10
rank=2 orbits=117 profiles=76
rank=3 orbits=1593 profiles=556
rank=4 orbits=10334 profiles=2948
rank=5 orbits=11219 profiles=6767
```

There are 301 mixed-sign unoriented ternary normals and 23,273 feasible orbits. The generator checks primitive coordinates before extension, bounding every intermediate extension coordinate by 180. Normalized exterior coordinates determine the spaces exactly; coordinate permutations define the orbits. Profiles pack cardinality counts into base-64 digits with no carries.

`verify_geometry.cpp` uses contained-normal sets and nullspace projections to check orbit distinctness and closure under every feasible extension. It prints `PASS: no duplicate orbits; every feasible successor has a listed orbit` and writes `geometry_results.json`: 3,628,254 total candidate extensions and 717,239 feasible ones, counting repetitions. `validate_geometry.py` separately checks every maximal minor, ternary generation, positive kernel point and profile, and writes `geometry_record_results.json`. Expected maxima of absolute primitive exterior coordinates by rank are `[1,1,2,4,12,30]`; all 23,273 records and 6,767 profiles pass.

Expected certificate replay results:

| Summands | Profiles | Cells | Largest stored upper bound | Endpoint cutoff |
| --- | --- | --- | --- | --- |
| 4 | 23 | 329 | -27813/1000000 | 1/4096 |
| 5 | 216 | 3757 | -14009/1000000 | 1/16384 |
| 6 | 6767 | 130172 | -1/20000 | 1/32768 |

## Regenerate the four- and five-summand profiles and certificates

```sh
nice -n 10 c++ -O2 -std=c++17 enumerate.cpp -o enumerate
nice -n 10 ./enumerate 4 16 enumerated_4.txt
nice -n 10 ./enumerate 5 55 enumerated_5.txt
nice -n 10 python3 -B enumerate_python.py
nice -n 10 python3 -B certify.py
```

Expected enumeration output is `m=4 bound=16 vectors=3876 profiles=23` and `m=5 bound=55 vectors=5006386 profiles=216`. The independent polynomial-product program prints `POLYNOMIAL ENUMERATION PASS m=4: 3876 vectors, 23 profiles` and `POLYNOMIAL ENUMERATION PASS m=5: 5006386 vectors, 216 profiles`, checking equality of the full profile sets. The certificate generator prints `CERTIFICATE PASS` summaries with 329 and 3,757 cells and the corresponding exact table bounds.

## Certificate and summary details

`replay.py` uses direct polar-form restriction and 18 logarithm terms. It prints `REPLAY PASS` and `ENDPOINT PASS` for four and five summands. The exact endpoint positive/negative constants are respectively `9/16`, `4285536255/4294967296`, and `175/512`, `72004829826842649/72057594037927936`.

`replay_six.py` uses monomial expansion, affine substitution and Bernstein conversion, with 18 logarithm terms. It checks every partition and stored bound, reconstructs every certificate witness, and prints a passing JSON summary plus `ENDPOINT PASS m=6`: positive upper `675/2048`, negative lower `9440410169765069258743/9444732965739290427392`.

`verify_six.py` independently uses normalized Bernstein differences, de Casteljau restriction, 22 logarithm terms and outward fibre rounding at 2^-42. It verifies the same 130,172 cells and every compact witness against the full profile set. Expected final fields include `pass: true`, `profiles: 6767`, `cells: 130172`, `largest_stored_upper: -1/20000`, `largest_recomputed_upper: -110169287/2199023255552`, and `compact_max_coordinate: 36`. Each singleton leading-coefficient sum is six. The certificate witnesses need not be the compact witnesses; equality of their full profiles is checked.

## Independent re-checks

These programs repeat computations that the paper's proofs take from the programs above; the paper does not rely on them.

- `enumerate_python.py` recomputes the four- and five-summand profiles from the polynomial prod_i (1 + X^((m+1) v_i + 1)): the coefficient of X^((m+1) s + k) is N_{s,k}. It prints the same vector and profile counts as `enumerate.cpp` and checks that the full profile sets, not only their sizes, are equal.
- `verify_geometry.cpp` represents each relation space by its set of contained ternary normals, computed from a nullspace basis, checks that the listed permutation orbits are distinct and closed under every feasible extension (3,628,254 candidate extensions from ranks zero through four, 717,239 of them feasible, counting repetitions), and checks that for four and five summands the orbit enumeration and the bounded integer enumeration give the same profile sets.
- The certificate generators use de Casteljau subdivision and 14 logarithm terms; the four- and five-summand replayer evaluates the polar form directly with 18 terms, and `verify_six.py` uses 22 terms (details above).

`summarize.py` prints the four-summand representative table, all certificate and orbit totals, agreement of the complete four/five profile sets, the compact coordinate bound, closure totals, endpoint-box costs and the three counterexample identities. In dimension seven its comparison-preserving bound is 896, its box count is 94162890744574080, its normal count is 966 and its permutation count is 5040. The stronger classical order-six determinant estimate is 160; the convenient divisibility estimate is 192. Rational logarithm lower bounds prove positive signs for all three counterexamples.

`check_latex.py` checks braces, nested environments, references, citations, macro names, math delimiters and duplicate labels/keys, and runs five deliberate rejection probes in memory. It defaults to the article source. A passing static check does not establish a PDF build or visual layout.

## Regenerate the bounded enumerations and certificates

```sh
nice -n 10 c++ -O2 -std=c++17 enumerate.cpp -o enumerate
nice -n 10 ./enumerate 4 16 enumerated_4.txt
nice -n 10 ./enumerate 5 55 enumerated_5.txt
nice -n 10 python3 -B enumerate_python.py
nice -n 10 python3 -B certify.py
nice -n 10 python3 -B certify_six.py --output-dir generated_six
```

The bounded enumerator prints `m=4 bound=16 vectors=3876 profiles=23` and `m=5 bound=55 vectors=5006386 profiles=216`. `enumerate_python.py` uses integer polynomial multiplication and prints `POLYNOMIAL ENUMERATION PASS` with those same counts, checking full set equality. Both complete bounded enumeration files are supplied. Their first line is `m bound vectors profiles`; each remaining line contains weights, support size and sorted base-64 count vectors.

`certify.py` replaces the four- and five-summand certificate JSON files and prints `CERTIFICATE PASS` with the table above. The six-summand generator reads the local orbit profiles and writes a new directory, refusing to append to an existing certificate. It uses local `bernstein_bounds.py`, 14 logarithm terms, outward fibre rounding at 2^-40, bisection and final upward rounding to millionths. Its final JSON summary has `complete: true`, `profiles_total: 6767`, `cells: 130172`, `largest_upper_bound: -1/20000`. Regenerated certificate records match the supplied records exactly; only timing metadata varies. The provided `certificate_six.jsonl` is approximately 10 MB, so no xz unpacking is needed.

Saved `.out`, `.resources` and result JSON files record the executed checks. Fractions and integer comparisons establish every proof inequality; floating-point values are used only for timings. The programs neither use multiple workers nor start background services.

## Encodings and subdivision identities

For four and five summands a count vector is encoded as sum_k N_{s,k} 64^k. No carry occurs because a fibre contains at most 32 subsets. Sorting these encodings preserves the multiset of count vectors. For six summands each count is at most binomial(6,k), hence at most 20; its base-64 profile encoding fits in an unsigned 64-bit integer. Before extension the largest absolute primitive exterior coordinate is 30, and an unnormalized extension coordinate is a sum of at most six such values, so its absolute value is at most 180.

The polar-form restriction used by `replay.py` averages the Bernstein coefficient c_k over a total of d independent Bernoulli variables, with d-j parameters equal to alpha and j equal to beta. Multiaffinity gives the coefficient at index j on [alpha,beta]. The monomial-to-Bernstein conversion used by `replay_six.py` is

```text
t^j = sum_{k=j}^d [binomial(k,j)/binomial(d,j)] binomial(d,k) t^k (1-t)^(d-k).
```

It follows from binomial(d,k) binomial(k,j) = binomial(d,j) binomial(d-j,k-j) and the binomial theorem.
