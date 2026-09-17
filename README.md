# BetaGeometry

An olympiad-geometry prover: a symbolic engine that proves what follows mechanically,
paired with a language model that proposes the auxiliary points it cannot invent on its
own.

The split is the point. Every step of every proof comes from the symbolic engine, so a
proof is correct whether or not the model's suggestion was any good; the model only ever
answers "which extra point should I draw?", and a bad answer costs time, not soundness.

## Layout

| | |
|---|---|
| `relations.py` | Predicates (`cong`, `eqangle`, `para`, ...) and the numeric primitives |
| `problem.py` | A problem: points, assumptions, goals, derived relations |
| `constructions.py` | `Canva` — the constructions that can draw a new point |
| `ar.py` | Algebraic reasoning: linear tables over segments and angles |
| `dd_ar.py` | The deductive engine — rules plus AR |
| `rabbits.py` | The construction vocabulary, the parser, and the proposer interface |
| `gemma.py` | Gemma-backed proposer, and the shared prompt format |
| `search.py` | The solve loop: deduce, propose a point, deduce again |
| `data_generation.py` | Mining `(premises, goal) -> construction` training examples |
| `read_in_geogebra_file.py` | `.ggb` diagram parsing |
| `read_in_relations.py` | `.txt` statement parsing |

Entry points: `solve.py`, `generate_batch.py`, `consolidate_data.py`,
`split_seeds.py`, `finetune.py`. Cluster scripts live in `oscar/` — see
[`oscar/README.md`](oscar/README.md).

## Install

Work inside a virtual environment so the project's dependencies stay separate from your
system Python — `torch` in particular is large and version-sensitive.

```bash
python3 -m venv .venv
source .venv/bin/activate          # Windows: .venv\Scripts\activate

python -m pip install -e .         # solving and data generation
python -m pip install -e '.[llm]'  # adds torch/transformers/peft, for Gemma
python -m pip install -e '.[dev]'  # pytest
```

Everything below assumes the environment is active — `source .venv/bin/activate` in each
new shell. `deactivate` leaves it. To start over, delete `.venv/` and recreate it.

If you would rather not activate it, call the interpreter directly:

```bash
.venv/bin/python solve.py problem1
.venv/bin/python -m pytest
```

## Solve a problem

A problem named `foo` needs `geogebra_files/foo.ggb` (the diagram) and
`text_files/foo.txt` (the statement).

```bash
python solve.py problem1                 # symbolic engine alone
python solve.py problem1 --plot
python solve.py usamo_2023_p1 --random-rabbits          # random auxiliary points
python solve.py usamo_2023_p1 --gemma gemma-finetuned-geometry
```

Statements use the predicate DSL, assumptions then goal:

```text
circle O B R D; cong O B O Y; col B O D; eqangle B O R Y O D; ? contri2 R B D Y D B
```

## Generate training data

Each example is a case where an auxiliary point was the key step: a fact that does not
mention the new point, but needed it in the derivation.

```bash
python data_generation.py --seed 0 --rounds 4     # one configuration
python generate_batch.py --seeds 0-999            # many, one file per seed
python consolidate_data.py --data-dir data        # merge, dropping duplicates
```

```json
{
  "input_text": "col X7 X3 X5; ? eqangle X1 X3 X5 X1 X3 X7",
  "output_text": "incenter(X5, X1, X7)"
}
```

At scale this belongs on a cluster — see [`oscar/README.md`](oscar/README.md).

## Fine-tune Gemma

```bash
python finetune.py --data data/training_data.json --output-dir gemma-finetuned-geometry
```

LoRA on `google/gemma-3-1b-pt` (gated: accept the licence and `huggingface-cli login`
first). Training and inference build their prompts through the same
`gemma.build_prompt`, so the two cannot drift apart.

## Tests

```bash
source .venv/bin/activate
pytest
```

The suite checks two independent things. **Completeness**: the problems that were
solvable stay solvable. **Soundness**: every relation the engine derives is true of the
diagram it came from, and every relation a construction emits is true of the figure it
just drew. The second is what catches a transposed index or a miscomputed construction —
the symbolic layer cannot detect those on its own, because it will happily reason from a
false premise.

Angles are compared as directed angles modulo pi throughout, matching the engine. Using
undirected angles instead reports every inscribed-angle relation as false.
