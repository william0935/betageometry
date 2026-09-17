"""Gemma-backed rabbit proposer.

Ported from the fine-tuning work in the sibling BetaGeometry repo and reworked so that
training and inference share one definition of the prompt format -- previously the
tokenisation in the training script and the prompt built at inference time were written
out separately, and any drift between them silently degrades the model.

`torch`/`transformers`/`peft` are imported lazily inside `GemmaProposer.load` so the rest
of the package (and the whole test suite) runs without them installed.
"""

from typing import List, Optional, Sequence

from rabbits import CONSTRUCTION_ARITY, RabbitProposer
from relations import Point


# Base checkpoint. `-pt` is the pretrained (non instruction-tuned) variant, which is what
# you want when fine-tuning on a rigid DSL like this one.
DEFAULT_MODEL_ID = "google/gemma-3-1b-pt"

# The prompt contract. Training and inference MUST build the prompt identically, so both
# go through `build_prompt` below rather than formatting their own strings.
SEPARATOR = "\n### Construction:\n"
INSTRUCTION = (
    "Propose one auxiliary construction that helps prove the goal. "
    "Answer with a single construction call and nothing else. "
    "Available constructions: " + ", ".join(sorted(CONSTRUCTION_ARITY)) + "."
)


def build_prompt(problem_text: str) -> str:
    """The exact text the model is conditioned on, for both training and inference."""
    return f"{problem_text} {INSTRUCTION}{SEPARATOR}"


def build_example(problem_text: str, construction: str, eos: str) -> str:
    """One full training sequence: prompt followed by the target construction."""
    return f"{build_prompt(problem_text)}{construction}{eos}"


class GemmaProposer:
    """Proposes constructions with a (fine-tuned) Gemma checkpoint.

    Example
    -------
        proposer = GemmaProposer(adapter_dir="gemma-finetuned-geometry")
        proposer.load()
        proposer.propose("cong A B A C; ? perp A D B C", points, k=4)
    """

    def __init__(self,
                 model_id: str = DEFAULT_MODEL_ID,
                 adapter_dir: Optional[str] = None,
                 max_new_tokens: int = 32,
                 temperature: float = 0.8,
                 device: Optional[str] = None):
        self.model_id = model_id
        self.adapter_dir = adapter_dir
        self.max_new_tokens = max_new_tokens
        self.temperature = temperature
        self.device = device
        self.model = None
        self.tokenizer = None

    def load(self):
        """Load tokenizer and weights. Safe to call more than once."""
        if self.model is not None:
            return self

        try:
            import torch
            from transformers import AutoModelForCausalLM, AutoTokenizer
        except ImportError as exc:  # pragma: no cover - depends on optional extras
            raise ImportError(
                "GemmaProposer needs the 'llm' extra: pip install -e '.[llm]'"
            ) from exc

        self.tokenizer = AutoTokenizer.from_pretrained(self.model_id)
        if self.tokenizer.pad_token is None:
            self.tokenizer.pad_token = self.tokenizer.eos_token

        self.model = AutoModelForCausalLM.from_pretrained(
            self.model_id,
            dtype=torch.bfloat16,
            device_map=self.device or "auto",
        )

        if self.adapter_dir:
            from peft import PeftModel
            self.model = PeftModel.from_pretrained(self.model, self.adapter_dir)

        self.model.eval()
        return self

    def propose(self, problem_text: str, points: Sequence[Point], k: int = 1) -> List[str]:
        """Sample up to `k` candidate constructions, most likely first.

        Candidates are returned as raw text; the caller parses and validates them with
        `rabbits.parse_construction`, which rejects anything outside the vocabulary.
        """
        import torch

        if self.model is None:
            self.load()

        prompt = build_prompt(problem_text)
        inputs = self.tokenizer(prompt, return_tensors="pt").to(self.model.device)

        with torch.inference_mode():
            outputs = self.model.generate(
                **inputs,
                max_new_tokens=self.max_new_tokens,
                do_sample=k > 1,
                temperature=self.temperature if k > 1 else None,
                num_return_sequences=k,
                eos_token_id=self.tokenizer.eos_token_id,
                pad_token_id=self.tokenizer.pad_token_id,
            )

        prompt_len = inputs["input_ids"].shape[-1]
        seen = set()
        candidates = []
        for sequence in outputs:
            text = self.tokenizer.decode(sequence[prompt_len:], skip_special_tokens=True)
            # The model is trained to stop after one call, but sampling can run on;
            # keep only the first line.
            text = text.strip().splitlines()[0].strip() if text.strip() else ""
            if text and text not in seen:
                seen.add(text)
                candidates.append(text)
        return candidates
