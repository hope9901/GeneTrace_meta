# How GeneTrace uses Meta

## Submission description

GeneTrace uses **Meta Model API with `muse-spark-1.3` as the investigation planner and evidence interpreter**. It is a custom application agent, rather than a prebuilt Meta agent product. At each step, the model receives the research question, accumulated evidence, currently available checks and the guarded provisional judgment. It returns a JSON tool selection, rationale and evidence-cited hypotheses. GeneTrace validates that response, executes an allowlisted analysis tool with fixed inputs, records reproducibility metadata, and returns the new evidence for the next decision. A final model call synthesizes provisional interpretations and limitations.

## What the model actually receives and returns

Input: question about nimT reliability, saved same-species population context, newly computed annotation/alignment results as they become available, allowed tool names and descriptions, partial-input caveats and the current structured judgment. The app does not send raw FASTQ files, the whole reference genome, arbitrary filesystem contents or the API key in prompt messages.

Output schema (application-validated JSON, not native API function calling):

```json
{
  "next_tool": "check_neighbor_regions",
  "reason": "Candidate read disagreement needs local-control context.",
  "hypotheses": [
    {"text": "Provisional interpretation, not a causal diagnosis.", "evidence_ids": ["E003"]}
  ]
}
```

This snippet illustrates the schema. It is not presented as the verbatim full model response. Real responses are preserved in `demo/verified_execution_report.json`.

## Verified real execution

Report `f64d1b2f2a5c4761b1d023d3bd125796`, model `muse-spark-1.3`, five validated model responses, four fresh tools, approximately 38.3 seconds in this run. Timings are a run observation, not a performance guarantee.

| Call | Model choice / purpose | Application action |
|---|---|---|
| 1 | `inspect_annotation`: inspect the recorded gene structure before interpretation | Parse the distributed original GFF and return overlapping features |
| 2 | `align_candidate_reads`: obtain a new candidate read result | Run minimap2, SAMtools and regional analysis on 756 reference-selected original read pairs |
| 3 | `inspect_read_provenance`: examine source linkage and read structure | Audit the newly aligned subset and disclose unverified assembly-version linkage |
| 4 | `check_neighbor_regions`: test whether disagreement is local | Calculate nimR and yeaO from the new alignment |
| 5 | Final synthesis (`next_tool: null`) | Validate cited evidence IDs and preserve provisional interpretations |

The fourth decision rationale in the recorded response was:

> Candidate nimT read disagreement (E003) needs local-control context; testing neighboring genes distinguishes locus-specific signal from broader regional assembly/read or dataset-linkage effects before any causal claim.

This is an investigation rationale, not proof that the check can distinguish those causal explanations conclusively. Both neighbors also showed disagreement; the final result kept the cause unresolved.

## Division of responsibility

| Meta model | GeneTrace code and established analysis tools |
|---|---|
| Selects among currently exposed checks and explains the choice | Defines allowed checks, required initial steps, conditional availability and fixed tool inputs |
| Revises hypotheses after evidence arrives | Executes minimap2 / SAMtools / GFF parsing / regional and provenance calculations |
| Produces evidence-cited provisional interpretations | Validates JSON, tool names and citation membership; records commands, versions and hashes |
| Explains competing possibilities and limits | Produces the deterministic structured judgment and prevents a confirmed biological-loss claim |

The initial annotation and alignment stages each expose only one required action. Later tool availability depends on observed disagreement; required follow-up cannot be skipped. The application chooses the search space. The model does not generate or execute shell commands, discover arbitrary tools, retrain a genome model, or replace alignment algorithms. Citation validation checks membership, not factual correctness of the model's prose.

## AI versus fixed policy

The fixed-policy arm uses the same data, tools and minimum-check rules without a model. Both produced equal deterministic measurements and can detect the neighboring discrepancy. Meta adds selection rationales and cited synthesis; this experiment does not establish improved accuracy, faster expert judgment or unique AI discovery. See [EXECUTION_COMPARISON.md](EXECUTION_COMPARISON.md).

## Code and artifacts

- `genetrace/agent.py`: official API request, prompt policy, JSON parsing, tool/citation validation.
- `genetrace/execution.py`: allowed tools, conditional follow-up, fixed executors and reassessment loop.
- `genetrace/demo.py`: local investigation endpoint and stored reports.
- `demo/verified_execution_report.json`: real model decisions, evidence and execution metadata.
- `demo/live_execution_verification.json`: compact verification record.

The HTTP interface used is `/v1/chat/completions` on `https://api.meta.ai`. API credentials stay in local `.env` and the authorization header. They are excluded from exported reports and the public repository.
