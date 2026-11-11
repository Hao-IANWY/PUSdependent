# Downstream Analysis Tools

Standalone scripts and R tools for the downstream analysis of called Ψ sites,
including statistical identification of writer-dependent and
temperature-responsive sites and coordinate annotation.

## Scripts and R Tools

- `call_dependent.py`: Wrapper that drives the differential Ψ-site callers.
- `call_dependent_betabinomial.R`: Beta-binomial GLMM (glmmTMB) test for writer-dependent Ψ sites.
- `call_dependent_methylKit.R`: Differential Ψ-site analysis using the methylKit package.
- `call_temperature_response_betabinomial.R`: Beta-binomial GLMM test for temperature-responsive Ψ sites.
- `convert_to_transcript_coord.py`: Converts genomic coordinates to transcript coordinates.
- `extract_siteinfo.py`: Extracts per-site annotations for called Ψ sites.

### tRNA Annotation
- `tRNA_Sprinzl_annoation/`: Annotation databases and scripts for mapping tRNA coordinates to canonical Sprinzl positions.
