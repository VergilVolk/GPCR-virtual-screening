# PACER-ContextKernel development contract

Frozen before the first model run on the US20260055116 named subset.

## Scientific question

Can an explicitly context-conditioned similarity kernel transfer structure-function information across human/rat and pERK/GTPgammaS without leaking another endpoint from the held-out molecule?

## Data and split

- Data: 25 unique named active moieties with ordinal A/B/C/D labels.
- Contexts: human-pERK, rat-pERK, human-GTPgammaS.
- Outer unit: molecule. Every endpoint row of the held-out molecule is removed from training.
- Inner selection: leave-one-training-molecule-out; no outer labels enter hyperparameter selection.

## Frozen methods

1. EndpointOnly-KRR: uses only the same endpoint.
2. PooledNoContext-KRR: pools all endpoint labels and ignores context.
3. PACER-ContextKRR: Tanimoto chemical kernel multiplied by one of four frozen context kernels: shared, same-endpoint, shared-plus-species/readout, or weak-shared-plus-species/readout.

Ridge alpha is selected from 0.01, 0.1, 1, 10 inside each outer fold. Ordinal bins are mapped to A=3, B=2, C=1, D=0. This mapping affects fitting but not the primary pairwise ordering metric.

## Endpoints and decision

- Primary: macro ordinal concordance across contexts.
- Secondary: per-context ordinal concordance and Spearman.
- Development pass: PACER-ContextKRR must improve macro concordance over both baselines and have no context below 0.50 concordance.

This is post-source-qualification method development, not blind external validation or SOTA evidence. A positive result must be tested on another untouched campaign.
