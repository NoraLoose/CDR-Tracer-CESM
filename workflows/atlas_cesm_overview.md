# atlas.py and cesm.py — Overview

## cesm.py

Low-level CESM case factory. Its job is to create, configure, and build individual CESM/SMYLE model cases — it handles cloning a reference case, writing namelists, setting CDR forcing (OAE, DOR, ERW, antitracers, beta/eta sensitivity fields), and submitting the build to the queue. It is essentially the interface between Python and the CESM XML/namelist system.

## atlas.py

High-level orchestration layer that sits on top of `cesm.py`. It manages entire ensembles of cases — batching and bundling SLURM job submissions, tracking case status, generating build scripts via Jinja templates, and running analysis/validation notebooks via Papermill. The central `global_irf_map` class organizes the full set of IRF (Impulse Response Function) simulation metadata for a global CDR atlas experiment (OAE, DOR, ERW, ANTITRACER vintages).

---

In short: `cesm.py` knows how to set up one CESM case; `atlas.py` knows how to run thousands of them and analyze the results.
