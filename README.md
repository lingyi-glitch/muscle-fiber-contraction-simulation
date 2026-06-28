# Muscle Fiber Contraction Simulation

MATLAB implementation of a biophysical skeletal muscle fiber model coupling calcium-mediated excitation, cross-bridge cycling, the First Reaction Method, and a one-dimensional finite element formulation of thin and thick filaments.

This repository accompanies the article:

> Jiaxiang Xu and Bin Chen. Effects of stimulation frequencies on energy efficiency of a muscle fiber during contraction. Biophysical Journal, 2026. https://doi.org/10.1016/j.bpj.2026.06.001

## Repository Contents

- `TaskList.m` - batch simulation entry point.
- `main.m` - runs one loading condition and writes the raw trajectory.
- `MonteCarlo.m` - stochastic cross-bridge state transitions using the First Reaction Method.
- `InitSKDD.m`, `MotorSKDD.m` - finite element stiffness assembly.
- `Ca2CaT.m`, `openratio.m`, `CalciumDynamics.m` - calcium activation dynamics and predicted binding-site activation.
- `datacreat.m`, `Dataread.m` - geometry input generation and reading.
- `examples/run_example.m` - short smoke-test simulation.
- `docs/` - parameter notes and output-column definitions.

## Requirements

- MATLAB R2023b or later is recommended.
- No MATLAB toolbox is required for the current code path.

The code was prepared and checked on MATLAB R2023b.

## Quick Start

From the repository root, run:

```matlab
addpath('examples');
run_example
```

The example creates a geometry file, runs a short single-case simulation, and writes output under:

```text
examples/output/example_short_run/
```

This example is only a smoke test for installation and code execution. Use a longer, carefully configured batch run for paper-level analysis.

For the batch script used as a template, run:

```matlab
TaskList
```

`TaskList.m` is intentionally written as a configurable driver. Adjust frequency, load, repeat number, and stop time there before launching large production simulations.

## Output

Each `Force=...` output file contains one row per stochastic event. See `docs/OUTPUTS.md` for the column definitions.

Generated simulation outputs, logs, and MATLAB temporary files are ignored by Git by default.

## Citation

If you use this code, please cite the paper above and this repository. GitHub can read `CITATION.cff` and show a "Cite this repository" button.

## License

This project is released under the MIT License. See `LICENSE`.
