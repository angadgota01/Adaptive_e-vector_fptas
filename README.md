# Adaptive ε-Vector FPTAS for the Multiobjective Shortest Path Problem

An interactive web demo for the **Adaptive ε-Vector FPTAS** — a fully polynomial-time approximation scheme for finding near-Pareto-optimal shortest paths across multiple cost objectives.

Open `adaptive-evector-fptas.html` in any modern browser. No server, build step, or dependencies required.

---

## Features

**Algorithm Overview**
A four-step visual breakdown of the approach: problem formulation, ε-approximation grid construction, adaptive per-objective scaling, and label-setting expansion with dominance pruning.

**Interactive Graph Builder**
- Click the canvas to place nodes
- Draw directed edges between nodes with auto-generated dual costs
- Set source and target nodes
- Load an 8-node preset graph instantly

**FPTAS Execution**
- Tune ε (epsilon), number of objectives k, and adaptive scaling mode via sliders
- Run the algorithm and watch approximate Pareto-optimal paths highlight live
- Execution log displays labels explored, paths found, and simulated runtime

**Pareto Front Visualization**
- Scatter plot comparing the exact vs. approximate Pareto front
- (1+ε) envelope shaded in real time as ε is adjusted
- Dominance staircase curve overlaid on the plot

**Complexity Reference**
A comparison table of this algorithm against Exact MOSP, Parametric Dijkstra, and Uniform ε-FPTAS, with the formal time and solution-size bounds displayed.
