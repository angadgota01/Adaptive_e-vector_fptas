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

**Algorithm Used**
ALGORITHM ComputeAdaptiveEpsilonVector(G, k, ε_total)
    //Computes a per-objective approximation tolerance vector for the Adaptive ε-Vector FPTAS
    //Input: Directed graph G = (V, E) with edge cost vectors of dimension k,
    //         number of objectives k,
    //         user-supplied total approximation budget ε_total ∈ (0, 1]
    //Output: Epsilon vector ε[1..k] where each ε[i] is the approximation tolerance for objective i

    for i ← 1 to k do
        min_cost[i] ← ∞
        max_cost[i] ← 0

    for each edge e ∈ E do
        for i ← 1 to k do
            if cost(e, i) > 0 then
                min_cost[i] ← min(min_cost[i], cost(e, i))
            max_cost[i] ← max(max_cost[i], cost(e, i))
            // cost(e,i) is the cost of edge e for objective i

    for i ← 1 to k do
        if min_cost[i] = ∞ then
            R[i] ← 0
        else
            R[i] ← log(max_cost[i] / min_cost[i])
    // R[i] indicates the logarithmic span of costs on objective i.
    // If R[i] is small, it means the costs for objective i are close to each other
    // If R[i] is large, it means the costs for objective i are spread out

    R_total ← 0
    for i ← 1 to k do
        R_total ← R_total + R[i]

    // If R_total is 0, it means all costs are the same
    // so we can set epsilon to epsilon_total for all objectives
    if R_total = 0 then
        for i ← 1 to k do
            ε[i] ← ε_total
        return ε[1..k]

    // Calculate weights inversely proportional to R[i] 
    for i ← 1 to k do
        w[i] ← 1 / (R[i] + 1)
    w_total ← 0
    for i ← 1 to k do
        w_total ← w_total + w[i]

    for i ← 1 to k do
        ε[i] ← ε_total × (w[i] / w_total)
    for i ← 1 to k do
        ε[i] ← max(ε[i], 0.01) // Minimum 1% error margin
        ε[i] ← min(ε[i], 1.0) // Maximum 100% error margin
    return ε[1..k]

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
