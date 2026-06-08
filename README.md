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

**Algorithms Used**
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

ALGORITHM SelectBestBalancedSolution(P, k)
    //Selects the single best solution from the approximate Pareto front
    //such that no objective deviates far from its individual optimum,
    //and large gains across objectives are preferred over small uniform losses
    //Input: Set P of (path, C[1..k]) pairs — the approximate Pareto front,
    //         number of objectives k
    //Output: Single best (path, C[1..k]) pair from P
    //--- Step 1: Compute the Ideal Vector ---
    //The ideal vector holds the best achievable value for each objective
    //independently, as if each objective were optimised in isolation
    for i ← 1 to k do
        ideal[i] ← ∞
    for each (path, C) ∈ P do
        for i ← 1 to k do
            if C[i] < ideal[i] then
                ideal[i] ← C[i]
    //--- Step 2: Compute the Nadir Vector ---
    //The nadir vector holds the worst value for each objective across the Pareto front
    //Used to normalise deviations to a common [0, 1] scale
    for i ← 1 to k do
        nadir[i] ← 0
    for each (path, C) ∈ P do
        for i ← 1 to k do
            if C[i] > nadir[i] then
                nadir[i] ← C[i]
    //--- Step 3: Compute Normalised Deviations for Each Solution ---
    //Normalise each objective's cost to [0, 1] relative to the ideal–nadir range
    //A value of 0 means this solution achieves the ideal on objective i
    //A value of 1 means this solution is at the worst point on objective i
    for each (path, C) ∈ P do
        for i ← 1 to k do
            if nadir[i] − ideal[i] = 0 then
                D(path, i) ← 0
            else
                D(path, i) ← (C[i] − ideal[i]) / (nadir[i] − ideal[i])
    //--- Step 4: Compute the Maximum Deviation per Solution ---
    //The maximum deviation is the worst any single objective deviates from ideal
    //This penalises solutions that are catastrophically bad on even one objective
    for each (path, C) ∈ P do
        MaxDev(path) ← 0
        for i ← 1 to k do
            if D(path, i) > MaxDev(path) then
                MaxDev(path) ← D(path, i)
    //--- Step 5: Compute the Average Deviation per Solution ---
    //The average deviation rewards solutions that are collectively close to ideal
    //across all objectives, not just on the worst one
    for each (path, C) ∈ P do
        SumDev ← 0
        for i ← 1 to k do
            SumDev ← SumDev + D(path, i)
        AvgDev(path) ← SumDev / k
    //--- Step 6: Compute the Gain-Weighted Balance Score ---
    //For each solution, compute how much improvement it offers across all
    //objectives relative to solutions that are better on the worst objective
    //A solution scores well if a small worsening on one objective
    //is compensated by large improvements on the remaining objectives
    for each (path_a, C_a) ∈ P do
        GainScore(path_a) ← 0
        for each (path_b, C_b) ∈ P do
            if path_a = path_b then continue
            TotalGain ← 0
            TotalLoss ← 0
            for i ← 1 to k do
                delta ← D(path_b, i) − D(path_a, i)
                if delta > 0 then
                    TotalGain ← TotalGain + delta
                else
                    TotalLoss ← TotalLoss + |delta|
            if TotalLoss = 0 then
                GainScore(path_a) ← GainScore(path_a) + TotalGain
            else
                GainScore(path_a) ← GainScore(path_a) + (TotalGain / TotalLoss)
    //--- Step 7: Compute the Final Combined Score ---
    //Combine MaxDev, AvgDev, and GainScore into one scalar
    //Lower MaxDev and AvgDev are better (penalise imbalance and distance from ideal)
    //Higher GainScore is better (reward large cross-objective gains)
    //α and β control how heavily to penalise worst-case vs average deviation
    α ← 0.5
    β ← 0.3
    γ ← 0.2
    for each (path, C) ∈ P do
        Score(path) ← α × MaxDev(path) + β × AvgDev(path) − γ × NormaliseGain(GainScore(path))
    //--- Normalise GainScore to [0,1] across all solutions ---
    MaxGain ← 0
    for each (path, C) ∈ P do
        if GainScore(path) > MaxGain then
            MaxGain ← GainScore(path)
    for each (path, C) ∈ P do
        if MaxGain > 0 then
            NormGain(path) ← GainScore(path) / MaxGain
        else
            NormGain(path) ← 0
        Score(path) ← α × MaxDev(path) + β × AvgDev(path) − γ × NormGain(path)
    //--- Step 8: Select the Solution with the Lowest Score ---
    best ← null
    bestScore ← ∞
    for each (path, C) ∈ P do
        if Score(path) < bestScore then
            bestScore ← Score(path)
            best ← (path, C)
    return best


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
