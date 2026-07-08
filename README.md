# Adaptive ε-Vector FPTAS for the Multiobjective Shortest Path Problem

An interactive web demo for the **Adaptive ε-Vector FPTAS** — a fully polynomial-time approximation scheme for finding near-Pareto-optimal shortest paths across multiple cost objectives.

Open `adaptive-evector-fptas.html` in any modern browser. No server, build step, or dependencies required.

---

## Features

**Algorithm Overview**
A four-step visual breakdown of the approach: problem formulation, ε-approximation grid construction, adaptive per-objective scaling, and label-setting expansion with dominance pruning.

**Interactive Graph Builder**
- Click the canvas to place nodes
- Draw directed edges with configurable multi-objective costs (2–4 objectives)
- Edit or delete nodes and edges directly from the canvas
- Set source and target nodes interactively
- Automatically prevents cycles and self-loops to ensure the graph remains valid
- Enforces non-negative edge costs required by the Adaptive ε-Vector FPTAS
- Load a built-in preset graph (Small, Medium, Dynamically Generated Large graph)
- Reset the graph at any time
- Import and export graph configurations as JSON files

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


ALGORITHM AdaptiveEVectorFPTAS(G, s, t, k, ε[1..k])
    //Finds all approximate Pareto-optimal s–t paths using adaptive per-objective bucketing
    //Input: Directed graph G = (V, E) with k-dimensional edge cost vectors,
    //         source node s, target node t, number of objectives k,
    //         epsilon vector ε[1..k] from ComputeAdaptiveEpsilonVector
    //Output: Set P of (path, C[1..k]) pairs where every returned path P̂ satisfies
    //         C(P̂, i) ≤ (1 + ε[i]) · C(P*, i) for all i, for some exact Pareto-optimal P*
    function B(v, i)
        if v = 0 then return 0
        return ⌊ log(v) / log(1 + ε[i]) ⌋
    function BucketDominates(C_a[1..k], C_b[1..k])
        for i ← 1 to k do
            if B(C_a[i], i) > B(C_b[i], i) then return false
        return true
    function TrueDominates(C_a[1..k], C_b[1..k])
        strict ← false
        for i ← 1 to k do
            if C_a[i] > C_b[i] then return false
            if C_a[i] < C_b[i] then strict ← true
        return strict
    //--- Initialisation ---
    for each node v ∈ V do
        Settled[v] ← ∅
    for i ← 1 to k do
        C_init[i] ← 0
    L_init ← (s, C_init[1..k], [ ])
    PQ ← empty min-heap ordered by C[1]
    Insert(PQ, L_init)
    P ← ∅
    //--- Main label-setting loop ---
    while PQ ≠ ∅ do
        (v, C[1..k], path) ← ExtractMin(PQ)
        //--- Discard if dominated by anything already settled at v ---
        skip ← false
        for each C_s ∈ Settled[v] do
            if BucketDominates(C_s, C) and C_s ≠ C then
                skip ← true
                break
        if skip then continue
        //--- Accept label: add to settled set at v ---
        Settled[v] ← Settled[v] ∪ {C}
        //--- Remove entries in Settled[v] now dominated by C ---
        for each C_s ∈ Settled[v] do
            if BucketDominates(C, C_s) and C_s ≠ C then
                Settled[v] ← Settled[v] \ {C_s}
        //--- If target reached, collect solution and do not expand ---
        if v = t then
            P ← P ∪ {(path, C[1..k])}
            continue
        //--- Expand: generate new labels along all outgoing edges ---
        for each edge e = (v, u) ∈ E do
            for i ← 1 to k do
                C_new[i] ← C[i] + cost(e, i)
            //--- Pre-insertion pruning: discard if dominated at u ---
            prune ← false
            for each C_s ∈ Settled[u] do
                if BucketDominates(C_s, C_new) then
                    prune ← true
                    break
            if not prune then
                L_new ← (u, C_new[1..k], path ∪ {e})
                Insert(PQ, L_new)
    //--- Post-processing: remove true-cost dominance from collected solutions ---
    for each (path_a, C_a) ∈ P do
        for each (path_b, C_b) ∈ P do
            if C_a ≠ C_b and TrueDominates(C_a, C_b) then
                P ← P \ {(path_b, C_b)}
    return P


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
- Configure the total approximation budget ε_total using a slider
- Select between 2 and 4 optimization objectives
- Automatically computes the adaptive ε-vector from the edge-cost distribution
- Executes the Adaptive ε-Vector FPTAS using bucket-based dominance pruning
- Displays execution statistics including labels explored, Pareto paths found, and runtime
- Highlights all approximate Pareto-optimal paths directly on the graph
- Ranks the resulting Pareto solutions using the Gain-Weighted Balance Score and identifies the best balanced solution

**Solution Visualization**
- Displays all approximate Pareto-optimal routes ranked by overall balance score
- Highlights the highest-ranked balanced solution on the graph
- Shows the complete node sequence and accumulated objective costs for every solution
- Reports the computed adaptive ε-vector used during execution

**Complexity Reference**
A comparison table of this algorithm against Exact MOSP, Parametric Dijkstra, and Uniform ε-FPTAS, with the formal time and solution-size bounds displayed. (Future)
