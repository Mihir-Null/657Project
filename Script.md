1\. Magic State Distillation Overview

As we know, Fault-Tolerant Computing is key to producing reliable quantum systems to perform calculations.

I'm going to mention a few terms here that will be covered in the following slides, so don't worry if there is anything unfamiliar to you.

In many stabilizer codes, all Clifford gates-like Hadamard, CNOT, and phase-can be implemented transversally, meaning errors don't spread within a block and that the qubits are acted upon separately.

But the T gate is different. If we apply T directly to encoded data, error can spread across the block, which breaks fault tolerance.

So the key idea of magic state distillation is not applying the T gate directly to data, but simulating its effect using only Clifford operations and measurements."

2\. Motivation

Here's the general motivation of why we care about Magic State Distillation:

Obviously, universal quantum computation requires a non-Clifford gate like T- Clifford gates alone are not enough.

Second, magic state distillation lets us take many imperfect, noisy copies of a magic state and refine them into fewer but much higher-quality ones.

And third, the entire process uses only measurements and Clifford operations, which we already know how to implement fault-tolerantly.

So in the end, we get the net effect of a T gate without ever applying a physical T gate directly to data.

3\. Transversal Gates

So, transversal can vary a little in meaning depending on the context but it usually refers to:

A gate that acts independently on corresponding qubits across multiple code blocks.  

Mathematically, it factors as a tensor product of small gates - one per qubit position.

The key idea is that no part of the transversal operation ever touches two qubits in the same block.

That's what makes it fault-tolerant: a single-qubit error cannot spread into a multi-qubit error inside the code.

4\. Example: Transversal CNOT

Here's a classic example: transversal CNOT. Each physical qubit in the control block interacts only with the corresponding qubit in the target block.

If one physical CNOT gate is faulty, it can only affect one pair of qubits - not the entire block.

This containment of errors is exactly what we want for fault-tolerant quantum computing.

5\. Fault-Tolerant Properties of Transversal Gates

There are two key properties that make transversal gates valuable. First, they avoid error blow-up inside a block-single faults remain localized - as explained in the previous slide

Second, when an error passes through a transversal gate, it stays an error. In other words, errors are just mapped onto each other, so bit-flip or phase-flip error may change type, but it never turns into a messy multi-qubit error.

A simple example is just the Hadamard gate, since conjugation maps X to Z and Z to X.

6\. Constraints: Eastin-Knill Theorem

There is, however, a fundamental limitation given by the Eastin-Knill theorem, which tells us that no quantum error-correcting code can implement a universal gate set using only transversal gates.

In practice, this means that transversal logical gates are restricted to a finite level of the Clifford hierarchy. So, non-Clifford gates like T cannot be implemented transversally in standard codes, which is why we need magic states.

7\. CSS codes

CSS codes are built from two nested classical linear codes, and the difference in their dimensions tells us how many logical qubits are encoded.

The most important feature is that X and Z errors are handled independently. One classical code detects bit-flip errors, and the other detects phase-flip errors.

This separation of error types makes CSS codes fit for magic state distillation

8\. Self-Dual Codes

A CSS code is self-dual when its underlying classical code equals its dual.

The dual of a code is just the set of vectors that annihilate it - its parity checks. So if a code is self-dual, it means it is its own annihilator.

In a CSS code, that makes the X- and Z-stabilizers come from the same structure. Since Hadamard swaps X and Z, that symmetry is exactly what allows something like the Hadamard to be transversal.

9\. Doubled QR Codes

First, a self-dual, doubly-even CSS code, where _doubly-even_ means every stabilizer has weight divisible by 4. This gives strong X/Z symmetry and transversal Clifford gates.

Second, a QR-derived doubly-even CSS code, where quadratic residue (QR) structure gives high distance at small block sizes.

We then apply a doubling map to combine these two codes.

The result is a weakly triply-even code, meaning certain key stabilizers have weight divisible by 8. This is exactly the condition needed for safe transversal TTT-gate distillation.

These codes achieve transversal T, high distance, and low overhead, which makes them very attractive for practical magic state distillation.

10\. Codes Diagram

So, here is the code diagram of the code families we used in our implementation, which is taken from the Jain and Albert paper. From left to right are the Doubly-Even Codes, the weak triply-even, and the triorthogonal codes.

11\. Bravyi-Haah Magic State Distillation

Finaly, The Bravyi-Haah protocol is the standard magic state distillation framework. It uses triorthogonal - in our case also weakly triply-even - CSS codes to distill high-fidelity T states.

We start with many noisy copies of the T state, each with some physical error rate p. The protocol applies Clifford operations and Pauli measurements on the encoded blocks, and we get a smaller number of magic states whose error rate is suppressed

p to the k, where k depends on the code's distance.

# Motivation and Approach

## Methodology / Bravyi–Haah Protocol for TE* and Triorthogonal Codes

* Remind audience: use same Bravyi–Haah distillation circuit, just w/ different codes (TE* / QR-based vs standard triorthogonal).
* Inputs:
  * The X-stabilizer matrix (H_X) (rows = X checks).
  * The logical Z operator as a binary vector (z_{\text{log}}).
  * Aphysical T-state error rate (p).
  * Noise model: i.i.d. (Z)-type noise on each injected (|T\rangle).
* Error patterns are binary vectors (e \in {0,1}^n) describing which injected (|T\rangle) are faulty.

* Two key probs per code:
  * Success probability (s(p) = \Pr[H_X e^T = 0]): fraction of blocks that pass syndrome checks.
  * Output error probability (p_{\text{out}}(p) = \Pr[z_{\text{log}}\cdot e = 1 \mid \text{accepted}]): chance distilled magic state is still bad given we accepted.
* don't derive closed forms, estimate tnumerically for each candidate QR-based TE code.

---

## Hypothesis

* expect QR-based TE codes to have better finite block len performance than:
  * Standard Bravyi–Haah triorthogonal code families.
  * Generic doubled self-dual constructions without QR structure.

* Why
  * QR structure gives high distance at relatively small block size n.
  * Weight distribution of QR codes leads to low-weight X-checks that still satisfy the mod-8 (weakly triply-even) divisibility constraints.
  * should improve error-suppression exponent (\alpha) for small (n).
* asymptotically meh:
  * Theory still suggests MSD exponent (\gamma \to 2) at large (n), same as Bravyi–Haah.
  * key question - do we win in the “practical” regime* (say (n \lesssim 30)–100)?

---

## Distillation Yield

* Define yield (Y(p)) as “magic states out per magic state in”:
  [
  Y(p) = \frac{k \cdot s(p)}{n}
  ]
  where:

  * (n) = number of input T states per block,
  * (k) = number of logical magic states produced per block,
  * (s(p)) = success probability.
* For QR TE* codes we study:

  * Typically (k = 1) (\Rightarrow) (Y(p) = s(p)/n).
* Intuition to explain:

  * Smaller (n) but still decent distance ⇒ even if (s(p)) is modest, (s(p)/n) can be quite good.
  * Compare against:

    * Standard Bravyi–Haah triorthogonal families with (n = 3k + 8).
    * Previously used doubled self-dual codes.
* Emphasize: yield is the practical cost metric—how many raw T-states we burn per good one.

---

## Scaling With Code Length

* Summarize Jain–Albert scaling result:

  * Code distance grows roughly as
    [
    d(n) \sim \Theta(\sqrt{n})
    ]
    for both TE* and associated triorthogonal families.
* For a BH-style protocol with (Z) noise:

  * Logical error behaves like
    [
    p_{\text{out}}(p) \sim C, p^{\alpha},\quad \alpha \approx d_Z
    ]
    where (d_Z) = minimum weight undetected Z-logical.
  * So (\alpha(n) \sim \Theta(\sqrt{n})) as well.
* Asymptotic exponent of magic state distillation cost:

  * Standard argument gives (\gamma_n = \log_{\alpha}(n/k) \to 2) as (n \to \infty).
* Key message:

  * Asymptotically do not beat Bravyi–Haah’s scaling.
  * Finite may still be better, which is what actual devices care about.

---

# Simulation

## Algorithm (Monte Carlo)

* For each code and each chosen physical error rate (p):

  1. Load code: (H_X), (z_{\text{log}}), block length (n), and any metadata (distance, type).
  2. Sample many error patterns (e \sim \text{Bernoulli}(p)^n):

     * Each bit of (e) says “this T-state suffered a Z error” or not.
  3. Check acceptance condition:

     * Compute syndrome (H_X e^T) over GF(2).
     * If syndrome is zero ⇒ block accepted, otherwise toss.
  4. For accepted blocks, evaluate log parity (z_{\text{log}} \cdot e):
     * 0 ⇒ distilled state is good.
     * 1 ⇒ distilled state is faulty.
* From these samples, est:

  * Success rate (s(p) = #\text{accepted} / N).
  * Conditional logical error (p_{\text{out}}(p) = #\text{“bad & accepted”} / #\text{accepted}).
  * Yield (Y(p) = s(p)/n) (since (k = 1)).
* Note explicitly: this is a Monte Carlo simulation of the BH protocol under i.i.d. Z noise.

---

## Results

* Explain axes before talking:

  * x-axis: physical T-state error rate (p) (log scale).
  * y-axis**: yield (Y(p)) (good magic states per input T-state).
* Describe what each curve corresponds to
* Qualitative observations from the plot:
  * At realistic physical error rates (e.g. (10^{-3})–(10^{-2})), QR-based TE* codes show higher yield than the baseline.
  * As (p) decreases, both protocols improve, but TE* curve stays above the Bravyi–Haah curve for all simulated points.
  * For very small (p), yields begin to flatten since rejection events become rare.
* Emphasize: even modest separation in yield translates to significant savings in T-coun* over large algorithms.

---

### Analysis / Contextualizing Results

* Link back to the hypothesis:

  * Data is consistent with: QR-based TE* codes provide better finite-size overhead for the error rates and sizes we can simulate.
* Interpretation:

  * High dist at smol (n) gives strong error suppression w/out huge blocks.
  * Div / TE* structure allows transversal diagonal gates while keeping checks relatively low-weight.
* Limitations of current simulations:

  * Only small block lengths explored so far, exp cost of scaling codes.
  * Only simple i.i.d. Z noise model—no correlated noise, no biased noise yet.
* Takeaway:

  * look promising for near-term architecture w/small low error n.
  * To fully validate, we’d need larger-n simulations and perhaps analytic approximations for (s(p)) and (p_{\text{out}}(p)).

---

# Conclusion

## Hypothesis vs Results

* Re-state original hypothesis in one sentence:

  * QR-based TE* codes should outperform Bravyi–Haah triorthogonal families for realistic, small-to-medium block sizes, even though the asymptotic exponent (\gamma) is the same.
* Summarize what we actually saw:

  * Simulation results support that TE* codes achieve higher yield in the tested regime.
  * Error suppression is consistent with expected exponents based on (d_Z).
* Strengths of the construction:

  * Good distance at low (n) from quadratic residue structure.
  * Transversal diagonal gates guaranteed by weak triply-even / divisibility conditions.
  * Compatible with standard Bravyi–Haah protocol—no exotic circuits needed.
* Limitations / caveats:

  * Asymptotically, (\gamma \to 2), so we don’t beat the best-known scaling.
  * Most codes we used have (k = 1) logical magic state, which limits asymptotic rate.
  * Simulations so far are numerical and finite-sample, not closed-form.

---

# Potential Future Work

# What’s Next?

* Multi-logical TE* / triorthogonal codes:

  * Construct codes w/ (k > 1) while preserving TE* / weak triply-even properties.
  * Could improve rate and asymptotic cost.
* Sharper analysis of (d_Z) and logical structure:

  * Characterize Z-logical operators more precisely.
  * better analytic bounds on (\alpha) and thus on (p_{\text{out}}(p)).
* Closed-form yield formulas for QR-based families:

  * Move from pure simulation to semi-analytic expressions for (s(p)), (p_{\text{out}}(p)).
* Generalize beyond qubits:

  * Explore qutrit or higher-dimensional analogues of doubled QR / TE* constructions.
  * Connect to diagonal Clifford hierarchy conditions for prime-dimensional qudits.
* Code switching and noise models:

  * Investigate code-switching protocols between TE* for > magic state perf.
  * Study performance under biased noise or more realistic device noise models.
* Concluding

  * QR-based TE* codes look like a practical upgrade path for magic state factories, especially in the finite-size, near-term regime.
