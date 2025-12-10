import numpy as np
import itertools
import matplotlib.pyplot as plt
from dataclasses import dataclass
from typing import Dict, Tuple, List

# ============================================================
# 0. GF(2) linear algebra helpers
#    (used to go from C1, C2 -> stabilizers + logicals)
# ============================================================

def gf2_rref(A: np.ndarray) -> Tuple[np.ndarray, List[int]]:
    """
    Row-reduced echelon form over GF(2).
    Returns (R, pivot_cols) where R is the RREF of A (mod 2).
    """
    A = (A.copy().astype(np.uint8) & 1)
    m, n = A.shape
    row = 0
    pivot_cols: List[int] = []

    for col in range(n):
        # Find a pivot row at or below current row
        pivot_row = None
        for r in range(row, m):
            if A[r, col] == 1:
                pivot_row = r
                break
        if pivot_row is None:
            continue

        # Swap to bring pivot to current row
        if pivot_row != row:
            A[[row, pivot_row]] = A[[pivot_row, row]]

        pivot_cols.append(col)

        # Eliminate in all other rows
        for r in range(m):
            if r != row and A[r, col] == 1:
                A[r, :] ^= A[row, :]

        row += 1
        if row == m:
            break

    return A, pivot_cols


def gf2_row_space_basis(A: np.ndarray) -> np.ndarray:
    """
    Return a row-basis of the row space of A over GF(2),
    as an array of independent rows.
    """
    if A.size == 0:
        return np.zeros((0, A.shape[1]), dtype=np.uint8)
    R, _ = gf2_rref(A)
    basis_rows = []
    for r in range(R.shape[0]):
        if np.any(R[r, :]):
            basis_rows.append(R[r, :].copy())
    if len(basis_rows) == 0:
        return np.zeros((0, A.shape[1]), dtype=np.uint8)
    return np.stack(basis_rows, axis=0).astype(np.uint8)


def gf2_nullspace(A: np.ndarray) -> np.ndarray:
    """
    Nullspace of A over GF(2): { v | A v^T = 0 (mod 2) }.
    Returns basis vectors as rows.
    """
    A = (A.copy().astype(np.uint8) & 1)
    m, n = A.shape
    R = A.copy()
    pivot_col_for_row = [-1] * m
    row = 0

    # Compute an RREF-like form to identify pivot columns
    for col in range(n):
        pivot_row = None
        for r in range(row, m):
            if R[r, col] == 1:
                pivot_row = r
                break
        if pivot_row is None:
            continue
        if pivot_row != row:
            R[[row, pivot_row]] = R[[pivot_row, row]]
        pivot_col_for_row[row] = col
        # Eliminate in all other rows
        for r in range(m):
            if r != row and R[r, col] == 1:
                R[r, :] ^= R[row, :]
        row += 1
        if row == m:
            break

    pivot_cols = {c for c in pivot_col_for_row if c != -1}
    free_cols = [j for j in range(n) if j not in pivot_cols]

    basis = []
    for free in free_cols:
        v = np.zeros(n, dtype=np.uint8)
        v[free] = 1
        # For each pivot row, solve for pivot entries
        for i in range(row):
            pc = pivot_col_for_row[i]
            if pc == -1:
                continue
            if R[i, free] == 1:
                v[pc] ^= 1
        basis.append(v)

    if len(basis) == 0:
        return np.zeros((0, n), dtype=np.uint8)
    return np.stack(basis, axis=0).astype(np.uint8)


def gf2_is_in_span(v: np.ndarray, basis: np.ndarray) -> bool:
    """
    Check if v is in span(basis) over GF(2).
    Assumes 'basis' rows are independent (as from gf2_row_space_basis).
    """
    v = (v.astype(np.uint8) & 1)
    if basis.size == 0:
        return not np.any(v)
    B = np.vstack([basis, v])
    R, _ = gf2_rref(B)
    # Rank is number of non-zero rows
    rank_before = basis.shape[0]
    rank_after = sum(1 for r in range(R.shape[0]) if np.any(R[r, :]))
    return rank_after == rank_before


# ============================================================
# 1. CSS(X, C2; Z, C1⊥) construction (Jain–Albert style)
#    Input: C1, C2 classical codes with C2 ⊂ C1 ⊂ C2⊥
# ============================================================

@dataclass
class CSSCode:
    n: int
    k: int
    HX: np.ndarray        # X-stabilizers (rows)
    HZ: np.ndarray        # Z-stabilizers (rows)
    X_logicals: np.ndarray  # shape (k, n)
    Z_logicals: np.ndarray  # shape (k, n)


def build_css_T_from_C1_C2(C1: np.ndarray, C2: np.ndarray) -> CSSCode:
    """
    Given classical binary generator matrices C1, C2 with
      C2 ⊂ C1 ⊂ C2⊥  (CSS-T recipe in Jain–Albert),
    and C1 stacked as [C1,(1); C2] (their eq. (1) style),
    build the corresponding CSS code:

      CSS(X, C2; Z, C1⊥)

    and extract:
      - HX: X-stabilizers (span of C2),
      - HZ: Z-stabilizers (span of C1⊥),
      - X_logicals: basis for C1 / C2,
      - Z_logicals: basis for C2⊥ / C1⊥.

    Assumes k = dim(C1) - dim(C2) is small (often k=1 in the paper).
    """

    C1 = (C1.astype(np.uint8) & 1)
    C2 = (C2.astype(np.uint8) & 1)
    n = C1.shape[1]

    # Row bases for C1 and C2
    C1_basis = gf2_row_space_basis(C1)
    C2_basis = gf2_row_space_basis(C2)

    dim_C1 = C1_basis.shape[0]
    dim_C2 = C2_basis.shape[0]
    k = dim_C1 - dim_C2
    if k <= 0:
        raise ValueError("C1 must properly contain C2 to have k>0 logical qubits.")

    # --- X side: stabilizers + logical X ---
    # HX: basis for C2
    HX = C2_basis

    # Logical X: rows in C1 that extend the span of C2
    X_logicals = []
    span_C2_plus_X = HX.copy()
    for row in C1_basis:
        if not gf2_is_in_span(row, span_C2_plus_X):
            X_logicals.append(row.copy())
            # Extend the span
            span_C2_plus_X = gf2_row_space_basis(
                np.vstack([span_C2_plus_X, row])
            )
            if len(X_logicals) == k:
                break
    if len(X_logicals) != k:
        raise RuntimeError("Could not find k independent logical X rows from C1/C2.")
    X_logicals = np.stack(X_logicals, axis=0).astype(np.uint8)

    # --- Z side: stabilizers + logical Z ---
    # Z-stabilizers: basis for C1⊥
    HZ = gf2_nullspace(C1_basis)

    # Z-logicals: complement of HZ inside C2⊥
    C2_perp_basis = gf2_nullspace(C2_basis)  # basis of C2⊥
    # Start a basis from HZ (these are stabilizers)
    Z_stab_basis = gf2_row_space_basis(HZ)
    Z_logicals = []
    span_Z_total = Z_stab_basis.copy()

    for row in C2_perp_basis:
        if not gf2_is_in_span(row, span_Z_total):
            Z_logicals.append(row.copy())
            span_Z_total = gf2_row_space_basis(
                np.vstack([span_Z_total, row])
            )
            if len(Z_logicals) == k:
                break

    if len(Z_logicals) != k:
        raise RuntimeError("Could not find k independent logical Z rows from C2⊥/C1⊥.")

    Z_logicals = np.stack(Z_logicals, axis=0).astype(np.uint8)

    return CSSCode(
        n=n,
        k=k,
        HX=HX,
        HZ=HZ,
        X_logicals=X_logicals,
        Z_logicals=Z_logicals,
    )


# ============================================================
# 2. Distillation code wrapper and algorithms
# ============================================================

@dataclass
class DistillationCode:
    name: str
    HX: np.ndarray       # X-check matrix (stabilizers)
    z_log: np.ndarray    # logical Z (for the magic state)
    k: int = 1           # number of logical qubits (usually 1)


def leading_exponent(HX: np.ndarray, z_log: np.ndarray, w_max: int) -> Tuple[int, int]:
    """
    Find smallest weight w_min of an undetected harmful Z-error:
      - undetected: HX * e^T = 0 mod 2
      - harmful:    z_log · e = 1 mod 2

    Returns:
      w_min, A_wmin (multiplicity of such errors at that weight)
    """
    HX = HX.astype(np.uint8)
    z_log = z_log.astype(np.uint8)
    r, n = HX.shape

    w_min = None
    A_wmin = 0

    for w in range(1, w_max + 1):
        found_at_this_weight = False
        for positions in itertools.combinations(range(n), w):
            e = np.zeros(n, dtype=np.uint8)
            e[list(positions)] = 1

            syndrome = (HX @ e) & 1
            if np.any(syndrome):
                continue  # detected

            if (z_log @ e) & 1:
                # harmful logical flip
                found_at_this_weight = True
                if w_min is None:
                    w_min = w
                    A_wmin = 1
                elif w == w_min:
                    A_wmin += 1

        if found_at_this_weight:
            break

    return w_min, A_wmin


def sample_block(
    HX: np.ndarray,
    z_log: np.ndarray,
    p: float,
    N: int,
    k: int = 1,
    rng: np.random.Generator = None,
) -> Tuple[float, float, float]:
    """
    Monte Carlo estimate of:
      s_hat:     success probability (block accepted)
      p_out_hat: conditional logical error prob, given acceptance
      Y_hat:     yield per input magic state = k * s_hat / n

    Noise model: i.i.d. Z errors on each qubit with rate p.
    """
    if rng is None:
        rng = np.random.default_rng()

    HX = HX.astype(np.uint8)
    z_log = z_log.astype(np.uint8)
    r, n = HX.shape

    accepted = 0
    accepted_and_bad = 0

    for _ in range(N):
        e = rng.binomial(1, p, size=n).astype(np.uint8)

        syndrome = (HX @ e) & 1
        if np.any(syndrome):
            continue  # rejected

        accepted += 1
        bad = (z_log @ e) & 1
        if bad == 1:
            accepted_and_bad += 1

    s_hat = accepted / N
    p_out_hat = (accepted_and_bad / accepted) if accepted > 0 else 0.0
    Y_hat = (k * s_hat) / n
    return s_hat, p_out_hat, Y_hat


# ============================================================
# 3. Example: plug in C1, C2 -> CSS -> H_X, logical Z
#    Replace this with your actual QR / TE* / triorthogonal codes.
# ============================================================

def get_example_codes() -> Dict[str, DistillationCode]:
    """
    Define one or more example codes by providing classical C1, C2.
    In your actual workflow, you'll:
      - construct C1, C2 from your QR / TE* / triorthogonal recipes
        (as in Jain–Albert),
      - call build_css_T_from_C1_C2(C1, C2),
      - then wrap HX and Z_logicals[0] into DistillationCode.

    Below is a tiny toy example using a simple [[3,1,1]]-style CSS
    (not a good distillation code, purely for structure demonstration).
    """

    # --- Toy example classical codes ---
    # C2 ⊂ C1 ⊂ C2⊥, k = dim(C1) - dim(C2) = 1
    # Here we just hand-pick something over 3 bits:

    # C2: span{ 111 }  (X-stabilizer space)
    C2 = np.array([[1, 1, 1]], dtype=np.uint8)

    # C1: span{ 111, 100 }  (extends C2 by one extra generator)
    C1 = np.array([
        [1, 0, 0],
        [1, 1, 1],
    ], dtype=np.uint8)

    css = build_css_T_from_C1_C2(C1, C2)
    # css.HX: X-stabilizers
    # css.Z_logicals[0]: logical Z

    toy_code = DistillationCode(
        name="Toy CSS (n=3)",
        HX=css.HX,
        z_log=css.Z_logicals[0],
        k=css.k
    )

    # In your actual use, something like:
    # C1_QR, C2_QR = <build from functions.py / Magma / Jain-Albert recipe>
    # css_QR = build_css_T_from_C1_C2(C1_QR, C2_QR)
    # code_QR = DistillationCode(
    #     name="QR-based TE* [[n,1,d]]",
    #     HX=css_QR.HX,
    #     z_log=css_QR.Z_logicals[0],
    #     k=1
    # )

    return {
        toy_code.name: toy_code,
        # "QR-based TE* (n=..., d=...)": code_QR,
        # etc.
    }


# ============================================================
# 4. Driver: run distillation experiments and plot
# ============================================================

def run_distillation_experiments(
    codes: Dict[str, DistillationCode],
    p_values: List[float],
    N: int = 100_000,
    estimate_alpha: bool = True,
    w_max: int = 10,
):
    rng = np.random.default_rng()
    results = {}

    for code_name, code in codes.items():
        HX, z_log, k = code.HX, code.z_log, code.k
        r, n = HX.shape

        print(f"=== Code: {code_name} (n={n}, k={k}) ===")

        if estimate_alpha:
            w_min, A_wmin = leading_exponent(HX, z_log, w_max=w_max)
            print(f"  Leading undetected harmful weight w_min = {w_min}, A_wmin = {A_wmin}")
        else:
            w_min, A_wmin = None, None

        s_list = []
        pout_list = []
        Y_list = []

        for p in p_values:
            s_hat, p_out_hat, Y_hat = sample_block(HX, z_log, p, N=N, k=k, rng=rng)
            s_list.append(s_hat)
            pout_list.append(p_out_hat)
            Y_list.append(Y_hat)
            print(f"  p={p:.2e}  s≈{s_hat:.3e}  p_out≈{p_out_hat:.3e}  Y≈{Y_hat:.3e}")

        # Empirical slope alpha_emp from log-log fit
        log_p = np.log(np.array(p_values))
        log_pout = np.log(np.array(pout_list) + 1e-20)
        slope, intercept = np.polyfit(log_p, log_pout, 1)
        alpha_emp = slope
        print(f"  Empirical alpha (slope of log p_out vs log p) ≈ {alpha_emp:.3f}\n")

        results[code_name] = {
            "p": np.array(p_values),
            "s": np.array(s_list),
            "p_out": np.array(pout_list),
            "Y": np.array(Y_list),
            "w_min": w_min,
            "A_wmin": A_wmin,
            "alpha_emp": alpha_emp,
        }

    # --------------------------
    # Plot p_out vs p (log-log)
    # --------------------------
    plt.figure()
    for code_name, data in results.items():
        p = data["p"]
        p_out = data["p_out"]
        plt.loglog(p, p_out, marker="o", label=code_name)

    plt.xlabel("Physical error rate p")
    plt.ylabel("Output error p_out(p)")
    plt.title("Magic State Distillation Error Map")
    plt.legend()
    plt.grid(True, which="both", linestyle=":")

    # --------------------------
    # Plot yield vs p
    # --------------------------
    plt.figure()
    for code_name, data in results.items():
        p = data["p"]
        Y = data["Y"]
        plt.semilogx(p, Y, marker="o", label=code_name)

    plt.xlabel("Physical error rate p")
    plt.ylabel("Yield Y(p) = k s(p) / n")
    plt.title("Magic State Distillation Yield")
    plt.legend()
    plt.grid(True, which="both", linestyle=":")

    plt.show()

    return results


# ============================================================
# 5. Main
# ============================================================

if __name__ == "__main__":
    # Physical error grid (tune as needed)
    p_values = [1e-2, 5e-3, 1e-3, 5e-4, 1e-4]

    # Get your codes (replace toy with your QR / TE* ones)
    codes = get_example_codes()

    # Run and plot
    results = run_distillation_experiments(
        codes=codes,
        p_values=p_values,
        N=50_000,        # increase for smoother curves
        estimate_alpha=True,
        w_max=15         # set >= distance if you know it
    )
