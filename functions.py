import galois
import sympy
import numpy as np
from sympy import factorint
import numpy as np
import numba
from numba import njit, prange
import multiprocessing as mp
from functools import partial
from sympy import primerange
from functools import reduce

def QR3primes(N: int) -> list[int]:
    return [p for p in primerange(0, N + 1) if p % 12 in (1, 11)]

def QR2primes(N: int) -> list[int]:
    return [p for p in primerange(0, N + 1) if p % 8 in (1, 7)]

def min_multiplicative_order(l, p):
    N = p - 1
    m = N
    factors = factorint(N)

    for q in factors:
        while m % q == 0 and pow(l, m // q, p) == 1:
            m //= q
    return m

def primitive_pth_root(p: int, l: int):
    if p == l:
        raise ValueError("Need p != l")
    if not galois.is_prime(p) or not galois.is_prime(l):
        raise ValueError("p and l must be prime")

    #minimal m 
    m = min_multiplicative_order(l, p)
    # m = p - 1
    if m < 78:
        GF = galois.GF(l**m) # here we could also just use p-1 if minimal m becomes a bottleneck
        g = GF.primitive_element
    else:
        print("starting poly search")
        g = galois.irreducible_poly(l, m)
        print("found")
    alpha = g ** ((l**m - 1) // p)

    # sanity-check
    assert alpha ** p == 1 and alpha != 1
    return alpha, m, GF

def polynom(block_length: int, version: str, dimension = 3, loud = True):
    modulus = block_length % 12
    if (not (modulus == 1 or modulus == 11)):
        print("Invalid block length!")

    GFp = galois.GF(block_length)
    GFl = galois.GF(dimension)
    primitive = GFp.primitive_element # more than one may exist this just grabs one
    if loud:
        print(primitive)
    quadratic_generator = primitive ** 2
    residues = { int(quadratic_generator**k) for k in range((block_length - 1) // 2) }
    field_elements = set(range(block_length))
    nonresidues = field_elements - residues - set([0])

    if loud:
        print(f"residues are {residues}")
        print(f"nonresidues are {nonresidues}")

    root_unity, min_multiplicative_order, GF_ext = primitive_pth_root(block_length,dimension)
    if loud:
        print(f"extension field used = {dimension}^{min_multiplicative_order}")
        print(f"root of unity used: {root_unity}")

    match version.startswith("Q"):
        case True:
            #print("Creating code Q")
            roots = [root_unity**q for q in residues]
            if loud:
                print(f"quadratic roots {roots}")
        case _:
            #print("Creating code N")
            roots = [root_unity**n for n in nonresidues]
            if loud:
                print(f"nonquadratic roots {roots}")
    
    match version[1:]:
        case "bar":
            #print("bar")
            generator = galois.Poly.Roots((roots + [GF_ext(1)]), field=GF_ext)
        case _:
            generator = galois.Poly.Roots(roots, field=GF_ext)
    
    coeffs = generator.coefficients(block_length)
    coeffs = np.array(coeffs,dtype=np.int8)
    if loud:
        print(f"coefficients are: {coeffs}")

    return coeffs

def generator_matrix(block_length, version, dimension = 3, loud = True):
    coeffs = polynom(block_length, version, dimension, loud=loud)
    G = np.zeros((block_length, block_length), dtype=np.int8)
    for i in range(block_length):
        G[i] = np.roll(coeffs,i)
    if loud:
        print(G)
    return G

# This is in the macwilliams and sloane convention where the row space is the codespace

# also check why self duality needed

def extended_generator_matrix(block_length, base = "Qbar", dimension = 3, loud = True):
    Qbar = generator_matrix(block_length, base, loud=loud)
    if loud:
        print("Warning: current implementation does not account for value of y in corner digit")
    final = -block_length % dimension
    Qbar0 = np.hstack((
            Qbar,
            np.zeros((Qbar.shape[0], 1), dtype=np.int8) 
    ))
    
    extra_row = np.ones((1, Qbar0.shape[1]), dtype=np.int8)
    extra_row[0, -1] = final
    Qbarext = np.vstack((
        Qbar0,
        extra_row
    ))
    return Qbarext

@njit(parallel=True, fastmath=True)
def delta_divisible(M_int8: np.ndarray, delta: int) -> bool:
    n, m = M_int8.shape
    M = M_int8.astype(np.int64)      # promotes once inside JIT

    # captured outer-scope var
    violated = False

    for i in prange(n):
        if violated:
            # some other thread already found a failure
            continue
        for j in range(i, n):
            s = 0
            for k in range(m):
                s += M[i, k] * M[j, k]
            if (s) % delta != 0:
                violated = True
                break
        if violated:
            break

    return not violated

def plain_delta_divisible_np(M_int8, delta, loud = False):
    M64 = M_int8.astype(np.int64, copy=False)
    gram = (M64 @ M64.T) % 3
    if loud:
        print(gram)
    return not np.any(gram)

def delta_divisible_np(M_int8, delta, loud = False):
    M64 = M_int8.astype(np.int64, copy=False)
    M64 = (M64 ** 2) % 3
    gram = (M64 @ M64.T) % 3
    if loud:
        print(gram)
    return not np.any(gram)

# def delta_divisible_np(M_int8, delta):
#     gen_divisible = reduce((lambda acc, x: acc and (np.count_nonzero(x) % delta) == 0), M_int8, True) # has some redundancy for cyclic generators
#     if gen_divisible:
#         return delta_divisible_prod(M_int8, delta)
#     else:
#         return False

def min_distance_polynom_brute(generator):
    block_length = len(generator)
    best    = np.count_nonzero(generator)
    target  = int(np.sqrt(block_length))

    def _brute_helper(idx, cw):
        nonlocal best
        w = np.count_nonzero(cw)
        if w < best:
            best = w
            if best <= target:
                raise StopIteration
        if idx == block_length:
            return
        _brute_helper(idx + 1, cw)
        _brute_helper(idx + 1, (cw + np.roll(generator,1)) % 3)
        _brute_helper(idx + 1, (cw + 2*np.roll(generator,1)) % 3)

    try:
        _brute_helper(0, generator)
    except StopIteration:
        pass
    return best

def min_distance_brute(generator):
    num_gen ,block_length = generator.shape
    best    = np.count_nonzero(generator[0])
    target  = int(np.sqrt(block_length))

    def _brute_helper(idx, cw):
        nonlocal best
        nonlocal generator
        w = np.count_nonzero(cw)
        if w < best:
            best = w
            if best <= target:
                raise StopIteration
        if idx >= num_gen:
            return
        _brute_helper(idx + 1, cw)
        _brute_helper(idx + 1, (cw + generator[idx]) % 3)
        _brute_helper(idx + 1, (cw + 2*generator[idx]) % 3)

    try:
        _brute_helper(1, generator)
    except StopIteration:
        pass
    return best

import numpy as np
import multiprocessing as mp
from functools import partial

__all__ = [
    "cyclic_generator_matrix",
    "min_distance_backtracking",
]

def cyclic_generator_matrix(g_coeffs: np.ndarray, n: int | None = None) -> np.ndarray:
    g_coeffs = np.asarray(g_coeffs, dtype=np.int8)
    if n is None:
        n = len(g_coeffs)
    deg = len(g_coeffs) - 1
    k = n - deg
    if k <= 0:
        raise ValueError("Generator polynomial degree must be < code length n.")
    row0 = np.zeros(n, dtype=np.int8)
    row0[: len(g_coeffs)] = g_coeffs % 3

    G = np.empty((k, n), dtype=np.int8)
    for i in range(k):
        G[i] = np.roll(row0, i)
    return G

class _DFSState:
    __slots__ = ("G", "n", "k", "codeword", "best")

    def __init__(self, G: np.ndarray, best: int):
        self.G = G
        self.n, self.k = G.shape[1], G.shape[0]
        self.codeword = np.zeros(self.n, dtype=np.int8)
        self.best = best
    def _dfs(self, row_idx: int):
        if row_idx == self.k:
            wt = int(np.count_nonzero(self.codeword))
            if 0 < wt < self.best.value:
                self.best.value = wt
            return
        if np.count_nonzero(self.codeword) >= self.best.value:
            return

        row = self.G[row_idx]
        saved = self.codeword.copy()
        for a in (0, 1, 2):
            if a:
                np.add(self.codeword, a * row, out=self.codeword, casting="safe")
                np.mod(self.codeword, 3, out=self.codeword)
            self._dfs(row_idx + 1)
            self.codeword[:] = saved
    def __call__(self, prefix: tuple[int, ...]):
        for i, a in enumerate(prefix):
            if a:
                np.add(self.codeword, a * self.G[i], out=self.codeword, casting="safe")
                np.mod(self.codeword, 3, out=self.codeword)
        self._dfs(len(prefix))
        return self.best.value


def _generate_prefixes(k: int, prefix_len: int):
    if prefix_len == 0:
        yield ()
        return
    digits = range(3)
    total = 3 ** prefix_len
    for i in range(total):
        x = i
        pref = []
        for _ in range(prefix_len):
            x, r = divmod(x, 3)
            pref.append(r)
        yield tuple(reversed(pref))


def min_distance_backtracking(G: np.ndarray, *, prefix_len: int | None = None, jobs: int | None = None) -> int:
    k, n = G.shape[0], G.shape[1]
    if prefix_len is None:
        prefix_len = min(3, k // 2)
    if not (0 <= prefix_len <= k):
        raise ValueError("prefix_len must be in [0,k].")

    with mp.Manager() as mgr:
        best = mgr.Value("i", n + 1)
        state_factory = partial(_DFSState, G, best)
        prefixes = list(_generate_prefixes(k, prefix_len))

        if jobs is None:
            jobs = mp.cpu_count()
        with mp.Pool(jobs, initializer=lambda: None) as pool:
            list(pool.imap_unordered(state_factory, prefixes))
        return best.value
    
    g = np.array(args.poly, dtype=np.int8) % 3
    G = cyclic_generator_matrix(g, n=args.n)
    t0 = time.perf_counter()
    d = min_distance_backtracking(G, prefix_len=args.prefix, jobs=args.jobs)
    dt = time.perf_counter() - t0
    print(f"[n={G.shape[1]}, k={G.shape[0]}]  d_min = {d}  (elapsed: {dt:.2f} s)")

def triorthogonal(M_int8: np.ndarray, loud=True):
    # M64 = M_int8.astype(np.int64, copy=False)
    # M64 = (M64 ** 2) % 3
    if loud:
        print(M64)
    print("Checking second triorthogonality condition not done yet")
    return delta_divisible_np(M64, 3, loud=loud)

def balanced(generator):
    def balance(generator):
        n1 = (generator == 1).sum()
        n2 = (generator == 2).sum()
        return (n1 == n2)
    
    if generator.ndim == 1:
        return balance(generator)
    else:
        for basis in generator:
            if not balance(basis):
                return False
    
    return True

    

    


