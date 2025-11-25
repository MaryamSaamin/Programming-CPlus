# fit_export_icn.py
import numpy as np
from scipy.stats import norm
from numpy.linalg import lstsq
import math
import sys

# ---------- User parameters ----------
# Central rational degrees (P: m, Q: n) -- Q has implicit leading 1 (b0 = 1)
m = 6
n = 6

# Tail rational degrees (C: p, D: q), D has implicit leading 1
p = 6
q = 6

# Join point (should match header)
x_low = 0.02425
x_high = 1.0 - x_low

# Sample counts
N_central = 4000    # central region linear samples (denser helps)
N_tail_log = 2000   # tail region log samples per side

# regularization (ridge)
lam = 1e-12

# weight parameters (upweight near join)
join_scale = 0.002  # width around x_low to upweight (absolute)
join_weight = 100.0

# ---------- Helper functions ----------
def build_central_system(xs, z_star, m, n, weights=None):
    # xs: array of x in central domain
    u = xs - 0.5
    r = u*u
    K = len(xs)
    A = np.zeros((K, (m+1) + n))  # a0..am, b1..bn
    # columns for a coefficients: u * r^i
    for i in range(m+1):
        A[:, i] = u * (r**i)
    # columns for b coefficients: -z_star * r^j  (j=1..n)
    for j in range(1, n+1):
        A[:, m + j - 0] = -z_star * (r**j)
    if weights is not None:
        W = np.sqrt(weights)
        A = A * W[:,None]
        y = (z_star * W)
    else:
        y = z_star
    # Solve (A^T A + lam I) theta = A^T y
    if lam > 0:
        # normal equations with ridge for stability
        ATA = A.T @ A
        rhs = A.T @ y
        ATA += lam * np.eye(ATA.shape[0])
        theta = np.linalg.solve(ATA, rhs)
    else:
        theta, *_ = lstsq(A, y, rcond=None)
    a = theta[:m+1]
    b_tail = theta[m+1:]
    # b array should be [1.0, b1, b2, ...]
    b = np.concatenate(([1.0], b_tail))
    return a, b

def build_tail_system(xs_tail, z_star_tail, p, q, left=True, weights=None):
    # xs_tail: array of x in tail region (very small) - for left tail use x small, for right mirror
    # transform to t = sqrt(-2 log m) where m = min(x,1-x) = x if left
    mvals = np.minimum(xs_tail, 1.0 - xs_tail)
    t = np.sqrt(-2.0 * np.log(mvals))
    K = len(xs_tail)
    A = np.zeros((K, (p+1) + q)) # c0..cp, d1..dq
    for i in range(p+1):
        A[:, i] = (t**i)
    for j in range(1, q+1):
        A[:, p + j - 0] = -z_star_tail * (t**j)
    if weights is not None:
        W = np.sqrt(weights)
        A = A * W[:,None]
        y = z_star_tail * W
    else:
        y = z_star_tail
    if lam > 0:
        ATA = A.T @ A
        rhs = A.T @ y
        ATA += lam * np.eye(ATA.shape[0])
        theta = np.linalg.solve(ATA, rhs)
    else:
        theta, *_ = lstsq(A, y, rcond=None)
    c = theta[:p+1]
    d_tail = theta[p+1:]
    d = np.concatenate(([1.0], d_tail))
    return c, d

def fmt_cpp_array(name, arr, per_line=4):
    s = f"constexpr double {name}[] = {{\n"
    for i,v in enumerate(arr):
        s += f"    {v: .17e}"
        if i+1 < len(arr):
            s += ","
        if (i+1) % per_line == 0:
            s += "\n"
        else:
            s += " "
    s += "\n};\n"
    return s

# ---------- Build nodes ----------
# central: linear grid between x_low and x_high (but include a fraction of tails too)
xs_central = np.linspace(x_low*0.5, 1.0 - x_low*0.5, N_central)
z_central = norm.ppf(xs_central)

# weights: upweight points near join
weights_c = np.ones_like(xs_central)
# gaussian bump near x_low and x_high
weights_c += join_weight * np.exp(-((xs_central - x_low)/join_scale)**2)
weights_c += join_weight * np.exp(-((xs_central - (1.0-x_low))/join_scale)**2)

# ---------- Fit central ----------
a_coeffs, b_coeffs = build_central_system(xs_central, z_central, m, n, weights=weights_c)

print("// ---- central coefficients (P: a0..am, Q: b0..bn where b0=1) ----")
print(fmt_cpp_array("a_central", a_coeffs))
# print only b1..bn for Horner in header (we will append +1.0 in code)
print(fmt_cpp_array("b_central_tail", b_coeffs[1:]))  # b1..bn

# ---------- Tail fitting (left tail) ----------
# choose log-spaced nodes from 1e-16 to x_low*0.9 (left tail)
xs_left = np.unique(np.concatenate([
    np.logspace(-16, math.log10(max(1e-16, x_low*0.9)), N_tail_log),
    np.linspace(max(1e-16, x_low*0.0001), x_low*0.9, 500)
]))
xs_left = xs_left[xs_left < x_low]
z_left = norm.ppf(xs_left)

# weights for tail: upweight extreme small and near-join
weights_left = np.ones_like(xs_left)
weights_left += 100.0 * np.exp(-((xs_left - x_low)/ (x_low*0.1))**2)
weights_left += 200.0 * (np.log10(xs_left + 1e-300) < -8)  # extra weight in extreme tail

c_left, d_left = build_tail_system(xs_left, z_left, p, q, left=True, weights=weights_left)

print("// ---- left tail coefficients (C: c0..cp, D: d0..dq where d0=1) ----")
print(fmt_cpp_array("c_tail", c_left))
print(fmt_cpp_array("d_tail_tail", d_left[1:]))  # d1..dq

# ---------- Right tail (mirror by symmetry) ----------
# For right tail we can mirror left coefficients: g(1-x) = -g(x)
print("// Note: right-tail coefficients can be the same as left (use sign on output).")

# ---------- Done ----------
print("// Paste a_central into central_rational numerator (a0..am).")
print("// Paste b_central_tail (b1..bn) into denominator coefficients (and use +1.0 tail).")
print("// Paste c_tail into tail numerator (c0..cp) and d_tail_tail (d1..dq) into tail denominator.")
