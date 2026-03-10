/// Chebyshev lowpass prototype g-values and bandpass element calculation.
///
/// Reference: Matthaei, Young & Jones – "Microwave Filters, Impedance-Matching
/// Networks, and Coupling Structures" (1980).

/// Computes the Chebyshev lowpass prototype g-values for an N-th order filter.
///
/// # Arguments
/// * `n`          – Filter order (number of reactive elements).
/// * `ripple_db`  – Passband ripple in dB (e.g. 0.1, 0.5, 3.0).
///
/// # Returns
/// A `Vec<f64>` of length `n + 2` containing [g0, g1, …, gN, g_{N+1}].
/// g0 = 1 always; g_{N+1} is 1 for odd N, or coth²(β/4) for even N.
///
/// # Panics
/// Panics if `n == 0` or `ripple_db <= 0.0`.
pub fn chebyshev_g_values(n: usize, ripple_db: f64) -> Vec<f64> {
    assert!(n >= 1, "Filter order must be at least 1");
    assert!(ripple_db > 0.0, "Ripple must be positive (dB)");

    // ── Step 1: ripple parameter (not used directly but confirms valid input) ──
    let _epsilon = (10f64.powf(ripple_db / 10.0) - 1.0).sqrt();

    // ── Step 2: β and γ ───────────────────────────────────────────────────────
    // β = ln( coth( Lr·ln10 / 40 ) )
    let x = ripple_db * 10f64.ln() / 40.0; // Lr·ln10/40
    let beta = {
        let coth_x = x.cosh() / x.sinh();  // coth(x)
        coth_x.ln()
    };
    let gamma = (beta / (2.0 * n as f64)).sinh();

    // ── Step 3-4: pole angles and auxiliary sequences ─────────────────────────
    let mut a = vec![0.0f64; n + 1]; // 1-indexed: a[k] for k = 1..=n
    let mut b = vec![0.0f64; n + 1]; // 1-indexed: b[k] for k = 1..=n

    for k in 1..=n {
        let theta_k = (2 * k - 1) as f64 * std::f64::consts::PI / (2.0 * n as f64);
        a[k] = theta_k.sin();
        b[k] = gamma * gamma + theta_k.sin().powi(2);
    }

    // ── Step 5: g-values ──────────────────────────────────────────────────────
    let mut g = vec![0.0f64; n + 2]; // indices 0..=n+1

    g[0] = 1.0;
    g[1] = 2.0 * a[1] / gamma;

    for k in 2..=n {
        g[k] = (4.0 * a[k - 1] * a[k]) / (b[k - 1] * g[k - 1]);
    }

    // Termination element
    g[n + 1] = if n % 2 == 1 {
        1.0
    } else {
        let coth_b4 = (beta / 4.0).cosh() / (beta / 4.0).sinh();
        coth_b4.powi(2)
    };

    g
}

/// Bandpass element values derived from Chebyshev lowpass prototype g-values.
#[derive(Debug, Clone)]
pub struct BandpassElement {
    /// Element index (1-based)
    pub index: usize,
    /// Series arm inductance (H)  — from series resonator
    pub l_series: f64,
    /// Series arm capacitance (F) — from series resonator
    pub c_series: f64,
    /// Shunt arm inductance (H)   — from shunt resonator
    pub l_shunt: f64,
    /// Shunt arm capacitance (F)  — from shunt resonator
    pub c_shunt: f64,
}

/// Compute bandpass filter LC element values from Chebyshev g-values.
///
/// # Arguments
/// * `n`           – Filter order.
/// * `ripple_db`   – Passband ripple in dB.
/// * `f_low_hz`    – Lower -3 dB (passband edge) frequency in Hz.
/// * `f_high_hz`   – Upper -3 dB (passband edge) frequency in Hz.
///
/// # Returns
/// A tuple of:
/// * `Vec<f64>`            — the raw g-values [g0 … g_{N+1}]
/// * `Vec<BandpassElement>` — LC element values for each stage
pub fn chebyshev_bandpass_elements(
    n: usize,
    ripple_db: f64,
    f_low_hz: f64,
    f_high_hz: f64,
) -> (Vec<f64>, Vec<BandpassElement>) {
    assert!(f_low_hz > 0.0 && f_high_hz > f_low_hz,
        "Require 0 < f_low < f_high");

    let two_pi = 2.0 * std::f64::consts::PI;
    let w1 = two_pi * f_low_hz;
    let w2 = two_pi * f_high_hz;
    let w0 = (w1 * w2).sqrt();          // geometric centre frequency
    let delta = (w2 - w1) / w0;         // fractional bandwidth

    let g = chebyshev_g_values(n, ripple_db);

    let elements: Vec<BandpassElement> = (1..=n)
        .map(|k| {
            let gk = g[k];
            BandpassElement {
                index: k,
                // Series resonator  (from series LP element)
                l_series: gk / (delta * w0),
                c_series: delta / (w0 * gk),
                // Shunt resonator   (from shunt LP element)
                l_shunt: delta / (w0 * gk),
                c_shunt: gk / (delta * w0),
            }
        })
        .collect();

    (g, elements)
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────
#[cfg(test)]
mod tests {
    use super::*;

    /// Classic N=3, 0.5 dB ripple Chebyshev prototype.
    /// Known values: g1≈1.5963, g2≈1.0967, g3≈1.5963, g4=1.0
    #[test]
    fn test_g_values_order3_half_db() {
        let g = chebyshev_g_values(3, 0.5);
        assert!((g[0] - 1.0).abs() < 1e-9, "g0 must be 1");
        assert!((g[1] - 1.5963).abs() < 1e-3, "g1 ≈ 1.5963, got {}", g[1]);
        assert!((g[2] - 1.0967).abs() < 1e-3, "g2 ≈ 1.0967, got {}", g[2]);
        assert!((g[3] - 1.5963).abs() < 1e-3, "g3 ≈ 1.5963, got {}", g[3]);
        assert!((g[4] - 1.0).abs()    < 1e-9, "g4 = 1 for odd N");
    }

    /// N=2, 0.5 dB ripple: g4 should be coth²(β/4), not 1.
    #[test]
    fn test_g_values_order2_even_termination() {
        let g = chebyshev_g_values(2, 0.5);
        assert!((g[0] - 1.0).abs() < 1e-9);
        // For even N, g_{N+1} ≠ 1
        assert!((g[3] - 1.0).abs() > 0.01,
            "Even-order termination should differ from 1, got {}", g[3]);
    }

    #[test]
    fn test_bandpass_elements_count() {
        let (g, elems) = chebyshev_bandpass_elements(5, 0.1, 1e9, 1.1e9);
        assert_eq!(g.len(), 7);      // g0..g6
        assert_eq!(elems.len(), 5);  // one entry per reactive stage
    }

    #[test]
    fn test_series_shunt_resonance() {
        // Both resonators in each element must resonate at ω0
        let f_low = 900e6_f64;
        let f_high = 1100e6_f64;
        let w0 = 2.0 * std::f64::consts::PI * (f_low * f_high).sqrt();
        let (_, elems) = chebyshev_bandpass_elements(3, 0.5, f_low, f_high);
        for e in &elems {
            let w_series = 1.0 / (e.l_series * e.c_series).sqrt();
            let w_shunt  = 1.0 / (e.l_shunt  * e.c_shunt ).sqrt();
            assert!((w_series - w0).abs() / w0 < 1e-9,
                "Series resonator off: {} vs {}", w_series, w0);
            assert!((w_shunt - w0).abs() / w0 < 1e-9,
                "Shunt resonator off:  {} vs {}", w_shunt,  w0);
        }
    }
}

