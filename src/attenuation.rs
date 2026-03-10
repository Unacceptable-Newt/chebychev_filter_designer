use crate::complex::*;
use crate::cheb_calc::*;

// ─────────────────────────────────────────────────────────────────────────────
// ABCD-matrix attenuation
// ─────────────────────────────────────────────────────────────────────────────

/// Compute attenuation (dB) of the bandpass filter at a given frequency.
///
/// Each g_k corresponds to ONE resonator in the ladder (alternating):
///   odd k  ->  series LC tank  (Z = jwL + 1/jwC)
///   even k ->  shunt  LC tank  (Y = jwC + 1/jwL)
///
/// The stored L/C values are normalised to a 1-ohm system; they are
/// scaled by Z0 here before building the ABCD matrices.
///
/// # Arguments
/// * `elements`  - Slice of `BandpassElement` from `chebyshev_bandpass_elements`.
/// * `freq_hz`   - Query frequency in Hz.
/// * `z0`        - Reference (source/load) impedance in Ohms (typically 50.0).
///
/// # Returns
/// Attenuation in dB (>= 0).  0 dB = perfect transmission.
pub fn attenuation_db(elements: &[BandpassElement], freq_hz: f64, z0: f64) -> f64 {
    let w   = 2.0 * std::f64::consts::PI * freq_hz;
    let z0c = C::re(z0);
    let mut total = identity();

    for e in elements {
        // Z0-denormalise: inductances scale by Z0, capacitances scale by 1/Z0
        let m: Mat2 = if e.index % 2 == 1 {
            // Odd stage -> series LC tank
            let l     = z0 * e.l_series;
            let c_val = e.c_series / z0;
            let z     = C::im(w * l) + C::re(1.0) / C::im(w * c_val);
            [[C::re(1.0), z      ],
             [C::re(0.0), C::re(1.0)]]
        } else {
            // Even stage -> shunt LC tank
            let l     = z0 * e.l_shunt;
            let c_val = e.c_shunt / z0;
            let y     = C::im(w * c_val) + C::re(1.0) / C::im(w * l);
            [[C::re(1.0), C::re(0.0)],
             [y,          C::re(1.0)]]
        };
        total = mat_mul(total, m);
    }

    // ABCD -> S21 = 2 / (A + B/Z0 + C*Z0 + D)
    let [a, b] = total[0];
    let [c, d] = total[1];
    let denom  = a + b / z0c + c * z0c + d;
    let s21    = C::re(2.0) / denom;
    let mag    = s21.abs();
    if mag <= 0.0 { f64::INFINITY } else { -20.0 * mag.log10() }
}

/// Sweep attenuation across a frequency range.
///
/// # Arguments
/// * `elements`  – Filter element values.
/// * `f_start`   – Start frequency in Hz.
/// * `f_stop`    – Stop frequency in Hz.
/// * `n_points`  – Number of evenly-spaced frequency points.
/// * `z0`        – Reference impedance in Ω.
///
/// # Returns
/// Vec of `(frequency_hz, attenuation_db)` pairs.
pub fn sweep_attenuation(
    elements: &[BandpassElement],
    f_start: f64,
    f_stop: f64,
    n_points: usize,
    z0: f64,
) -> Vec<(f64, f64)> {
    assert!(n_points >= 2);
    let step = (f_stop - f_start) / (n_points - 1) as f64;
    (0..n_points)
        .map(|i| {
            let f = f_start + i as f64 * step;
            (f, attenuation_db(elements, f, z0))
        })
        .collect()
}

// ─────────────────────────────────────────────────────────────────────────────
// Additional tests
// ─────────────────────────────────────────────────────────────────────────────
#[cfg(test)]
mod attenuation_tests {
    use super::*;

    fn make_filter() -> (Vec<f64>, Vec<BandpassElement>, f64, f64) {
        let (g, elems) = chebyshev_bandpass_elements(5, 0.5, 900e6, 1100e6);
        (g, elems, 900e6, 1100e6)
    }

    /// At centre frequency (1 GHz) attenuation should be near 0 dB (≤ ripple).
    #[test]
    fn test_centre_frequency_low_loss() {
        let (_, elems, f_low, f_high) = make_filter();
        let f_centre = (f_low * f_high).sqrt();
        let att = attenuation_db(&elems, f_centre, 50.0);
        assert!(att < 1.0, "Centre-freq loss should be < 1 dB, got {att:.3} dB");
    }

    /// Well outside the band attenuation should be high.
    #[test]
    fn test_stopband_high_attenuation() {
        let (_, elems, _, _) = make_filter();
        let att = attenuation_db(&elems, 500e6, 50.0);
        assert!(att > 20.0, "500 MHz should be well attenuated, got {att:.1} dB");
    }

    /// Sweep returns the right number of points.
    #[test]
    fn test_sweep_length() {
        let (_, elems, _, _) = make_filter();
        let pts = sweep_attenuation(&elems, 500e6, 1500e6, 101, 50.0);
        assert_eq!(pts.len(), 101);
    }
}

