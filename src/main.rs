mod cheb_calc;
mod complex;
mod attenuation;


use cheb_calc::*;
use attenuation::*;
// ─────────────────────────────────────────────────────────────────────────────
// Demo main
// ─────────────────────────────────────────────────────────────────────────────
fn main() {
    let n         = 5;
    let ripple_db = 0.5;
    let f_low     = 900e6_f64;   // 900 MHz
    let f_high    = 1100e6_f64;  // 1.1 GHz

    println!("═══════════════════════════════════════════════════════");
    println!(" Chebyshev Bandpass Filter  |  N={n}  |  {ripple_db} dB ripple");
    println!(" Passband: {:.0} MHz – {:.0} MHz",
             f_low / 1e6, f_high / 1e6);
    println!("═══════════════════════════════════════════════════════");

    let (g, elements) = chebyshev_bandpass_elements(n, ripple_db, f_low, f_high);

    println!("\nPrototype g-values:");
    for (i, gv) in g.iter().enumerate() {
        println!("  g[{i}] = {gv:.6}");
    }

    println!("\nBandpass LC elements (series & shunt per stage):");
    println!("{:<6} {:>14} {:>14} {:>14} {:>14}",
             "Stage", "L_series (nH)", "C_series (pF)", "L_shunt (nH)", "C_shunt (pF)");
    println!("{}", "─".repeat(66));
    for e in &elements {
        println!("{:<6} {:>14.4} {:>14.4} {:>14.4} {:>14.4}",
                 e.index,
                 e.l_series * 1e9,
                 e.c_series * 1e12,
                 e.l_shunt  * 1e9,
                 e.c_shunt  * 1e12);
    }
    println!();

    // ── Attenuation sweep ─────────────────────────────────────────────────
    println!("Attenuation sweep (50 Ω system):");
    println!("{:>12}  {:>12}", "Freq (MHz)", "Atten (dB)");
    println!("{}", "─".repeat(27));
    let sweep = sweep_attenuation(&elements, 500e6, 1500e6, 21, 50.0);
    for (f, att) in &sweep {
        let marker = if *att < 1.0 { " ← passband" } else { "" };
        println!("{:>12.1}  {:>12.3}{}", f / 1e6, att, marker);
    }
    println!();
}
