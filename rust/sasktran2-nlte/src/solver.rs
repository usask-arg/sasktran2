//! Steady-state excited-state populations for a column of atmosphere.

use crate::linalg::solve_linear_system;
use crate::mechanism::{Coefficient, Mechanism, Species};
use crate::prelude::*;

const BOLTZMANN_J_PER_K: f64 = 1.380649e-23;

/// Inputs on a common vertical grid (any number of levels, solved independently).
pub struct Column<'a> {
    pub temperature_k: ArrayView1<'a, f64>,
    /// Used only to derive the total number density "M" when a mechanism needs
    /// it and it is not given in `densities_m3`.
    pub pressure_pa: Option<ArrayView1<'a, f64>>,
    /// Number densities [m^-3] of the mechanism's background species.
    pub densities_m3: HashMap<String, ArrayView1<'a, f64>>,
    /// The mechanism's rate inputs [s^-1], e.g. photolysis rates.
    pub rates_per_s: HashMap<String, ArrayView1<'a, f64>>,
}

/// Results on the column grid; the leading axis follows the mechanism's
/// state and process order.
pub struct Solution {
    pub state_density_m3: Array2<f64>,
    /// Event rate of every process [m^-3 s^-1]. For radiative processes this is
    /// the photon volume emission rate.
    pub process_rate_m3_s: Array2<f64>,
    pub production_m3_s: Array2<f64>,
    pub loss_m3_s: Array2<f64>,
    /// max |production - loss| / max production at each level.
    pub relative_residual: Array1<f64>,
}

/// Solves `production(n) = loss(n)` for the state populations at every level.
///
/// Mechanisms are validated to be linear in the state populations, so each
/// level is a single dense linear solve.
pub fn solve_steady_state(mechanism: &Mechanism, column: &Column) -> Result<Solution> {
    let n_levels = column.temperature_k.len();
    let n_states = mechanism.states().len();
    let n_processes = mechanism.processes().len();

    if n_states == 0 {
        return Err(anyhow!("Mechanism '{}' has no states", mechanism.name));
    }
    let background = background_profiles(mechanism, column, n_levels)?;
    let rate_inputs = rate_input_profiles(mechanism, column, n_levels)?;

    let mut solution = Solution {
        state_density_m3: Array2::zeros((n_states, n_levels)),
        process_rate_m3_s: Array2::zeros((n_processes, n_levels)),
        production_m3_s: Array2::zeros((n_states, n_levels)),
        loss_m3_s: Array2::zeros((n_states, n_levels)),
        relative_residual: Array1::zeros(n_levels),
    };

    for level in 0..n_levels {
        let temperature = column.temperature_k[level];
        if !(temperature.is_finite() && temperature > 0.0) {
            return Err(anyhow!(
                "Temperature must be positive and finite, got {temperature} at level {level}"
            ));
        }

        // Event rate per unit state population (or absolute rate when the
        // process has no state reactant).
        let coefficients: Vec<f64> = mechanism
            .processes()
            .iter()
            .map(|process| {
                let mut c = match &process.coefficient {
                    Coefficient::Rate(law) => law.evaluate(temperature),
                    Coefficient::Input(index) => rate_inputs[*index][level],
                    Coefficient::EinsteinA(a) => *a,
                };
                for reactant in &process.reactants {
                    if let Species::Background(index) = reactant {
                        c *= background[*index][level];
                    }
                }
                c
            })
            .collect();

        let mut matrix = Array2::<f64>::zeros((n_states, n_states));
        let mut rhs = Array2::<f64>::zeros((n_states, 1));
        for (process, &c) in mechanism.processes().iter().zip(&coefficients) {
            match process.state_reactant() {
                Some(source) => {
                    matrix[[source, source]] -= c;
                    for_each_state_product(process, |product, fraction| {
                        matrix[[product, source]] += fraction * c;
                    });
                }
                None => for_each_state_product(process, |product, fraction| {
                    rhs[[product, 0]] -= fraction * c;
                }),
            }
        }

        let densities = solve_linear_system(&matrix, &rhs).map_err(|err| {
            anyhow!(
                "Steady-state solve failed at level {level}: {err}. States with no loss \
                 process make the system singular: {}",
                states_without_loss(mechanism, &matrix).join(", ")
            )
        })?;

        for state in 0..n_states {
            solution.state_density_m3[[state, level]] = densities[[state, 0]];
        }
        for (p, (process, &c)) in mechanism.processes().iter().zip(&coefficients).enumerate() {
            let rate = match process.state_reactant() {
                Some(source) => c * densities[[source, 0]],
                None => c,
            };
            solution.process_rate_m3_s[[p, level]] = rate;
            if let Some(source) = process.state_reactant() {
                solution.loss_m3_s[[source, level]] += rate;
            }
            for_each_state_product(process, |product, fraction| {
                solution.production_m3_s[[product, level]] += fraction * rate;
            });
        }

        let production = solution.production_m3_s.column(level);
        let loss = solution.loss_m3_s.column(level);
        let scale = production.iter().fold(0.0_f64, |m, &v| m.max(v.abs()));
        let imbalance = production
            .iter()
            .zip(loss.iter())
            .fold(0.0_f64, |m, (p, l)| m.max((p - l).abs()));
        solution.relative_residual[level] = if scale > 0.0 { imbalance / scale } else { 0.0 };
    }

    Ok(solution)
}

fn for_each_state_product(process: &crate::mechanism::Process, mut f: impl FnMut(usize, f64)) {
    for channel in &process.channels {
        for product in &channel.products {
            if let Species::State(index) = product {
                f(*index, channel.fraction);
            }
        }
    }
}

fn background_profiles(
    mechanism: &Mechanism,
    column: &Column,
    n_levels: usize,
) -> Result<Vec<Array1<f64>>> {
    let mut missing = Vec::new();
    let mut profiles = Vec::with_capacity(mechanism.background().len());
    for name in mechanism.background() {
        let profile = match (
            column.densities_m3.get(name),
            name.as_str(),
            &column.pressure_pa,
        ) {
            (Some(profile), _, _) => profile.to_owned(),
            (None, "M", Some(pressure)) => {
                check_len("pressure", pressure.len(), n_levels)?;
                Zip::from(pressure)
                    .and(&column.temperature_k)
                    .map_collect(|p, t| p / (BOLTZMANN_J_PER_K * t))
            }
            (None, _, _) => {
                missing.push(name.clone());
                continue;
            }
        };
        check_len(&format!("{name} density"), profile.len(), n_levels)?;
        if let Some(bad) = profile.iter().find(|d| !(d.is_finite() && **d >= 0.0)) {
            return Err(anyhow!(
                "{name} density must be non-negative and finite, got {bad}"
            ));
        }
        profiles.push(profile);
    }
    if !missing.is_empty() {
        return Err(anyhow!(
            "Mechanism '{}' needs background densities that were not given: {}",
            mechanism.name,
            missing.join(", ")
        ));
    }
    Ok(profiles)
}

fn rate_input_profiles<'a>(
    mechanism: &Mechanism,
    column: &Column<'a>,
    n_levels: usize,
) -> Result<Vec<ArrayView1<'a, f64>>> {
    let missing: Vec<&str> = mechanism
        .rate_inputs()
        .iter()
        .filter(|name| !column.rates_per_s.contains_key(*name))
        .map(String::as_str)
        .collect();
    if !missing.is_empty() {
        return Err(anyhow!(
            "Mechanism '{}' needs rate inputs that were not given: {}",
            mechanism.name,
            missing.join(", ")
        ));
    }
    mechanism
        .rate_inputs()
        .iter()
        .map(|name| {
            let profile = column.rates_per_s[name];
            check_len(name, profile.len(), n_levels)?;
            if let Some(bad) = profile.iter().find(|r| !(r.is_finite() && **r >= 0.0)) {
                return Err(anyhow!(
                    "Rate input {name} must be non-negative and finite, got {bad}"
                ));
            }
            Ok(profile)
        })
        .collect()
}

fn check_len(name: &str, len: usize, n_levels: usize) -> Result<()> {
    if len != n_levels {
        return Err(anyhow!(
            "{name} has {len} levels but temperature has {n_levels}"
        ));
    }
    Ok(())
}

fn states_without_loss(mechanism: &Mechanism, matrix: &Array2<f64>) -> Vec<String> {
    mechanism
        .states()
        .iter()
        .enumerate()
        .filter(|(i, _)| matrix[[*i, *i]] == 0.0)
        .map(|(_, state)| state.id.clone())
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::array;

    const CHAIN: &str = r#"
[mechanism]
name = "chain"
version = "0"
background = ["O2", "O(3P)", "M"]

[references]
ref = "test"

[[state]]
id = "O(1D)"

[[state]]
id = "O2(b)"

[[photo_reaction]]
id = "src"
reactant = "O2"
products = ["O(1D)", "O(3P)"]
rate_input = "J"
reference = "ref"

[[reaction]]
id = "o1d_o2"
reactants = ["O(1D)", "O2"]
rate = { law = "constant", value = 1.0e-6, units = "cm3 s-1" }
channels = [{ yield = 0.25, products = ["O2(b)", "O(3P)"] }]
reference = "ref"

[[reaction]]
id = "o2b_m"
reactants = ["O2(b)", "M"]
products = ["O2", "M"]
rate = { law = "constant", value = 0.0, units = "m3 s-1" }
reference = "ref"

[[transition]]
id = "a_band"
upper = "O2(b)"
lower = "O2"
einstein_a_s = 2.0
reference = "ref"
"#;

    fn column<'a>(
        temperature: &'a Array1<f64>,
        o2: &'a Array1<f64>,
        o3p: &'a Array1<f64>,
        j: &'a Array1<f64>,
        pressure: &'a Array1<f64>,
    ) -> Column<'a> {
        Column {
            temperature_k: temperature.view(),
            pressure_pa: Some(pressure.view()),
            densities_m3: HashMap::from([
                ("O2".to_string(), o2.view()),
                ("O(3P)".to_string(), o3p.view()),
            ]),
            rates_per_s: HashMap::from([("J".to_string(), j.view())]),
        }
    }

    #[test]
    fn two_state_chain_matches_closed_form() {
        let mechanism = Mechanism::from_toml_str(CHAIN).unwrap();
        let (t, o2, o3p, j, p) = (
            array![200.0, 250.0],
            array![1.0e12, 2.0e12],
            array![0.0, 0.0],
            array![1.0, 3.0],
            array![1.0, 1.0],
        );

        let solution = solve_steady_state(&mechanism, &column(&t, &o2, &o3p, &j, &p)).unwrap();

        for level in 0..2 {
            // O(1D): production J [O2], loss k [O2] with k = 1e-12 m3 s-1.
            let o1d = j[level] * o2[level] / (1.0e-12 * o2[level]);
            // O2(b): 25% of O(1D) quenching, lost only by emission (A = 2).
            let o2b = 0.25 * 1.0e-12 * o2[level] * o1d / 2.0;
            assert!((solution.state_density_m3[[0, level]] - o1d).abs() < 1e-12 * o1d);
            assert!((solution.state_density_m3[[1, level]] - o2b).abs() < 1e-12 * o2b);

            // The radiative process rate is the photon volume emission rate.
            let a_band = mechanism
                .processes()
                .iter()
                .position(|p| p.id == "a_band")
                .unwrap();
            let ver = solution.process_rate_m3_s[[a_band, level]];
            assert!((ver - 2.0 * o2b).abs() < 1e-12 * ver);
            assert!(solution.relative_residual[level] < 1e-12);
        }
    }

    #[test]
    fn budgets_balance_at_steady_state() {
        let mechanism = Mechanism::from_toml_str(CHAIN).unwrap();
        let (t, o2, o3p, j, p) = (
            array![200.0],
            array![1.0e12],
            array![0.0],
            array![1.0],
            array![1.0],
        );
        let solution = solve_steady_state(&mechanism, &column(&t, &o2, &o3p, &j, &p)).unwrap();
        for state in 0..2 {
            let production = solution.production_m3_s[[state, 0]];
            let loss = solution.loss_m3_s[[state, 0]];
            assert!((production - loss).abs() < 1e-12 * production);
        }
    }

    #[test]
    fn total_density_is_derived_from_pressure() {
        let text = CHAIN.replace(
            "value = 0.0, units = \"m3 s-1\"",
            "value = 1.0e-25, units = \"m3 s-1\"",
        );
        let mechanism = Mechanism::from_toml_str(&text).unwrap();
        let (t, o2, o3p, j, p) = (
            array![200.0],
            array![1.0e12],
            array![0.0],
            array![1.0],
            array![10.0],
        );
        let solution = solve_steady_state(&mechanism, &column(&t, &o2, &o3p, &j, &p)).unwrap();

        let m = 10.0 / (BOLTZMANN_J_PER_K * 200.0);
        let o1d = 1.0 / 1.0e-12;
        let o2b = 0.25 * 1.0e-12 * 1.0e12 * o1d / (2.0 + 1.0e-25 * m);
        assert!((solution.state_density_m3[[1, 0]] - o2b).abs() < 1e-12 * o2b);
    }

    #[test]
    fn reaction_between_background_species_is_a_source_in_si_units() {
        let text = r#"
[mechanism]
name = "source"
version = "0"
background = ["O3", "O(3P)", "O2"]

[references]
ref = "test"

[[state]]
id = "O2(X, v=1)"

[[reaction]]
id = "o3_o"
reactants = ["O3", "O(3P)"]
products = ["O2(X, v=1)", "O2"]
rate = { law = "constant", value = 1.0e-12, units = "cm3 s-1" }
reference = "ref"

[[transition]]
id = "decay"
upper = "O2(X, v=1)"
lower = "O2"
einstein_a_s = 1.0
reference = "ref"
"#;
        let mechanism = Mechanism::from_toml_str(text).unwrap();
        let (t, o3, o3p, o2) = (array![200.0], array![2.0e12], array![3.0e12], array![0.0]);
        let column = Column {
            temperature_k: t.view(),
            pressure_pa: None,
            densities_m3: HashMap::from([
                ("O3".to_string(), o3.view()),
                ("O(3P)".to_string(), o3p.view()),
                ("O2".to_string(), o2.view()),
            ]),
            rates_per_s: HashMap::new(),
        };

        let solution = solve_steady_state(&mechanism, &column).unwrap();

        // 1e-18 m3 s-1 * 2e12 m-3 * 3e12 m-3, lost at 1 s-1.
        assert!((solution.state_density_m3[[0, 0]] - 6.0e6).abs() < 1e-12 * 6.0e6);
    }

    #[test]
    fn missing_inputs_are_listed() {
        let mechanism = Mechanism::from_toml_str(CHAIN).unwrap();
        let t = array![200.0];
        let column = Column {
            temperature_k: t.view(),
            pressure_pa: None,
            densities_m3: HashMap::new(),
            rates_per_s: HashMap::new(),
        };
        let err = solve_steady_state(&mechanism, &column)
            .err()
            .unwrap()
            .to_string();
        assert!(err.contains("O2, O(3P), M"), "{err}");
    }

    #[test]
    fn state_without_loss_is_reported() {
        let text = CHAIN.replace("einstein_a_s = 2.0", "einstein_a_s = 0.0");
        let mechanism = Mechanism::from_toml_str(&text).unwrap();
        let (t, o2, o3p, j, p) = (
            array![200.0],
            array![1.0e12],
            array![0.0],
            array![1.0],
            array![0.0],
        );
        let err = solve_steady_state(&mechanism, &column(&t, &o2, &o3p, &j, &p))
            .err()
            .unwrap()
            .to_string();
        assert!(err.contains("O2(b)"), "{err}");
    }

    #[test]
    fn negative_rate_input_is_an_error() {
        let mechanism = Mechanism::from_toml_str(CHAIN).unwrap();
        let (t, o2, o3p, j, p) = (
            array![200.0],
            array![1.0e12],
            array![0.0],
            array![-1.0],
            array![1.0],
        );
        assert!(solve_steady_state(&mechanism, &column(&t, &o2, &o3p, &j, &p)).is_err());
    }
}
