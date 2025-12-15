use crate::core::population::PopulationMap;
use color_eyre::{eyre::bail, eyre::eyre, Result};
use ndarray::{Array1, Array2, Axis};

/// Helper to stack depth vectors
pub fn stack_depths(depths: Vec<Vec<u32>>, sample_names: Vec<String>) -> Result<Array2<u32>> {
    if depths.is_empty() {
        bail!("No sample depths provided");
    }

    if depths.len() != sample_names.len() {
        bail!(
            "Mismatch: {} depth vectors but {} sample names",
            depths.len(),
            sample_names.len()
        );
    }

    let num_positions = depths[0].len();
    let num_samples = depths.len();

    // Validate all samples have same number of positions
    for (i, sample_depths) in depths.iter().enumerate() {
        if sample_depths.len() != num_positions {
            bail!(
                "Sample '{}' has {} positions, expected {}",
                sample_names[i],
                sample_depths.len(),
                num_positions
            );
        }
    }
    
    // Directly build array in correct layout (positions, samples)
    let mut array = Array2::<u32>::zeros((num_positions, num_samples));
    for (sample_idx, sample_depths) in depths.into_iter().enumerate() {
        // Assign entire column at once
        array.column_mut(sample_idx).assign(&Array1::from(sample_depths));
    }
    
    Ok(array)
}

pub fn build_pop_membership(sample_names: &[String], pop_map: &PopulationMap) -> Array2<bool> {
    let num_samples = sample_names.len();
    let num_pops = pop_map.num_populations();

    let mut pop_membership = Array2::<bool>::from_elem((num_samples, num_pops), false);
    for (sample_idx, name) in sample_names.iter().enumerate() {
        if let Some((pop_idx, _)) = pop_map.lookup(name) {
            pop_membership[[sample_idx, pop_idx]] = true;
        }
    }
    pop_membership
}

/// Depth values across multiple samples at contiguous positions
pub struct MultisampleDepthArray {
    /// Depth values with shape (positions, samples)
    pub depths: Array2<u32>,
    /// Sample names corresponding to columns in depths array
    pub sample_names: Vec<String>,
}

impl MultisampleDepthArray {
    /// Build from per-sample depth vectors

    pub fn from_samples(depths: Vec<Vec<u32>>, sample_names: Vec<String>) -> Result<Self> {
        let depths = stack_depths(depths, sample_names.clone())?;

        Ok(Self {
            depths,
            sample_names,
        })
    }
    /// Apply per-sample depth thresholds
    ///
    /// Returns a boolean mask where true means the sample passes thresholds at that position
    pub fn apply_sample_thresholds(&self, min_depth: f64, max_depth: f64) -> Array2<bool> {
        self.depths.mapv(|d| {
            let d = d as f64;
            d >= min_depth && d <= max_depth
        })
    }

    /// Count callable samples per population using site-level filters

    pub fn count_callable_per_population(
        &self,
        callable_mask: &Array2<bool>,
        pop_membership: &Array2<bool>,
        min_proportion: f64,
        mean_depth_range: (f64, f64),
    ) -> Result<Array2<u16>> {
        let min_samples = (self.num_samples() as f64 * min_proportion).ceil() as usize;
        let site_passes = self.apply_site_filters(callable_mask, min_samples, mean_depth_range);
        let callable_u16 = callable_mask.mapv(|x| x as u16);

        let pop_membership_u16 = pop_membership.mapv(|x| x as u16);
        let counts = callable_u16.dot(&pop_membership_u16);

        let site_mask = site_passes.insert_axis(Axis(1)).mapv(|x| x as u16);
        let result = &counts * &site_mask;

        Ok(result)
    }

    /// Apply site-level filters to determine which positions pass thresholds
    fn apply_site_filters(
        &self,
        callable_mask: &Array2<bool>,
        min_samples: usize,
        mean_depth_range: (f64, f64),
    ) -> Array1<bool> {
        let mean_depths = self.depths.mapv(|x| x as f64).mean_axis(Axis(1)).unwrap();

        let mean_depth_ok =
            mean_depths.mapv(|mean| mean >= mean_depth_range.0 && mean <= mean_depth_range.1);

        let callable_counts = callable_mask.mapv(|x| x as usize).sum_axis(Axis(1));
        let proportion_ok = callable_counts.mapv(|count| count >= min_samples);

        &mean_depth_ok & &proportion_ok
    }
    /// Number of positions
    pub fn num_positions(&self) -> usize {
        self.depths.nrows()
    }

    /// Number of samples
    pub fn num_samples(&self) -> usize {
        self.depths.ncols()
    }
}
