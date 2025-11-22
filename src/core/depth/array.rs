use crate::core::population::PopulationMap;
use color_eyre::{eyre, eyre::bail, Result};
use ndarray::{Array1, Array2, Axis};

/// Depth values across multiple samples at contiguous positions
pub struct MultisampleDepthArray {
    /// Depth values with shape (positions, samples)
    pub depths: Array2<u32>,
    /// Sample names corresponding to columns in depths array
    pub sample_names: Vec<String>,
    /// Population membership matrix with shape (samples, populations)
    /// Element [i, j] is true if sample i belongs to population j
    pub pop_membership: Array2<bool>,
}

impl MultisampleDepthArray {
    /// Build from per-sample depth vectors
    ///
    /// # Arguments
    /// * `depths` - Vec of depth vectors, one per sample (samples-major layout)
    /// * `sample_names` - Names corresponding to each sample
    /// * `pop_map` - Population map to assign samples to populations
    pub fn from_samples(
        depths: Vec<Vec<u32>>,
        sample_names: Vec<String>,
        pop_map: &PopulationMap,
    ) -> Result<Self> {
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
        let num_samples = sample_names.len();
        let num_pops = pop_map.num_populations();

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

        // Build population membership matrix
        let mut pop_membership = Array2::<bool>::from_elem((num_samples, num_pops), false);
        for (sample_idx, name) in sample_names.iter().enumerate() {
            let (pop_idx, _) = pop_map
                .lookup(name)
                .ok_or_else(|| eyre::eyre!("Sample '{}' not found in population map", name))?;
            pop_membership[[sample_idx, pop_idx]] = true;
        }

        // Convert to Array1s and stack along axis 1 to get (positions, samples)
        let arrays: Vec<Array1<u32>> = depths.into_iter().map(Array1::from_vec).collect();

        let views: Vec<_> = arrays.iter().map(|a| a.view()).collect();
        let depths = ndarray::stack(Axis(1), &views)?;

        Ok(Self {
            depths,
            sample_names,
            pop_membership,
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
    ///
    /// # Arguments
    /// * `callable_mask` - Per-sample callable mask from `apply_sample_thresholds`
    /// * `min_proportion` - Minimum proportion of samples that must be callable
    /// * `mean_depth_range` - (min, max) mean depth range for site to be callable
    ///
    /// Returns Array2<u8> with shape (positions, populations) where each value
    /// is the number of callable samples in that population at that position
    pub fn count_callable_per_population(
        &self,
        callable_mask: &Array2<bool>,
        min_proportion: f64,
        mean_depth_range: (f64, f64),
    ) -> Result<Array2<u8>> {
        let min_samples = (self.num_samples() as f64 * min_proportion).ceil() as usize;
        let site_passes = self.apply_site_filters(callable_mask, min_samples, mean_depth_range);
        let callable_u8 = callable_mask.mapv(|x| x as u8);

        // Matrix multiply: (positions, samples) × (samples, populations) = (positions, populations)
        let pop_membership_u8 = self.pop_membership.mapv(|x| x as u8);
        let counts = callable_u8.dot(&pop_membership_u8);

        // Reshape to (positions, 1) and convert to u8 for multiplication
        let site_mask = site_passes.insert_axis(Axis(1)).mapv(|x| x as u8);

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
