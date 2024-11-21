use super::callable::*;
use anyhow::{anyhow, bail, Context, Result};
use crossbeam::channel::{Receiver, Sender};
use fnv::FnvHashMap;
use noodles::tabix;
use noodles::{
    bgzf,
    core::{region, Region},
    vcf::{
        self,
        variant::record::{
            samples::{keys::key, series::Value, Series},
            AlternateBases,
        },
    },
};
use rayon::{iter::ParallelBridge, prelude::*};
use std::collections::HashSet;
use std::num::NonZeroUsize;
use std::{fs::File, io, io::BufReader, path::Path, thread};
pub type ProcessedRecord = (String, usize, usize, Option<Vec<[u32; 2]>>); //chrom idx, window_idx, position, allele counts

#[derive(Debug)]
pub struct PopulationData {
    pub refs: u32,
    pub alts: u32,
    pub within_diffs: u32,
    pub within_comps: u32,
}

#[derive(Debug)]
pub struct WindowedData {
    pub chrom: String,
    pub begin: u32,
    pub end: u32,
    pub populations: Vec<PopulationData>,
    pub dxy_comps: Option<Vec<u32>>,
    pub dxy_diffs: Option<Vec<u32>>,
    pub sites_skipped: HashSet<u32>,
    pub sites: Option<HashSet<u32>>, // for specifiying specific sites in this window we are only interersted in
    pub ploidy: u32,
}

impl WindowedData {
    pub fn new(
        chrom: &str,
        begin: u32,
        end: u32,
        num_pops: usize,
        sites: Option<HashSet<u32>>,
        ploidy: u32,
    ) -> Self {
        let chrom = chrom.to_string();
        let mut populations = Vec::with_capacity(num_pops);
        let mut sites_skipped: HashSet<u32> = HashSet::new();
        for _ in 0..num_pops {
            populations.push(PopulationData {
                refs: 0,
                alts: 0,
                within_diffs: 0,
                within_comps: 0,
            });
        }

        let num_pop_combs = if num_pops > 1 {
            num_pops * (num_pops - 1) / 2
        } else {
            0
        };

        let dxy = if num_pop_combs > 0 {
            Some(vec![0; num_pop_combs]) // Initialize with `num_pop_combs` zeros
        } else {
            None
        };
        let dxy_comps = dxy.clone();
        let dxy_diffs = dxy.clone();
        WindowedData {
            chrom,
            begin,
            end,
            populations,
            dxy_comps,
            dxy_diffs,
            sites_skipped,
            sites,
            ploidy,
        }
    }
    pub fn get_pair_index(&self, pop1: usize, pop2: usize) -> usize {
        if pop1 < pop2 {
            pop1 * (self.populations.len() - 1) - (pop1 * (pop1 + 1) / 2) + pop2 - 1
        } else {
            pop2 * (self.populations.len() - 1) - (pop2 * (pop2 + 1) / 2) + pop1 - 1
        }
    }
    pub fn update_dxy_diffs(&mut self, pop1: usize, pop2: usize, diffs: u32, comps: u32) {
        let pair_index = self.get_pair_index(pop1, pop2);

        if let Some(dxy_diffs) = &mut self.dxy_diffs {
            dxy_diffs[pair_index] += diffs;
        }

        if let Some(dxy_comps) = &mut self.dxy_comps {
            dxy_comps[pair_index] += comps;
        }
    }
    pub fn get_population_mut(&mut self, population_idx: usize) -> Option<&mut PopulationData> {
        if let Some(population) = self.populations.get_mut(population_idx) {
            Some(population)
        } else {
            None
        }
    }
    pub fn update_population(&mut self, population_idx: usize, values: [u32; 2]) {
        if let Some(population) = self.get_population_mut(population_idx) {
            log::trace!("Update pop values: {:?}", values);
            let gts = values[0] + values[1];
            let diffs = values[0] * values[1];
            let comps = (gts * (gts - 1)) / 2;
            population.refs += values[0];
            population.alts += values[1];
            population.within_diffs += diffs;
            population.within_comps += comps;
        } else {
            panic!("Invalid population index: {}", population_idx);
        }
    }
    pub fn query_d4(&mut self, d4_reader: &mut D4CallableSites) -> Result<()> {
        if let Some(ref sites) = self.sites {
            let mut updates = Vec::new();

            for &site in sites {
                let (within_comps, dxy_comps) = d4_reader.query(
                    &self.chrom,
                    site,
                    site + 1, // `query` expects an exclusive end, so use `site + 1`
                    self.ploidy,
                    &self.sites_skipped.clone(),
                )?;
                updates.push((within_comps, dxy_comps));
            }

            for (within_comps, dxy_comps) in updates {
                self.update_comparisons(&within_comps, &dxy_comps);
            }
        } else {
            // Adjust `begin` from 1-based to 0-based
            let adjusted_begin = self.begin - 1;
            let (within_comps, dxy_comps) = d4_reader.query(
                &self.chrom,
                adjusted_begin,
                self.end,
                self.ploidy,
                &self.sites_skipped.clone(),
            )?;
            self.update_comparisons(&within_comps, &dxy_comps);
        }
        Ok(())
    }

    fn update_comparisons(&mut self, within_comps: &[u32], dxy_comps: &Option<Vec<u32>>) {
        for (pop_idx, within_comp) in within_comps.iter().enumerate() {
            if let Some(population) = self.get_population_mut(pop_idx) {
                population.within_comps = *within_comp;
            }
        }

        if let Some(dxy_comps) = dxy_comps {
            if let Some(ref mut self_dxy_comps) = self.dxy_comps {
                for (comp_idx, dxy_comp) in dxy_comps.iter().enumerate() {
                    self_dxy_comps[comp_idx] = *dxy_comp;
                }
            }
        }
    }
    pub fn calc_pi(&mut self) -> Result<Vec<f32>> {
        let mut res = vec![0.0; self.populations.len()];
        for (pop_idx, pop) in self.populations.iter().enumerate() {
            let pi = pop.within_diffs as f32 / pop.within_comps as f32;
            res[pop_idx] = pi;
        }
        Ok(res)
    }
    pub fn calc_dxy(&self) -> Result<Vec<f32>> {
        if let (Some(dxy_diffs), Some(dxy_comps)) = (&self.dxy_diffs, &self.dxy_comps) {
            let mut res = vec![0.0; dxy_diffs.len()];
            for (idx, (diffs, comps)) in dxy_diffs.iter().zip(dxy_comps.iter()).enumerate() {
                res[idx] = *diffs as f32 / *comps as f32
            }
            Ok(res)
        } else {
            bail!("No dxy comps or diffs!")
        }
    }
}
