use super::callable::D4CallableSites;
use crate::stat::callable;
use crate::utils::{count_combinations, PopulationMapping};
use anyhow::{bail, Context, Result};
use bstr::ByteSlice;
use fnv::FnvHashMap;
use indicatif::ProgressBar;
use log::{info, trace};
use noodles::core::Region;
use noodles::tabix;
use noodles::vcf;
use noodles::vcf::variant::record::AlternateBases;
use std::cmp::Ordering;
use std::collections::{HashSet, VecDeque};
use std::num::NonZeroUsize;
use std::ops::Bound;
use std::path::Path;
use std::sync::{Arc, Mutex};
use std::thread;
use std::time::Instant;

pub trait RegionExt {
    fn start_as_u32(&self) -> u32;
    fn end_as_u32(&self) -> u32;
    fn name_as_str(&self) -> &str;
}

impl RegionExt for Region {
    fn start_as_u32(&self) -> u32 {
        match self.start() {
            Bound::Included(pos) => pos.get().try_into().unwrap(),
            Bound::Excluded(_) => panic!("Region starts should always be Included."),
            Bound::Unbounded => panic!("Region starts should always be bounded."),
        }
    }

    fn end_as_u32(&self) -> u32 {
        match self.end() {
            Bound::Included(pos) => pos.get().try_into().unwrap(),
            Bound::Excluded(_) => panic!("Region ends should always be Included."),
            Bound::Unbounded => panic!("Region ends should always be bounded."),
        }
    }

    fn name_as_str(&self) -> &str {
        self.name()
            .to_str()
            .expect("Region name is not valid UTF-8")
    }
}

#[derive(Debug, Default)]
pub struct Population {
    pub name: String,
    pub refs: u32,
    pub alts: u32,
    pub within_diffs: u32,
    pub within_comps: u32,
}

impl Population {
    pub fn new(name: String) -> Self {
        Self {
            name,
            ..Default::default()
        }
    }
}

#[derive(Debug)]
pub struct Window {
    pub region: Region,
    pub populations: Vec<Population>,
    pub dxy_comps: Option<Vec<u32>>,
    pub dxy_diffs: Option<Vec<u32>>,
    pub sites: Option<HashSet<u32>>, // for specifiying specific sites in this window we are only interersted in
    pub ploidy: u32,
}

impl Window {
    pub fn from_regions(
        regions: Vec<Region>,
        population_info: &PopulationMapping,
        sites: Option<Vec<HashSet<u32>>>,
        ploidy: u32,
    ) -> VecDeque<Self> {
        let sites_iter = match sites {
            Some(sites_vec) => sites_vec.into_iter().map(Some).collect::<Vec<_>>(),
            None => vec![None; regions.len()],
        };

        regions
            .into_iter()
            .zip(sites_iter.into_iter()) // Pair regions with corresponding sites (or None)
            .map(|(region, site)| {
                let populations: Vec<Population> = population_info
                    .get_popname_refs()
                    .unwrap_or(vec!["pop1"])
                    .iter()
                    .map(|name| Population::new(name.to_string()))
                    .collect();
                let population_combinations =
                    count_combinations(population_info.num_populations as u32, 2);

                let (dxy_comps, dxy_diffs) = if population_combinations > 0 {
                    (
                        Some(Vec::<u32>::with_capacity(population_combinations as usize)),
                        Some(Vec::<u32>::with_capacity(population_combinations as usize)),
                    )
                } else {
                    (None, None)
                };

                Self {
                    region,
                    populations,
                    dxy_comps,
                    dxy_diffs,
                    sites: site, // This can be Some(HashSet) or None
                    ploidy,
                }
            })
            .collect() // Collect the iterator into a VecDeque<Window>
    }

    pub fn get_region_info(&self) -> (&str, u32, u32) {
        let chrom = self.region.name_as_str();
        let start = self.region.start_as_u32();
        let end = self.region.end_as_u32();

        (chrom, start, end)
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
    pub fn get_population_mut(&mut self, population_idx: usize) -> Option<&mut Population> {
        if let Some(population) = self.populations.get_mut(population_idx) {
            Some(population)
        } else {
            None
        }
    }
    pub fn update_population(&mut self, population_idx: usize, values: [u32; 2]) {
        if let Some(population) = self.get_population_mut(population_idx) {
            let gts = values[0] + values[1];
            let diffs = values[0] * values[1];
            let comps = count_combinations(gts, 2);
            population.refs += values[0];
            population.alts += values[1];
            population.within_diffs += diffs;
            population.within_comps += comps;
        } else {
            panic!("Invalid population index: {}", population_idx);
        }
    }

    pub fn update_counts(&mut self, counts: Vec<[u32; 2]>) {
        if counts.len() >= 2 {
            for pop1 in 0..counts.len() {
                let pop1_vals = counts[pop1];
                self.update_population(pop1, pop1_vals);
                for pop2 in (pop1 + 1)..counts.len() {
                    let pop2_vals = counts[pop2];
                    self.update_population(pop2, pop2_vals);
                    let diffs = (pop1_vals[0] * pop2_vals[1]) + (pop1_vals[1] * pop2_vals[0]);
                    let comps = (pop1_vals[0] + pop2_vals[1]) * (pop1_vals[1] + pop2_vals[0]);
                    self.update_dxy_diffs(pop1, pop2, diffs, comps);
                }
            }
        } else {
            self.update_population(0, counts[0]); // only 1 population
        }
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
}

impl Ord for Window {
    fn cmp(&self, other: &Self) -> Ordering {
        // Compare by region name first
        match self.region.name().cmp(other.region.name()) {
            Ordering::Equal => {
                // Compare by region start
                let self_start = self.region.start_as_u32();
                let other_start = other.region.start_as_u32();

                self_start.cmp(&other_start)
            }
            other => other,
        }
    }
}

impl PartialOrd for Window {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl PartialEq for Window {
    fn eq(&self, other: &Self) -> bool {
        self.region.name() == other.region.name() && self.region.start() == other.region.start()
    }
}

impl Eq for Window {}

pub fn process_windows<P: AsRef<Path>>(
    vcf_path: P,
    d4_path: Option<P>,
    roh_path: Option<P>,
    worker_count: NonZeroUsize,
    windows: VecDeque<Window>,
    progress_bar: Option<ProgressBar>,
    population_info: PopulationMapping,
) -> Result<VecDeque<Window>> {
    let num_windows = windows.len();
    let res = Arc::new(Mutex::new(VecDeque::<Window>::with_capacity(num_windows)));
    let work_queue = Arc::new(Mutex::new(windows));

    let bar = if let Some(bar) = progress_bar {
        bar.set_length(num_windows as u64);
        Some(bar)
    } else {
        None
    };
    let start_time = Instant::now();
    info!("Processing {} regions...", num_windows);
    thread::scope(|scope| {
        for i in 0..worker_count.get() {
            trace!("Spawned worker {}", i);
            let work_queue = Arc::clone(&work_queue);
            let res = Arc::clone(&res);
            let vcf_path = vcf_path.as_ref();
            let d4_path = d4_path.as_ref().map(|p| p.as_ref().to_path_buf());
            let roh_path = roh_path.as_ref().map(|p| p.as_ref().to_path_buf());
            let bar = bar.clone();
            let pop_info = population_info.clone();

            scope.spawn(move || -> Result<()> {
                trace!("Worker {} starting...", i);

                let mut vcf_reader =
                    vcf::io::indexed_reader::Builder::default().build_from_path(vcf_path)?;

                let header = vcf_reader.read_header()?;

                // actual work loop
                while let Some(mut window) = {
                    let mut queue = work_queue.lock().unwrap();
                    queue.pop_front()
                } {
                    trace!("Worker {} aquired window", i);

                    let sample_map = &pop_info.samp_idx_to_pop_idx;
                    let population_names = pop_info.get_popname_refs();
                    let sample_name_to_samp_idx = &pop_info.sample_name_to_sample_idx;

                    let mut d4_reader = if let Some(ref path) = d4_path {
                        Some(D4CallableSites::from_file(&population_names, path)?)
                    } else {
                        None
                    };

                    let roh_tabix_query = if let Some(ref path) = roh_path {
                        let mut tabix_records = vec![];
                        let mut reader =
                            tabix::io::indexed_reader::Builder::default().build_from_path(path)?;
                        let query = reader.query(&window.region)?;
                        for result in query {
                            let record = result?;
                            let fields: Vec<&str> = record.as_ref().split("\t").collect();
                            if let (Some(start), Some(end), Some(sample)) =
                                (fields.get(1), fields.get(2), fields.get(3))
                            {
                                let mut start: u32 = start.parse()?;
                                let end: u32 = end.parse()?;
                                start += 1;
                                tabix_records
                                    .push((sample.to_string(), std::ops::Range { start, end }));
                            }
                        }
                        Some(tabix_records)
                    } else {
                        None
                    };

                    let query = vcf_reader.query(&header, &window.region)?;
                    let mut sites_skipped: HashSet<u32> = HashSet::new();

                    for result in query {
                        let record = result?;
                        let start = record.variant_start().unwrap().unwrap().get();

                        let samples_in_roh: HashSet<usize> =
                            if let Some(ref tabix_query) = roh_tabix_query {
                                tabix_query
                                    .iter()
                                    .filter_map(|(sample, interval)| {
                                        if interval.contains(&(start as u32)) {
                                            sample_name_to_samp_idx.get(sample).copied()
                                        // Get index and copy value if it exists
                                        } else {
                                            None
                                        }
                                    })
                                    .collect() // Collect results into a HashSet
                            } else {
                                HashSet::with_capacity(0)
                            };

                        let should_process_site = match &window.sites {
                            Some(sites) => sites.contains(&(start as u32)),
                            None => true,
                        };

                        if should_process_site && record.alternate_bases().len() <= 1 {
                            let num_pops = pop_info.num_populations;
                            let mut counts_vec: Vec<[u32; 2]> = vec![[0; 2]; num_pops];
                            super::alleles::count_alleles(
                                &record,
                                &header,
                                &sample_map,
                                &mut counts_vec,
                                samples_in_roh,
                            )?;
                            trace!("Counts: {:?}", counts_vec);
                            window.update_counts(counts_vec);
                        } else {
                            // Skip the site
                            sites_skipped.insert(start as u32);
                        }
                    }
                    if let Some(reader) = d4_reader.as_mut() {
                        let query_d4_time = Instant::now();
                        let (chrom, window_begin, window_end) = &window.get_region_info();

                        match reader.query(
                            chrom,
                            *window_begin-1,
                            *window_end,
                            window.ploidy,
                            &sites_skipped,
                        ) {
                            Ok((within, dxy)) => {
                                trace!(
                                    "Worker {} queried callable counts d4 in: {:#?}",
                                    i,
                                    query_d4_time.elapsed()
                                );
                                window.update_comparisons(&within, &dxy);
                            }
                            Err(e) => {
                                bail!(
                                    "Worker {} encountered an error querying callable counts: {}",
                                    i,
                                    e
                                );
                            }
                        }
                    }
                    {
                        let mut result_queue = res.lock().unwrap();
                        result_queue.push_back(window);
                    }

                    if let Some(ref bar) = bar {
                        bar.inc(1);
                    }
                }
                Ok(())
            });
        }
    });
    info!("Finished processing in {:#?}", start_time.elapsed());
    let sort_time = Instant::now();

    {
        let mut result_queue = res.lock().unwrap();
        result_queue.make_contiguous().sort();
    }
    info!("Sorted results in {:#?}", sort_time.elapsed());
    let result = Ok(std::mem::take(&mut *res.lock().unwrap()));
    result
}
