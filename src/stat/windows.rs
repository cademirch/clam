use std::cmp::Ordering;
use std::collections::{HashSet, VecDeque};
use std::num::NonZeroUsize;
use std::ops::Bound;
use std::path::Path;
use std::sync::{Arc, Mutex};
use std::thread;
use std::time::Instant;

use anyhow::{bail, Result};
use bstr::ByteSlice;
use indicatif::ProgressBar;
use log::{info, trace};
use noodles::core::Region;
use noodles::vcf::variant::record::AlternateBases;
use noodles::{tabix, vcf};

use super::alleles::VCFData;
use super::build_vcf_reader;
use super::callable::*;
use super::output::{DxyRecord, FstRecord, HetRecord, PiRecord};
use crate::utils::{count_combinations, PopulationMapping};

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
    pub sample_names: Vec<String>,
    pub count_het_sites: Vec<u32>,
    pub bases_in_roh: Vec<u32>,
}

impl Population {
    pub fn new(name: String, sample_names: Vec<String>) -> Self {
        let count_het_sites: Vec<u32> = vec![0; sample_names.len()];
        let bases_in_roh = count_het_sites.clone();

        Self {
            name,
            sample_names,
            count_het_sites,
            bases_in_roh,
            ..Default::default()
        }
    }

    pub fn to_het_records(
        &self,
        chrom: &str,
        start: u32,
        end: u32,
        callable_bases: u32,
    ) -> Vec<HetRecord> {
        let mut records = Vec::with_capacity(self.sample_names.len());

        for i in 0..self.sample_names.len() {
            let record = HetRecord::new(
                chrom,
                start,
                end,
                callable_bases,
                self.bases_in_roh[i],
                &self.sample_names[i],
                self.count_het_sites[i],
            );
            records.push(record);
        }

        records
    }

    pub fn name(&self) -> &str {
        &self.name
    }
    pub fn calc_pi(&self) -> f32 {
        if self.within_comps > 0 {
            self.within_diffs as f32 / self.within_comps as f32
        } else {
            f32::NAN
        }
    }

    pub fn to_pi_record(&self, chrom: &str, begin: u32, end: u32) -> PiRecord {
        let pi = self.calc_pi();
        PiRecord::new(
            &self.name,
            chrom,
            begin,
            end,
            pi,
            self.within_comps,
            self.within_diffs,
        )
    }
}

#[derive(Debug)]
pub struct Window {
    pub region: Region,
    pub populations: Vec<Population>,
    pub dxy_stats: DxyStats,
    pub fst_stats: FstStats,
    pub sites: HashSet<u32>,
    pub ploidy: u32,
    pub callable_sites: u32,
}

#[derive(Debug, Default)]
pub struct DxyStats {
    pub comparisons: Vec<u32>,
    pub differences: Vec<u32>,
}

impl DxyStats {
    pub fn dxy(&self, pair_index: usize) -> f32 {
        if pair_index >= self.differences.len() || pair_index >= self.comparisons.len() {
            panic!("pair_index out of bounds");
        }

        let differences = self.differences[pair_index];
        let comparisons = self.comparisons[pair_index];

        if comparisons > 0 {
            differences as f32 / comparisons as f32
        } else {
            f32::NAN
        }
    }
}
#[derive(Debug, Default)]
pub struct FstStats {
    numerator: Vec<f32>,
    denominator: Vec<f32>,
}

impl FstStats {
    pub fn fst(&self, pop_pair_idx: usize) -> f32 {
        if self.denominator[pop_pair_idx] > 0.0 {
            self.numerator[pop_pair_idx] / self.denominator[pop_pair_idx]
        } else {
            f32::NAN
        }
    }
}

impl Window {
    pub fn from_regions(
        regions: Vec<Region>,
        population_info: &PopulationMapping,
        sites: Option<Vec<HashSet<u32>>>,
        ploidy: u32,
    ) -> VecDeque<Self> {
        let sites_iter = sites.unwrap_or_else(|| vec![HashSet::new(); regions.len()]);
        let dxy_size = if population_info.num_populations() >= 2 {
            count_combinations(population_info.num_populations() as u32, 2) as usize
        } else {
            0
        };

        regions
            .into_iter()
            .zip(sites_iter)
            .map(|(region, sites)| {
                let populations = population_info
                    .get_popname_refs()
                    .iter()
                    .enumerate()
                    .map(|(idx, name)| {
                        Population::new(
                            name.to_string(),
                            population_info
                                .get_sample_refs(idx)
                                .iter()
                                .map(ToString::to_string)
                                .collect(),
                        )
                    })
                    .collect();

                Self {
                    region,
                    populations,
                    dxy_stats: DxyStats {
                        comparisons: vec![0; dxy_size],
                        differences: vec![0; dxy_size],
                    },
                    fst_stats: FstStats {
                        numerator: vec![0.0; dxy_size],
                        denominator: vec![0.0; dxy_size],
                    },
                    sites,
                    ploidy,
                    callable_sites: 0,
                }
            })
            .collect()
    }

    pub fn to_fst_records(&self) -> Vec<FstRecord> {
        let mut records = Vec::new();
        let (chrom, start, end) = self.get_region_info();
        let pop_names: Vec<&str> = self.populations.iter().map(|x| x.name()).collect();

        // Iterate over population pairs
        for i in 0..pop_names.len() {
            for j in (i + 1)..pop_names.len() {
                let idx = self.get_pair_index(i, j);

                let fst = self.fst_stats.fst(idx);
                records.push(FstRecord::new(
                    pop_names[i],
                    pop_names[j],
                    chrom,
                    start,
                    end,
                    fst,
                ));
            }
        }

        records
    }

    pub fn to_dxy_records(&self) -> Vec<DxyRecord> {
        let mut records = Vec::new();
        let (chrom, start, end) = self.get_region_info();
        let pop_names: Vec<&str> = self.populations.iter().map(|x| x.name()).collect();
        for i in 0..pop_names.len() {
            for j in (i + 1)..pop_names.len() {
                let idx = self.get_pair_index(i, j);
                let dxy = self.dxy_stats.dxy(idx);

                records.push(DxyRecord::new(
                    pop_names[i],
                    pop_names[j],
                    chrom,
                    start,
                    end,
                    dxy,
                    self.dxy_stats.comparisons[idx],
                    self.dxy_stats.differences[idx],
                ));
            }
        }

        records
    }

    pub fn to_pi_records(&self) -> Vec<PiRecord> {
        let (chrom, start, end) = self.get_region_info();
        self.populations
            .iter()
            .map(|x| x.to_pi_record(chrom, start, end))
            .collect()
    }

    pub fn to_het_records(&self) -> Vec<HetRecord> {
        let (chrom, start, end) = self.get_region_info();

        self.populations
            .iter()
            .flat_map(|pop| pop.to_het_records(chrom, start, end, self.callable_sites))
            .collect()
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

    pub fn hudson_fst(&mut self, within1: f32, within2: f32, between: f32, pop_pair_idx: usize) {
        // Calculate the average within-population heterozygosity
        let average_within = (within1 + within2) / 2.0;

        // If any value is NaN or invalid, skip this position
        if average_within.is_nan() || between.is_nan() || between <= 0.0 {
            return;
        }

        // Calculate and accumulate numerator and denominator separately
        let numerator = between - average_within;
        let denominator = between;
        trace!(
            "window: {:?}, avg_within: {}, btwn: {}, num: {}, den: {}",
            self.get_region_info(),
            average_within,
            between,
            numerator,
            denominator
        );

        self.fst_stats.numerator[pop_pair_idx] += numerator;
        self.fst_stats.denominator[pop_pair_idx] += denominator;
    }

    pub fn update_dxy_diffs(&mut self, pop1: usize, pop2: usize, diffs: u32, comps: u32) {
        let idx = self.get_pair_index(pop1, pop2);
        self.dxy_stats.differences[idx] += diffs;
        self.dxy_stats.comparisons[idx] += comps;
    }

    pub fn get_population_mut(&mut self, population_idx: usize) -> Option<&mut Population> {
        self.populations.get_mut(population_idx)
    }

    pub fn update_population(&mut self, population_idx: usize, values: [u32; 2]) -> f32 {
        let population = self
            .populations
            .get_mut(population_idx)
            .unwrap_or_else(|| panic!("Invalid population index: {}", population_idx));

        let gts = values[0] + values[1];
        let diffs = values[0] * values[1];
        let comps = count_combinations(gts, 2);
        trace!("update population {}, vals: {:?}, diffs: {} comps: {}, withindiffs: {}, withincomps: {}", population_idx, values, diffs, comps, population.within_diffs, population.within_comps);

        population.refs += values[0];
        population.alts += values[1];
        population.within_diffs += diffs;
        population.within_comps += comps;

        if comps > 0 {
            diffs as f32 / comps as f32
        } else {
            f32::NAN
        }
    }

    pub fn update_allele_counts(&mut self, counts: Vec<[u32; 2]>) {
        trace!("window: {:?} counts: {:?}", self.get_region_info(), counts);
        if counts.len() >= 2 {
            // First update all populations once
            let within_values: Vec<f32> = (0..counts.len())
                .map(|pop_idx| self.update_population(pop_idx, counts[pop_idx]))
                .collect();

            // Then do FST calculations
            for pop1 in 0..counts.len() {
                let pop1_vals = counts[pop1];
                for pop2 in (pop1 + 1)..counts.len() {
                    let pop2_vals = counts[pop2];

                    let pop1_total = pop1_vals[0] + pop1_vals[1]; // total alleles in pop1
                    let pop2_total = pop2_vals[0] + pop2_vals[1]; // total alleles in pop2

                    // Number of differences when comparing between populations
                    let diffs = (pop1_vals[0] * pop2_vals[1]) + (pop1_vals[1] * pop2_vals[0]);

                    // Total number of possible comparisons between populations
                    let comps = pop1_total * pop2_total;
                    let between = if comps > 0 {
                        diffs as f32 / comps as f32
                    } else {
                        f32::NAN
                    };
                    self.update_dxy_diffs(pop1, pop2, diffs, comps);
                    let pair_idx = self.get_pair_index(pop1, pop2);
                    self.hudson_fst(within_values[pop1], within_values[pop2], between, pair_idx);
                }
            }
        } else {
            self.update_population(0, counts[0]);
        }
    }

    pub fn update_het_counts(&mut self, het_counts: Vec<Vec<u32>>) {
        self.populations
            .iter_mut()
            .zip(het_counts)
            .for_each(|(population, counts)| {
                population
                    .count_het_sites
                    .iter_mut()
                    .zip(counts)
                    .for_each(|(site_count, count)| *site_count += count);
            });
    }
    fn count_alleles(
        &mut self,
        record: &vcf::Record,
        header: &vcf::Header,
        population_info: &PopulationMapping,
        samples_in_roh: HashSet<String>,
    ) -> Result<()> {
        let mut counts = VCFData::new(population_info.get_sample_counts_per_population());

        let result = super::alleles::count_alleles(
            record,
            header,
            population_info,
            &mut counts,
            samples_in_roh,
        );
        match result {
            Ok(_) => {}
            Err(e) => log::error!("count_alleles failed: {:?}", e),
        }

        self.update_allele_counts(counts.allele_counts);
        self.update_het_counts(counts.het_counts);

        Ok(())
    }

    fn update_comparisons_d4(
        &mut self,
        within_comps: &[u32],
        dxy_comps: &[u32],
        callable_sites: u32,
    ) {
        for (pop_idx, within_comp) in within_comps.iter().enumerate() {
            if let Some(population) = self.get_population_mut(pop_idx) {
                population.within_comps += *within_comp;
            }
        }
        self.callable_sites = callable_sites;
        self.dxy_stats
            .comparisons
            .iter_mut()
            .zip(dxy_comps)
            .for_each(|(existing, new)| *existing += new);
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
    callable_path: Option<P>,
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
            let callable_path = callable_path.as_ref().map(|p| p.as_ref().to_path_buf());
            let roh_path = roh_path.as_ref().map(|p| p.as_ref().to_path_buf());
            let bar = bar.clone();
            let pop_info = population_info.clone();

            scope.spawn(move || -> Result<()> {
                trace!("Worker {} starting...", i);
                let (mut vcf_reader, header) = build_vcf_reader(vcf_path)?;

                while let Some(mut window) = {
                    let mut queue = work_queue.lock().unwrap();
                    queue.pop_front()
                } {
                    trace!("Worker {} aquired window", i);

                    let mut callable_sites = if let Some(ref path) = callable_path {
                        let extension = path.extension().and_then(|ext| ext.to_str());
                        match extension {
                            Some("d4") => Some(CallableSites::D4(D4CallableSites::from_file(
                                pop_info.get_popname_refs(),
                                &path,
                            )?)),
                            Some("gz") | Some("bed.gz") => {
                                Some(CallableSites::Bed(BedCallableSites::from_file(
                                    path,
                                    pop_info.get_sample_counts_per_population(),
                                )?))
                            }
                            _ => bail!(
                                "Unsupported callable sites file format: {}",
                                &path.display()
                            ),
                        }
                    } else {
                        None
                    };

                    let roh_tabix_query = if let Some(ref path) = roh_path {
                        trace!("Roh path: {}", &path.display());
                        let mut tabix_records = vec![];
                        let mut reader =
                            tabix::io::indexed_reader::Builder::default().build_from_path(path)?;
                        trace!("roh region: {}", &window.region);
                        let query = reader.query(&window.region)?;
                        for result in query {
                            let record = result?;
                            let fields: Vec<&str> = record.as_ref().split("\t").collect();
                            trace!("roh fields: {:?}", fields);
                            if let (Some(start), Some(end), Some(sample)) =
                                (fields.get(1), fields.get(2), fields.get(3))
                            {
                                let start: u32 = start.parse()?;
                                let end: u32 = end.parse()?;

                                trace!(
                                    "Roh: region:{}, start:{}, end:{}, sample:{}",
                                    &window.region,
                                    start,
                                    end,
                                    sample
                                );

                                tabix_records.push((sample.to_string(), start..=end));
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

                        let samples_in_roh: HashSet<String> =
                            if let Some(ref tabix_query) = roh_tabix_query {
                                trace!("Tabix results: {:?}", tabix_query);
                                tabix_query
                                    .iter()
                                    .filter_map(|(sample, interval)| {
                                        if interval.contains(&(start as u32)) {
                                            Some(sample.to_owned())
                                        // Get index and copy value if it exists
                                        } else {
                                            None
                                        }
                                    })
                                    .collect() // Collect results into a HashSet
                            } else {
                                HashSet::with_capacity(0)
                            };
                        trace!("Samples in roh: {:?}", samples_in_roh);
                        let should_process_site =
                            window.sites.is_empty() || window.sites.contains(&(start as u32));

                        if should_process_site && record.alternate_bases().len() <= 1 {
                            window.count_alleles(&record, &header, &pop_info, samples_in_roh)?;
                            sites_skipped.insert(start as u32); // skip this site in callable
                        } else {
                            sites_skipped.insert(start as u32);
                        }
                    }

                    if let Some(reader) = callable_sites.as_mut() {
                        let query_d4_time = Instant::now();
                        let (chrom, window_begin, window_end) = &window.get_region_info();

                        match reader.query(
                            chrom,
                            *window_begin,
                            *window_end,
                            window.ploidy,
                            &sites_skipped,
                        ) {
                            Ok(QueryResult(within, dxy, callable_sites)) => {
                                trace!(
                                    "Worker {} queried callable counts d4 in: {:#?}",
                                    i,
                                    query_d4_time.elapsed()
                                );
                                window.update_comparisons_d4(&within, &dxy, callable_sites);
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

                    if let Some(roh) = roh_tabix_query {
                        for (sample, range) in roh {
                            let (pop_idx, sample_idx) = pop_info.lookup_sample_name(&sample)?;
                            window.get_population_mut(pop_idx).unwrap().bases_in_roh[sample_idx] =
                                range.end() - range.start()
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
