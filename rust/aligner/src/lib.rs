//! Implementations of the CutAdapt aligners (Martin et al. 2012).
//! 
//! An `Aligner` operates on a pair of (possibly differently sized) byte arrays (`[u8]`) and 
//! searches for an optimal alignment between then. The two arrays are called `reference` and
//! `query`. An `Aligner` is intended to perform multiple queries against a reference, and so
//! `reference` is a required field when creating an `Aligner`. 
//! 
//! An optimal alignment fulfills all of these criteria:
//! 
//! - its error_rate is at most `max_error_rate`
//! - among those alignments with `error_rate <= max_error_rate`, the optimal alignment 
//!   contains a maximal number of matches (i.e. there is no alignment with more matches)
//! - among those alignments with the same numbers of matches (if there are multiple), the 
//!   optimal alignment has the minimal number of errors
//! - among the remaning alignments (if there are multiple), the optimal alignment is the one 
//!   with the leftmost starting position
//! 
//! The comparison between array elements is performed by a (user-configurable) function that 
//! returns an unsized integer (`usize`) cost. The default cost function returns `0` for matches 
//! and `1` for mismatches. The cost function is not told whether a mismatch is a substitution, 
//! insertion, or deletion; instead, different multipliers can be configured for each type of 
//! mismatch.
//! 
//! There is no special handling for ambiguous elements. Ambiguity is achieved by translating the
//! character sequences into binary representation, where an ambiguous character is a bitwise-or
//! of one-hot encoded canonical characters (e.g. `A=0b01, C=0b10, M=A|C == 0b11`. In addition,
//! the scoring function uses bitwise-and to determine if a pair of elements represent intersecting
//! sets of canonical characters, and thus match. See `mod utils` for implementations of common
//! translations and scoring functions.
//! 
//! ## `AnchoredAligner`
//! 
//! An `AnchoredAligner` is one that lines up the pair of arrays at one of their ends, and searches 
//! for an ungapped match starting from that end. `PrefixAligner` and `SuffixAligner` are 
//! `AnchoredAligner`s that match at the starts and ends of a pair, respectively. An 
//! `AnchoredAligner` uses simple pairwise comparison of array elements to search for the optimal 
//! alignment.
//! 
//! ## `DPAligner`
//! 
//! A `DPAligner` uses a dynamic programming algorithm to search for the optimal alignment
//! between a pair of arrays given a set of constraints. Global, local, and semi-global alignments
//! are possible. Gaps are allowed.
//! 
//! ## Todo
//! 
//! * Parameterize the data container for cost and positional information - for short reads,
//!   u8 is sufficient, but for long reads u16 is necessary

/// An optimal alignment between reference and query.
pub struct OptimalAlignment {
    pub ref_start: u16,
    pub ref_stop: u16,
    pub query_start: u16,
    pub query_stop: u16,
    pub num_matches: u16,
    pub total_cost: u16,
}

#[derive(Debug)]
pub trait Aligner {
    /// Locates one sequence within another by computing an optimal alignment between `reference` 
    /// and `query`.
    /// 
    /// # Examples
    ///
    /// An optimal semiglobal alignment of query 'SISSI' to reference 'MISSISSIPPI' looks like 
    /// this:
    /// 
    /// ```
    /// ---SISSI---
    /// MISSISSIPPI
    /// ```
    /// 
    /// The result is:
    /// 
    /// ```
    /// OptimalAlignment {
    ///     ref_start: 3,
    ///     ref_stop: 8,
    ///     query_start: 0,
    ///     query_stop: 5,
    ///     num_matches: 5,
    ///     cost: 0,
    /// }
    /// ```
    /// 
    /// The aligned parts are `reference[ref_start:ref_stop]` and `query[query_start:query_stop]`. 
    /// The error rate is: `errors / length`, where length is `ref_stop - ref_start == 5`.
    fn locate(&self, query: &[u8]) -> Optional<OptimalAlignment>;
}

/// A matching function, used to score a pair of elements.
pub trait MatchFunc {
    fn score(&self, a: u8, b: u8) -> usize;
}

/// The trait Matchfunc is implemented for Fn(u8, u8) -> usize so that Scoring can be 
/// instantiated using closures and custom user defined functions.
impl<F> MatchFunc for F
where
    F: Fn(u8, u8) -> usize,
{
    fn score(&self, a: u8, b: u8) -> usize {
        (self)(a, b)
    }
}

/// An exact matching function.
pub struct ExactMatch {
    pub match_cost: usize,
    pub mismatch_cost: usize,
}

impl MatchFunc for ExactMatch {
    fn score(&self, a: u8, b: u8) -> usize {
        if a == b {
            self.match_cost
        } else {
            self.mismatch_cost
        }
    }
}

/// A matching function for sequences with ambiguous elements
pub struct AmbiguousMatch {
    pub match_cost: usize,
    pub mismatch_cost: usize,
}

impl MatchFunc for AmbiguousMatch {
    fn score(&self, a: u8, b: u8) -> usize {
        if a & b != 0 {
            self.match_cost
        } else {
            self.mismatch_cost
        }
    }
}

lazy_static! {
    static ref EXACT_MATCH: ExactMatch = {
        ExactMatch {
            match_cost: 0,
            mismatch_cost: 1,
        }
    }

    static ref AMBIGUOUS_MATCH: AmbiguousMatch = {
        AmbiguousMatch {
            match_cost: 0,
            mismatch_cost: 1,
        }
    }
}

bitflags! {
    #[derive(Default)]
    pub struct AlignerFlags: u8 {
        const START_WITHIN_SEQ1 = 0b00001;
        const START_WITHIN_SEQ2 = 0b00010;
        const STOP_WITHIN_SEQ1 = 0b00100;
        const STOP_WITHIN_SEQ2 = 0b01000;
        const SEMIGLOBAL = 0b10000;
    }

    impl Default for AlignerFlags {
        fn default() -> Self { AlignerFlags::SEMIGLOBAL }
    }
}



/// Aligner options.
#[derive(Default)]
struct AlignerOpts {
    /// The maximum allowed error rate, as a fraction of the number of bases in the alignment
    pub max_error_rate: f64,
    /// A bitwise-or of of options that determine the cost for gaps at the starts and ends of the 
    /// reference and query strings. To allow skipping of a prefix of string1 at no cost, set the
    /// START_WITHIN_SEQ1 flag, and to allow skipping of a prefix of string2 at no cost, set the 
    /// START_WITHIN_SEQ2 flag. If both are set, a prefix of string1 or of string2 is skipped,
    /// but never both. Similarly, set STOP_WITHIN_SEQ1 and/or STOP_WITHIN_SEQ2 to allow skipping 
    /// of suffixes of string1 and/or string2 (again, when both flags are set, only one of the 
    /// suffixes is skipped). If all flags are set, this results in standard semiglobal alignment.
    pub flags: AlignerFlags,
    /// The minimum number of overlapping bases required between the reference and query
    pub min_overlap: u16,
    /// The cost for an indel; default is 1; ignored by anchored aligners
    pub indel_cost: u16,
    /// Whether wildcards are allowed in the reference sequence
    pub wildcard_ref: bool, 
    /// Whether wildcards are allowed in the query sequence
    pub wildcard_query: bool,
}

/// Anchored aligners perform simple pairwise comparison of bases starting at a fixed point, i.e. 
/// they do not use dynamic programming. Anchored aligners do not allow indels.
mod anchored {
    pub enum Anchor {
        Prefix,
        Suffix,
    }

    pub struct AnchoredAligner {
        pub reference: &[u8],
        max_mismatches: u16,
        min_overlap: u16,
    }
    
    impl AnchoredAligner {
        /// Creates a new `AnchoredAligner`.
        /// 
        /// `wildcard_symbols` is an optional set of values that may appear in the reference and
        /// that are interepreted as wildcards. Wildcards are ignored when determining the
        /// "effective" length of the sequence (for the purpose of computing the maxiumum number
        /// of errors allowed based on `max_error_rate`). A reference that consists only of
        /// wildcards is not allowed and will cause an error to be returned.
        fn new(
            reference: &[u8],
            anchor: Anchor,
            max_error_rate: f64, 
            min_overlap: u16,
            wildcard_values: Option<&HashSet<u8>>,
        ) -> Result<Self, &str> {
            let effective_len = match ignore_symbols {
                Some(x) => reference.len() - reference
                    .iter()
                    .filter(|s| s in wildcard_values.unwrap())
                    .count(),
                None => reference.len(),
            }
            if effective_len == 0 {
                return Err("Cannot have only wildcards in the sequence")
            }
            Ok(
                AnchoredAligner {
                    reference: reference,
                    max_mismatches: max_error_rate * effective_len,
                    min_overlap: min_overlap,
                }
            )
        }
    }
}

/// Aligners that use dynamic programming.
mod dp {
    /// A cell in the dynamic programming matrix.
    struct DPCell {
        /// The cost of the alignment
        cost: u16,
        /// The number of matches in this alignment 
        matches: u16,
        /// Where the alignment originated; negative for positions within seq1, positive for 
        /// positions within seq2
        origin: i16,
    }

    /// A candidate alignment between reference and query strings.
    struct DPAlignment {
        /// Where the alignment originated; negative for positions within seq1, positive for 
        /// positions within seq2
        origin: i16,
        /// Position at which the alignment stops in the reference string
        ref_stop: u16,
        /// Position at which the alignment stops in the query string
        query_stop: u16,
        /// The number of matches in this alignment 
        matches: u16,
        /// The cost of the alignment
        cost: u16,
    }

    impl DPAlignment {
        fn to_optimal(&self) -> OptimalAlignment {
            if self.origin >= 0 {
                let ref_start = 0
                let query_start = self.origin
            } else {
                let ref_start = -self.origin
                let query_start = 0
            }

            OptimalAlignment {
                ref_start: ref_start,
                ref_stop: self.ref_stop,
                query_start = query_start,
                query_stop = self.query_stop,
                matches = self.matches,
                cost = self.cost,
            }
        }
    }
}
