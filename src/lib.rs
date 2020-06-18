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
//! - among the remaining alignments (if there are multiple), the optimal alignment is the one
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
//! of one-hot encoded canonical characters (e.g. `A=0b0001, C=0b0010, M=A|C == 0b0011`. In 
//! addition, the scoring function uses bitwise-and to determine if a pair of elements represent 
//! intersecting sets of canonical characters, and thus match. See `mod utils` for implementations 
//! of common translations and scoring functions.
//! 
//! ## `AnchoredAligner`
//! 
//! An `AnchoredAligner` is one that lines up the pair of arrays at one of their ends, and searches 
//! for an ungapped match starting from that end. An `AnchoredAligner` uses simple pairwise 
//! comparison of array elements to search for the optimal alignment.
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

use std::cmp;

/// An optimal alignment between reference and query.
pub struct OptimalAlignment {
    pub reference_start: u16,
    pub reference_stop: u16,
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
    ///     reference_start: 3,
    ///     reference_stop: 8,
    ///     query_start: 0,
    ///     query_stop: 5,
    ///     num_matches: 5,
    ///     cost: 0,
    /// }
    /// ```
    /// 
    /// The aligned parts are `reference[reference_start:reference_stop]` and `query
    /// [query_start:query_stop]`. The error rate is: `errors / length`, where length is 
    /// `reference_stop - reference_start == 5`.
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

/// An exact matching function - two elements match if they are equal, otherwise mismatch.
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

/// A matching function for sequences with ambiguous elements - two elements match if they
/// intersect (i.e. share at least one bit in common), otherwise mismatch
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

/// Default matching functions
lazy_static! {
    static ref DEFAULT_EXACT_MATCH: ExactMatch = {
        ExactMatch {
            match_cost: 0,
            mismatch_cost: 1,
        }
    }

    static ref DEFAULT_AMBIGUOUS_MATCH: AmbiguousMatch = {
        AmbiguousMatch {
            match_cost: 0,
            mismatch_cost: 1,
        }
    }
}

pub struct SequenceTranslator {
    pub reference_wildcards: bool,
    pub query_wildcards: bool
}

impl SequenceTranslator {
    /// Translates the reference sequence.
    /// 
    /// Returns tuple (encoded_reference, effective_length)
    fn translate_reference(&self, seq: &str): (&[u8], usize) {
        if self.reference_wildcards || self.query_wildcards {
            
        } else {

        }

        let effective_length = match ignore_symbols {
            Some(x) => reference.len() - reference
                .iter()
                .filter(|s| s in wildcard_values.unwrap())
                .count(),
            None => reference.len(),
        }
        if effective_length == 0 {
            return Err("Cannot have only wildcards in the sequence")
        }
    }

    /// Translates a query sequence
    fn translate_query(&self, seq: &str): &[u8] {

    }
}

/// Anchored aligners perform simple pairwise comparison of bases starting at a fixed point, i.e. 
/// they do not use dynamic programming. Anchored aligners do not allow indels.
mod anchored {
    /// Enumeration of anchor points
    pub enum Anchor {
        /// Anchor at the beginning of the sequences
        Prefix,
        /// Anchor at the end of the sequences
        Suffix,
    }

    pub struct AnchoredAligner {
        /// The reference sequence
        pub reference: &[u8],
        /// Translator for query sequences
        sequence_translator: &SequenceTranslator,
        /// Maximum number of mismatches allowed
        max_mismatches: u16,
        /// Minimum number of overlapping positions required
        min_overlap: u16,
    }
    
    impl AnchoredAligner {
        /// Creates a new `AnchoredAligner`.
        /// 
        /// Arguments:
        /// 
        /// * `reference` - The reference sequence
        /// * `anchor` - The end of the sequences to anchor
        /// * `max_error_rate` - The maximum allowed error rate, as a fraction of the number of 
        ///    bases in the alignment
        /// * `min_overlap` - The minimum number of overlapping bases required between the 
        ///    reference and query
        /// * `wildcard_values` - an optional set of values that may appear in the reference and
        ///    that are interepreted as wildcards. Wildcards are ignored when determining the
        ///    "effective" length of the sequence (for the purpose of computing the maxiumum number
        ///    of errors allowed based on `max_error_rate`). A reference that consists only of
        ///    wildcards is not allowed and will cause an error to be returned.
        fn new(
            reference: &str,
            sequence_translator: &SequenceTranslator,
            anchor: Anchor,
            max_error_rate: f64, 
            min_overlap: u16,
            wildcard_values: Option<&HashSet<u8>>,
        ) -> Result<Self, &str> {
            let (translated_reference, effective_length) = 
                sequence_translator.translate_reference(reference)
            if not (0 <= max_error_rate <= 1.):
                raise ValueError("max_error_rate must be between 0 and 1")
            Ok(
                AnchoredAligner {
                    reference: translated_reference,
                    sequence_translator: sequence_translator,
                    max_mismatches: max_error_rate * effective_length,
                    min_overlap: min_overlap,
                }
            )
        }
    }
}

/// Aligners that use dynamic programming.
mod dp {
    /// Enumeration of flags that control how alignment is performed
    bitflags! {
        #[derive(Default)]
        pub struct DPAlignerFlags: u8 {
            const START_WITHIN_REFERENCE = 0b0001;
            const START_WITHIN_QUERY = 0b0010;
            const STOP_WITHIN_REFERENCE = 0b0100;
            const STOP_WITHIN_QUERY = 0b1000;
            const SEMIGLOBAL = 0b1111;
        }
    }

    impl Default for DPAlignerFlags {
        fn default() -> Self { AlignerFlags::SEMIGLOBAL }
    }

    /// Aligner options.
    #[derive(Default)]
    struct CostCoefficients {
        /// The cost for substitutions
        pub substitution: u16,
        /// The cost for insertions
        pub insertion: u16,
        /// The cost for deletions
        pub deletion: u16,
    }

    impl Default for CostCoefficients {
        fn default() -> Self {
            substitution: 1,
            insertion: 1,
            deletion: 1,
        }
    }

    /// A cell in the dynamic programming matrix.
    #[derive(Clone)]
    struct DPCell {
        /// The cost of the alignment
        cost: u16,
        /// The number of matches in this alignment 
        matches: u16,
        /// Where the alignment originated; negative for positions within seq1, positive for 
        /// positions within seq2
        origin: i16,
    }

    /// Define constructors for the four possible combinations of START_WITHIN_REFERENCE and
    /// START_WITHIN_QUERY (where a false value implies being anchored to the left side of the
    /// DP matrix)
    impl DPCell {
        fn new_both_anchored(column: u16, min_column: u16, costs: CostCoefficients) -> Self {
            DPCell {
                cost: cmp::max(columns, min_column) * costs.insertion
                matches: 0,
                origin: 0
            }
        }

        fn new_query_anchored(column: u16, min_column: u16, costs: CostCoefficients) -> Self {
            DPCell {
                cost: min_column * costs.insertion,
                matches = 0,
                origin = cmp::min(0, min_column - column)
            }
        }

        fn new_reference_anchored(column: u16, min_column: u16, costs: CostCoefficients) -> Self {
            DPCell {
                cost: column * costs.insertion,
                matches: 0,
                origin: cmp::max(0, min_column - column)
            }
        }

        fn new_neither_anchored(column: u16, min_column: u16, costs: CostCoefficients) -> Self {
            DPCell {
                cost: min(column, min_column) * costs.insertion,
                matches = 0,
                origin = min_column - column
            }
        }

        fn update(mut &self, (new_cost, new_matches, new_origin): (u16, u16, u16)) {
            self.cost = new_cost
            self.matches = new_matches
            self.origin = new_origin
        }
    }
    
    /// A candidate alignment between reference and query strings.
    struct DPAlignment {
        /// Where the alignment originated; negative for positions within seq1, positive for 
        /// positions within seq2
        origin: i16,
        /// Position at which the alignment stops in the reference string
        reference_stop: u16,
        /// Position at which the alignment stops in the query string
        query_stop: u16,
        /// The number of matches in this alignment 
        matches: u16,
        /// The cost of the alignment
        cost: u16,
    }

    impl DPAlignment {
        fn update(&mut self, &cell: DPCell) {
            self.origin = cell.origin
            self.matches = cell.matches
            self.cost = cell.cost
        }

        fn to_optimal(&self) -> OptimalAlignment {
            if self.origin >= 0 {
                let reference_start = 0
                let query_start = self.origin
            } else {
                let reference_start = -self.origin
                let query_start = 0
            }

            OptimalAlignment {
                reference_start: reference_start,
                reference_stop: self.reference_stop,
                query_start = query_start,
                query_stop = self.query_stop,
                matches = self.matches,
                cost = self.cost,
            }
        }
    }

    pub struct DPAligner {
        columns: Vector[]
    }

    impl DPAligner {
        /// 
        /// Arguments:
        /// 
        /// * `reference` - The reference sequence
        /// * `sequence_translator` - Translates between character and byte sequences
        /// * `flags` - A bitwise-or of of options that determine the cost for gaps at the starts 
        ///   and ends of the reference and query strings. To allow skipping of a prefix of the 
        ///   reference at no cost, set the START_WITHIN_REFERENCE flag, and to allow skipping of a 
        ///   prefix of the query at no cost, set the START_WITHIN_QUERY flag. If both are set, a 
        ///   prefix of reference or of query is skipped, but never both. Similarly, set 
        ///   STOP_WITHIN_REFERENCE and/or STOP_WITHIN_QUERY to allow skipping of suffixes of 
        ///   reference and/or query (again, when both flags are set, only one of the suffixes is 
        ///   skipped). If all flags are set, this results in standard semiglobal alignment.
        /// * `max_error_rate` - The maximum allowed error rate, as a fraction of the number of 
        ///    bases in the alignment
        /// * `min_overlap` - The minimum number of overlapping bases required between the 
        ///    reference and query
        /// * `wildcard_values` - an optional set of values that may appear in the reference and
        ///    that are interepreted as wildcards. Wildcards are ignored when determining the
        ///    "effective" length of the sequence (for the purpose of computing the maxiumum number
        ///    of errors allowed based on `max_error_rate`). A reference that consists only of
        ///    wildcards is not allowed and will cause an error to be returned.
        fn new(
            reference: &str,
            sequence_translator: &SequenceTranslator,
            pub flags: AlignerFlags,
            max_error_rate: f64, 
            min_overlap: u16,
            wildcard_values: Option<&HashSet<u8>>,
        ) -> Result<Self, &str> {
            let (translated_reference, effective_length) = 
                sequence_translator.translate_reference(reference);
            
            if not (0 <= max_error_rate <= 1.) {
                raise ValueError("max_error_rate must be between 0 and 1");
            }
            
            let start_in_reference: bool = flags & DPAlignerFlags::START_WITHIN_REFERENCE;
            let start_in_query: bool = flags & DPAlignerFlags::START_WITHIN_QUERY;
            let column_factory = match (start_in_reference, start_in_query) {
                (false, false) => DPCell::new_both_anchored,
                (true, false) => DPCell::new_query_anchored,
                (false, true) => DPCell::new_reference_anchored,
                (true, true) => DPCell::new_neither_anchored
            };

            Ok(
                DPAligner {
                    
                }
            )
        }
    }

    impl Aligner for DPAligner {
        fn locate(&self, query: &[u8]) -> Optional<OptimalAlignment> {
            let query_length = query.len()
            let start_in_reference: bool = flags & DPAlignerFlags::START_WITHIN_REFERENCE
            let start_in_query: bool = flags & DPAlignerFlags::START_WITHIN_QUERY;
            let stop_in_query: bool = flags & DPAlignerFlags::STOP_WITHIN_QUERY;
            
            // Determine largest and smallest columns we need to compute
            let min_column = if stop_in_query { 
                0
            } else { 
                cmp::max(0, query_length - self.reference-length - self.max_mismatches)
            }
            let max_column = if start_in_query {
                query_length
            } else {
                min(query_length, self.reference_length + self.max_mismatches)
            }

            // Create the DP matrix
            // TODO: make columns thread-local and pre-allocate the vector, then just update the 
            // values each time
            let columns: Vec<DPCell> = (0..(reference.length + 1))
                .iter()
                .map(|i| column_factory(i, min_column, self.costs))
                .collect()

            // The index of the last column
            let mut last_column = if start_in_reference {
                self.reference_length;
            } else {
                // Use Ukkonen's trick: compute the index of the last cell that is at most k
                cmp::min(self.reference_length, self.max_mismatches + 1);
            }

            let mut best_alignment = DPAlignment {
                reference_stop: self.reference_length,
                query_stop: self.query_length,
                cost: self.reference_length + self.query_length,
                origin: 0,
                matches: 0
            }

            for query_idx in (min_column + 1)..(max_column + 1) {
                let mut diagonal_cell = columns[0].clone()

                if start_in_query {
                    columns[0].origin = query_idx;
                } else {
                    columns[0].cost = query_idx * self.costs.insertion;
                }

                for reference_idx in 1..(last_column + 1) {
                    let score = self
                        .match_func
                        .score(self.reference[reference_idx - 1], query[query_idx - 1]));
                    
                    let mut reference_cell = columns[reference_idx]

                    let new_values = if score == 0 {
                        (diagonal_cell.cost, diagonal_cell.origin, diagonal_cell.matches + 1);
                    } else {
                        let cost_diagonal = diagonal_cell.cost + self.costs.mismatch
                        let cost_insertion = reference_cell.cost + self.costs.insertion
                        let prev_cell = columns[reference_idx - 1]
                        let cost_deletion = prev_cell.cost + self.costs.deletion

                        if cost_diagonal <= cost_insertion && cost_diagonal <= cost_deletion {
                            // mismatch
                            (cost_diagonal, diagonal_cell.origin, diagonal_cell.matches);
                        } else if cost_insertion <= cost_deletion {
                            // insertion
                            (cost_insertion, prev_cell.origin, prev_cell.matches);
                        } else {
                            // deletion
                            (cost_deletion, reference_cell.origin, reference_cell.matches);
                        }
                    };
                    
                    diagonal_cell = reference_cell.clone();

                    reference_cell.update(new_values);
                }

                while last_column >= 0 && columns[last_column].cost > self.max_errors {
                    last -= 1;
                }

                if last < self.reference_length {
                    last += 1;
                } else if stop_in_query {
                    // Found a match
                    let match_column = columns[self.reference_length];
                    let reference_alignment_length = self.reference_length 
                        + cmp::min(match_column.origin, 0);
                    
                    let effective_length = reference_alignment_length;
                    
                    if reference_alignment_length >= self.min_overlap
                        && match_column.cost <= effective_length * self.max_error_rate
                        && (match_column.matches > best_alignment.matches
                            || (match_column.matches == best_alignment.matches
                                && cost < best_alignment.cost)) {
                        // update best match
                        best_alignment.update(match_column);

                        if match_column.cost == 0 && match_column.matches == self.reference_length {
                            break;
                        }
                    }
                }
            }

            if max_column == self.query_length {
                let stop_in_reference: bool = flags & DPAlignerFlags.STOP_WITHIN_REFERENCE
                let first_column = if stop_in_reference {
                    0
                } else {
                    self.reference_length
                }
                for reference_idx in first_column..(self.reference_length + 1) {
                    
                }
            }
        }
    }
}
