extern create ascii;
#[macro_use] extern crate bitflags;

use ascii::{AsciiChar, AsciiStr};
use std::collections::HashMap;

trait Sequence {
    /// Computes the effective length of the sequence, which is the length minus the number of
    /// 'N' characters.
    fn effective_len(&self) -> int;

    /// Translate this sequence into an array of values using a translation table.
    fn translate(&self, table: HashMap<AsciiChar, T>) -> &ArrayVec<T>;
}

impl Sequence for AsciiStr {
    fn effective_len(&self) -> usize {
        self.chars().filter(|c| c != AsciiChar::N && c != AsciiChar::n).count()
    }

    fn translate(&self, table: HashMap<AsciiChar, T>) -> ArrayVec<T> {
        self.map(|x| table[x]).collect::<ArrayVec<T; self.len()>>();
    }
}

/// Compare position `i` in `self` with position `j` in `other`. Returns `true` if the
    /// sets of possible characters at the two positions intersect.
    fn compare_positions(&self, i: u16, other: Sequence, j: u16) -> bool;
/// A `Sequence` that allows ASCII characters. Comparison is performed between sequences by 
/// checking for equality of the characters at corresponding positions.
mod ascii {
    struct AsciiSequence {
        pub seq: &AsciiStr
    }

    impl Sequence for AsciiSequence {
        fn compare(&self, i: u16, other: Sequence, j: u16) -> bool {
            return self.seq[i] == other.seq[j]
        }
    }
}

/// Wildcards may be used in the reference and/or query string. IUPAC wildcard characters can be 
/// allowed in the reference and/or query by setting the appropriate flags in `AlignerOpts`. If 
/// any of the flags is set, all non-IUPAC characters in the sequences compare as 'not equal'.
mod wildcards {
    bitflags! {
        pub struct Bases: u8 {
            const A = 0b0001;
            const C = 0b0010;
            const G = 0b0100;
            const T = 0b1000;  // also used for 'U'
        }
    }

    // char -> int translation
    // Borrowed code from:
    // * https://stackoverflow.com/questions/27582739/how-do-i-create-a-hashmap-literal
    // * https://users.rust-lang.org/t/how-to-do-pythons-str-translate-in-rus/16229/6

    macro_rules! char_int_map(
        { $($key:expr => $value:expr),+ } => {
            {
                let mut m = HashMap::new();
                $(
                    m.insert($key.to_ascii_uppercase(), $value);
                    m.insert($key.to_ascii_lowercase(), $value)
                )+
                m
            }
        };
    );

    lazy_static! {
        static ref ACGT_TABLE: HashMap<AsciiChar, Bases> = {
            char_int_map!{
                AsciiChar::A => Bases::A,
                AsciiChar::C => Bases::C,
                AsciiChar::G => Bases::G,
                AsciiChar::T => Bases::T,
                AsciiChar::U => Bases::T,
            }
        };

        static ref IUPAC_TABLE: HashMap<AsciiChar, Bases> = {
            char_int_map!{
                AsciiChar::A => Bases::A,
                AsciiChar::C => Bases::C,
                AsciiChar::G => Bases::G,
                AsciiChar::T => Bases::T,
                AsciiChar::U => Bases::T,
                AsciiChar::M => Bases::A | Bases::C,
                AsciiChar::R => Bases::A | Bases::G,
                AsciiChar::W => Bases::A | Bases::T,
                AsciiChar::S => Bases::C | Bases::G,
                AsciiChar::Y => Bases::C | Bases::T,
                AsciiChar::K => Bases::G | Bases::T,
                AsciiChar::V => Bases::A | Bases::C | Bases::G,
                AsciiChar::H => Bases::A | Bases::C | Bases::T,
                AsciiChar::D => Bases::A | Bases::G | Bases::T,
                AsciiChar::B => Bases::C | Bases::G | Bases::T,
                AsciiChar::N => Bases::A | Bases::C | Bases::G | Bases::T,
            }
        };
    }

    struct WildcardSequence {
        /// The original sequence
        pub seq: &AsciiStr,
        /// The translated sequence
        translated: Box<ArrayVec<u8>>
    }
    
    impl BaseComparator for WildcardSequence {
        fn new(seq: &AsciiStr, table: HashMap<AsciiChar, Bases>) {
            
            WildcardSequence { seq: seq, translated: translated } 
        }

        fn compare(&self, i: u16, other: WildcardSequence, j: u16) -> bool {
            return (self.translated[i] & other.translated[j]) != 0
        }
    }
}
