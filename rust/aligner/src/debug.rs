/// Defines a data structure (DPMatrix) that stores a dynamic programming matrix and can print it
/// nicely for the user. For debugging purposes (debug_assertions must be true).
#[cfg(debug_assertions)]
mod dpmatrix {
    extern crate ndarray;
    use ndarray::Array;

    /// A dynamic programming matrix.
    #[derive(Debug)]
    struct DPMatrix {
        /// The matrix itself; uses a signed data container (e.g. i16) so that a negative value
        /// can be used to represent uninitialized cells.
        matrix: Array2<i16>,
        /// The reference used to create the matrix
        reference: &'static [u8],
        /// The query used to create the matrix
        query: &'static [u8],
    }

    impl DPMatrix {
        /// Builds a dynamic programming matrix from a reference and query string.
        /// 
        /// # Arguments
        /// 
        /// * reference: The reference string
        /// * query: The query string
        fn new(reference: &'static [u8], query: &'static [u8]) -> DPMatrix {
            let m = reference.len()
            let n = query.len()

            DPMatrix {
                matrix: Array::from_element((m + 1, n + 1), -1),
                reference: reference,
                query: query
            }
        }

        fn set_entry(&self, i: usize, j: usize, cost: i16) {
            self.matrix.view_mut()[(i, j)] = cost
        }
    }

    impl fmt::Debug for DPMatrix {
        // TODO: translate python to rust
        // rows = ["     " + " ".join(c.rjust(2) for c in self.query)]
        // for c, row in zip(" " + self.reference, self._rows):
        //     r = c + " " + " ".join("  " if v is None else f"{v:2d}" for v in row)
        //     rows.append(r)
        // return "\n".join(rows)
        fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
            write!(f, "matrix")
        }
    }
}
