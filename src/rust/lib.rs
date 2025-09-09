#[cfg(feature = "python")]
use pyo3::prelude::*;
#[cfg(feature = "python")]
use pyo3::types::PyType;
#[cfg(feature = "python")]
use std::path::PathBuf;

pub mod aminoacid;
pub mod config;
pub mod encoding;
pub mod errors;
pub mod index;
pub mod io;
pub mod iterators;
pub mod kmer;
pub mod metrics;
pub mod search;
pub mod signature;
pub mod types;

#[cfg(test)]
mod tests;

// Re-export main types for easier access
pub use aminoacid::AminoAcidAmbiguity;
pub use index::ProteomeIndex;
pub use signature::SEED;

#[cfg(feature = "python")]
#[pyclass]
#[derive(Clone, Copy)]
pub enum PyProteinEncoding {
    Raw,
    Dayhoff,
    HP,
}

#[cfg(feature = "python")]
#[pymethods]
impl PyProteinEncoding {
    #[classmethod]
    fn raw(_cls: Py<PyType>) -> Self {
        Self::Raw
    }

    #[classmethod]
    fn dayhoff(_cls: Py<PyType>) -> Self {
        Self::Dayhoff
    }

    #[classmethod]
    fn hp(_cls: Py<PyType>) -> Self {
        Self::HP
    }
}

#[cfg(feature = "python")]
impl From<PyProteinEncoding> for String {
    fn from(py_encoding: PyProteinEncoding) -> Self {
        match py_encoding {
            PyProteinEncoding::Raw => "protein".to_string(),
            PyProteinEncoding::Dayhoff => "dayhoff".to_string(),
            PyProteinEncoding::HP => "hp".to_string(),
        }
    }
}

#[cfg(feature = "python")]
/// Formats the sum of two numbers as string.
#[pyfunction]
fn sum_as_string(a: usize, b: usize) -> PyResult<String> {
    Ok((a + b).to_string())
}

#[cfg(feature = "python")]
#[pyclass]
pub struct PyProteomeIndex {
    // Not currently used, but will be used in the future
    #[allow(dead_code)]
    index: crate::index::ProteomeIndex,
}

#[cfg(feature = "python")]
#[pymethods]
impl PyProteomeIndex {
    #[new]
    fn new(ksize: u32, scaled: u32, moltype: &str, db_path: &str) -> PyResult<Self> {
        let index =
            crate::index::ProteomeIndex::new(PathBuf::from(db_path), ksize, scaled, moltype, false)
                .map_err(|e| PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(e.to_string()))?;

        Ok(PyProteomeIndex { index })
    }
}

#[cfg(feature = "python")]
/// A Python module implemented in Rust.
#[pymodule]
fn kmerseek(_py: Python<'_>, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(sum_as_string, m)?)?;
    m.add_class::<PyProteinEncoding>()?;
    m.add_class::<PyProteomeIndex>()?;
    Ok(())
}
