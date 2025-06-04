use pyo3::prelude::*;
use pyo3::types::PyType;
mod index;

use crate::index::ProteomeIndex;
use std::path::PathBuf;

#[pyclass]
#[derive(Clone, Copy)]
pub enum PyProteinEncoding {
    Raw,
    Dayhoff,
    HP,
}

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

impl From<PyProteinEncoding> for String {
    fn from(py_encoding: PyProteinEncoding) -> Self {
        match py_encoding {
            PyProteinEncoding::Raw => "protein".to_string(),
            PyProteinEncoding::Dayhoff => "dayhoff".to_string(),
            PyProteinEncoding::HP => "hp".to_string(),
        }
    }
}

/// Formats the sum of two numbers as string.
#[pyfunction]
fn sum_as_string(a: usize, b: usize) -> PyResult<String> {
    Ok((a + b).to_string())
}

/// A Python module implemented in Rust.
#[pymodule]
fn kmerseek(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(sum_as_string, m)?)?;
    m.add_class::<PyProteinEncoding>()?;
    m.add_class::<PyProteomeIndex>()?;
    Ok(())
}

#[pyclass]
struct PyProteomeIndex {
    inner: ProteomeIndex,
}

#[pymethods]
impl PyProteomeIndex {
    #[new]
    fn new(ksize: u32, scaled: u32, moltype: &str, encoding: PyProteinEncoding) -> Self {
        Self {
            inner: ProteomeIndex::new(ksize, scaled, moltype, &String::from(encoding)),
        }
    }

    fn process_protein_files(&self, files: Vec<String>) -> PyResult<()> {
        // Convert strings to PathBuf
        let paths: Vec<PathBuf> = files.into_iter().map(PathBuf::from).collect();

        self.inner
            .process_protein_files(&paths)
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(e.to_string()))?;
        Ok(())
    }

    fn get_kmers(&self, hash: u64) -> Option<(String, String)> {
        self.inner.get_kmers(hash)
    }

    fn add_sequence(&mut self, seq: &str, force: bool) {
        self.inner.add_sequence(seq, force);
    }
}
