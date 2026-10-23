//! A tiny "Python number" type used where the original code may hold either an int or a
//! float and the textual output depends on which one it is.

use crate::py::fmt::float_repr;

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum PyNum {
    Int(i64),
    Float(f64),
}

impl PyNum {
    pub fn f(self) -> f64 {
        match self {
            PyNum::Int(i) => i as f64,
            PyNum::Float(x) => x,
        }
    }
    /// Python `str()`
    pub fn py_str(self) -> String {
        match self {
            PyNum::Int(i) => i.to_string(),
            PyNum::Float(x) => float_repr(x),
        }
    }
    /// Python `int()`
    pub fn int(self) -> i64 {
        match self {
            PyNum::Int(i) => i,
            PyNum::Float(x) => x.trunc() as i64,
        }
    }
    pub fn is_float(self) -> bool {
        matches!(self, PyNum::Float(_))
    }
}

impl std::fmt::Display for PyNum {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.py_str())
    }
}

/// Python `a + b` on PyNum
pub fn add(a: PyNum, b: PyNum) -> PyNum {
    match (a, b) {
        (PyNum::Int(x), PyNum::Int(y)) => PyNum::Int(x + y),
        _ => PyNum::Float(a.f() + b.f()),
    }
}

/// A Python value that is written with `str()` into xTea's text files.
#[derive(Clone, Debug, PartialEq)]
pub enum PyVal {
    Int(i64),
    Float(f64),
    Str(String),
}

impl PyVal {
    pub fn py_str(&self) -> String {
        match self {
            PyVal::Int(i) => i.to_string(),
            PyVal::Float(x) => float_repr(*x),
            PyVal::Str(s) => s.clone(),
        }
    }
    /// Python `int(str(v))`-like conversion used when values are re-read.
    pub fn as_i64(&self) -> i64 {
        match self {
            PyVal::Int(i) => *i,
            PyVal::Float(x) => x.trunc() as i64,
            PyVal::Str(s) => s.trim().parse::<i64>().unwrap_or(0),
        }
    }
    pub fn as_f64(&self) -> f64 {
        match self {
            PyVal::Int(i) => *i as f64,
            PyVal::Float(x) => *x,
            PyVal::Str(s) => s.trim().parse::<f64>().unwrap_or(0.0),
        }
    }
}

impl From<i64> for PyVal {
    fn from(v: i64) -> Self {
        PyVal::Int(v)
    }
}
impl From<f64> for PyVal {
    fn from(v: f64) -> Self {
        PyVal::Float(v)
    }
}
impl From<String> for PyVal {
    fn from(v: String) -> Self {
        PyVal::Str(v)
    }
}
impl From<&str> for PyVal {
    fn from(v: &str) -> Self {
        PyVal::Str(v.to_string())
    }
}
