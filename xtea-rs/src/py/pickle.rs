//! A small pickle reader (protocols 0-5, without out-of-band buffers) with support for
//! joblib's `NumpyArrayWrapper`, i.e. the format `joblib.dump` uses for numpy arrays.
//!
//! It is only meant to load the objects xTea needs from Python-generated files (the
//! deep-forest genotyping model).  Python objects are not instantiated: class instances
//! become [`Object`] records holding the class path, the constructor arguments and the
//! state passed to `__setstate__`, which is enough to read their attributes.

use anyhow::{anyhow, bail, Context, Result};
use std::cell::RefCell;
use std::collections::HashMap;
use std::rc::Rc;

/// A numpy array read from a pickle (data kept as raw bytes).
#[derive(Debug, Clone)]
pub struct NdArray {
    /// numpy dtype kind: 'i', 'u', 'f' or 'b'.
    pub kind: char,
    pub itemsize: usize,
    pub little_endian: bool,
    pub shape: Vec<usize>,
    /// true if `data` is stored in Fortran (column-major) order.
    pub fortran: bool,
    pub data: Vec<u8>,
}

impl NdArray {
    pub fn len(&self) -> usize {
        self.shape.iter().product()
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Offset (in items) of the element with C-order flat index `i`.
    fn offset_of(&self, i: usize) -> usize {
        if !self.fortran || self.shape.len() <= 1 {
            return i;
        }
        // decompose the C-order index, then recompose in Fortran order
        let mut rem = i;
        let mut idx = vec![0usize; self.shape.len()];
        for d in (0..self.shape.len()).rev() {
            idx[d] = rem % self.shape[d];
            rem /= self.shape[d];
        }
        let mut off = 0usize;
        let mut stride = 1usize;
        for d in 0..self.shape.len() {
            off += idx[d] * stride;
            stride *= self.shape[d];
        }
        off
    }

    fn raw(&self, i: usize) -> &[u8] {
        let o = self.offset_of(i) * self.itemsize;
        &self.data[o..o + self.itemsize]
    }

    fn item_bytes<const N: usize>(&self, i: usize) -> [u8; N] {
        let mut b = [0u8; N];
        b.copy_from_slice(self.raw(i));
        if !self.little_endian {
            b.reverse();
        }
        b
    }

    /// Element `i` (C order) as i64; integer and boolean arrays only.
    pub fn get_i64(&self, i: usize) -> Result<i64> {
        Ok(match (self.kind, self.itemsize) {
            ('i', 1) => i8::from_le_bytes(self.item_bytes::<1>(i)) as i64,
            ('i', 2) => i16::from_le_bytes(self.item_bytes::<2>(i)) as i64,
            ('i', 4) => i32::from_le_bytes(self.item_bytes::<4>(i)) as i64,
            ('i', 8) => i64::from_le_bytes(self.item_bytes::<8>(i)),
            ('u', 1) | ('b', 1) => self.item_bytes::<1>(i)[0] as i64,
            ('u', 2) => u16::from_le_bytes(self.item_bytes::<2>(i)) as i64,
            ('u', 4) => u32::from_le_bytes(self.item_bytes::<4>(i)) as i64,
            ('u', 8) => u64::from_le_bytes(self.item_bytes::<8>(i)) as i64,
            (k, s) => bail!("unsupported integer dtype {}{}", k, s),
        })
    }

    /// Element `i` (C order) as f64 (integers are converted).
    pub fn get_f64(&self, i: usize) -> Result<f64> {
        Ok(match (self.kind, self.itemsize) {
            ('f', 4) => f32::from_le_bytes(self.item_bytes::<4>(i)) as f64,
            ('f', 8) => f64::from_le_bytes(self.item_bytes::<8>(i)),
            _ => self.get_i64(i)? as f64,
        })
    }

    pub fn to_i64_vec(&self) -> Result<Vec<i64>> {
        (0..self.len()).map(|i| self.get_i64(i)).collect()
    }

    pub fn to_f64_vec(&self) -> Result<Vec<f64>> {
        (0..self.len()).map(|i| self.get_f64(i)).collect()
    }
}

/// An instance of a Python class (or the result of calling a global).
#[derive(Debug)]
pub struct Object {
    pub module: String,
    pub name: String,
    pub args: Vec<Value>,
    pub state: Option<Value>,
    /// set when the object is (or wraps) a numpy array
    pub array: Option<Rc<NdArray>>,
}

#[derive(Debug, Clone)]
pub enum Value {
    None,
    Bool(bool),
    Int(i64),
    Float(f64),
    Str(Rc<str>),
    Bytes(Rc<Vec<u8>>),
    Tuple(Rc<Vec<Value>>),
    List(Rc<RefCell<Vec<Value>>>),
    Dict(Rc<RefCell<Vec<(Value, Value)>>>),
    Set(Rc<RefCell<Vec<Value>>>),
    Global(Rc<(String, String)>),
    Object(Rc<RefCell<Object>>),
    Array(Rc<NdArray>),
}

impl Value {
    pub fn type_name(&self) -> String {
        match self {
            Value::None => "None".into(),
            Value::Bool(_) => "bool".into(),
            Value::Int(_) => "int".into(),
            Value::Float(_) => "float".into(),
            Value::Str(_) => "str".into(),
            Value::Bytes(_) => "bytes".into(),
            Value::Tuple(_) => "tuple".into(),
            Value::List(_) => "list".into(),
            Value::Dict(_) => "dict".into(),
            Value::Set(_) => "set".into(),
            Value::Global(g) => format!("global {}.{}", g.0, g.1),
            Value::Object(o) => {
                let o = o.borrow();
                format!("{}.{}", o.module, o.name)
            }
            Value::Array(_) => "ndarray".into(),
        }
    }

    /// `module.name` of an object's class.
    pub fn class_path(&self) -> Option<String> {
        match self {
            Value::Object(o) => {
                let o = o.borrow();
                Some(format!("{}.{}", o.module, o.name))
            }
            _ => None,
        }
    }

    fn key_matches(k: &Value, key: &str) -> bool {
        matches!(k, Value::Str(s) if &**s == key)
    }

    /// `d[key]` for a dict with string keys.
    pub fn get_item(&self, key: &str) -> Option<Value> {
        match self {
            Value::Dict(d) => d
                .borrow()
                .iter()
                .find(|(k, _)| Self::key_matches(k, key))
                .map(|(_, v)| v.clone()),
            _ => None,
        }
    }

    /// Attribute lookup on an object restored through `__setstate__` / `__dict__`.
    pub fn getattr(&self, name: &str) -> Option<Value> {
        match self {
            Value::Object(o) => {
                let o = o.borrow();
                match &o.state {
                    Some(st @ Value::Dict(_)) => st.get_item(name),
                    Some(Value::Tuple(t)) if t.len() == 2 => {
                        t[0].get_item(name).or_else(|| t[1].get_item(name))
                    }
                    _ => None,
                }
            }
            Value::Dict(_) => self.get_item(name),
            _ => None,
        }
    }

    pub fn req_attr(&self, name: &str) -> Result<Value> {
        self.getattr(name)
            .ok_or_else(|| anyhow!("pickled object {} has no attribute '{}'", self.type_name(), name))
    }

    pub fn as_i64(&self) -> Result<i64> {
        match self {
            Value::Int(v) => Ok(*v),
            Value::Bool(b) => Ok(*b as i64),
            Value::Array(a) if a.len() == 1 => a.get_i64(0),
            Value::Object(o) => match &o.borrow().array {
                Some(a) if a.len() == 1 => a.get_i64(0),
                _ => bail!("expected an integer, found {}", self.type_name()),
            },
            _ => bail!("expected an integer, found {}", self.type_name()),
        }
    }

    pub fn as_f64(&self) -> Result<f64> {
        match self {
            Value::Float(v) => Ok(*v),
            Value::Int(v) => Ok(*v as f64),
            Value::Bool(b) => Ok(*b as i64 as f64),
            _ => bail!("expected a float, found {}", self.type_name()),
        }
    }

    /// Python truthiness for the simple types used in model parameters.
    pub fn truthy(&self) -> bool {
        match self {
            Value::None => false,
            Value::Bool(b) => *b,
            Value::Int(v) => *v != 0,
            Value::Float(v) => *v != 0.0,
            Value::Str(s) => !s.is_empty(),
            Value::Bytes(b) => !b.is_empty(),
            Value::Tuple(t) => !t.is_empty(),
            Value::List(l) => !l.borrow().is_empty(),
            Value::Dict(d) => !d.borrow().is_empty(),
            Value::Set(s) => !s.borrow().is_empty(),
            _ => true,
        }
    }

    pub fn as_str(&self) -> Result<String> {
        match self {
            Value::Str(s) => Ok(s.to_string()),
            _ => bail!("expected a str, found {}", self.type_name()),
        }
    }

    /// Items of a list or tuple.
    pub fn as_list(&self) -> Result<Vec<Value>> {
        match self {
            Value::List(l) => Ok(l.borrow().clone()),
            Value::Tuple(t) => Ok((**t).clone()),
            _ => bail!("expected a list, found {}", self.type_name()),
        }
    }

    pub fn as_array(&self) -> Result<Rc<NdArray>> {
        match self {
            Value::Array(a) => Ok(a.clone()),
            Value::Object(o) => o
                .borrow()
                .array
                .clone()
                .ok_or_else(|| anyhow!("expected a numpy array, found {}", self.type_name())),
            _ => bail!("expected a numpy array, found {}", self.type_name()),
        }
    }
}

struct Unpickler<'a> {
    data: &'a [u8],
    pos: usize,
    stack: Vec<Value>,
    marks: Vec<usize>,
    memo: HashMap<usize, Value>,
}

fn dtype_of(v: &Value) -> Result<(char, usize, bool)> {
    let (descr, order) = match v {
        Value::Object(o) => {
            let o = o.borrow();
            if o.name != "dtype" {
                bail!("expected a numpy dtype, found {}.{}", o.module, o.name);
            }
            let descr = o
                .args
                .first()
                .ok_or_else(|| anyhow!("numpy dtype without arguments"))?
                .as_str()?;
            let order = match &o.state {
                Some(Value::Tuple(t)) if t.len() > 1 => match &t[1] {
                    Value::Str(s) => s.to_string(),
                    _ => "=".to_string(),
                },
                _ => "=".to_string(),
            };
            (descr, order)
        }
        Value::Str(s) => (s.to_string(), "=".to_string()),
        _ => bail!("expected a numpy dtype, found {}", v.type_name()),
    };
    let mut chars = descr.chars();
    let mut kind = chars.next().ok_or_else(|| anyhow!("empty dtype"))?;
    let mut rest: String = chars.collect();
    let mut order = order;
    if kind == '<' || kind == '>' || kind == '|' || kind == '=' {
        order = kind.to_string();
        let mut c2 = rest.chars();
        kind = c2.next().ok_or_else(|| anyhow!("bad dtype {}", descr))?;
        rest = c2.collect();
    }
    if kind == 'O' {
        bail!("numpy object arrays are not supported");
    }
    if !matches!(kind, 'i' | 'u' | 'f' | 'b') {
        bail!("unsupported numpy dtype '{}'", descr);
    }
    let itemsize: usize = if rest.is_empty() && kind == 'b' {
        1
    } else {
        rest.parse().map_err(|_| anyhow!("unsupported numpy dtype '{}'", descr))?
    };
    let little = order != ">";
    Ok((kind, itemsize, little))
}

impl<'a> Unpickler<'a> {
    fn read(&mut self, n: usize) -> Result<&'a [u8]> {
        if self.pos + n > self.data.len() {
            bail!("pickle data was truncated");
        }
        let s = &self.data[self.pos..self.pos + n];
        self.pos += n;
        Ok(s)
    }

    fn u8(&mut self) -> Result<u8> {
        Ok(self.read(1)?[0])
    }

    fn u16(&mut self) -> Result<u16> {
        let b = self.read(2)?;
        Ok(u16::from_le_bytes([b[0], b[1]]))
    }

    fn u32(&mut self) -> Result<u32> {
        let b = self.read(4)?;
        Ok(u32::from_le_bytes([b[0], b[1], b[2], b[3]]))
    }

    fn i32(&mut self) -> Result<i32> {
        Ok(self.u32()? as i32)
    }

    fn u64(&mut self) -> Result<u64> {
        let b = self.read(8)?;
        let mut a = [0u8; 8];
        a.copy_from_slice(b);
        Ok(u64::from_le_bytes(a))
    }

    fn line(&mut self) -> Result<&'a str> {
        let start = self.pos;
        while self.pos < self.data.len() && self.data[self.pos] != b'\n' {
            self.pos += 1;
        }
        if self.pos >= self.data.len() {
            bail!("pickle data was truncated");
        }
        let s = std::str::from_utf8(&self.data[start..self.pos])
            .map_err(|_| anyhow!("invalid text line in pickle"))?;
        self.pos += 1;
        Ok(s)
    }

    fn pop(&mut self) -> Result<Value> {
        if let Some(&m) = self.marks.last() {
            if self.stack.len() <= m {
                bail!("unpickling stack underflow");
            }
        }
        self.stack.pop().ok_or_else(|| anyhow!("unpickling stack underflow"))
    }

    fn top(&self) -> Result<&Value> {
        self.stack.last().ok_or_else(|| anyhow!("unpickling stack underflow"))
    }

    fn pop_mark(&mut self) -> Result<Vec<Value>> {
        let m = self.marks.pop().ok_or_else(|| anyhow!("could not find MARK"))?;
        Ok(self.stack.split_off(m))
    }

    fn put_memo(&mut self, idx: usize) -> Result<()> {
        let v = self.top()?.clone();
        self.memo.insert(idx, v);
        Ok(())
    }

    fn get_memo(&mut self, idx: usize) -> Result<()> {
        let v = self
            .memo
            .get(&idx)
            .cloned()
            .ok_or_else(|| anyhow!("memo value {} not found", idx))?;
        self.stack.push(v);
        Ok(())
    }

    fn str_value(b: &[u8]) -> Result<Value> {
        let s = std::str::from_utf8(b).map_err(|_| anyhow!("invalid utf-8 string in pickle"))?;
        Ok(Value::Str(Rc::from(s)))
    }

    fn long_from_le(b: &[u8]) -> Result<Value> {
        if b.is_empty() {
            return Ok(Value::Int(0));
        }
        if b.len() > 8 {
            // accept only values that still fit into an i64
            let neg = b[b.len() - 1] & 0x80 != 0;
            let fill = if neg { 0xff } else { 0x00 };
            if b[8..].iter().any(|&x| x != fill) || ((b[7] & 0x80 != 0) != neg) {
                bail!("integer in pickle does not fit into 64 bits");
            }
            return Self::long_from_le(&b[..8]);
        }
        let neg = b[b.len() - 1] & 0x80 != 0;
        let mut a = if neg { [0xffu8; 8] } else { [0u8; 8] };
        a[..b.len()].copy_from_slice(b);
        Ok(Value::Int(i64::from_le_bytes(a)))
    }

    fn make_object(module: &str, name: &str, args: Vec<Value>) -> Value {
        Value::Object(Rc::new(RefCell::new(Object {
            module: module.to_string(),
            name: name.to_string(),
            args,
            state: None,
            array: None,
        })))
    }

    /// `callable(*args)` for the globals xTea's pickles use.
    fn call(&mut self, func: Value, args: Vec<Value>) -> Result<Value> {
        let (module, name) = match &func {
            Value::Global(g) => (g.0.clone(), g.1.clone()),
            Value::Object(o) => {
                let o = o.borrow();
                (o.module.clone(), o.name.clone())
            }
            _ => bail!("cannot call a {} while unpickling", func.type_name()),
        };
        let numpy_core = module == "numpy.core.multiarray" || module == "numpy._core.multiarray";
        Ok(match (module.as_str(), name.as_str()) {
            ("collections", "OrderedDict") | ("builtins", "dict") | ("__builtin__", "dict") => {
                let mut items = Vec::new();
                if let Some(a) = args.first() {
                    match a {
                        Value::Dict(d) => items = d.borrow().clone(),
                        Value::List(_) | Value::Tuple(_) => {
                            for p in a.as_list()? {
                                let kv = p.as_list()?;
                                if kv.len() == 2 {
                                    items.push((kv[0].clone(), kv[1].clone()));
                                }
                            }
                        }
                        _ => {}
                    }
                }
                Value::Dict(Rc::new(RefCell::new(items)))
            }
            ("builtins", "set") | ("builtins", "frozenset") | ("__builtin__", "set")
            | ("__builtin__", "frozenset") => {
                let items = match args.first() {
                    Some(a) => a.as_list()?,
                    None => Vec::new(),
                };
                Value::Set(Rc::new(RefCell::new(items)))
            }
            ("builtins", "list") | ("__builtin__", "list") => {
                let items = match args.first() {
                    Some(a) => a.as_list()?,
                    None => Vec::new(),
                };
                Value::List(Rc::new(RefCell::new(items)))
            }
            ("builtins", "tuple") | ("__builtin__", "tuple") => {
                let items = match args.first() {
                    Some(a) => a.as_list()?,
                    None => Vec::new(),
                };
                Value::Tuple(Rc::new(items))
            }
            ("builtins", "bytearray") | ("__builtin__", "bytearray") => match args.first() {
                Some(Value::Bytes(b)) => Value::Bytes(b.clone()),
                _ => Value::Bytes(Rc::new(Vec::new())),
            },
            ("_codecs", "encode") => {
                // protocol 2 encodes bytes as _codecs.encode(str, 'latin1')
                let s = args.first().ok_or_else(|| anyhow!("_codecs.encode without data"))?.as_str()?;
                Value::Bytes(Rc::new(s.chars().map(|c| c as u32 as u8).collect()))
            }
            ("copyreg", "_reconstructor") | ("copy_reg", "_reconstructor") => {
                let cls = args.first().cloned().unwrap_or(Value::None);
                match cls {
                    Value::Global(g) => Self::make_object(&g.0, &g.1, Vec::new()),
                    _ => bail!("unsupported copyreg._reconstructor call"),
                }
            }
            (_, "scalar") if numpy_core => {
                // numpy scalar: (dtype, raw bytes)
                if args.len() < 2 {
                    bail!("bad numpy scalar in pickle");
                }
                let (kind, itemsize, little) = dtype_of(&args[0])?;
                let raw = match &args[1] {
                    Value::Bytes(b) => (**b).clone(),
                    Value::Str(s) => s.chars().map(|c| c as u32 as u8).collect(),
                    _ => bail!("bad numpy scalar data in pickle"),
                };
                let arr = NdArray { kind, itemsize, little_endian: little, shape: vec![], fortran: false, data: raw };
                if arr.data.len() < itemsize {
                    bail!("bad numpy scalar data in pickle");
                }
                match kind {
                    'f' => Value::Float(arr.get_f64(0)?),
                    'b' => Value::Bool(arr.get_i64(0)? != 0),
                    _ => Value::Int(arr.get_i64(0)?),
                }
            }
            (_, "_reconstruct") if numpy_core => {
                // ndarray placeholder, filled by the following BUILD
                Self::make_object("numpy", "ndarray", Vec::new())
            }
            _ => Self::make_object(&module, &name, args),
        })
    }

    /// Fill a `numpy.core.multiarray._reconstruct` placeholder from its state
    /// `(version, shape, dtype, is_fortran, rawdata)`.
    fn build_ndarray(state: &Value) -> Result<NdArray> {
        let st = state.as_list()?;
        if st.len() < 5 {
            bail!("unsupported numpy array state in pickle");
        }
        let shape: Vec<usize> = st[1].as_list()?.iter().map(|v| v.as_i64().map(|x| x as usize)).collect::<Result<_>>()?;
        let (kind, itemsize, little) = dtype_of(&st[2])?;
        let fortran = st[3].truthy();
        let data = match &st[4] {
            Value::Bytes(b) => (**b).clone(),
            Value::Str(s) => s.chars().map(|c| c as u32 as u8).collect(),
            _ => bail!("numpy object arrays are not supported"),
        };
        let arr = NdArray { kind, itemsize, little_endian: little, shape, fortran, data };
        if arr.data.len() < arr.len() * itemsize {
            bail!("numpy array data in pickle is truncated");
        }
        Ok(arr)
    }

    /// Read the raw array that follows a joblib `NumpyArrayWrapper` in the stream.
    fn read_joblib_array(&mut self, wrapper: &Value) -> Result<NdArray> {
        let shape: Vec<usize> = wrapper
            .req_attr("shape")?
            .as_list()?
            .iter()
            .map(|v| v.as_i64().map(|x| x as usize))
            .collect::<Result<_>>()?;
        let order = wrapper.getattr("order").map(|v| v.as_str()).transpose()?.unwrap_or_else(|| "C".to_string());
        let (kind, itemsize, little) = dtype_of(&wrapper.req_attr("dtype")?)?;
        if let Some(al) = wrapper.getattr("numpy_array_alignment_bytes") {
            if !matches!(al, Value::None) {
                let pad = self.u8()? as usize;
                if pad != 0 {
                    self.read(pad)?;
                }
            }
        }
        let count: usize = shape.iter().product();
        let data = self.read(count * itemsize)?.to_vec();
        Ok(NdArray { kind, itemsize, little_endian: little, shape, fortran: order == "F", data })
    }

    fn build(&mut self) -> Result<()> {
        let state = self.pop()?;
        let obj = self.top()?.clone();
        match &obj {
            Value::Object(o) => {
                let (module, name) = {
                    let ob = o.borrow();
                    (ob.module.clone(), ob.name.clone())
                };
                if module == "numpy" && name == "ndarray" && o.borrow().array.is_none() && o.borrow().state.is_none() {
                    let arr = Self::build_ndarray(&state)?;
                    let mut ob = o.borrow_mut();
                    ob.array = Some(Rc::new(arr));
                    ob.state = Some(state);
                    return Ok(());
                }
                o.borrow_mut().state = Some(state);
                if module == "joblib.numpy_pickle" && name == "NumpyArrayWrapper" {
                    let arr = Rc::new(self.read_joblib_array(&obj)?);
                    o.borrow_mut().array = Some(arr.clone());
                    // joblib replaces the wrapper by the array on the stack
                    self.stack.pop();
                    self.stack.push(Value::Array(arr));
                } else if module == "joblib.numpy_pickle" && name == "NDArrayWrapper" {
                    bail!("joblib files written by joblib <= 0.9 are not supported");
                }
            }
            Value::Dict(d) => {
                // default __setstate__ for dict subclasses: update the instance dict
                if let Value::Dict(sd) = &state {
                    let items = sd.borrow().clone();
                    d.borrow_mut().extend(items);
                }
            }
            _ => {}
        }
        Ok(())
    }

    fn run(&mut self) -> Result<Value> {
        loop {
            let op = self.u8()?;
            match op {
                0x80 => {
                    self.u8()?; // PROTO
                }
                0x95 => {
                    self.u64()?; // FRAME
                }
                b'.' => return self.pop(),
                b'(' => self.marks.push(self.stack.len()),
                b'0' => {
                    if self.stack.is_empty() || self.marks.last() == Some(&self.stack.len()) {
                        self.pop_mark()?;
                    } else {
                        self.pop()?;
                    }
                }
                b'1' => {
                    self.pop_mark()?;
                }
                b'2' => {
                    let v = self.top()?.clone();
                    self.stack.push(v);
                }
                b'N' => self.stack.push(Value::None),
                0x88 => self.stack.push(Value::Bool(true)),
                0x89 => self.stack.push(Value::Bool(false)),
                b'I' => {
                    let l = self.line()?.trim();
                    let v = match l {
                        "01" => Value::Bool(true),
                        "00" => Value::Bool(false),
                        _ => Value::Int(l.parse().with_context(|| format!("bad INT '{}'", l))?),
                    };
                    self.stack.push(v);
                }
                b'J' => {
                    let v = self.i32()?;
                    self.stack.push(Value::Int(v as i64));
                }
                b'K' => {
                    let v = self.u8()?;
                    self.stack.push(Value::Int(v as i64));
                }
                b'M' => {
                    let v = self.u16()?;
                    self.stack.push(Value::Int(v as i64));
                }
                b'L' => {
                    let l = self.line()?.trim().trim_end_matches('L');
                    self.stack.push(Value::Int(l.parse().with_context(|| format!("bad LONG '{}'", l))?));
                }
                0x8a => {
                    let n = self.u8()? as usize;
                    let b = self.read(n)?;
                    let v = Self::long_from_le(b)?;
                    self.stack.push(v);
                }
                0x8b => {
                    let n = self.i32()?;
                    if n < 0 {
                        bail!("negative LONG4 size");
                    }
                    let b = self.read(n as usize)?;
                    let v = Self::long_from_le(b)?;
                    self.stack.push(v);
                }
                b'F' => {
                    let l = self.line()?.trim();
                    self.stack.push(Value::Float(l.parse().with_context(|| format!("bad FLOAT '{}'", l))?));
                }
                b'G' => {
                    let b = self.read(8)?;
                    let mut a = [0u8; 8];
                    a.copy_from_slice(b);
                    self.stack.push(Value::Float(f64::from_be_bytes(a)));
                }
                b'S' => {
                    let l = self.line()?.trim();
                    let inner = if l.len() >= 2
                        && ((l.starts_with('\'') && l.ends_with('\'')) || (l.starts_with('"') && l.ends_with('"')))
                    {
                        &l[1..l.len() - 1]
                    } else {
                        bail!("the STRING opcode argument must be quoted");
                    };
                    self.stack.push(Value::Str(Rc::from(inner)));
                }
                b'T' => {
                    let n = self.i32()?;
                    if n < 0 {
                        bail!("negative BINSTRING length");
                    }
                    let b = self.read(n as usize)?;
                    let s: String = b.iter().map(|&c| c as char).collect();
                    self.stack.push(Value::Str(Rc::from(s.as_str())));
                }
                b'U' => {
                    let n = self.u8()? as usize;
                    let b = self.read(n)?;
                    let s: String = b.iter().map(|&c| c as char).collect();
                    self.stack.push(Value::Str(Rc::from(s.as_str())));
                }
                b'V' => {
                    let l = self.line()?;
                    self.stack.push(Value::Str(Rc::from(l)));
                }
                b'X' => {
                    let n = self.u32()? as usize;
                    let b = self.read(n)?;
                    self.stack.push(Self::str_value(b)?);
                }
                0x8c => {
                    let n = self.u8()? as usize;
                    let b = self.read(n)?;
                    self.stack.push(Self::str_value(b)?);
                }
                0x8d => {
                    let n = self.u64()? as usize;
                    let b = self.read(n)?;
                    self.stack.push(Self::str_value(b)?);
                }
                b'B' => {
                    let n = self.u32()? as usize;
                    let b = self.read(n)?.to_vec();
                    self.stack.push(Value::Bytes(Rc::new(b)));
                }
                b'C' => {
                    let n = self.u8()? as usize;
                    let b = self.read(n)?.to_vec();
                    self.stack.push(Value::Bytes(Rc::new(b)));
                }
                0x8e | 0x96 => {
                    let n = self.u64()? as usize;
                    let b = self.read(n)?.to_vec();
                    self.stack.push(Value::Bytes(Rc::new(b)));
                }
                b')' => self.stack.push(Value::Tuple(Rc::new(Vec::new()))),
                b't' => {
                    let items = self.pop_mark()?;
                    self.stack.push(Value::Tuple(Rc::new(items)));
                }
                0x85 => {
                    let a = self.pop()?;
                    self.stack.push(Value::Tuple(Rc::new(vec![a])));
                }
                0x86 => {
                    let b = self.pop()?;
                    let a = self.pop()?;
                    self.stack.push(Value::Tuple(Rc::new(vec![a, b])));
                }
                0x87 => {
                    let c = self.pop()?;
                    let b = self.pop()?;
                    let a = self.pop()?;
                    self.stack.push(Value::Tuple(Rc::new(vec![a, b, c])));
                }
                b']' => self.stack.push(Value::List(Rc::new(RefCell::new(Vec::new())))),
                b'l' => {
                    let items = self.pop_mark()?;
                    self.stack.push(Value::List(Rc::new(RefCell::new(items))));
                }
                b'a' => {
                    let v = self.pop()?;
                    match self.top()? {
                        Value::List(l) => l.borrow_mut().push(v),
                        other => bail!("APPEND to a {}", other.type_name()),
                    }
                }
                b'e' => {
                    let items = self.pop_mark()?;
                    match self.top()? {
                        Value::List(l) => l.borrow_mut().extend(items),
                        other => bail!("APPENDS to a {}", other.type_name()),
                    }
                }
                b'}' => self.stack.push(Value::Dict(Rc::new(RefCell::new(Vec::new())))),
                b'd' => {
                    let items = self.pop_mark()?;
                    let pairs = items.chunks(2).filter(|c| c.len() == 2).map(|c| (c[0].clone(), c[1].clone())).collect();
                    self.stack.push(Value::Dict(Rc::new(RefCell::new(pairs))));
                }
                b's' => {
                    let v = self.pop()?;
                    let k = self.pop()?;
                    Self::dict_set(self.top()?, k, v)?;
                }
                b'u' => {
                    let items = self.pop_mark()?;
                    let target = self.top()?.clone();
                    for c in items.chunks(2) {
                        if c.len() == 2 {
                            Self::dict_set(&target, c[0].clone(), c[1].clone())?;
                        }
                    }
                }
                0x8f => self.stack.push(Value::Set(Rc::new(RefCell::new(Vec::new())))),
                0x90 => {
                    let items = self.pop_mark()?;
                    match self.top()? {
                        Value::Set(s) => s.borrow_mut().extend(items),
                        other => bail!("ADDITEMS to a {}", other.type_name()),
                    }
                }
                0x91 => {
                    let items = self.pop_mark()?;
                    self.stack.push(Value::Set(Rc::new(RefCell::new(items))));
                }
                b'c' => {
                    let module = self.line()?.to_string();
                    let name = self.line()?.to_string();
                    self.stack.push(Value::Global(Rc::new((module, name))));
                }
                0x93 => {
                    let name = self.pop()?.as_str()?;
                    let module = self.pop()?.as_str()?;
                    self.stack.push(Value::Global(Rc::new((module, name))));
                }
                b'R' => {
                    let args = self.pop()?.as_list()?;
                    let func = self.pop()?;
                    let v = self.call(func, args)?;
                    self.stack.push(v);
                }
                0x81 => {
                    let args = self.pop()?.as_list()?;
                    let cls = self.pop()?;
                    let v = self.call(cls, args)?;
                    self.stack.push(v);
                }
                0x92 => {
                    let _kwargs = self.pop()?;
                    let args = self.pop()?.as_list()?;
                    let cls = self.pop()?;
                    let v = self.call(cls, args)?;
                    self.stack.push(v);
                }
                b'o' => {
                    let mut items = self.pop_mark()?;
                    if items.is_empty() {
                        bail!("OBJ without a class");
                    }
                    let cls = items.remove(0);
                    let v = self.call(cls, items)?;
                    self.stack.push(v);
                }
                b'i' => {
                    let module = self.line()?.to_string();
                    let name = self.line()?.to_string();
                    let items = self.pop_mark()?;
                    let v = self.call(Value::Global(Rc::new((module, name))), items)?;
                    self.stack.push(v);
                }
                b'b' => self.build()?,
                b'p' => {
                    let idx: usize = self.line()?.trim().parse().map_err(|_| anyhow!("bad PUT index"))?;
                    self.put_memo(idx)?;
                }
                b'q' => {
                    let idx = self.u8()? as usize;
                    self.put_memo(idx)?;
                }
                b'r' => {
                    let idx = self.u32()? as usize;
                    self.put_memo(idx)?;
                }
                0x94 => {
                    let idx = self.memo.len();
                    self.put_memo(idx)?;
                }
                b'g' => {
                    let idx: usize = self.line()?.trim().parse().map_err(|_| anyhow!("bad GET index"))?;
                    self.get_memo(idx)?;
                }
                b'h' => {
                    let idx = self.u8()? as usize;
                    self.get_memo(idx)?;
                }
                b'j' => {
                    let idx = self.u32()? as usize;
                    self.get_memo(idx)?;
                }
                other => bail!("unsupported pickle opcode 0x{:02x} at offset {}", other, self.pos - 1),
            }
        }
    }

    fn dict_set(target: &Value, k: Value, v: Value) -> Result<()> {
        match target {
            Value::Dict(d) => {
                let mut d = d.borrow_mut();
                // keep insertion order, replace existing string/int keys
                let pos = d.iter().position(|(ek, _)| match (ek, &k) {
                    (Value::Str(a), Value::Str(b)) => a == b,
                    (Value::Int(a), Value::Int(b)) => a == b,
                    _ => false,
                });
                match pos {
                    Some(p) => d[p].1 = v,
                    None => d.push((k, v)),
                }
                Ok(())
            }
            other => bail!("SETITEM on a {}", other.type_name()),
        }
    }
}

/// Unpickle an in-memory buffer (joblib arrays included).
pub fn loads(data: &[u8]) -> Result<Value> {
    let mut u = Unpickler { data, pos: 0, stack: Vec::new(), marks: Vec::new(), memo: HashMap::new() };
    u.run()
}

/// Unpickle a file written with `pickle.dump` or (uncompressed) `joblib.dump`.
pub fn load_file(path: &str) -> Result<Value> {
    let data = std::fs::read(path).with_context(|| format!("cannot read {}", path))?;
    if data.len() >= 2 && data[0] == 0x1f && data[1] == 0x8b {
        // joblib with compress=True (zlib/gzip): decompress first
        use std::io::Read;
        let mut d = flate2::read::MultiGzDecoder::new(&data[..]);
        let mut out = Vec::new();
        d.read_to_end(&mut out).with_context(|| format!("cannot decompress {}", path))?;
        return loads(&out).with_context(|| format!("cannot unpickle {}", path));
    }
    if data.len() >= 2 && data[0] == 0x78 {
        use std::io::Read;
        let mut d = flate2::read::ZlibDecoder::new(&data[..]);
        let mut out = Vec::new();
        if d.read_to_end(&mut out).is_ok() {
            return loads(&out).with_context(|| format!("cannot unpickle {}", path));
        }
    }
    loads(&data).with_context(|| format!("cannot unpickle {}", path))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn basic_protocol3() {
        // pickle.dumps({'a': [1, 2.5, 'x'], 'b': (True, None), 'c': -5, 'd': 70000}, protocol=3)
        let data: &[u8] = b"\x80\x03}q\x00(X\x01\x00\x00\x00aq\x01]q\x02(K\x01G@\x04\x00\x00\x00\x00\x00\x00X\x01\x00\x00\x00xq\x03eX\x01\x00\x00\x00bq\x04\x88N\x86q\x05X\x01\x00\x00\x00cq\x06J\xfb\xff\xff\xffX\x01\x00\x00\x00dq\x07Jp\x11\x01\x00u.";
        let v = loads(data).unwrap();
        let a = v.get_item("a").unwrap().as_list().unwrap();
        assert_eq!(a[0].as_i64().unwrap(), 1);
        assert_eq!(a[1].as_f64().unwrap(), 2.5);
        assert_eq!(a[2].as_str().unwrap(), "x");
        let b = v.get_item("b").unwrap().as_list().unwrap();
        assert!(b[0].truthy());
        assert!(matches!(b[1], Value::None));
        assert_eq!(v.get_item("c").unwrap().as_i64().unwrap(), -5);
        assert_eq!(v.get_item("d").unwrap().as_i64().unwrap(), 70000);
    }

    #[test]
    fn basic_protocol4() {
        // the same dict plus 'e': 2**40, protocol=4 (FRAME, MEMOIZE, SHORT_BINUNICODE, LONG1)
        let data: &[u8] = b"\x80\x04\x95B\x00\x00\x00\x00\x00\x00\x00}\x94(\x8c\x01a\x94]\x94(K\x01G@\x04\x00\x00\x00\x00\x00\x00\x8c\x01x\x94e\x8c\x01b\x94\x88N\x86\x94\x8c\x01c\x94J\xfb\xff\xff\xff\x8c\x01d\x94Jp\x11\x01\x00\x8c\x01e\x94\x8a\x06\x00\x00\x00\x00\x00\x01u.";
        let v = loads(data).unwrap();
        assert_eq!(v.get_item("a").unwrap().as_list().unwrap()[2].as_str().unwrap(), "x");
        assert_eq!(v.get_item("d").unwrap().as_i64().unwrap(), 70000);
        assert_eq!(v.get_item("e").unwrap().as_i64().unwrap(), 1i64 << 40);
    }

    #[test]
    fn long1_values() {
        assert!(matches!(Unpickler::long_from_le(&[0xff]).unwrap(), Value::Int(-1)));
        assert!(matches!(Unpickler::long_from_le(&[0x00, 0x80, 0x00]).unwrap(), Value::Int(32768)));
        assert!(matches!(
            Unpickler::long_from_le(&[0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0x7f, 0x00]).unwrap(),
            Value::Int(i64::MAX)
        ));
    }
}
