//! Process-wide Python-compatible random generator (the original uses the global `random`).

use crate::py::random::{new_rng, PyRandom};
use std::sync::{Mutex, OnceLock};

static RNG: OnceLock<Mutex<PyRandom>> = OnceLock::new();
static SEED: OnceLock<Option<u64>> = OnceLock::new();

/// Fix the seed (call before the first random number is drawn).  `None` keeps the default
/// behaviour: `XTEA_SEED` from the environment if set, otherwise OS entropy (like Python).
pub fn set_seed(seed: Option<u64>) {
    let _ = SEED.set(seed);
}

pub fn with_rng<T>(f: impl FnOnce(&mut PyRandom) -> T) -> T {
    let m = RNG.get_or_init(|| {
        let seed = SEED.get().cloned().flatten();
        Mutex::new(new_rng(seed))
    });
    let mut g = m.lock().unwrap();
    f(&mut g)
}
