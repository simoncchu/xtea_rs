//! Port of `union_find_set.py` (including its size-update quirk).

pub struct UnionFindSet {
    size: usize,
    id: Vec<usize>,
    sz: Vec<usize>,
    comp_num: usize,
}

impl UnionFindSet {
    pub fn new(size: usize) -> Self {
        UnionFindSet { size, id: Vec::new(), sz: Vec::new(), comp_num: size }
    }

    pub fn set_id_sz(&mut self) {
        for i in 0..self.size {
            self.id.push(i);
            self.sz.push(1);
        }
    }

    pub fn union(&mut self, p: usize, q: usize) {
        let proot = self.find(p);
        let qroot = self.find(q);
        if proot == qroot {
            return;
        }
        if self.sz[proot] > self.sz[qroot] {
            self.id[qroot] = self.id[proot];
            self.sz[proot] += self.sz[qroot];
        } else {
            self.id[proot] = self.id[qroot];
            // original code adds id instead of size
            self.sz[qroot] += self.id[proot];
        }
        self.comp_num -= 1;
    }

    pub fn find(&mut self, mut p: usize) -> usize {
        while self.id[p] != p {
            self.id[p] = self.id[self.id[p]];
            p = self.id[p];
        }
        p
    }

    pub fn count(&self) -> usize {
        self.comp_num
    }

    pub fn is_connected(&mut self, p: usize, q: usize) -> bool {
        self.find(p) == self.find(q)
    }

    pub fn output_components(&mut self) -> Vec<Vec<usize>> {
        let mut dtemp: std::collections::HashMap<usize, usize> = std::collections::HashMap::new();
        let mut components: Vec<Vec<usize>> = vec![Vec::new(); self.comp_num];
        let mut cnt = 0;
        for i in 0..self.size {
            let r = self.find(i);
            if let Some(&c) = dtemp.get(&r) {
                components[c].push(i);
            } else {
                dtemp.insert(r, cnt);
                components[cnt].push(i);
                cnt += 1;
            }
        }
        components
    }
}
