use ff::PrimeField;
use std::cmp::max;

pub trait Hasher<F: PrimeField>: Clone {
    fn new() -> Self;
    fn hash(&self, input: Vec<F>) -> F;
}

pub struct MerkleTree<F: PrimeField, H: Hasher<F>> {
    pub layers: Vec<Vec<F>>, // to root
    hasher: H,
}

#[derive(Debug)]
pub struct MerkleProof<F: PrimeField> {
    pub root: F,
    pub leaf: F,
    pub siblings: Vec<F>,
}

impl<F: PrimeField> MerkleProof<F> {
    pub fn verify(&self, hasher: &impl Hasher<F>) -> bool {
        let mut current_hash = self.leaf;
        for sibling in &self.siblings {
            current_hash = hasher.hash(vec![current_hash, *sibling]);
        }

        current_hash == self.root
    }
}

impl<F: PrimeField, H: Hasher<F>> MerkleTree<F, H> {
    pub fn new(hasher: H) -> Self {
        Self {
            layers: vec![],
            hasher,
        }
    }

    fn compute_layers(&mut self, leaves: Vec<F>) -> Vec<F> {
        self.layers.push(leaves.clone());

        if leaves.len() == 1 {
            return leaves;
        }

        let mut parent_nodes = vec![];
        for i in (0..leaves.len()).step_by(2) {
            if leaves.len() - 1 == i {
                // Sibling.
                let parent = self.hasher.hash(vec![leaves[i], leaves[i]]);
                parent_nodes.push(parent);
            } else {
                let left = leaves[i];
                let right = leaves[i + 1];
                let parent = self.hasher.hash(vec![left, right]);
                parent_nodes.push(parent);
            }
        }
        self.compute_layers(parent_nodes)
    }

    pub fn commit(&mut self, leaves: Vec<F>) -> F {
        let n = leaves.len();
        assert!(n & (n - 1) == 0);

        let mut leaves = {
            if leaves.len() % 2 == 0 {
                leaves
            } else {
                let mut leaves = leaves;
                leaves.push(F::zero());
                leaves
            }
        };

        self.layers.push(leaves.clone());

        while leaves.len() != 1 {
            let mut layer = vec![];
            for i in (0..leaves.len()).step_by(2) {
                let left = leaves[i];
                let right = leaves[i + 1];
                let parent = self.hasher.hash(vec![left, right]);
                layer.push(parent);
            }
            self.layers.push(layer.clone());
            leaves = layer;
        }

        leaves[0]
    }

    pub fn open(&self, leaf: F) -> MerkleProof<F> {
        let mut siblings = vec![];

        let leaf_index = self.layers[0].iter().position(|&x| x == leaf).unwrap();
        let mut sibling_indices = vec![];

        for i in 0..(self.layers.len() - 1) {
            if i == 0 {
                sibling_indices.push(if leaf_index % 2 == 0 {
                    leaf_index + 1
                } else {
                    leaf_index - 1
                });
            } else {
                let current_index = sibling_indices[i - 1] / 2;
                sibling_indices.push(if current_index % 2 == 0 {
                    if (current_index == self.layers[i].len() - 1) {
                        current_index
                    } else {
                        current_index + 1
                    }
                } else {
                    current_index - 1
                });
            }
        }

        for (i, index) in sibling_indices.iter().enumerate() {
            siblings.push(self.layers[i][*index]);
        }

        MerkleProof {
            root: self.layers.last().unwrap()[0],
            leaf,
            siblings,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use pasta_curves::Fp;

    #[derive(Clone)]
    struct PoseidonHasher {}

    impl Hasher<Fp> for PoseidonHasher {
        fn new() -> Self {
            Self {}
        }

        fn hash(&self, inputs: Vec<Fp>) -> Fp {
            inputs[0] + inputs[1]
        }
    }

    #[test]
    fn test_tree() {
        let hasher = PoseidonHasher::new();
        let mut tree = MerkleTree::new(hasher.clone());
        let leaves = vec![
            Fp::from(1),
            Fp::from(2),
            Fp::from(3),
            Fp::from(4),
            Fp::from(5),
            Fp::from(6),
        ];
        let root = tree.commit(leaves.clone());
        for i in 0..leaves.len() {
            let proof = tree.open(leaves[i]);
            assert!(proof.verify(&hasher));
        }
    }
}
