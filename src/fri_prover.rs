use std::marker::PhantomData;

use crate::tree::Hasher;
use crate::tree::{MerkleProof, MerkleTree};
use crate::utils::sample_indices;
use merlin::Transcript;
use pasta_curves::{arithmetic::FieldExt, group::ff::PrimeField};

pub struct ProofStream<F: PrimeField<Repr = [u8; 32]> + FieldExt, H: Hasher<F>> {
    trees: Vec<MerkleTree<F, H>>,
    codewords: Vec<Vec<F>>,
}

pub struct FRIProver<F: PrimeField<Repr = [u8; 32]> + FieldExt, H: Hasher<F>> {
    domain: Vec<F>,
    domain_length: usize,
    expansion_factor: usize, // (i.e. expansion factor) (info bits) / (total bits)
    num_colinearity_checks: usize,
    omega: F,
    _marker: PhantomData<H>,
}

#[derive(Debug)]
pub struct FRIProof<F: PrimeField<Repr = [u8; 32]> + FieldExt> {
    pub reduced_codeword: Vec<F>,
    pub queries: Vec<LayerProof<F>>,
}

#[derive(Debug)]
pub struct LayerProof<F: PrimeField<Repr = [u8; 32]> + FieldExt> {
    pub openings: Vec<(MerkleProof<F>, MerkleProof<F>, MerkleProof<F>)>,
}

impl<F: PrimeField<Repr = [u8; 32]> + FieldExt, H: Hasher<F>> FRIProver<F, H> {
    fn construct(
        omega: F,
        domain_length: usize,
        expansion_factor: usize,
        num_colinearity_checks: usize,
    ) -> Self {
        let mut domain = vec![];
        for i in 0..domain_length {
            domain.push(omega.pow(&[i as u64, 0, 0, 0]));
        }

        Self {
            domain,
            domain_length,
            expansion_factor,
            num_colinearity_checks,
            omega,
            _marker: PhantomData,
        }
    }

    fn num_rounds(&self) -> usize {
        // TODO: Enable num_colinearity_check and rounds configuration
        /*
        let mut colinearity_checks_count = 0;
        let mut n_rounds = 0;
        while (4 * colinearity_checks_count) < self.domain_length {
            colinearity_checks_count += self.num_colinearity_checks;
            n_rounds += 1;
        }
         */
        // TODO Choose a secure num_rounds!
        ((self.domain_length as f64).log2() as usize) - 3 // this `3` is just random
    }

    fn fold(&self, codeword: Vec<F>, domain: Vec<F>, alpha: F) -> Vec<F> {
        assert!(codeword.len() == domain.len());
        let two_inv = F::from(2).invert().unwrap();
        let one = F::from(1);

        let n = domain.len();

        let mut folded_codeword = vec![];
        for i in 0..(n / 2) {
            // f*(w^2i) = 1/2 * ((1 + alpha * w^-i) * f(w^i) + (1 - alpha * w^-i) * f(-w^i))
            // w^(n/2) = -1
            // -w^i = domain[i + n/2]

            //  let omega_pow_minus_i = self.omega.pow(&[i as u64, 0, 0, 0]).invert().unwrap();
            let omega_pow_minus_i = domain[n - 1 - i];

            let f_star_eval = two_inv
                * ((one + alpha * omega_pow_minus_i) * codeword[i]
                    + (one - alpha * omega_pow_minus_i) * codeword[i + (n / 2)]);
            folded_codeword.push(f_star_eval);
        }

        folded_codeword
    }

    fn commit(
        &self,
        codeword: Vec<F>,
        transcript: &mut Transcript,
    ) -> (Vec<Vec<F>>, Vec<MerkleTree<F, H>>) {
        // Commit to the codeword
        let hasher = H::new();

        let mut domain = self.domain.clone();

        let mut codewords = vec![codeword.clone()];
        let mut trees = vec![];

        for i in 0..self.num_rounds() {
            let current_codeword = codewords[i].clone();

            let mut tree = MerkleTree::new(hasher.clone());
            let root = tree.commit(current_codeword.clone());

            transcript.append_message(b"root", &root.to_repr());
            trees.push(tree);

            let mut alpha = [0u8; 64];
            transcript.challenge_bytes(b"alpha", &mut alpha);
            let alpha = F::from_bytes_wide(&alpha);

            let next_codeword = self.fold(current_codeword, domain.to_vec(), alpha);
            let mut domain_unique = vec![];
            domain.iter().map(|x| x.square()).for_each(|x| {
                if !domain_unique.contains(&x) {
                    domain_unique.push(x);
                }
            });
            domain = domain_unique;

            codewords.push(next_codeword)
        }

        (codewords, trees)
    }

    fn query(
        &self,
        codewords: Vec<Vec<F>>,
        trees: Vec<MerkleTree<F, H>>,
        indices: Vec<usize>,
    ) -> Vec<LayerProof<F>> {
        // A domain: w^i
        // B domain: w^{n/2 + i}
        // C domain: w^{2i}

        assert!(indices.len() == self.num_colinearity_checks);
        let mut indices = indices;

        let mut queries = vec![];

        for (i, codeword) in codewords.iter().enumerate() {
            // Skip the last codeword since the verifier's gonna check it directly
            if i == codewords.len() - 1 {
                continue;
            }

            // Halve the range of the indices
            indices = indices
                .iter()
                .map(|index| {
                    if codeword.len() == 1 {
                        0
                    } else {
                        index % (codeword.len() / 2)
                    }
                })
                .collect::<Vec<usize>>();

            let a_indices = indices.clone();
            let b_indices = indices
                .iter()
                .map(|index| (codeword.len() / 2) + index)
                .collect::<Vec<usize>>();
            let c_indices = indices.clone();

            let mut openings = vec![];
            for j in 0..self.num_colinearity_checks {
                let a_y = codeword[a_indices[j]];
                let a_y_proof = trees[i].open(a_y);
                let b_y = codeword[b_indices[j]];
                let b_y_proof = trees[i].open(b_y);

                let c_y = codeword[c_indices[j]];
                let c_y_proof = trees[i].open(c_y);

                openings.push((a_y_proof, b_y_proof, c_y_proof));
            }

            queries.push(LayerProof { openings })
        }

        queries
    }

    fn prove(&self, codeword: Vec<F>) -> FRIProof<F> {
        let mut transcript =
            Transcript::new(b"Fast Reed-Solomon Interactive Oracle Proof of Proximity");

        let (codewords, trees) = self.commit(codeword, &mut transcript);

        let indices = sample_indices(
            self.num_colinearity_checks,
            codewords[0].len(),                   // Length of the initial codeword
            codewords[codewords.len() - 2].len(), // Length of the reduced codeword
            &mut transcript,
        );

        let queries = self.query(codewords.clone(), trees, indices);

        FRIProof {
            reduced_codeword: codewords[codewords.len() - 1].clone(),
            queries,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fft::fft;
    use crate::fri_verifier::FriVerifier;
    use crate::unipoly::UniPoly;
    use pasta_curves::group::ff::PrimeField;
    use pasta_curves::Fp;

    #[derive(Clone)]
    struct PoseidonHasher {}

    impl Hasher<Fp> for PoseidonHasher {
        fn new() -> Self {
            Self {}
        }

        fn hash(&self, inputs: Vec<Fp>) -> Fp {
            inputs[0] + inputs[1] + Fp::one()
        }
    }

    #[test]
    fn test_prove() {
        let poly_degree = 1024;

        let mut coeffs = vec![];
        for i in 0..(poly_degree + 1) {
            coeffs.push(Fp::from(i as u64));
        }

        let poly = UniPoly::new(coeffs);

        let root_of_unity = Fp::root_of_unity();

        let expansion_factor = 2;
        let num_colinearity_checks = 2;

        let subgroup_order = (poly.degree() * expansion_factor).next_power_of_two();

        // Generator for the subgroup with order _subgroup_order_ in the field
        let omega = root_of_unity.pow(&[
            2u32.pow(32 - ((subgroup_order as f64).log2() as u32)) as u64,
            0,
            0,
            0,
        ]);

        let mut domain = vec![];
        for i in 0..subgroup_order {
            domain.push(omega.pow(&[i as u64, 0, 0, 0]));
        }

        let mut coeffs_expanded = poly.coeffs.clone();
        coeffs_expanded.resize(domain.len(), Fp::zero());

        let evals = fft(coeffs_expanded, domain.clone());

        let prover = FRIProver::<Fp, PoseidonHasher>::construct(
            omega,
            subgroup_order,
            expansion_factor,
            num_colinearity_checks,
        );

        println!("d {:?}", poly.degree());
        println!("domain_len {}", subgroup_order);
        println!("rounds {}", prover.num_rounds());

        let proof = prover.prove(evals);

        // C of the first round is the polynomial we're committing to.
        let poly_commitment = proof.queries[0].openings[0].2.root;

        let verifier = FriVerifier::<Fp, PoseidonHasher>::construct(
            omega,
            subgroup_order,
            expansion_factor,
            num_colinearity_checks,
        );

        let hasher = PoseidonHasher::new();
        verifier.verify(proof, hasher, poly_commitment);
    }
}
