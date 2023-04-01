use std::marker::PhantomData;

use crate::fri_prover::FRIProof;
use crate::tree::Hasher;
use crate::unipoly::UniPoly;
use crate::utils::sample_indices;
use merlin::Transcript;
use pasta_curves::arithmetic::FieldExt;
use pasta_curves::group::ff::PrimeField;

pub struct FriVerifier<F: PrimeField<Repr = [u8; 32]> + FieldExt, H: Hasher<F>> {
    domain: Vec<F>,
    expansion_factor: usize, // (i.e. expansion factor) (info bits) / (total bits)
    domain_length: usize,
    num_colinearity_checks: usize,
    _marker: PhantomData<H>,
}

impl<F: PrimeField<Repr = [u8; 32]> + FieldExt, H: Hasher<F>> FriVerifier<F, H> {
    pub fn construct(
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
            _marker: PhantomData,
        }
    }

    pub fn verify(&self, proof: FRIProof<F>, hasher: H, com: F) {
        let mut transcript =
            Transcript::new(b"Fast Reed-Solomon Interactive Oracle Proof of Proximity");

        let final_codeword = proof.reduced_codeword;

        // TODO: Precompute this
        let mut domain_reduced = self.domain.clone();
        for _ in 0..proof.queries.len() {
            let mut domain_unique = vec![];
            domain_reduced.iter().map(|x| x.square()).for_each(|x| {
                if !domain_unique.contains(&x) {
                    domain_unique.push(x);
                }
            });
            domain_reduced = domain_unique;
        }

        let interpolant = UniPoly::interpolate(domain_reduced.clone(), final_codeword.clone());

        /*
        let eval1 = interpolant.eval(domain_reduced[0]);
        let eval2 = interpolant.eval(domain_reduced[1]);
        println!("eval1 {:?}", eval1);
        println!("eval1 {:?}", eval2);
        println!("exected evals {:?}", final_codeword);
         */

        let degree = if final_codeword.len() == 1 {
            0
        } else {
            final_codeword.len() / self.expansion_factor
        };

        println!("degree {}", degree);
        println!("interpolant.degree() {}", interpolant.degree());
        assert!(interpolant.degree() == degree);

        let mut indices = sample_indices(
            self.num_colinearity_checks,
            self.domain_length,
            self.domain_length / 2usize.pow(proof.queries.len() as u32),
            &mut transcript,
        );

        for (i, layer) in proof.queries.iter().enumerate() {
            assert!(
                layer.openings.len() == self.num_colinearity_checks,
                "Invalid number of colinearity checks"
            );

            // Halve the range of the indices
            indices = indices
                .iter()
                .map(|index| index % (self.domain_length / 2 >> (i + 1)))
                .collect::<Vec<usize>>();

            let a_indices = indices.clone();
            let b_indices = indices
                .iter()
                .map(|index| (self.domain_length / 2 >> (i + 1)) + index)
                .collect::<Vec<usize>>();
            let c_indices = indices.clone();

            // Colinearity checks
            for (j, (a, b, c)) in layer.openings.iter().enumerate() {
                let a_y = a.leaf;
                let b_y = b.leaf;
                let c_y = c.leaf;

                let a_x = self.domain[a_indices[j]];
                let b_x = self.domain[b_indices[j]];
                let c_x = self.domain[c_indices[j]];

                // Check that (a_x, a_y), (b_x, b_y) and (c_x, c_y) are on a straight line.
                let coeff = (a_y - b_y) * (a_x - b_x).invert().unwrap();
                let intercept = a_y - coeff * a_x;
                assert!(c_y == coeff * c_x + intercept);

                // Check Merkle proofs
                a.verify(&hasher);
                b.verify(&hasher);
                c.verify(&hasher);

                // Check that the root is correct
                assert_eq!(
                    a.root, b.root,
                    "Roots of the two Merkle proofs are not equal"
                );

                if i == 0 && j == 0 {
                    assert_eq!(c.root, com, "c.root != com");
                } else if j > 0 {
                    assert_eq!(
                        layer.openings[j - 1].0.root,
                        c.root,
                        "a_prev.root != c.root"
                    );

                    assert_eq!(
                        layer.openings[j - 1].1.root,
                        c.root,
                        "b_prev.root != c.root"
                    );
                }
            }
        }
    }
}
