use criterion::{black_box, criterion_group, criterion_main, Criterion};
use fri::{FriProver, Transcript, UniPoly};
use pasta_curves::Fp;

fn criterion_benchmark(c: &mut Criterion) {
    let poly_degree = 2u32.pow(13u32);

    let mut coeffs = vec![];
    for i in 0..(poly_degree + 1) {
        coeffs.push(Fp::from(i as u64));
    }

    let poly = &UniPoly::new(coeffs);
    let prover = FriProver::<Fp>::new(poly.degree());

    let transcript = &mut Transcript::new(b"bench_fri");

    c.bench_function("bench_fri", |b| {
        b.iter(|| prover.prove_degree(black_box(poly), black_box(transcript)))
    });
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
