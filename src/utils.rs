use merlin::Transcript;
use pasta_curves::arithmetic::FieldExt;

pub fn hash_two<F>(values: &[F; 2]) -> F
where
    F: FieldExt<Repr = [u8; 32]>,
{
    let mut bytes = vec![];
    bytes.extend_from_slice(&values[0].to_repr());
    bytes.resize(32, 0);
    bytes.extend_from_slice(&values[1].to_repr());
    bytes.resize(32, 0);

    let mut bytes_8 = vec![];
    for i in 0..(bytes.len() / 8) {
        let mut acc: u64 = 0;
        for j in 0..8 {
            acc += (bytes[8 * i + j] as u64) << j;
        }
        bytes_8.push(acc);
    }

    bytes_8.resize(25, 0);

    println!("bytes_8: {:?}", bytes_8.len());

    keccak::f1600(&mut bytes_8.as_slice().try_into().unwrap());

    let mut bytes = bytes_8[0..16]
        .iter()
        .flat_map(|x| x.to_le_bytes().to_vec())
        .collect::<Vec<u8>>();
    bytes.reverse();

    F::from_bytes_wide(&bytes[0..64].try_into().unwrap())
}

fn sample_index(random_bytes: [u8; 64], size: usize) -> usize {
    let mut acc: u64 = 0;
    for b in random_bytes {
        acc = acc << 8 ^ (b as u64);
    }

    (acc % (size as u64)) as usize
}

pub fn sample_indices(
    num_indices: usize,
    max_index: usize,
    reduced_max_index: usize,
    transcript: &mut Transcript,
) -> Vec<usize> {
    assert!(num_indices <= 2 * reduced_max_index, "not enough entropy!");
    assert!(num_indices <= reduced_max_index);

    let mut indices = vec![];

    let mut reduced_indices = vec![];

    let mut counter = 0;

    while indices.len() < num_indices {
        let mut random_bytes = [0u8; 64];
        transcript.append_u64(b"counter", counter);
        transcript.challenge_bytes(b"index", &mut random_bytes);
        let index = sample_index(random_bytes, max_index);
        let reduced_index = index % reduced_max_index;

        counter += 1;
        if !reduced_indices.contains(&reduced_index) {
            reduced_indices.push(reduced_index);
            indices.push(index);
        };
    }

    indices
}
