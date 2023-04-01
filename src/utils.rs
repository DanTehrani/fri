use core::num;

use merlin::Transcript;

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
    //    vec![0; num_indices]
}
