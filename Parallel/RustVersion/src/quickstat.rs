

use rayon::iter::{IndexedParallelIterator, IntoParallelRefMutIterator, ParallelIterator};
use rayon::prelude::*;

const MIN_CHUNK_SIZE: usize = 1000;

/// Do a quickstat on a slice of indices.
pub fn quickstat_index<F: Fn(usize, usize) -> bool>(indices: &mut [usize], goal: usize, lt: F) {
  let mut s = 0;
  let mut e = indices.len();
  while s + 1 < e {
    let pivot = fastrand::usize(s..e);
    indices.swap(s, pivot);
    let mut low = s + 1;
    let mut high = e - 1;
    while low <= high {
        if lt(indices[low], indices[s]) {
            low += 1;
        } else {
            indices.swap(low, high);
            high -= 1;
        }
    }
    indices.swap(s, high);
    if high < goal {
        s = high + 1;
    } else if high > goal {
        e = high;
    } else {
        s = e;
    }
  }
}

/// Do a quickstat on a slice of indices in parallel. Buffer size should match the indices. It
/// is needed because the approach taken uses something like a counting sort and the values
/// can't be ordered in place. It is passed in so we don't have to allocate a new buffer for
/// every call.
pub fn quickstat_index_par<F: Fn(usize, usize) -> bool + Sync>(indices: &mut [usize], buffer: &mut [usize], goal: usize, lt: F) {
  let mut s = 0;
  let mut e = indices.len();
  while s + 1 < e {
    let chunk_size = (e - s) / num_cpus::get();
    let pivot = fastrand::usize(s..e);
    indices.swap(s, pivot);
    let high_index = if chunk_size > MIN_CHUNK_SIZE {
        // Count how many are below/above
        let cnts: Vec<(&[usize], usize)> = indices[s+1..e].par_chunks(chunk_size).map( 
            |inds| {
                (inds, inds.into_iter().fold(0, |cnt, ind| { cnt + if lt(*ind, indices[s]) {1} else {0} } ))
        }).collect();
        // Calculate indices for destination
        let mut cnts_with_buffer = Vec::new();
        let mut buffer_slice = &mut buffer[s..e];
        let mut num_less = 0;
        for i in 0..cnts.len() {
            let (low, rest) = buffer_slice.split_at_mut(cnts[i].1);
            let (middle, high) = rest.split_at_mut(rest.len() - (cnts[i].0.len() - cnts[i].1));
            cnts_with_buffer.push((cnts[i].0, low, high));
            buffer_slice = middle;
            num_less += cnts[i].1;
        }
        // Move values
        cnts_with_buffer.par_iter_mut().for_each(|(inds, low_buffer, high_buffer)| {
            let mut low = 0;
            let mut high = 0;
            for i in 0..inds.len() {
                if lt(inds[i], indices[s]) {
                    low_buffer[low] = inds[i];
                    low += 1;
                } else {
                    high_buffer[high] = inds[i];
                    high += 1;
                }
            }
        });
        buffer[s + num_less] = indices[s];
        indices.par_iter_mut().zip(buffer.par_iter_mut()).for_each(|(i, b)| *i = *b);
        s + num_less
    } else {
        let mut low = s + 1;
        let mut high = e - 1;
        while low <= high {
            if lt(indices[low], indices[s])
            {
                low += 1;
            } else {
                indices.swap(low, high);
                high -= 1;
            }
        }
        indices.swap(s, high);
        high
    };
    
    if high_index < goal {
        s = high_index + 1;
    } else if high_index > goal {
        e = high_index;
    } else {
        s = e;
    }
  }
}

#[cfg(test)]
mod tests {
    use crate::quickstat::{quickstat_index, quickstat_index_par};

    #[test]
    fn small_test() {
        let vals = vec!(2.3, 9.8, 3.1, 1.6, 6.7, 7.8, 8.6);
        let mut indices = vec!(0, 1, 2, 3, 4, 5, 6);
        quickstat_index(&mut indices, 3, |i1, i2| vals[i1] < vals[i2]);
        assert_eq!(indices[3], 4);
    }

    #[test]
    fn random_test() {
        let n = 100000;
        let mut vals = vec!();
        let mut indices = vec!();
        for i in 0..n {
            vals.push(fastrand::f64());
            indices.push(i);
        }
        for _ in 0..200 {
            let goal = fastrand::usize(0..n);
            quickstat_index(&mut indices, goal, |i1, i2| vals[i1] < vals[i2]);
            for i in 0..n {
                if i < goal {
                    assert!(vals[indices[i]] < vals[indices[goal]])
                } else if i >= goal {
                    assert!(vals[indices[i]] >= vals[indices[goal]])
                } 
            }
        }
    }

    #[test]
    fn random_slice_test() {
        let n = 100000;
        let mut vals = vec!();
        let mut indices = vec!();
        for i in 0..n {
            vals.push(fastrand::f64());
            indices.push(i);
        }
        for _ in 0..200 {
            let start = fastrand::usize(0..n / 4);
            let end = start + fastrand::usize(0..3 * n / 4);
            let goal = fastrand::usize(start..end);
            for i in 0..n {
                vals[i] = fastrand::f64();
                indices[i] = i;
            }
            quickstat_index(&mut indices[start..end], goal - start, |i1, i2| vals[i1] < vals[i2]);
            for i in 0..n {
                if i < start {
                    assert_eq!(indices[i], i);
                } else if i < goal {
                    assert!(vals[indices[i]] < vals[indices[goal]])
                } else if i < end {
                    assert!(vals[indices[i]] >= vals[indices[goal]])
                } else {
                    assert_eq!(indices[i], i);
                }
            }
        }
    }

    #[test]
    fn random_test_par() {
        let n = 1000000;
        let mut vals = vec!();
        let mut indices = vec!();
        let mut buffer: Vec<usize> = vec!();
        for i in 0..n {
            vals.push(fastrand::f64());
            indices.push(i);
            buffer.push(0);
        }
        for _ in 0..20 {
            let goal = fastrand::usize(0..n);
            quickstat_index_par(&mut indices, &mut buffer, goal, |i1, i2| vals[i1] < vals[i2]);
            for i in 0..n {
                if i < goal {
                    assert!(vals[indices[i]] < vals[indices[goal]])
                } else if i >= goal {
                    assert!(vals[indices[i]] >= vals[indices[goal]])
                } 
            }
        }
    }

    #[test]
    fn random_slice_test_par() {
        let n = 1000000;
        let mut vals = vec!();
        let mut indices = vec!();
        let mut buffer: Vec<usize> = vec!();
        for i in 0..n {
            vals.push(fastrand::f64());
            indices.push(i);
            buffer.push(0);
        }
        for _ in 0..20 {
            let start = fastrand::usize(0..n / 4);
            let end = start + fastrand::usize(0..3 * n / 4);
            let goal = fastrand::usize(start..end);
            for i in 0..n {
                vals[i] = fastrand::f64();
                indices[i] = i;
            }
            quickstat_index_par(&mut indices[start..end], &mut buffer, goal - start, |i1, i2| vals[i1] < vals[i2]);
            for i in 0..n {
                if i < start {
                    assert_eq!(indices[i], i);
                } else if i < goal {
                    assert!(vals[indices[i]] < vals[indices[goal]])
                } else if i < end {
                    assert!(vals[indices[i]] >= vals[indices[goal]])
                } else {
                    assert_eq!(indices[i], i);
                }
            }
        }
    }

}