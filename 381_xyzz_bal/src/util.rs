// Copyright Supranational LLC
// Licensed under the Apache License, Version 2.0, see LICENSE for details.
// SPDX-License-Identifier: Apache-2.0

use std::ops::AddAssign;

use blst::{blst_fr, blst_fr_from_uint64, blst_p1, blst_p1_affine, blst_p1_generator, blst_p1_mult, blst_p1_to_affine, blst_scalar, blst_scalar_from_fr};
use rand::SeedableRng;
use rand_chacha::ChaCha20Rng;

use ark_ec::{AffineRepr, CurveGroup, Group};
use ark_std::UniformRand;

use ark_ff::prelude::*;
use ark_std::vec::Vec;

pub fn generate_points_scalars(
    len: usize,
    batch_size: usize,
) -> (Vec<blst_p1_affine>, Vec<blst_fr>) {
    let rand_gen: usize = 1 << 11;
    let mut rng = ChaCha20Rng::from_entropy();

    let mut points = (0..rand_gen)
        .map(|_| {
            let val: [u64; 4] = [
                rand::random(),
                rand::random(),
                rand::random(),
                rand::random(),
            ];
            let mut ret = blst_fr::default();
            unsafe { blst_fr_from_uint64(&mut ret, val.as_ptr()) };
            let mut scalar = blst_scalar::default();
            unsafe { blst_scalar_from_fr(&mut scalar, &ret) };

            let mut point = unsafe { *blst_p1_generator() };
            unsafe { blst_p1_mult(&mut point, &point, scalar.b.as_ptr(), 255) };
            
            let mut affine = blst_p1_affine::default();
            unsafe { blst_p1_to_affine(&mut affine, &point) };
            affine
        })
        .collect::<Vec<_>>();

    // Sprinkle in some infinity points
    //    points[3] = G::zero();
    while points.len() < len {
        points.append(&mut points.clone());
    }

    let points = &points[0..len];

    let scalars = (0..len * batch_size)
        .map(|_| {
            let val: [u64; 4] = [
                rand::random(),
                rand::random(),
                rand::random(),
                rand::random(),
            ];
            let mut ret = blst_fr::default();
            unsafe { blst_fr_from_uint64(&mut ret, val.as_ptr()) };
            ret
        })
        .collect::<Vec<_>>();

    (points, scalars)
}

#[cfg(feature = "parallel")]
use rayon::prelude::*;

pub struct VariableBaseMSM2;

impl VariableBaseMSM2 {
    pub fn multi_scalar_mul<G: AffineRepr>(
        bases: &[G],
        scalars: &[<G::ScalarField as PrimeField>::BigInt],
    ) -> G::Group {
        let size = ark_std::cmp::min(bases.len(), scalars.len());
        let scalars = &scalars[..size];
        let bases = &bases[..size];
        let scalars_and_bases_iter =
            scalars.iter().zip(bases).filter(|(s, _)| !s.is_zero());

        // let c = if size < 32 {
        //     3
        // } else {
        //     super::ln_without_floats(size) + 2
        // };
        let c = 21;

        let num_bits =
            <G::ScalarField as PrimeField>::MODULUS_BIT_SIZE as usize;
        let fr_one = G::ScalarField::one().into();

        let zero = G::Group::zero();
        let window_starts: Vec<_> = (0..num_bits).step_by(c).collect();

        // Each window is of size `c`.
        // We divide up the bits 0..num_bits into windows of size `c`, and
        // in parallel process each such window.
        let window_sums: Vec<_> = ark_std::cfg_into_iter!(window_starts)
            .map(|w_start| {
                let mut res = zero;
                // let mut count = 0;

                // We don't need the "zero" bucket, so we only have 2^c - 1 buckets.
                let mut buckets = vec![zero; (1 << c) - 1];
                // This clone is cheap, because the iterator contains just a
                // pointer and an index into the original vectors.
                scalars_and_bases_iter.clone().for_each(|(&scalar, base)| {
                    if scalar == fr_one {
                        // We only process unit scalars once in the first window.
                        if w_start == 0 {
                            res.add_assign(base);
                        }
                    } else {
                        let mut scalar = scalar;

                        // We right-shift by w_start, thus getting rid of the
                        // lower bits.
                        scalar.divn(w_start as u32);

                        // We mod the remaining bits by 2^{window size}, thus taking `c` bits.
                        let scalar = scalar.as_ref()[0] % (1 << c);

                        // If the scalar is non-zero, we update the corresponding
                        // bucket.
                        // (Recall that `buckets` doesn't have a zero bucket.)
                        // if w_start != 252{
                        if scalar != 0 {
                            buckets[(scalar - 1) as usize].add_assign(base);
                        }
                        // }
                    }
                });

                // Compute sum_{i in 0..num_buckets} (sum_{j in i..num_buckets} bucket[j])
                // This is computed below for b buckets, using 2b curve additions.
                //
                // We could first normalize `buckets` and then use mixed-addition
                // here, but that's slower for the kinds of groups we care about
                // (Short Weierstrass curves and Twisted Edwards curves).
                // In the case of Short Weierstrass curves,
                // mixed addition saves ~4 field multiplications per addition.
                // However normalization (with the inversion batched) takes ~6
                // field multiplications per element,
                // hence batch normalization is a slowdown.

                // `running_sum` = sum_{j in i..num_buckets} bucket[j],
                // where we iterate backward from i = num_buckets to 0.
                let mut running_sum = G::Group::zero();
                buckets.into_iter().rev().for_each(|b| {
                    running_sum += &b;
                    res += &running_sum;
                });
                res
            })
            .collect();

        // We store the sum for the lowest window.
        let lowest = *window_sums.first().unwrap();

        // let mm = window_sums[1..]
        // .iter()
        // .rev()
        // .fold(zero, |mut total, sum_i| {
        //     total += sum_i;
        //     // for _ in 0..c {
        //     //     total.double_in_place();
        //     // }
        //     total
        // });

        // We're traversing windows from high to low.
        lowest
            + &window_sums[1..]
                .iter()
                .rev()
                .fold(zero, |mut total, sum_i| {
                    total += sum_i;
                    for _ in 0..c {
                        total.double_in_place();
                    }
                    total
                })

        // if mm == zero{
        //     panic!("zero!");
        // }

        // mm
    }
}
