// Copyright Supranational LLC
// Licensed under the Apache License, Version 2.0, see LICENSE for details.
// SPDX-License-Identifier: Apache-2.0

use std::{ptr, str::FromStr};

use blst::*;
use wlc_msm_bal::*;

#[test]
fn msm_correctness() {
    let test_npow = std::env::var("TEST_NPOW").unwrap_or("19".to_string());
    let npoints_npow = i32::from_str(&test_npow).unwrap();

    let batches = 4;
    let (points, scalars) =
        util::generate_points_scalars(1usize << npoints_npow, batches);

    let mut context = multi_scalar_mult_init(points.as_slice());
    let msm_results =
        multi_scalar_mult(&mut context, points.as_slice(), scalars.as_slice());

    let scalars = scalars.iter().map(|fr| {
        let mut scalar = blst_scalar::default();
        unsafe { blst_scalar_from_fr(&mut scalar, &fr) };
        scalar
    }).collect::<Vec<_>>();

    let mut scratch = [0u8; unsafe { blst_p1s_mult_pippenger_scratch_sizeof(points.len()) }];
    for b in 0..batches {
        let start = b * points.len();
        let end = (b + 1) * points.len();

        let mut point = blst_p1::default();
        let points = [points.as_ptr(), ptr::null()];
        let scalars = [scalars[start].b.as_ptr(), ptr::null()];
        unsafe {
            blst_p1s_mult_pippenger(
                &mut point,
                points.as_ptr(),
                points.len(),
                scalars,
                255,
                &mut scratch,
            )
        };

        let arkworks_result = <G1Projective as VariableBaseMSM>::msm(
            points.as_slice(),
            &scalars[start..end],
        )
        .unwrap()
        .into_affine();

        assert_eq!(msm_results[b].into_affine(), arkworks_result);
    }
}
