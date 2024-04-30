// Copyright Supranational LLC
// Licensed under the Apache License, Version 2.0, see LICENSE for details.
// SPDX-License-Identifier: Apache-2.0

use std::os::raw::c_void;
use ark_bls12_381::{Fr, G1Affine};
use ark_ec::AffineRepr;
use ark_ff::PrimeField;
use ark_std::Zero;

#[allow(unused_imports)]
use blst::*;

sppark::cuda_error!();

pub mod util;

#[repr(C)]
pub struct MultiScalarMultContext {
    context: *mut c_void,
}

#[cfg_attr(feature = "quiet", allow(improper_ctypes))]
extern "C" {
    fn mult_pippenger_faster_init(
        context: *mut MultiScalarMultContext,
        points_with_infinity: *const blst_p1_affine,
        npoints: usize,
    ) -> cuda::Error;

    // fn mult_pippenger_init(
    //     context: *mut MultiScalarMultContext,
    //     points_with_infinity: *const G1Affine,
    //     npoints: usize,
    //     ffi_affine_sz: usize,
    // ) -> cuda::Error;

    // fn mult_pippenger_faster_inf2();

    fn mult_pippenger_faster_inf(
        context: *mut MultiScalarMultContext,
        out: *mut u64,
        points_with_infinity: *const blst_p1_affine,
        npoints: usize,
        batch_size: usize,
        scalars: *const blst_fr,
    ) -> cuda::Error;
    
    // fn mult_pippenger_inf(
    //     context: *mut MultiScalarMultContext,
    //     out: *mut u64,
    //     points_with_infinity: *const G1Affine,
    //     npoints: usize,
    //     batch_size: usize,
    //     scalars: *const Fr,
    //     ffi_affine_sz: usize,
    // ) -> cuda::Error;
}

pub fn multi_scalar_mult_init(
    points: &[blst_p1_affine],
) -> MultiScalarMultContext {
    let mut ret = MultiScalarMultContext {
        context: std::ptr::null_mut(),
    };
        
    let err = unsafe {
        mult_pippenger_faster_init(
            &mut ret,
            points as *const _ as *const blst_p1_affine,
            points.len(),
        )
    };
    if err.code != 0 {
        panic!("{}", String::from(err));
    }

    ret
}
    
pub fn multi_scalar_mult(
    context: &mut MultiScalarMultContext,
    points: &[blst_p1_affine],
    scalars: &[blst_fr],
) -> Vec<blst_p1> {
    let npoints = points.len();
    if scalars.len() % npoints != 0 {
        panic!("length mismatch")
    }

    //let mut context = multi_scalar_mult_init(points);

    let batch_size = scalars.len() / npoints;
    let mut ret = vec![blst_p1::default(); batch_size];
    let err = unsafe {
        let result_ptr = 
            &mut *(&mut ret as *mut Vec<blst_p1>
                   as *mut Vec<u64>);
        
        // mult_pippenger_faster_inf2();

        mult_pippenger_faster_inf(
            context,
            result_ptr.as_mut_ptr(),
            points as *const _ as *const blst_p1_affine,
            npoints, batch_size,
            scalars as *const _ as *const blst_fr,
        )
    };
    if err.code != 0 {
        panic!("{}", String::from(err));
    }

    ret
}
