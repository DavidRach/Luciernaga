use extendr_api::prelude::*;
use nnls::nnls;

/// Verify the matrix type works
/// @export
/// @noRd
#[extendr]
fn describe_matrix(matrix: ArrayView2<f64>){
    println!("This R matrix has shape {:?}", matrix.dim())
}

/// Verify the vector type works
/// @export
/// @noRd
#[extendr]
fn describe_vector(vector: ArrayView1<f64>){
    println!("This R vector has length {:?}", vector.len())
}

/// Perform nnls in Rust, return object to R. 
/// @export
/// @noRd
#[extendr]
fn nnls_rust(matrix: ArrayView2<f64>, vector: ArrayView2<f64>) -> Robj {

    let mut results = Vec::new();

    for row in vector.axis_iter(Axis(0)) {
        let result = nnls(matrix, row);
        let (coefficients, _) = result;
        let coefficients_array: Array1<f64> = coefficients.to_owned();
        results.push(coefficients_array);
    }

    let results_array = Array2::from_shape_vec(
        (results.len(), results[0].len()), results.into_iter().flatten().collect()).unwrap();

    results_array.try_into().unwrap()
}

/// Perform nnls in Rust, return object to R. 
/// @export
/// @noRd
#[extendr]
fn nnls_rust2(matrix: ArrayView2<f64>, vector: ArrayView1<f64>) -> Robj {

    let result = nnls(matrix, vector);

    let (coefficients, _) = result;

    let coefficients_array: Array1<f64> = coefficients.to_owned();

    coefficients_array.try_into().unwrap()
}

// Macro to generate exports.
// This ensures exported functions are registered with R.
// See corresponding C code in `entrypoint.c`.
extendr_module! {
    mod Luciernaga;
    fn describe_matrix;
    fn describe_vector;
    fn nnls_rust;
    fn nnls_rust2;
}
