use anyhow::Result;
use cxx::CxxVector;
use gauss_quad::legendre::GaussLegendre;
use std::pin::Pin;

#[cxx::bridge(namespace = "sasktran2::rust::math")]
pub mod ffi {
    extern "Rust" {
        fn gauss_quad_nodes_weights(
            n: usize,
            nodes: Pin<&mut CxxVector<f64>>,
            weights: Pin<&mut CxxVector<f64>>,
        ) -> Result<()>;
    }
}

fn gauss_quad_nodes_weights(
    n: usize,
    mut nodes: Pin<&mut CxxVector<f64>>,
    mut weights: Pin<&mut CxxVector<f64>>,
) -> Result<()> {
    let quad = GaussLegendre::new(n)?;

    for (node, w) in quad.iter() {
        nodes.as_mut().push(*node);
        weights.as_mut().push(*w);
    }

    Ok(())
}
