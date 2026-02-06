use crate::optical::storage::OpticalQuantities;
use rebasis::grid::MappingMatrix;

pub fn reduce_optical(
    optical_quantities: &OpticalQuantities,
    mapping: &MappingMatrix,
) -> OpticalQuantities {
    OpticalQuantities {
        cross_section: mapping.dot(optical_quantities.cross_section.view()),
        ssa: mapping.dot(optical_quantities.ssa.view()),
        legendre: None, // todo: reduce legendre
        fortran_ordering: optical_quantities.fortran_ordering,
    }
}
