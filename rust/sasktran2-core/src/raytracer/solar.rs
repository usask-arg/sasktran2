use super::grid::GeometryKind;
use super::layer::Layer;
use super::vec3::Vec3;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SolarContext {
    pub sun_unit: Vec3,
    pub geometry: GeometryKind,
}

impl SolarContext {
    #[inline(always)]
    pub fn new(sun_unit: Vec3, geometry: GeometryKind) -> Self {
        Self {
            sun_unit: sun_unit.normalized(),
            geometry,
        }
    }
}

pub fn calculate_csz_saz(
    sun_unit: Vec3,
    position: Vec3,
    look_away: Vec3,
    geometry: GeometryKind,
) -> (f64, f64) {
    let local_up = match geometry {
        GeometryKind::Spherical => position.normalized(),
        GeometryKind::PlaneParallel => Vec3::UNIT_Z,
    };

    let csz = local_up.dot(sun_unit);
    let los_projected = (look_away - local_up * look_away.dot(local_up)).normalized();
    let sun_projected = (sun_unit - local_up * sun_unit.dot(local_up)).normalized();
    let y_axis = local_up.cross(sun_projected);
    let saz = y_axis
        .dot(los_projected)
        .atan2(sun_projected.dot(los_projected));

    (csz, saz)
}

pub fn add_solar_parameters(layer: &mut Layer, context: SolarContext) {
    let (cos_sza_entrance, saz_entrance) = calculate_csz_saz(
        context.sun_unit,
        layer.entrance.position,
        layer.average_look_away,
        context.geometry,
    );
    let (cos_sza_exit, saz_exit) = calculate_csz_saz(
        context.sun_unit,
        layer.exit.position,
        layer.average_look_away,
        context.geometry,
    );

    layer.cos_sza_entrance = cos_sza_entrance;
    layer.saz_entrance = saz_entrance;
    layer.cos_sza_exit = cos_sza_exit;
    layer.saz_exit = saz_exit;
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn plane_parallel_solar_zenith_uses_z_axis() {
        let (csz, _) = calculate_csz_saz(
            Vec3::UNIT_Z,
            Vec3::new(100.0, 200.0, 1.0),
            Vec3::UNIT_X,
            GeometryKind::PlaneParallel,
        );
        assert!((csz - 1.0).abs() < 1e-12);
    }
}
