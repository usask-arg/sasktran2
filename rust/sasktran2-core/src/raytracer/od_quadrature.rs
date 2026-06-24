use super::grid::{GeometryKind, InterpolationMethod};
use super::layer::Layer;

pub fn add_od_quadrature(
    layer: &mut Layer,
    geometry: GeometryKind,
    interpolation: InterpolationMethod,
) {
    let r0 = layer.entrance.position.norm();
    let r1 = layer.exit.position.norm();
    let dr = r1 - r0;
    let distance = layer.effective_distance();
    layer.average_look_away = (layer.exit.position - layer.entrance.position).normalized();

    if interpolation == InterpolationMethod::Lower {
        if r0 < r1 {
            layer.od_quad_start = distance;
            layer.od_quad_end = 0.0;
        } else {
            layer.od_quad_start = 0.0;
            layer.od_quad_end = distance;
        }
        layer.od_quad_start_fraction = 0.5;
        layer.od_quad_end_fraction = 0.5;
        return;
    }

    if dr.abs() < 0.001
        || geometry == GeometryKind::PlaneParallel
        || interpolation == InterpolationMethod::Shell
    {
        layer.od_quad_start = distance / 2.0;
        layer.od_quad_end = distance / 2.0;
        layer.od_quad_start_fraction = 0.5;
        layer.od_quad_end_fraction = 0.5;
        return;
    }

    let average_norm = layer.average_look_away.norm();
    let costheta0 = layer.entrance.position.dot(layer.average_look_away) / (r0 * average_norm);
    let costheta1 = layer.exit.position.dot(layer.average_look_away) / (r1 * average_norm);

    let t0 = r0 * costheta0;
    let t1 = r1 * costheta1;
    let rt = r0 * (1.0 - costheta0 * costheta0).max(0.0).sqrt();

    let (dt1, dt2) = if t1 >= t0 {
        let dt1 = t1 - t0;
        let dt2 = if rt.abs() < 10.0 {
            0.5 * (r1 * t1 - r0 * t0)
        } else {
            0.5 * ((r1 * t1 - r0 * t0) + rt * rt * ((r1 + t1) / (r0 + t0)).ln())
        };
        (dt1, dt2)
    } else {
        let dt1 = t0 - t1;
        let dt2 = if rt.abs() < 10.0 {
            0.5 * (r0 * t0 - r1 * t1)
        } else {
            0.5 * ((r0 * t0 - r1 * t1) + rt * rt * ((r0 + t0) / (r1 + t1)).ln())
        };
        (dt1, dt2)
    };

    layer.od_quad_start = (r1 * dt1 - dt2) / dr * layer.curvature_factor;
    layer.od_quad_end = -(r0 * dt1 - dt2) / dr * layer.curvature_factor;

    let denom = layer.od_quad_start + layer.od_quad_end;
    if denom == 0.0 || !denom.is_finite() {
        layer.od_quad_start_fraction = 0.5;
        layer.od_quad_end_fraction = 0.5;
    } else {
        layer.od_quad_start_fraction = layer.od_quad_start / denom;
        layer.od_quad_end_fraction = layer.od_quad_end / denom;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::raytracer::grid::{
        BoundarySet, CellId, GeometryKind, InterpolationMethod, InterpolationStencil,
    };
    use crate::raytracer::layer::{Layer, TraceEvent, TracePoint};
    use crate::raytracer::vec3::Vec3;

    fn point(position: Vec3) -> TracePoint {
        TracePoint {
            distance: 0.0,
            position,
            altitude: position.z,
            event: TraceEvent::boundary(BoundarySet::new()),
            interpolation: InterpolationStencil::new(),
            cell: Some(CellId::AltitudeLayer(0)),
            on_exact_vertical_boundary: true,
        }
    }

    #[test]
    fn plane_parallel_quadrature_splits_distance_evenly() {
        let mut layer = Layer::new(
            point(Vec3::new(0.0, 0.0, 1.0)),
            point(Vec3::new(0.0, 0.0, 5.0)),
            Some(CellId::AltitudeLayer(0)),
        );

        add_od_quadrature(
            &mut layer,
            GeometryKind::PlaneParallel,
            InterpolationMethod::Linear,
        );

        assert!((layer.od_quad_start - 2.0).abs() < 1e-12);
        assert!((layer.od_quad_end - 2.0).abs() < 1e-12);
    }
}
