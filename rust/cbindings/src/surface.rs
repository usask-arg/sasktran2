use crate::ffi;

pub struct Surface {
    pub surface: *mut ffi::Surface,
}

impl Surface {
    pub fn new(num_wavel: usize, num_stokes: usize) -> Self {
        Surface {
            surface: unsafe { ffi::sk_surface_create(num_wavel as i32, num_stokes as i32) },
        }
    }
}

impl Drop for Surface {
    fn drop(&mut self) {
        unsafe {
            ffi::sk_surface_destroy(self.surface);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_surface() {
        let num_wavel = 10;
        let num_stokes = 3;
        let mut surface = Surface::new(num_wavel, num_stokes);
    }
}
