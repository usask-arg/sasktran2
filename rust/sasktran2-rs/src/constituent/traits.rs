use crate::atmosphere::*;
use crate::prelude::*;

pub trait Constituent {
    fn add_to_atmosphere(&self, storage: &mut impl AtmosphereStorageAccess) -> Result<()>;

    fn register_derivatives(
        &self,
        storage: &mut impl AtmosphereStorageAccess,
        constituent_name: &str,
    ) -> Result<()>;
}
