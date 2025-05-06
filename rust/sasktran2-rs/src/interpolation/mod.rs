pub mod grid1d;
pub mod linear;

#[derive(Debug, Clone, Copy)]
pub enum OutOfBoundsMode {
    Zero,   // Sets to 0 if out of bounds
    Extend, // Extends the last value if out of bounds
}
