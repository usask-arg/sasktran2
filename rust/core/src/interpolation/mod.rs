pub mod linear;

pub enum OutOfBoundsMode {
    Zero,   // Sets to 0 if out of bounds
    Extend  // Extends the last value if out of bounds
}
