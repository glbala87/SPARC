//! Barcode detection and matching module

mod matcher;
mod whitelist;

pub use matcher::{BarcodeCorrector, BarcodeMatcher};
pub use whitelist::Whitelist;

/// Result of barcode matching
#[derive(Debug, Clone)]
pub enum BarcodeMatch {
    /// Exact match found
    Exact(String),
    /// Corrected match found (original, corrected, edit distance)
    Corrected(String, String, u32),
    /// No match found
    NoMatch(String),
}

impl BarcodeMatch {
    /// Get the final barcode (corrected if applicable)
    pub fn barcode(&self) -> Option<&str> {
        match self {
            BarcodeMatch::Exact(bc) => Some(bc),
            BarcodeMatch::Corrected(_, bc, _) => Some(bc),
            BarcodeMatch::NoMatch(_) => None,
        }
    }

    /// Check if a valid match was found
    pub fn is_valid(&self) -> bool {
        matches!(self, BarcodeMatch::Exact(_) | BarcodeMatch::Corrected(..))
    }
}
