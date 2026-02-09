//! Barcode matching and correction

use super::{BarcodeMatch, Whitelist};
use ahash::AHashMap;

/// Barcode matcher with exact matching
pub struct BarcodeMatcher {
    whitelist: Whitelist,
}

impl BarcodeMatcher {
    pub fn new(whitelist: Whitelist) -> Self {
        Self { whitelist }
    }

    /// Match a barcode exactly
    pub fn match_exact(&self, barcode: &str) -> BarcodeMatch {
        if self.whitelist.contains(barcode) {
            BarcodeMatch::Exact(barcode.to_string())
        } else {
            BarcodeMatch::NoMatch(barcode.to_string())
        }
    }

    /// Get the whitelist
    pub fn whitelist(&self) -> &Whitelist {
        &self.whitelist
    }
}

/// Barcode corrector with fuzzy matching using Hamming distance
pub struct BarcodeCorrector {
    whitelist: Whitelist,
    /// Maximum Hamming distance for correction
    max_distance: u32,
    /// Pre-computed index for 1-mismatch lookup
    mismatch_index: AHashMap<String, Vec<String>>,
}

impl BarcodeCorrector {
    /// Create a new barcode corrector
    pub fn new(whitelist: Whitelist, max_distance: u32) -> Self {
        let mismatch_index = if max_distance >= 1 {
            Self::build_mismatch_index(&whitelist)
        } else {
            AHashMap::new()
        };

        Self {
            whitelist,
            max_distance,
            mismatch_index,
        }
    }

    /// Build index for 1-mismatch lookup
    fn build_mismatch_index(whitelist: &Whitelist) -> AHashMap<String, Vec<String>> {
        let mut index: AHashMap<String, Vec<String>> = AHashMap::new();

        for barcode in whitelist.iter() {
            // Generate all 1-mismatch variants
            let chars: Vec<char> = barcode.chars().collect();
            for i in 0..chars.len() {
                for &base in &['A', 'C', 'G', 'T', 'N'] {
                    if chars[i] != base {
                        let mut variant = chars.clone();
                        variant[i] = base;
                        let variant_str: String = variant.into_iter().collect();
                        index
                            .entry(variant_str)
                            .or_default()
                            .push(barcode.clone());
                    }
                }
            }
        }

        index
    }

    /// Calculate Hamming distance between two sequences
    fn hamming_distance(a: &str, b: &str) -> u32 {
        if a.len() != b.len() {
            return u32::MAX;
        }
        a.chars().zip(b.chars()).filter(|(a, b)| a != b).count() as u32
    }

    /// Match a barcode with correction
    pub fn match_barcode(&self, barcode: &str) -> BarcodeMatch {
        // First try exact match
        if self.whitelist.contains(barcode) {
            return BarcodeMatch::Exact(barcode.to_string());
        }

        if self.max_distance == 0 {
            return BarcodeMatch::NoMatch(barcode.to_string());
        }

        // Try 1-mismatch lookup using index
        if let Some(candidates) = self.mismatch_index.get(barcode) {
            if candidates.len() == 1 {
                return BarcodeMatch::Corrected(barcode.to_string(), candidates[0].clone(), 1);
            }
            // Multiple candidates - ambiguous, no correction
            return BarcodeMatch::NoMatch(barcode.to_string());
        }

        // For higher distances, do brute force search
        if self.max_distance > 1 {
            let mut best_match: Option<(String, u32)> = None;
            let mut ambiguous = false;

            for wl_barcode in self.whitelist.iter() {
                let dist = Self::hamming_distance(barcode, wl_barcode);
                if dist <= self.max_distance {
                    match &best_match {
                        None => best_match = Some((wl_barcode.clone(), dist)),
                        Some((_, best_dist)) => {
                            if dist < *best_dist {
                                best_match = Some((wl_barcode.clone(), dist));
                                ambiguous = false;
                            } else if dist == *best_dist {
                                ambiguous = true;
                            }
                        }
                    }
                }
            }

            if let Some((corrected, dist)) = best_match {
                if !ambiguous {
                    return BarcodeMatch::Corrected(barcode.to_string(), corrected, dist);
                }
            }
        }

        BarcodeMatch::NoMatch(barcode.to_string())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_exact_match() {
        let barcodes = vec![
            "AAACCCAAGAAACACT".to_string(),
            "AAACCCAAGAAACCAT".to_string(),
        ];
        let whitelist = Whitelist::from_vec(barcodes).unwrap();
        let matcher = BarcodeMatcher::new(whitelist);

        assert!(matches!(
            matcher.match_exact("AAACCCAAGAAACACT"),
            BarcodeMatch::Exact(_)
        ));
        assert!(matches!(
            matcher.match_exact("AAACCCAAGAAACXXX"),
            BarcodeMatch::NoMatch(_)
        ));
    }

    #[test]
    fn test_barcode_correction() {
        let barcodes = vec!["AAACCCAAGAAACACT".to_string()];
        let whitelist = Whitelist::from_vec(barcodes).unwrap();
        let corrector = BarcodeCorrector::new(whitelist, 1);

        // Exact match
        let result = corrector.match_barcode("AAACCCAAGAAACACT");
        assert!(matches!(result, BarcodeMatch::Exact(_)));

        // 1 mismatch - should be corrected
        let result = corrector.match_barcode("TAACCCAAGAAACACT");
        assert!(matches!(result, BarcodeMatch::Corrected(_, _, 1)));

        // 2 mismatches - should not match with max_distance=1
        let result = corrector.match_barcode("TTACCCAAGAAACACT");
        assert!(matches!(result, BarcodeMatch::NoMatch(_)));
    }
}
