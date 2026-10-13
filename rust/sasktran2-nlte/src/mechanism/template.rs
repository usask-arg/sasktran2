//! `for_v` expansion and species-id canonicalisation.

use crate::prelude::*;

/// The vibrational levels an entry expands to: `[lo, hi]` inclusive, or a
/// single un-templated entry when `for_v` is absent.
pub(super) fn levels(id: &str, for_v: Option<[u32; 2]>) -> Result<Vec<Option<u32>>> {
    match for_v {
        None => Ok(vec![None]),
        Some([lo, hi]) if lo <= hi => Ok((lo..=hi).map(Some).collect()),
        Some([lo, hi]) => Err(anyhow!("'{id}': for_v range [{lo}, {hi}] is empty")),
    }
}

/// Substitutes `{v}`, `{v-1}`, `{v+2}`, ... with the level `v`.
///
/// Placeholders are an error in entries without `for_v`, as is a substitution
/// that would produce a negative level.
pub(super) fn substitute(text: &str, v: Option<u32>) -> Result<String> {
    let mut out = String::with_capacity(text.len());
    let mut rest = text;
    while let Some(start) = rest.find("{v") {
        out.push_str(&rest[..start]);
        let after = &rest[start + 2..];
        let end = after
            .find('}')
            .ok_or_else(|| anyhow!("unterminated placeholder in '{text}'"))?;
        let offset: i64 = match after[..end].trim() {
            "" => 0,
            delta => delta
                .replace(' ', "")
                .parse()
                .map_err(|_| anyhow!("invalid placeholder '{{v{}}}' in '{text}'", &after[..end]))?,
        };
        let v = v.ok_or_else(|| anyhow!("placeholder in '{text}' but the entry has no for_v"))?;
        let level = i64::from(v) + offset;
        if level < 0 {
            return Err(anyhow!("'{text}' gives a negative level for v={v}"));
        }
        out.push_str(&level.to_string());
        rest = &after[end + 1..];
    }
    out.push_str(rest);
    Ok(out)
}

/// Canonical species id: drops a `v=0` qualifier, so `O2(b, v=0)` is `O2(b)`
/// and `O2(v=0)` is `O2`. Whitespace inside the parentheses is normalised to
/// `", "` separators.
pub(super) fn canonical_species(id: &str) -> String {
    let id = id.trim();
    let Some(open) = id.find('(') else {
        return id.to_string();
    };
    if !id.ends_with(')') {
        return id.to_string();
    }
    let base = &id[..open];
    let items: Vec<&str> = id[open + 1..id.len() - 1]
        .split(',')
        .map(str::trim)
        .filter(|item| !item.is_empty() && item.replace(' ', "") != "v=0")
        .collect();
    if items.is_empty() {
        base.to_string()
    } else {
        format!("{base}({})", items.join(", "))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn substitutes_offsets() {
        assert_eq!(
            substitute("O2(X, v={v}) -> O2(X, v={v-1}) + J_{v+2}", Some(5)).unwrap(),
            "O2(X, v=5) -> O2(X, v=4) + J_7"
        );
    }

    #[test]
    fn placeholder_without_for_v_is_an_error() {
        assert!(substitute("O2(X, v={v})", None).is_err());
    }

    #[test]
    fn negative_level_is_an_error() {
        assert!(substitute("O2(X, v={v-2})", Some(1)).is_err());
    }

    #[test]
    fn plain_text_is_unchanged() {
        assert_eq!(substitute("O(1D)", None).unwrap(), "O(1D)");
    }

    #[test]
    fn empty_range_is_an_error() {
        assert!(levels("x", Some([3, 2])).is_err());
        assert_eq!(levels("x", Some([2, 3])).unwrap(), vec![Some(2), Some(3)]);
    }

    #[test]
    fn canonical_species_drops_ground_vibrational_level() {
        assert_eq!(canonical_species("O2(b, v=0)"), "O2(b)");
        assert_eq!(canonical_species("O2(v=0)"), "O2");
        assert_eq!(canonical_species("O2(b,v=1)"), "O2(b, v=1)");
        assert_eq!(canonical_species("O(1D)"), "O(1D)");
        assert_eq!(canonical_species(" N2 "), "N2");
    }
}
