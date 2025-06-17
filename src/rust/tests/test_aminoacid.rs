use crate::aminoacid::AminoAcidAmbiguity;

#[test]
fn test_valid_amino_acids() {
    let aa = AminoAcidAmbiguity::new();

    // Test standard amino acids
    for c in [
        'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V',
        'W', 'Y',
    ]
    .iter()
    {
        assert!(aa.is_valid_aa(*c));
    }

    // Test ambiguous codes
    assert!(aa.is_valid_aa('X'));
    assert!(aa.is_valid_aa('B'));
    assert!(aa.is_valid_aa('Z'));
    assert!(aa.is_valid_aa('J'));

    // Test invalid characters
    assert!(!aa.is_valid_aa('1'));
    assert!(!aa.is_valid_aa('$'));
}

#[test]
fn test_ambiguity_resolution() {
    let aa = AminoAcidAmbiguity::new();

    // Test standard amino acids (should return themselves)
    for c in [
        'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V',
        'W', 'Y',
    ]
    .iter()
    {
        assert_eq!(aa.resolve_ambiguity(*c), *c);
    }

    // Test ambiguous codes (should return one of their possible resolutions)
    let b_resolution = aa.resolve_ambiguity('B');
    assert!(vec!['D', 'N'].contains(&b_resolution));

    let z_resolution = aa.resolve_ambiguity('Z');
    assert!(vec!['E', 'Q'].contains(&z_resolution));

    let j_resolution = aa.resolve_ambiguity('J');
    assert!(vec!['I', 'L'].contains(&j_resolution));
}

#[test]
fn test_sequence_validation() {
    let aa = AminoAcidAmbiguity::new();

    // Test valid sequence
    assert!(aa.validate_sequence("ACDEFGHIKLMNPQRSTVWY").is_ok());
    assert!(aa.validate_sequence("ACDEFXBZJ").is_ok());

    // Test invalid sequence
    let result = aa.validate_sequence("ACDEF1GHIKL");
    assert!(result.is_err());
    assert!(result.unwrap_err().to_string().contains("Invalid amino acid '1'"));
}
