use anyhow::Result;
use std::fs::File;
use std::io::Write;
use tempfile::tempdir;

use crate::uniprot::UniProtSequence;

const TEST_BCL2_XML: &str =
    "tests/testdata/index/uniprotkb_gene_bcl2_AND_reviewed_true_2025_06_04.xml";

#[test]
fn test_bcl2_xml_parsing() -> Result<()> {
    // Parse the BCL2 test file
    let proteins = UniProtSequence::from_xml(TEST_BCL2_XML)?;

    // Basic validation
    assert!(
        !proteins.is_empty(),
        "Should have found at least one protein"
    );

    // Check the first protein
    let protein = &proteins[0];

    // BCL2 proteins should have an ID containing "Bcl-2"
    assert!(
        protein.id.contains("Bcl-2"),
        "Protein ID should contain Bcl-2"
    );

    // Should have a non-empty accession
    assert!(
        !protein.accession.is_empty(),
        "Protein should have an accession number"
    );

    // Should have a sequence
    assert!(
        !protein.sequence.is_empty(),
        "Protein should have a sequence"
    );

    // Should have some features
    assert!(!protein.features.is_empty(), "Protein should have features");

    // Print some debug information
    println!("Found {} proteins", proteins.len());
    println!("First protein:");
    println!("  ID: {}", protein.id);
    println!("  Accession: {}", protein.accession);
    println!("  Sequence length: {}", protein.sequence.len());
    println!("  Number of features: {}", protein.features.len());
    println!("\nFeature types:");
    for feature in &protein.features {
        println!("  {} ({})", feature.feature_type, feature.description);
    }

    Ok(())
}

#[test]
fn test_simple_xml_parsing() -> Result<()> {
    // Create a temporary directory for test files
    let dir = tempdir()?;
    let file_path = dir.path().join("test.xml");

    // Create a test XML file
    let xml_content = r#"<?xml version="1.0" encoding="UTF-8"?>
    <uniprot xmlns="http://uniprot.org/uniprot">
        <entry>
            <protein>
                <recommendedName>
                    <fullName>Apoptosis regulator Bcl-2</fullName>
                </recommendedName>
            </protein>
            <accession>P10415</accession>
            <sequence>MAHAGRTGYDNREIVMKYIHYKLSQRGYEWDAGDVGAAPPGAAPAPGIFSSQPGHTPHPAASRDPVARTSPLQTPAAPGAAAGPALSPVPPVVHLTLRQAGDDFSRRYRRDFAEMSSQLHLTPFTARGRFATVVEELFRDGVNWGRIVAFFEFGGVMCVESVNREMSPLVDNIALWMTEYLNRHLHTWIQDNGGWDAFVELYGPSMRPLFDFSWLSLKTLLSLALVGACITLGAYLGHK</sequence>
            <feature type="chain" description="Apoptosis regulator Bcl-2">
                <location>
                    <begin position="1"/>
                    <end position="239"/>
                </location>
            </feature>
            <feature type="domain" description="BH4">
                <location>
                    <begin position="10"/>
                    <end position="30"/>
                </location>
            </feature>
        </entry>
    </uniprot>"#;

    let mut file = File::create(&file_path)?;
    file.write_all(xml_content.as_bytes())?;

    // Parse the test file
    let proteins = UniProtSequence::from_xml(&file_path)?;

    assert_eq!(proteins.len(), 1);
    let protein = &proteins[0];
    assert_eq!(protein.id, "Apoptosis regulator Bcl-2");
    assert_eq!(protein.accession, "P10415");
    assert!(!protein.sequence.is_empty());
    assert_eq!(protein.features.len(), 2);

    let feature = &protein.features[0];
    assert_eq!(feature.feature_type, "chain");
    assert_eq!(feature.description, "Apoptosis regulator Bcl-2");

    Ok(())
}
