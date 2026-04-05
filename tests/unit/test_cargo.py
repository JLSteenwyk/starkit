from starkit.cargo import extract_cargo


def test_extract_cargo_basic(simple_record):
    captain_feature = simple_record.features[0]
    cargo = extract_cargo(simple_record, 4000, 20000, captain_feature)
    # Should find the 4 cargo genes but not the captain
    assert len(cargo) == 4
    assert all(c.gene_id != "captain_001" for c in cargo)


def test_extract_cargo_excludes_outside():
    """Cargo genes outside the boundary should be excluded."""
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.SeqFeature import SeqFeature, SimpleLocation

    record = SeqRecord(Seq("ATGC" * 25000), id="c1")
    # Captain at 1000-2000
    captain = SeqFeature(
        location=SimpleLocation(1000, 2000, strand=1),
        type="CDS",
        qualifiers={"protein_id": ["cap"], "translation": ["MAA"]},
    )
    # Cargo inside boundary
    inside = SeqFeature(
        location=SimpleLocation(3000, 4000, strand=1),
        type="CDS",
        qualifiers={"protein_id": ["inside"], "product": ["kinase"], "translation": ["MAA"]},
    )
    # Cargo outside boundary
    outside = SeqFeature(
        location=SimpleLocation(50000, 51000, strand=1),
        type="CDS",
        qualifiers={"protein_id": ["outside"], "product": ["pump"], "translation": ["MAA"]},
    )
    record.features = [captain, inside, outside]

    cargo = extract_cargo(record, 500, 10000, captain)
    assert len(cargo) == 1
    assert cargo[0].gene_id == "inside"
