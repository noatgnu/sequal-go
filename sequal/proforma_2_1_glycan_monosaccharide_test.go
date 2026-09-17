package sequal

import "testing"

func TestGlycanMonosaccharide_AllStandardNames(t *testing.T) {
	names := []string{
		"Hex", "HexNAc", "HexS", "HexP", "HexNAcS", "HexN", "HexNS",
		"dHex", "Fuc", "aHex", "en,aHex", "Neu", "NeuAc", "NeuGc",
		"Sug", "Tri", "Tet", "Pen", "Hep", "Oct", "Non", "Dec",
		"Kdo", "Kdn", "Sulfo", "Phospho",
	}

	for _, name := range names {
		t.Run(name, func(t *testing.T) {
			proformaString := "N[Glycan:" + name + "1]K"
			seq, err := FromProforma(proformaString)
			if err != nil {
				t.Fatalf("Failed to parse %s: %v", proformaString, err)
			}

			mod := seq.GetSeq()[0].GetMods()[0]
			pv := mod.GetModificationValue().GetPipeValues()[0]
			if !pv.IsValidGlycan() {
				t.Errorf("Expected %s to be a valid monosaccharide", name)
			}
		})
	}
}

func TestGlycanMonosaccharide_KdoKdnRoundTrip(t *testing.T) {
	tests := []string{
		"N[Glycan:Kdo1]K",
		"N[Glycan:Kdn1]K",
		"NEEYN[Glycan:Kdo1Kdn1Hex2HexNAc2]K",
	}

	for _, proforma := range tests {
		t.Run(proforma, func(t *testing.T) {
			seq, err := FromProforma(proforma)
			if err != nil {
				t.Fatalf("Failed to parse %s: %v", proforma, err)
			}

			output := seq.ToProforma()
			if _, err := FromProforma(output); err != nil {
				t.Fatalf("Failed round-trip for %s -> %s: %v", proforma, output, err)
			}
		})
	}
}

func TestGlycanMonosaccharide_ReferenceRepoCases(t *testing.T) {
	cases := []string{
		"Glycan:Hex2HexNAc",
		"Glycan:{C8H14N1O5:z+1}1Hex2",
		"Glycan:{C8H13[15N1]O5}1Hex2",
	}

	for _, mod := range cases {
		t.Run(mod, func(t *testing.T) {
			proformaString := "N[" + mod + "]K"
			seq, err := FromProforma(proformaString)
			if err != nil {
				t.Fatalf("Failed to parse %s: %v", proformaString, err)
			}

			pv := seq.GetSeq()[0].GetMods()[0].GetModificationValue().GetPipeValues()[0]
			if !pv.IsValidGlycan() {
				t.Errorf("Expected %s to be valid", mod)
			}
		})
	}
}
