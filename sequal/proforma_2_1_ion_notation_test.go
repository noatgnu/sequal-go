package sequal

import (
	"testing"
)

func TestIonNotation_ATypeIon(t *testing.T) {
	proforma := "PEPT[a-type-ion]IDE"
	seq, err := FromProforma(proforma)
	if err != nil {
		t.Fatalf("Failed to parse %s: %v", proforma, err)
	}

	if len(seq.GetSeq()) == 0 || len(seq.GetSeq()[3].GetMods()) == 0 {
		t.Fatal("Expected modification on position 3")
	}

	mod := seq.GetSeq()[3].GetMods()[0]
	if !mod.IsIonType() {
		t.Errorf("Expected IsIonType()=true for 'a-type-ion'")
	}
	if mod.GetValue() != "a-type-ion" {
		t.Errorf("Expected value 'a-type-ion', got '%s'", mod.GetValue())
	}
}

func TestIonNotation_BTypeIon(t *testing.T) {
	proforma := "PEPT[b-type-ion]IDE"
	seq, err := FromProforma(proforma)
	if err != nil {
		t.Fatalf("Failed to parse %s: %v", proforma, err)
	}

	if len(seq.GetSeq()) == 0 || len(seq.GetSeq()[3].GetMods()) == 0 {
		t.Fatal("Expected modification on position 3")
	}

	mod := seq.GetSeq()[3].GetMods()[0]
	if !mod.IsIonType() {
		t.Errorf("Expected IsIonType()=true for 'b-type-ion'")
	}
}

func TestIonNotation_AllIonTypes(t *testing.T) {
	ionTypes := []string{
		"a-type-ion",
		"b-type-ion",
		"c-type-ion",
		"x-type-ion",
		"y-type-ion",
		"z-type-ion",
	}

	for _, ionType := range ionTypes {
		t.Run(ionType, func(t *testing.T) {
			proforma := "PEPT[" + ionType + "]IDE"
			seq, err := FromProforma(proforma)
			if err != nil {
				t.Fatalf("Failed to parse %s: %v", proforma, err)
			}

			if len(seq.GetSeq()) == 0 || len(seq.GetSeq()[3].GetMods()) == 0 {
				t.Fatal("Expected modification on position 3")
			}

			mod := seq.GetSeq()[3].GetMods()[0]
			if !mod.IsIonType() {
				t.Errorf("Expected IsIonType()=true for '%s'", ionType)
			}
		})
	}
}

func TestIonNotation_UnimodIonReferences(t *testing.T) {
	tests := []struct {
		name        string
		unimodID    string
		shouldBeIon bool
	}{
		{"a-type-ion", "UNIMOD:140", true},
		{"b-type-ion", "UNIMOD:2132", true},
		{"c-type-ion", "UNIMOD:4", true},
		{"x-type-ion", "UNIMOD:24", true},
		{"y-type-ion", "UNIMOD:2133", true},
		{"z-type-ion", "UNIMOD:23", true},
		{"non-ion", "UNIMOD:21", false}, // Phospho
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			proforma := "PEPT[" + tt.unimodID + "]IDE"
			seq, err := FromProforma(proforma)
			if err != nil {
				t.Fatalf("Failed to parse %s: %v", proforma, err)
			}

			if len(seq.GetSeq()) == 0 || len(seq.GetSeq()[3].GetMods()) == 0 {
				t.Fatal("Expected modification on position 3")
			}

			mod := seq.GetSeq()[3].GetMods()[0]
			if mod.IsIonType() != tt.shouldBeIon {
				t.Errorf("Expected IsIonType()=%v for %s, got %v", tt.shouldBeIon, tt.unimodID, mod.IsIonType())
			}
		})
	}
}

func TestIonNotation_ShortUnimodReferences(t *testing.T) {
	tests := []struct {
		name        string
		unimodID    string
		shouldBeIon bool
	}{
		{"a-type-ion short", "U:140", true},
		{"b-type-ion short", "U:2132", true},
		{"c-type-ion short", "U:4", true},
		{"non-ion short", "U:21", false},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			proforma := "PEPT[" + tt.unimodID + "]IDE"
			seq, err := FromProforma(proforma)
			if err != nil {
				t.Fatalf("Failed to parse %s: %v", proforma, err)
			}

			if len(seq.GetSeq()) == 0 || len(seq.GetSeq()[3].GetMods()) == 0 {
				t.Fatal("Expected modification on position 3")
			}

			mod := seq.GetSeq()[3].GetMods()[0]
			if mod.IsIonType() != tt.shouldBeIon {
				t.Errorf("Expected IsIonType()=%v for %s, got %v", tt.shouldBeIon, tt.unimodID, mod.IsIonType())
			}
		})
	}
}

func TestIonNotation_CaseInsensitive(t *testing.T) {
	tests := []string{
		"a-TYPE-ION",
		"B-Type-Ion",
		"c-type-ION",
	}

	for _, ionType := range tests {
		t.Run(ionType, func(t *testing.T) {
			proforma := "PEPT[" + ionType + "]IDE"
			seq, err := FromProforma(proforma)
			if err != nil {
				t.Fatalf("Failed to parse %s: %v", proforma, err)
			}

			if len(seq.GetSeq()) == 0 || len(seq.GetSeq()[3].GetMods()) == 0 {
				t.Fatal("Expected modification on position 3")
			}

			mod := seq.GetSeq()[3].GetMods()[0]
			if !mod.IsIonType() {
				t.Errorf("Expected IsIonType()=true for '%s' (case insensitive)", ionType)
			}
		})
	}
}

func TestIonNotation_NonIonModifications(t *testing.T) {
	tests := []string{
		"Phospho",
		"Acetyl",
		"Oxidation",
		"+79.966",
		"UNIMOD:21", // Phospho
	}

	for _, modType := range tests {
		t.Run(modType, func(t *testing.T) {
			proforma := "PEPT[" + modType + "]IDE"
			seq, err := FromProforma(proforma)
			if err != nil {
				t.Fatalf("Failed to parse %s: %v", proforma, err)
			}

			if len(seq.GetSeq()) == 0 || len(seq.GetSeq()[3].GetMods()) == 0 {
				t.Fatal("Expected modification on position 3")
			}

			mod := seq.GetSeq()[3].GetMods()[0]
			if mod.IsIonType() {
				t.Errorf("Expected IsIonType()=false for '%s', got true", modType)
			}
		})
	}
}

func TestIonNotation_RoundTrip(t *testing.T) {
	tests := []string{
		"PEPT[a-type-ion]IDE",
		"PEPT[b-type-ion]IDE",
		"PEPT[c-type-ion]IDE",
		"PEPT[UNIMOD:140]IDE",
		"PEPT[UNIMOD:2132]IDE",
		"[b-type-ion]-PEPTIDE",
	}

	for _, proforma := range tests {
		t.Run(proforma, func(t *testing.T) {
			seq, err := FromProforma(proforma)
			if err != nil {
				t.Fatalf("Failed to parse %s: %v", proforma, err)
			}

			output := seq.ToProforma()
			seq2, err := FromProforma(output)
			if err != nil {
				t.Fatalf("Failed to re-parse %s: %v", output, err)
			}

			// Check that IsIonType is preserved
			var mod1, mod2 *Modification
			if len(seq.GetSeq()) > 3 && len(seq.GetSeq()[3].GetMods()) > 0 {
				mod1 = seq.GetSeq()[3].GetMods()[0]
				mod2 = seq2.GetSeq()[3].GetMods()[0]
			}

			if mod1 != nil && mod2 != nil {
				if mod1.IsIonType() != mod2.IsIonType() {
					t.Errorf("IsIonType not preserved after roundtrip: %v -> %v", mod1.IsIonType(), mod2.IsIonType())
				}
			}
		})
	}
}
