package sequal

import (
	"testing"
)

func TestIntegration_PlacementControlsWithIonNotation(t *testing.T) {
	tests := []struct {
		name           string
		proformaString string
		description    string
		shouldParse    bool
	}{
		{
			name:           "Global placement controls with ion type in sequence",
			proformaString: "<[Phospho|Position:S,T,Y|Limit:2]@S,T,Y>PEPT[a-type-ion]IDE",
			description:    "Combines global placement controls with ion notation",
			shouldParse:    true,
		},
		{
			name:           "CoMKP with multiple ion types",
			proformaString: "<[Oxidation|Position:M|CoMKP]@M>ME[b-type-ion]PT[c-type-ion]IDE",
			description:    "Global CoMKP with multiple different ion types",
			shouldParse:    true,
		},
		{
			name:           "CoMUP with Unimod ion reference",
			proformaString: "<[Acetyl|Position:K|CoMUP]@K>PEPT[UNIMOD:140]IDE",
			description:    "Global CoMUP with Unimod ion type reference",
			shouldParse:    true,
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			seq, err := FromProforma(tt.proformaString)
			if tt.shouldParse {
				if err != nil {
					t.Errorf("Failed to parse %s: %v", tt.proformaString, err)
					return
				}

				if len(seq.GetGlobalMods()) == 0 {
					t.Errorf("Expected global modification for %s", tt.proformaString)
					return
				}

				globalMod := seq.GetGlobalMods()[0]
				if globalMod.GetPositionConstraint() == nil {
					t.Errorf("Expected position constraint in global mod")
				}

				hasIonType := false
				for _, aa := range seq.GetSeq() {
					for _, mod := range aa.GetMods() {
						if mod.IsIonType() {
							hasIonType = true
							break
						}
					}
				}

				if !hasIonType {
					t.Errorf("Expected at least one ion type modification")
				}
			} else {
				if err == nil {
					t.Errorf("Expected error parsing %s", tt.proformaString)
				}
			}
		})
	}
}

func TestIntegration_ChargedFormulasWithIonTypes(t *testing.T) {
	tests := []struct {
		name           string
		proformaString string
		description    string
		shouldParse    bool
	}{
		{
			name:           "Charged formula with ion type",
			proformaString: "PE[a-type-ion]PT[Formula:Zn1:z+2]IDE",
			description:    "Ion type and charged formula in sequence",
			shouldParse:    true,
		},
		{
			name:           "Multiple charged formulas with placement controls",
			proformaString: "<[Phospho|Position:S|Limit:1]@S>PEPT[Formula:Cu1:z+1]IDE",
			description:    "Global placement controls with charged formula",
			shouldParse:    true,
		},
		{
			name:           "Charged formula with multiple ion types",
			proformaString: "PE[a-type-ion]PT[Formula:Cu1:z+1]ID[b-type-ion]E",
			description:    "Charged formula with multiple ion types",
			shouldParse:    true,
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			seq, err := FromProforma(tt.proformaString)
			if tt.shouldParse {
				if err != nil {
					t.Errorf("Failed to parse %s: %v", tt.proformaString, err)
					return
				}

				hasChargedFormula := false
				hasIonType := false

				for _, aa := range seq.GetSeq() {
					for _, mod := range aa.GetMods() {
						if mod.GetModificationValue() != nil {
							for _, pv := range mod.GetModificationValue().GetPipeValues() {
								if pv.GetCharge() != nil {
									hasChargedFormula = true
									break
								}
							}
						}
						if mod.IsIonType() {
							hasIonType = true
						}
					}
				}

				if tt.name == "Charged formula with ion type" {
					if !hasChargedFormula {
						t.Errorf("Expected charged formula")
					}
					if !hasIonType {
						t.Errorf("Expected ion type")
					}
				}

				if tt.name == "Multiple charged formulas with placement controls" {
					if !hasChargedFormula {
						t.Errorf("Expected charged formula")
					}
					if len(seq.GetGlobalMods()) == 0 {
						t.Errorf("Expected global modification")
					}
				}

				if tt.name == "Charged formula with multiple ion types" {
					if !hasChargedFormula || !hasIonType {
						t.Errorf("Expected charged formula and ion types")
					}
				}
			} else {
				if err == nil {
					t.Errorf("Expected error parsing %s", tt.proformaString)
				}
			}
		})
	}
}

func TestIntegration_ComplexGlobalModifications(t *testing.T) {
	tests := []struct {
		name           string
		proformaString string
		description    string
		shouldParse    bool
	}{
		{
			name:           "Multiple global mods with different placement controls",
			proformaString: "<[Phospho|Position:S,T,Y|Limit:2|CoMKP]@S,T,Y><[Oxidation|Position:M|CoMUP]@M>STMYSM",
			description:    "Two global mods with different placement control combinations",
			shouldParse:    true,
		},
		{
			name:           "Global isotope with global fixed mod",
			proformaString: "<13C><[Carbamidomethyl|Position:C]@C>PEPTCDE",
			description:    "Isotope label combined with fixed modification with position constraint",
			shouldParse:    true,
		},
		{
			name:           "Terminal global with placement controls",
			proformaString: "<[Acetyl|Position:K|Limit:1]@K>[Acetyl]-KKKK",
			description:    "Global placement controls with N-terminal modification",
			shouldParse:    true,
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			seq, err := FromProforma(tt.proformaString)
			if tt.shouldParse {
				if err != nil {
					t.Errorf("Failed to parse %s: %v", tt.proformaString, err)
					return
				}

				if tt.name == "Multiple global mods with different placement controls" {
					if len(seq.GetGlobalMods()) != 2 {
						t.Errorf("Expected 2 global modifications, got %d", len(seq.GetGlobalMods()))
						return
					}

					mod1 := seq.GetGlobalMods()[0]
					mod2 := seq.GetGlobalMods()[1]

					if !mod1.GetColocalizeKnown() {
						t.Errorf("Expected first mod to have CoMKP")
					}

					if !mod2.GetColocalizeUnknown() {
						t.Errorf("Expected second mod to have CoMUP")
					}
				}

				if tt.name == "Global isotope with global fixed mod" {
					if len(seq.GetGlobalMods()) != 2 {
						t.Errorf("Expected 2 global modifications, got %d", len(seq.GetGlobalMods()))
					}
				}
			} else {
				if err == nil {
					t.Errorf("Expected error parsing %s", tt.proformaString)
				}
			}
		})
	}
}

func TestIntegration_AllFeaturesCombined(t *testing.T) {
	tests := []struct {
		name           string
		proformaString string
		description    string
	}{
		{
			name:           "Maximum feature combination",
			proformaString: "<[Phospho|Position:S,T,Y|Limit:2|CoMKP]@S,T,Y>[Acetyl]-PE[a-type-ion]PT[Formula:Zn1:z+2]IDE",
			description:    "Placement controls + terminal mod + ion notation + charged formula",
		},
		{
			name:           "Complex real-world scenario",
			proformaString: "<13C><[Carbamidomethyl|Position:C]@C>[Acetyl]-ACDE[Phospho]FGH[b-type-ion]IKL[Formula:Cu1:z+1]M",
			description:    "Isotope + placement controls + terminal mod + phosphorylation + ion type + charged formula",
		},
		{
			name:           "Multiple placement controls with various modifications",
			proformaString: "<[Oxidation|Position:M|Limit:1|CoMUP]@M><[Phospho|Position:S,T,Y|CoMKP]@S,T,Y>MS[Phospho]TM[UNIMOD:140]Y",
			description:    "Two global mods with placement controls + local phospho + Unimod ion",
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			seq, err := FromProforma(tt.proformaString)
			if err != nil {
				t.Fatalf("Failed to parse %s: %v", tt.proformaString, err)
			}

			output := seq.ToProforma()
			seq2, err := FromProforma(output)
			if err != nil {
				t.Fatalf("Failed to re-parse %s: %v", output, err)
			}

			if len(seq.GetGlobalMods()) != len(seq2.GetGlobalMods()) {
				t.Errorf("Global mod count mismatch: %d != %d", len(seq.GetGlobalMods()), len(seq2.GetGlobalMods()))
			}

			for i := range seq.GetGlobalMods() {
				if i >= len(seq2.GetGlobalMods()) {
					break
				}
				mod1 := seq.GetGlobalMods()[i]
				mod2 := seq2.GetGlobalMods()[i]

				if mod1.GetColocalizeKnown() != mod2.GetColocalizeKnown() {
					t.Errorf("CoMKP mismatch at global mod %d", i)
				}
				if mod1.GetColocalizeUnknown() != mod2.GetColocalizeUnknown() {
					t.Errorf("CoMUP mismatch at global mod %d", i)
				}
			}

			totalIonTypes1 := 0
			totalIonTypes2 := 0
			for _, aa := range seq.GetSeq() {
				for _, mod := range aa.GetMods() {
					if mod.IsIonType() {
						totalIonTypes1++
					}
				}
			}
			for _, aa := range seq2.GetSeq() {
				for _, mod := range aa.GetMods() {
					if mod.IsIonType() {
						totalIonTypes2++
					}
				}
			}

			if totalIonTypes1 != totalIonTypes2 {
				t.Errorf("Ion type count mismatch: %d != %d", totalIonTypes1, totalIonTypes2)
			}
		})
	}
}

func TestIntegration_EdgeCases(t *testing.T) {
	tests := []struct {
		name           string
		proformaString string
		description    string
		shouldParse    bool
	}{
		{
			name:           "Empty position constraint with other controls",
			proformaString: "<[Phospho|Limit:2|CoMKP]@S,T,Y>PEPTIDES",
			description:    "Limit and CoMKP without Position tag",
			shouldParse:    true,
		},
		{
			name:           "Ion type at N-terminus with global mod",
			proformaString: "<[Oxidation|Position:M]@M>[a-type-ion]-PEPTIDE",
			description:    "Ion type as N-terminal modification with global placement control",
			shouldParse:    true,
		},
		{
			name:           "Multiple ion types in sequence",
			proformaString: "PE[a-type-ion]PT[b-type-ion]ID[c-type-ion]E",
			description:    "Multiple different ion types in same sequence",
			shouldParse:    true,
		},
		{
			name:           "Charged formula with terminal mod and placement controls",
			proformaString: "<[Phospho|Position:S|Limit:1]@S>[Acetyl]-PEPT[Formula:Zn1:z+2]IDES",
			description:    "All major 2.1 features combined",
			shouldParse:    true,
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			seq, err := FromProforma(tt.proformaString)
			if tt.shouldParse {
				if err != nil {
					t.Errorf("Failed to parse %s: %v", tt.proformaString, err)
					return
				}

				output := seq.ToProforma()
				_, err = FromProforma(output)
				if err != nil {
					t.Errorf("Failed round-trip for %s: %v", tt.proformaString, err)
				}
			} else {
				if err == nil {
					t.Errorf("Expected error parsing %s", tt.proformaString)
				}
			}
		})
	}
}

func TestIntegration_SerializationConsistency(t *testing.T) {
	tests := []string{
		"<[Phospho|Position:S,T,Y|Limit:2|CoMKP]@S,T,Y>PEPTIDES",
		"<[Oxidation|Position:M|CoMUP]@M>MMMM",
		"[Acetyl]-PEPT[Formula:Zn1:z+2]IDE",
		"PEPT[a-type-ion]IDE",
		"<13C><[Carbamidomethyl|Position:C]@C>PEPTCDE",
		"<[Phospho|Position:S,T,Y|Limit:2|CoMKP]@S,T,Y>[Acetyl]-PE[a-type-ion]PT[Formula:Zn1:z+2]IDE",
	}

	for _, proforma := range tests {
		t.Run(proforma, func(t *testing.T) {
			seq1, err := FromProforma(proforma)
			if err != nil {
				t.Fatalf("Failed to parse %s: %v", proforma, err)
			}

			output := seq1.ToProforma()
			seq2, err := FromProforma(output)
			if err != nil {
				t.Fatalf("Failed to re-parse %s: %v", output, err)
			}

			output2 := seq2.ToProforma()
			seq3, err := FromProforma(output2)
			if err != nil {
				t.Fatalf("Failed second round-trip for %s: %v", output2, err)
			}

			if len(seq1.GetGlobalMods()) != len(seq3.GetGlobalMods()) {
				t.Errorf("Global mod count changed after multiple round-trips")
			}

			if len(seq1.GetSeq()) != len(seq3.GetSeq()) {
				t.Errorf("Sequence length changed after multiple round-trips")
			}
		})
	}
}
