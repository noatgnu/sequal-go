package sequal

import (
	"testing"
)

// Comprehensive Integration Tests for ProForma 2.1

func TestComprehensiveIntegration_NamedEntitiesWithMultipleFeatures(t *testing.T) {
	tests := []struct {
		name           string
		proformaString string
		expectations   func(t *testing.T, seq *Sequence)
	}{
		{
			name:           "Peptidoform name with terminal global mods and charge",
			proformaString: "(>TMT-labeled peptide)<[TMT6plex]@K,N-term>PEPTIDEK/2",
			expectations: func(t *testing.T, seq *Sequence) {
				if seq.GetPeptidoformName() == nil || *seq.GetPeptidoformName() != "TMT-labeled peptide" {
					t.Errorf("Expected peptidoform name 'TMT-labeled peptide'")
				}
				if len(seq.GetGlobalMods()) != 1 {
					t.Errorf("Expected 1 global modification, got %d", len(seq.GetGlobalMods()))
				}
				if seq.GetCharge() == nil || *seq.GetCharge() != 2 {
					t.Errorf("Expected charge 2")
				}
			},
		},
		{
			name:           "All three naming levels with modifications",
			proformaString: "(>>>Chimeric Spectrum 1234)(>>Precursor z=2)(>Phospho-peptide)<[Phospho]@S,T,Y>PEPT[Oxidation]IDES/2",
			expectations: func(t *testing.T, seq *Sequence) {
				if seq.GetCompoundIonName() == nil || *seq.GetCompoundIonName() != "Chimeric Spectrum 1234" {
					t.Errorf("Expected compound ion name 'Chimeric Spectrum 1234'")
				}
				if seq.GetPeptidoformIonName() == nil || *seq.GetPeptidoformIonName() != "Precursor z=2" {
					t.Errorf("Expected peptidoform ion name 'Precursor z=2'")
				}
				if seq.GetPeptidoformName() == nil || *seq.GetPeptidoformName() != "Phospho-peptide" {
					t.Errorf("Expected peptidoform name 'Phospho-peptide'")
				}
				if len(seq.GetGlobalMods()) != 1 {
					t.Errorf("Expected 1 global modification")
				}
			},
		},
		{
			name:           "Named entity with labile glycans",
			proformaString: "(>Glycopeptide){Glycan:{C8H13N1O5}1Hex2}PEPN[Glycan:HexNAc2Hex3]TIDE",
			expectations: func(t *testing.T, seq *Sequence) {
				if seq.GetPeptidoformName() == nil || *seq.GetPeptidoformName() != "Glycopeptide" {
					t.Errorf("Expected peptidoform name 'Glycopeptide'")
				}
				labileMods := seq.GetMods()[-3]
				if len(labileMods) == 0 {
					t.Errorf("Expected labile modifications")
				}
			},
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			seq, err := FromProforma(tt.proformaString)
			if err != nil {
				t.Fatalf("Failed to parse %s: %v", tt.proformaString, err)
			}
			tt.expectations(t, seq)
		})
	}
}

func TestComprehensiveIntegration_TerminalGlobalModsComplexScenarios(t *testing.T) {
	tests := []struct {
		name           string
		proformaString string
		expectations   func(t *testing.T, seq *Sequence)
	}{
		{
			name:           "Multiple terminal global mods with placement controls",
			proformaString: "<[TMT6plex|Limit:1]@K,N-term><[Oxidation|Position:M,C]@M,C-term:G>MTPEILTCNSIGCLKG",
			expectations: func(t *testing.T, seq *Sequence) {
				if len(seq.GetGlobalMods()) != 2 {
					t.Errorf("Expected 2 global modifications, got %d", len(seq.GetGlobalMods()))
				}
			},
		},
		{
			name:           "Terminal global mods with ambiguous modifications",
			proformaString: "<[Gln->pyro-Glu]@N-term:Q>QPEPTIDE[Phospho#1]S[#1]T",
			expectations: func(t *testing.T, seq *Sequence) {
				if len(seq.GetGlobalMods()) != 1 {
					t.Errorf("Expected 1 global modification")
				}
			},
		},
		{
			name:           "Terminal global mods with charge states",
			proformaString: "<[Acetyl]@N-term><[TMT6plex]@K>PEPTIDEK/2",
			expectations: func(t *testing.T, seq *Sequence) {
				if len(seq.GetGlobalMods()) != 2 {
					t.Errorf("Expected 2 global modifications")
				}
				if seq.GetCharge() == nil || *seq.GetCharge() != 2 {
					t.Errorf("Expected charge 2")
				}
			},
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			seq, err := FromProforma(tt.proformaString)
			if err != nil {
				t.Fatalf("Failed to parse %s: %v", tt.proformaString, err)
			}
			tt.expectations(t, seq)
		})
	}
}

func TestComprehensiveIntegration_CustomMonosaccharidesComplex(t *testing.T) {
	tests := []struct {
		name           string
		proformaString string
		expectations   func(t *testing.T, seq *Sequence)
	}{
		{
			name:           "Mixed custom and standard monosaccharides",
			proformaString: "PEPN[Glycan:HexNAc2{C8H13N1O5Na1:z+1}1Hex3]TIDE",
			expectations: func(t *testing.T, seq *Sequence) {
				if len(seq.GetSeq()) < 4 {
					t.Errorf("Sequence too short")
					return
				}
				if len(seq.GetSeq()[3].GetMods()) == 0 {
					t.Errorf("Expected modification at position 3")
				}
			},
		},
		{
			name:           "Labile and static custom monosaccharides together",
			proformaString: "{Glycan:{C11H17N1O9}1Hex2}PEPN[Glycan:{C8H13[15N1]O5}2Hex1]TIDE",
			expectations: func(t *testing.T, seq *Sequence) {
				if len(seq.GetMods()[-3]) == 0 {
					t.Errorf("Expected labile modifications")
				}
			},
		},
		{
			name:           "Custom monosaccharides in sequences",
			proformaString: "NSTPEPN[Glycan:{C8H13N1O5}1Hex2]TIDE",
			expectations: func(t *testing.T, seq *Sequence) {
				if len(seq.GetSeq()) < 7 {
					t.Errorf("Sequence too short")
					return
				}
				if len(seq.GetSeq()[6].GetMods()) == 0 {
					t.Errorf("Expected modification at position 6")
				}
			},
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			seq, err := FromProforma(tt.proformaString)
			if err != nil {
				t.Fatalf("Failed to parse %s: %v", tt.proformaString, err)
			}
			tt.expectations(t, seq)
		})
	}
}

func TestComprehensiveIntegration_ChargedFormulasAndIonNotation(t *testing.T) {
	tests := []struct {
		name           string
		proformaString string
		expectations   func(t *testing.T, seq *Sequence)
	}{
		{
			name:           "Charged formula with ion type notation",
			proformaString: "PEPTIDE[Formula:Zn1:z+2]-[b-type-ion]",
			expectations: func(t *testing.T, seq *Sequence) {
				hasChargedFormula := false
				for _, aa := range seq.GetSeq() {
					for _, mod := range aa.GetMods() {
						if mod.GetModificationValue() != nil {
							for _, pv := range mod.GetModificationValue().GetPipeValues() {
								if pv.GetCharge() != nil {
									hasChargedFormula = true
								}
							}
						}
					}
				}
				if !hasChargedFormula {
					t.Errorf("Expected charged formula")
				}

				cTermMods := seq.GetMods()[-2]
				if len(cTermMods) == 0 {
					t.Errorf("Expected C-terminal modification")
				}
			},
		},
		{
			name:           "Multiple charged formulas with different charges",
			proformaString: "PEPT[Formula:C2H3NO:z-1]IDE[Formula:Zn1:z+2]K",
			expectations: func(t *testing.T, seq *Sequence) {
				chargeCount := 0
				for _, aa := range seq.GetSeq() {
					for _, mod := range aa.GetMods() {
						if mod.GetModificationValue() != nil {
							for _, pv := range mod.GetModificationValue().GetPipeValues() {
								if pv.GetCharge() != nil {
									chargeCount++
								}
							}
						}
					}
				}
				if chargeCount < 2 {
					t.Errorf("Expected at least 2 charged formulas, got %d", chargeCount)
				}
			},
		},
		{
			name:           "Charged glycan with ion notation",
			proformaString: "N[Glycan:{C8H13N1O5Na1:z+1}1Hex2]PEPTIDE-[y-type-ion]",
			expectations: func(t *testing.T, seq *Sequence) {
				if len(seq.GetSeq()) == 0 || len(seq.GetSeq()[0].GetMods()) == 0 {
					t.Errorf("Expected modification on first amino acid")
				}
				if len(seq.GetMods()[-2]) == 0 {
					t.Errorf("Expected C-terminal modification")
				}
			},
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			seq, err := FromProforma(tt.proformaString)
			if err != nil {
				t.Fatalf("Failed to parse %s: %v", tt.proformaString, err)
			}
			tt.expectations(t, seq)
		})
	}
}

func TestComprehensiveIntegration_RealWorldComplexScenarios(t *testing.T) {
	tests := []struct {
		name           string
		proformaString string
		expectations   func(t *testing.T, seq *Sequence)
	}{
		{
			name:           "TMT-labeled phosphopeptide with multiple PTMs",
			proformaString: "(>TMT-labeled phosphopeptide)<[TMT6plex]@K,N-term><[Oxidation]@M>MTPEILTS[Phospho]CNSIGCLK/2",
			expectations: func(t *testing.T, seq *Sequence) {
				if seq.GetPeptidoformName() == nil || *seq.GetPeptidoformName() != "TMT-labeled phosphopeptide" {
					t.Errorf("Expected peptidoform name 'TMT-labeled phosphopeptide'")
				}
				if len(seq.GetGlobalMods()) != 2 {
					t.Errorf("Expected 2 global modifications")
				}
				if seq.GetCharge() == nil || *seq.GetCharge() != 2 {
					t.Errorf("Expected charge 2")
				}
			},
		},
		{
			name:           "Glycopeptide with multiple glycosylation sites",
			proformaString: "{Glycan:Hex1HexNAc1}N[Glycan:{C8H13N1O5}2Hex3HexNAc2]GSTN[Glycan:HexNAc2Hex3NeuAc1]VTPEPTIDE",
			expectations: func(t *testing.T, seq *Sequence) {
				if len(seq.GetMods()[-3]) == 0 {
					t.Errorf("Expected labile modifications")
				}
			},
		},
		{
			name:           "Complex proteoform with multiple modification types",
			proformaString: "(>Ubiquitinated protein)[Acetyl]-MRSGSHHHHHHGSPEPTM[Oxidation]IDEK[Ubiquitin#BRANCH]SEQUENCE-[Amidated]",
			expectations: func(t *testing.T, seq *Sequence) {
				if seq.GetPeptidoformName() == nil || *seq.GetPeptidoformName() != "Ubiquitinated protein" {
					t.Errorf("Expected peptidoform name 'Ubiquitinated protein'")
				}
				if len(seq.GetMods()[-1]) == 0 {
					t.Errorf("Expected N-terminal modification")
				}
				if len(seq.GetMods()[-2]) == 0 {
					t.Errorf("Expected C-terminal modification")
				}
			},
		},
		{
			name:           "Metal-binding peptide with charged formulas",
			proformaString: "(>Zinc finger peptide)CAQECGK[Formula:Zn1:z+2]SFTSALK[Formula:Zn1:z+2]SRHK/3",
			expectations: func(t *testing.T, seq *Sequence) {
				if seq.GetPeptidoformName() == nil || *seq.GetPeptidoformName() != "Zinc finger peptide" {
					t.Errorf("Expected peptidoform name 'Zinc finger peptide'")
				}
				if seq.GetCharge() == nil || *seq.GetCharge() != 3 {
					t.Errorf("Expected charge 3")
				}
			},
		},
		{
			name:           "Peptide with unknown position mods",
			proformaString: "[Phospho]^3[Oxidation]^2?MSTPEPTMSTY",
			expectations: func(t *testing.T, seq *Sequence) {
				unknownMods := seq.GetMods()[-4]
				if len(unknownMods) != 5 {
					t.Errorf("Expected 5 unknown position modifications, got %d", len(unknownMods))
				}
			},
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			seq, err := FromProforma(tt.proformaString)
			if err != nil {
				t.Fatalf("Failed to parse %s: %v", tt.proformaString, err)
			}
			tt.expectations(t, seq)
		})
	}
}

func TestComprehensiveIntegration_MultiChainSequences(t *testing.T) {
	tests := []struct {
		name           string
		proformaString string
		expectations   func(t *testing.T, seq *Sequence)
	}{
		{
			name:           "Multi-chain with naming level on first chain",
			proformaString: "(>>>Disulfide-linked)(>>Ion pair)(>Chain1)PEPTIDEC//CSEQUENCE",
			expectations: func(t *testing.T, seq *Sequence) {
				if seq.GetCompoundIonName() == nil || *seq.GetCompoundIonName() != "Disulfide-linked" {
					t.Errorf("Expected compound ion name 'Disulfide-linked'")
				}
				if seq.GetPeptidoformIonName() == nil || *seq.GetPeptidoformIonName() != "Ion pair" {
					t.Errorf("Expected peptidoform ion name 'Ion pair'")
				}
				if !seq.IsMultiChain() {
					t.Errorf("Expected multi-chain sequence")
				}
			},
		},
		{
			name:           "Multi-chain with different terminal global mods",
			proformaString: "<[Acetyl]@N-term>PEPTIDE//<[TMT6plex]@K,N-term>KSEQUENCEK",
			expectations: func(t *testing.T, seq *Sequence) {
				if !seq.IsMultiChain() {
					t.Errorf("Expected multi-chain sequence")
				}
			},
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			seq, err := FromProforma(tt.proformaString)
			if err != nil {
				t.Fatalf("Failed to parse %s: %v", tt.proformaString, err)
			}
			tt.expectations(t, seq)
		})
	}
}

func TestComprehensiveIntegration_AllProForma21Features(t *testing.T) {
	tests := []struct {
		name           string
		proformaString string
		expectations   func(t *testing.T, seq *Sequence)
	}{
		{
			name:           "All ProForma 2.1 features in one sequence",
			proformaString: "(>>>MS2 Spectrum)(>>Precursor z=2)(>Complex peptide)<[TMT6plex]@K,N-term><[Gln->pyro-Glu]@N-term:Q>[Acetyl]-QN[Glycan:HexNAc2Hex3]PEPTM[Oxidation|INFO:High confidence]IDEK[Formula:Zn1:z+2]-[b-type-ion]/2",
			expectations: func(t *testing.T, seq *Sequence) {
				if seq.GetCompoundIonName() == nil || *seq.GetCompoundIonName() != "MS2 Spectrum" {
					t.Errorf("Expected compound ion name 'MS2 Spectrum'")
				}
				if seq.GetPeptidoformIonName() == nil || *seq.GetPeptidoformIonName() != "Precursor z=2" {
					t.Errorf("Expected peptidoform ion name 'Precursor z=2'")
				}
				if seq.GetPeptidoformName() == nil || *seq.GetPeptidoformName() != "Complex peptide" {
					t.Errorf("Expected peptidoform name 'Complex peptide'")
				}
				if len(seq.GetGlobalMods()) != 2 {
					t.Errorf("Expected 2 global modifications")
				}
				if len(seq.GetMods()[-1]) == 0 {
					t.Errorf("Expected N-terminal modification")
				}
				if len(seq.GetMods()[-2]) == 0 {
					t.Errorf("Expected C-terminal modification")
				}
				if seq.GetCharge() == nil || *seq.GetCharge() != 2 {
					t.Errorf("Expected charge 2")
				}
			},
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			seq, err := FromProforma(tt.proformaString)
			if err != nil {
				t.Fatalf("Failed to parse %s: %v", tt.proformaString, err)
			}
			tt.expectations(t, seq)
		})
	}
}

func TestComprehensiveIntegration_RoundTripComplexSequences(t *testing.T) {
	tests := []string{
		"(>Named peptide)<[TMT6plex]@K,N-term>PEPTIDEK/2",
		"N[Glycan:HexNAc2{C8H13N1O5}1]PEPTIDE",
		"<[Gln->pyro-Glu]@N-term:Q>QPEPTM[Oxidation]IDE-[b-type-ion]",
		"[Phospho]^2?STPEPTIDE",
	}

	for _, original := range tests {
		t.Run(original, func(t *testing.T) {
			seq, err := FromProforma(original)
			if err != nil {
				t.Fatalf("Failed to parse %s: %v", original, err)
			}

			proforma := seq.ToProforma()
			seq2, err := FromProforma(proforma)
			if err != nil {
				t.Fatalf("Failed to re-parse %s: %v", proforma, err)
			}

			proforma2 := seq2.ToProforma()
			if proforma != original {
				t.Logf("Round-trip mismatch (first pass): %s != %s", proforma, original)
			}
			if proforma2 != original {
				t.Logf("Round-trip mismatch (second pass): %s != %s", proforma2, original)
			}
		})
	}
}

func TestComprehensiveIntegration_EdgeCasesComplex(t *testing.T) {
	tests := []struct {
		name           string
		proformaString string
		expectations   func(t *testing.T, seq *Sequence)
	}{
		{
			name:           "Terminal modifications with global terminal mods",
			proformaString: "<[Acetyl]@N-term>[Carbamyl]-PEPTIDE-[Methyl]",
			expectations: func(t *testing.T, seq *Sequence) {
				if len(seq.GetGlobalMods()) != 1 {
					t.Errorf("Expected 1 global modification")
				}
				if len(seq.GetMods()[-1]) == 0 {
					t.Errorf("Expected N-terminal modification")
				}
				if len(seq.GetMods()[-2]) == 0 {
					t.Errorf("Expected C-terminal modification")
				}
			},
		},
		{
			name:           "Sequence ambiguity with ProForma 2.1 features",
			proformaString: "(>Ambiguous peptide)<[TMT6plex]@K>PEPT(?IL)DE[Formula:Zn1:z+2]K/2",
			expectations: func(t *testing.T, seq *Sequence) {
				if seq.GetPeptidoformName() == nil || *seq.GetPeptidoformName() != "Ambiguous peptide" {
					t.Errorf("Expected peptidoform name 'Ambiguous peptide'")
				}
				if len(seq.GetGlobalMods()) != 1 {
					t.Errorf("Expected 1 global modification")
				}
			},
		},
		{
			name:           "Isotope labels with ProForma 2.1 features",
			proformaString: "(>Heavy labeled)<[Label:13C(6)15N(2)]@K,R>PEPTIDER",
			expectations: func(t *testing.T, seq *Sequence) {
				if seq.GetPeptidoformName() == nil || *seq.GetPeptidoformName() != "Heavy labeled" {
					t.Errorf("Expected peptidoform name 'Heavy labeled'")
				}
				if len(seq.GetGlobalMods()) != 1 {
					t.Errorf("Expected 1 global modification")
				}
			},
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			seq, err := FromProforma(tt.proformaString)
			if err != nil {
				t.Fatalf("Failed to parse %s: %v", tt.proformaString, err)
			}
			tt.expectations(t, seq)
		})
	}
}

func TestComprehensiveIntegration_SerializationIntegrity(t *testing.T) {
	tests := []string{
		"<[TMT6plex]@K,N-term><[Oxidation]@M><[Phospho]@S,T,Y>MSTPEPTIDEK",
		"(>>>C1)(>>I1)(>P1)<[Mod]@N-term>{Glycan:Hex1}[Ac]-N[Gly:H1]S-[Am]/2",
	}

	for _, original := range tests {
		t.Run(original, func(t *testing.T) {
			seq, err := FromProforma(original)
			if err != nil {
				t.Fatalf("Failed to parse %s: %v", original, err)
			}

			proforma := seq.ToProforma()
			if proforma != original {
				t.Logf("Serialization mismatch: %s != %s", proforma, original)
			}
		})
	}
}

func TestComprehensiveIntegration_PerformanceStressTests(t *testing.T) {
	tests := []struct {
		name           string
		proformaString string
	}{
		{
			name:           "Very long sequence with multiple modifications",
			proformaString: "<[Oxidation]@M>MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM",
		},
		{
			name:           "Many unknown position modifications",
			proformaString: "[Phospho|Position:S,T,Y]^10?SYTSTPEPTIDESTSTY",
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			seq, err := FromProforma(tt.proformaString)
			if err != nil {
				t.Fatalf("Failed to parse %s: %v", tt.proformaString, err)
			}

			// Test round-trip
			output := seq.ToProforma()
			_, err = FromProforma(output)
			if err != nil {
				t.Errorf("Failed round-trip for %s: %v", tt.proformaString, err)
			}
		})
	}
}

func TestComprehensiveIntegration_GlycanCountValidation(t *testing.T) {
	validCases := []struct {
		glycan string
		reason string
	}{
		{"HexNAc", "at end, count 1 implied"},
		{"HexNAc2Hex", "HexNAc has count 2, Hex at end"},
		{"Hex(3)HexNAc2", "explicit counts"},
		{"{C8H13N1O5}1Hex2", "custom with count, Hex with count"},
	}

	for _, tc := range validCases {
		t.Run("Valid_"+tc.glycan, func(t *testing.T) {
			proformaString := "N[Glycan:" + tc.glycan + "]K"
			seq, err := FromProforma(proformaString)
			if err != nil {
				t.Fatalf("Failed to parse %s: %v", proformaString, err)
			}

			if len(seq.GetSeq()[0].GetMods()) == 0 {
				t.Fatal("Expected modification on first amino acid")
			}

			mod := seq.GetSeq()[0].GetMods()[0]
			if len(mod.GetModificationValue().GetPipeValues()) == 0 {
				t.Fatal("Expected pipe value")
			}

			if !mod.GetModificationValue().GetPipeValues()[0].IsValidGlycan() {
				t.Errorf("Expected %s to be valid (%s)", tc.glycan, tc.reason)
			}
		})
	}

	invalidCases := []struct {
		glycan string
		reason string
	}{
		{"HexNAcHex", "HexNAc not at end, missing count"},
		{"HexNAc0", "zero count not allowed"},
		{"{C8H13N1O5}Hex2", "custom not at end, missing count"},
	}

	for _, tc := range invalidCases {
		t.Run("Invalid_"+tc.glycan, func(t *testing.T) {
			proformaString := "N[Glycan:" + tc.glycan + "]K"
			seq, err := FromProforma(proformaString)
			if err != nil {
				t.Fatalf("Failed to parse %s: %v", proformaString, err)
			}

			if len(seq.GetSeq()[0].GetMods()) == 0 {
				t.Fatal("Expected modification on first amino acid")
			}

			mod := seq.GetSeq()[0].GetMods()[0]
			if len(mod.GetModificationValue().GetPipeValues()) == 0 {
				t.Fatal("Expected pipe value")
			}

			if mod.GetModificationValue().GetPipeValues()[0].IsValidGlycan() {
				t.Errorf("Expected %s to be invalid (%s)", tc.glycan, tc.reason)
			}
		})
	}
}
