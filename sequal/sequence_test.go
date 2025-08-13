package sequal

import (
	"testing"
)

func TestBasicSequenceParsing(t *testing.T) {
	t.Run("basic peptide with modification", func(t *testing.T) {
		proforma := "PEP[Phospho]TIDE"
		seq, err := FromProforma(proforma)
		if err != nil {
			t.Fatalf("Failed to parse ProForma '%s': %v", proforma, err)
		}

		if seq.ToStrippedString() != "PEPTIDE" {
			t.Errorf("Expected 'PEPTIDE', got '%s'", seq.ToStrippedString())
		}
		if len(seq.seq) > 2 && len(seq.seq[2].GetMods()) > 0 {
			if seq.seq[2].GetMods()[0].GetValue() != "Phospho" {
				t.Errorf("Expected 'Phospho', got '%s'", seq.seq[2].GetMods()[0].GetValue())
			}
		}
		if seq.ToProforma() != proforma {
			t.Errorf("Expected roundtrip to match, got '%s'", seq.ToProforma())
		}
	})

	t.Run("mass shift notation", func(t *testing.T) {
		proforma := "PEP[+79.966]TIDE"
		seq, err := FromProforma(proforma)
		if err != nil {
			t.Fatalf("Failed to parse ProForma '%s': %v", proforma, err)
		}

		if seq.ToStrippedString() != "PEPTIDE" {
			t.Errorf("Expected 'PEPTIDE', got '%s'", seq.ToStrippedString())
		}
		if len(seq.seq) > 2 && len(seq.seq[2].GetMods()) > 0 {
			mod := seq.seq[2].GetMods()[0]
			if mod.GetValue() != "+79.966" {
				t.Errorf("Expected '+79.966', got '%s'", mod.GetValue())
			}
			if mod.GetMass() != nil {
				if *mod.GetMass() < 79.965 || *mod.GetMass() > 79.967 {
					t.Errorf("Expected mass around 79.966, got %f", *mod.GetMass())
				}
			}
		}
		if seq.ToProforma() != "PEP[+79.966]TIDE" {
			t.Errorf("Expected 'PEP[+79.966]TIDE', got '%s'", seq.ToProforma())
		}
	})

	t.Run("multiple modifications", func(t *testing.T) {
		proforma := "PEPS[Phospho][Acetyl]TIDE"
		seq, err := FromProforma(proforma)
		if err != nil {
			t.Fatalf("Failed to parse ProForma '%s': %v", proforma, err)
		}

		if seq.ToStrippedString() != "PEPSTIDE" {
			t.Errorf("Expected 'PEPSTIDE', got '%s'", seq.ToStrippedString())
		}
		if len(seq.seq) > 3 && len(seq.seq[3].GetMods()) >= 2 {
			mods := seq.seq[3].GetMods()
			if mods[0].GetValue() != "Phospho" {
				t.Errorf("Expected first mod 'Phospho', got '%s'", mods[0].GetValue())
			}
			if mods[1].GetValue() != "Acetyl" {
				t.Errorf("Expected second mod 'Acetyl', got '%s'", mods[1].GetValue())
			}
		}
		if seq.ToProforma() != proforma {
			t.Errorf("Expected roundtrip to match, got '%s'", seq.ToProforma())
		}
	})

	t.Run("ambiguous modifications", func(t *testing.T) {
		proforma := "PEPS{Phospho}TIDE"
		seq, err := FromProforma(proforma)
		if err != nil {
			t.Fatalf("Failed to parse ProForma '%s': %v", proforma, err)
		}

		if seq.ToStrippedString() != "PEPSTIDE" {
			t.Errorf("Expected 'PEPSTIDE', got '%s'", seq.ToStrippedString())
		}
		if len(seq.seq) > 3 && len(seq.seq[3].GetMods()) > 0 {
			mod := seq.seq[3].GetMods()[0]
			if mod.GetValue() != "Phospho" {
				t.Errorf("Expected 'Phospho', got '%s'", mod.GetValue())
			}
			if mod.GetModType() != "ambiguous" {
				t.Errorf("Expected 'ambiguous', got '%s'", mod.GetModType())
			}
		}
		if seq.ToProforma() != proforma {
			t.Errorf("Expected roundtrip to match, got '%s'", seq.ToProforma())
		}
	})

	t.Run("complex sequence", func(t *testing.T) {
		proforma := "[Acetyl]-PEP[Phospho]T{Oxidation}IDE-[Amidated]"
		seq, err := FromProforma(proforma)
		if err != nil {
			t.Fatalf("Failed to parse ProForma '%s': %v", proforma, err)
		}

		if seq.ToStrippedString() != "PEPTIDE" {
			t.Errorf("Expected 'PEPTIDE', got '%s'", seq.ToStrippedString())
		}

		// Check N-terminal modification
		nTermMods := seq.mods[-1]
		if nTermMods == nil || len(nTermMods) == 0 {
			t.Error("Expected N-terminal modification")
		} else {
			if nTermMods[0].GetValue() != "Acetyl" {
				t.Errorf("Expected 'Acetyl', got '%s'", nTermMods[0].GetValue())
			}
		}

		// Check residue modifications
		if len(seq.seq) > 2 && len(seq.seq[2].GetMods()) > 0 {
			if seq.seq[2].GetMods()[0].GetValue() != "Phospho" {
				t.Errorf("Expected 'Phospho', got '%s'", seq.seq[2].GetMods()[0].GetValue())
			}
		}
		if len(seq.seq) > 3 && len(seq.seq[3].GetMods()) > 0 {
			mod := seq.seq[3].GetMods()[0]
			if mod.GetValue() != "Oxidation" {
				t.Errorf("Expected 'Oxidation', got '%s'", mod.GetValue())
			}
			if mod.GetModType() != "ambiguous" {
				t.Errorf("Expected 'ambiguous', got '%s'", mod.GetModType())
			}
		}

		// Check C-terminal modification
		cTermMods := seq.mods[-2]
		if cTermMods == nil || len(cTermMods) == 0 {
			t.Error("Expected C-terminal modification")
		} else {
			if cTermMods[0].GetValue() != "Amidated" {
				t.Errorf("Expected 'Amidated', got '%s'", cTermMods[0].GetValue())
			}
		}

		if seq.ToProforma() != proforma {
			t.Errorf("Expected roundtrip to match, got '%s'", seq.ToProforma())
		}
	})

	t.Run("negative mass shift", func(t *testing.T) {
		proforma := "PEP[-17.027]TIDE"
		seq, err := FromProforma(proforma)
		if err != nil {
			t.Fatalf("Failed to parse ProForma '%s': %v", proforma, err)
		}

		if seq.ToStrippedString() != "PEPTIDE" {
			t.Errorf("Expected 'PEPTIDE', got '%s'", seq.ToStrippedString())
		}
		if len(seq.seq) > 2 && len(seq.seq[2].GetMods()) > 0 {
			mod := seq.seq[2].GetMods()[0]
			if mod.GetValue() != "-17.027" {
				t.Errorf("Expected '-17.027', got '%s'", mod.GetValue())
			}
			if mod.GetMass() != nil {
				mass := *mod.GetMass()
				if mass > -17.026 || mass < -17.028 {
					t.Errorf("Expected mass around -17.027, got %f", mass)
				}
			}
		}
	})

	t.Run("gap notation", func(t *testing.T) {
		proforma := "RTAAX[+367.0537]WT"
		seq, err := FromProforma(proforma)
		if err != nil {
			t.Fatalf("Failed to parse ProForma '%s': %v", proforma, err)
		}

		// Check sequence contains 'X'
		if seq.ToStrippedString() != "RTAAXWT" {
			t.Errorf("Expected 'RTAAXWT', got '%s'", seq.ToStrippedString())
		}

		// Check gap modification
		if len(seq.seq) > 4 {
			if seq.seq[4].GetValue() != "X" {
				t.Errorf("Expected 'X' at position 4, got '%s'", seq.seq[4].GetValue())
			}
			if len(seq.seq[4].GetMods()) > 0 {
				mod := seq.seq[4].GetMods()[0]
				if mod.GetModType() != "gap" {
					t.Errorf("Expected 'gap' mod type, got '%s'", mod.GetModType())
				}
				if mod.GetValue() != "+367.0537" {
					t.Errorf("Expected '+367.0537', got '%s'", mod.GetValue())
				}

				// Verify mass of the gap
				if mod.GetMass() != nil {
					mass := *mod.GetMass()
					if mass < 367.053 || mass > 367.055 {
						t.Errorf("Expected mass around 367.0537, got %f", mass)
					}
				}
			}
		}

		// Verify roundtrip
		if seq.ToProforma() != proforma {
			t.Errorf("Expected roundtrip to match, got '%s'", seq.ToProforma())
		}

		// Test with a negative mass gap
		proforma2 := "PEPTX[-10.0]IDE"
		seq2, err := FromProforma(proforma2)
		if err != nil {
			t.Fatalf("Failed to parse ProForma '%s': %v", proforma2, err)
		}
		if seq2.ToStrippedString() != "PEPTXIDE" {
			t.Errorf("Expected 'PEPTXIDE', got '%s'", seq2.ToStrippedString())
		}
		if len(seq2.seq) > 4 && len(seq2.seq[4].GetMods()) > 0 {
			if seq2.seq[4].GetMods()[0].GetValue() != "-10.0" {
				t.Errorf("Expected '-10.0', got '%s'", seq2.seq[4].GetMods()[0].GetValue())
			}
		}
		expectedRoundtrip := "PEPTX[-10]IDE"
		if seq2.ToProforma() != expectedRoundtrip {
			t.Errorf("Expected '%s', got '%s'", expectedRoundtrip, seq2.ToProforma())
		}
	})
}

func TestModificationHandling(t *testing.T) {
	tests := []struct {
		name            string
		proforma        string
		expectedModPos  int
		expectedModName string
	}{
		{
			name:            "Phosphorylation",
			proforma:        "ELVIS[Phospho]K",
			expectedModPos:  4, // Position of S
			expectedModName: "Phospho",
		},
		{
			name:            "Mass shift",
			proforma:        "ELVIS[+79.966]K",
			expectedModPos:  4,
			expectedModName: "+79.966",
		},
		{
			name:            "Unimod modification",
			proforma:        "ELVIS[U:Phospho]K",
			expectedModPos:  4,
			expectedModName: "Phospho",
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			seq, err := FromProforma(tt.proforma)
			if err != nil {
				t.Fatalf("Failed to parse ProForma '%s': %v", tt.proforma, err)
			}

			if tt.expectedModPos < len(seq.seq) {
				aa := seq.seq[tt.expectedModPos]
				mods := aa.GetMods()
				if len(mods) == 0 {
					t.Errorf("Expected modification at position %d, but found none", tt.expectedModPos)
					return
				}

				mod := mods[0]
				modValue := mod.GetValue()
				if modValue != tt.expectedModName {
					t.Errorf("Expected modification '%s', got '%s'", tt.expectedModName, modValue)
				}
			}
		})
	}
}

func TestTerminalModifications(t *testing.T) {
	t.Run("single terminal modifications", func(t *testing.T) {
		proforma1 := "[Acetyl]-PEPTIDE-[Amidated]"
		seq1, err := FromProforma(proforma1)
		if err != nil {
			t.Fatalf("Failed to parse ProForma '%s': %v", proforma1, err)
		}

		if seq1.ToStrippedString() != "PEPTIDE" {
			t.Errorf("Expected 'PEPTIDE', got '%s'", seq1.ToStrippedString())
		}

		// Check N-terminal modification
		nTermMods1 := seq1.mods[-1]
		if nTermMods1 == nil || len(nTermMods1) == 0 {
			t.Error("Expected N-terminal modification")
		} else {
			if nTermMods1[0].GetValue() != "Acetyl" {
				t.Errorf("Expected 'Acetyl', got '%s'", nTermMods1[0].GetValue())
			}
		}

		// Check C-terminal modification
		cTermMods1 := seq1.mods[-2]
		if cTermMods1 == nil || len(cTermMods1) == 0 {
			t.Error("Expected C-terminal modification")
		} else {
			if cTermMods1[0].GetValue() != "Amidated" {
				t.Errorf("Expected 'Amidated', got '%s'", cTermMods1[0].GetValue())
			}
		}
	})

	t.Run("multiple terminal modifications", func(t *testing.T) {
		proforma2 := "[Acetyl][Methyl]-PEPTIDE-[Amidated][Phosphorylated]"
		seq2, err := FromProforma(proforma2)
		if err != nil {
			t.Fatalf("Failed to parse ProForma '%s': %v", proforma2, err)
		}

		if seq2.ToStrippedString() != "PEPTIDE" {
			t.Errorf("Expected 'PEPTIDE', got '%s'", seq2.ToStrippedString())
		}

		// Check multiple N-terminal modifications
		nTermMods2 := seq2.mods[-1]
		if nTermMods2 == nil || len(nTermMods2) != 2 {
			t.Errorf("Expected 2 N-terminal modifications, got %d", len(nTermMods2))
		} else {
			if nTermMods2[0].GetValue() != "Acetyl" {
				t.Errorf("Expected first N-term mod 'Acetyl', got '%s'", nTermMods2[0].GetValue())
			}
			if nTermMods2[1].GetValue() != "Methyl" {
				t.Errorf("Expected second N-term mod 'Methyl', got '%s'", nTermMods2[1].GetValue())
			}
		}

		// Check multiple C-terminal modifications
		cTermMods2 := seq2.mods[-2]
		if cTermMods2 == nil || len(cTermMods2) != 2 {
			t.Errorf("Expected 2 C-terminal modifications, got %d", len(cTermMods2))
		} else {
			if cTermMods2[0].GetValue() != "Amidated" {
				t.Errorf("Expected first C-term mod 'Amidated', got '%s'", cTermMods2[0].GetValue())
			}
			if cTermMods2[1].GetValue() != "Phosphorylated" {
				t.Errorf("Expected second C-term mod 'Phosphorylated', got '%s'", cTermMods2[1].GetValue())
			}
		}
	})

	t.Run("terminal modifications with hyphens in names", func(t *testing.T) {
		proforma3 := "[N-Terminal-Acetyl]-PEPTIDE-[C-Terminal-Amidation]"
		seq3, err := FromProforma(proforma3)
		if err != nil {
			t.Fatalf("Failed to parse ProForma '%s': %v", proforma3, err)
		}

		if seq3.ToStrippedString() != "PEPTIDE" {
			t.Errorf("Expected 'PEPTIDE', got '%s'", seq3.ToStrippedString())
		}

		// Check N-terminal modification with hyphen in name
		nTermMods3 := seq3.mods[-1]
		if nTermMods3 == nil || len(nTermMods3) == 0 {
			t.Error("Expected N-terminal modification")
		} else {
			if nTermMods3[0].GetValue() != "N-Terminal-Acetyl" {
				t.Errorf("Expected 'N-Terminal-Acetyl', got '%s'", nTermMods3[0].GetValue())
			}
		}

		// Check C-terminal modification with hyphen in name
		cTermMods3 := seq3.mods[-2]
		if cTermMods3 == nil || len(cTermMods3) == 0 {
			t.Error("Expected C-terminal modification")
		} else {
			if cTermMods3[0].GetValue() != "C-Terminal-Amidation" {
				t.Errorf("Expected 'C-Terminal-Amidation', got '%s'", cTermMods3[0].GetValue())
			}
		}
	})

	// Verify roundtrip for all cases
	testCases := []string{
		"[Acetyl]-PEPTIDE-[Amidated]",
		"[Acetyl][Methyl]-PEPTIDE-[Amidated][Phosphorylated]",
		"[N-Terminal-Acetyl]-PEPTIDE-[C-Terminal-Amidation]",
		"[N-Acetyl][alpha-amino]-PEPTIDE-[C-Terminal][beta-COOH]",
	}

	for _, proforma := range testCases {
		t.Run("roundtrip_"+proforma, func(t *testing.T) {
			seq, err := FromProforma(proforma)
			if err != nil {
				t.Fatalf("Failed to parse ProForma '%s': %v", proforma, err)
			}
			
			regenerated := seq.ToProforma()
			if regenerated != proforma {
				t.Errorf("Roundtrip failed. Expected '%s', got '%s'", proforma, regenerated)
			}
		})
	}
}

func TestChargeHandling(t *testing.T) {
	tests := []struct {
		name         string
		proforma     string
		expectedCharge *int
		expectedSpecies *string
	}{
		{
			name:         "Simple charge",
			proforma:     "PEPTIDE/2",
			expectedCharge: IntPtr(2),
			expectedSpecies: nil,
		},
		{
			name:         "Negative charge",
			proforma:     "PEPTIDE/-3",
			expectedCharge: IntPtr(-3),
			expectedSpecies: nil,
		},
		{
			name:         "Charge with ionic species",
			proforma:     "PEPTIDE/2[+Na+]",
			expectedCharge: IntPtr(2),
			expectedSpecies: StringPtr("+Na+"),
		},
		{
			name:         "Complex with modification and charge",
			proforma:     "ELVIS[Phospho]K/3",
			expectedCharge: IntPtr(3),
			expectedSpecies: nil,
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			seq, err := FromProforma(tt.proforma)
			if err != nil {
				t.Fatalf("Failed to parse ProForma '%s': %v", tt.proforma, err)
			}

			charge := seq.GetCharge()
			if tt.expectedCharge == nil {
				if charge != nil {
					t.Errorf("Expected no charge, got %d", *charge)
				}
			} else {
				if charge == nil {
					t.Errorf("Expected charge %d, got nil", *tt.expectedCharge)
				} else if *charge != *tt.expectedCharge {
					t.Errorf("Expected charge %d, got %d", *tt.expectedCharge, *charge)
				}
			}

			species := seq.GetIonicSpecies()
			if tt.expectedSpecies == nil {
				if species != nil {
					t.Errorf("Expected no ionic species, got %s", *species)
				}
			} else {
				if species == nil {
					t.Errorf("Expected ionic species %s, got nil", *tt.expectedSpecies)
				} else if *species != *tt.expectedSpecies {
					t.Errorf("Expected ionic species %s, got %s", *tt.expectedSpecies, *species)
				}
			}
		})
	}
}

func TestGlobalModifications(t *testing.T) {
	tests := []struct {
		name              string
		proforma          string
		expectedGlobalMods int
	}{
		{
			name:              "Fixed modification",
			proforma:          "<[Carbamidomethyl]@C>PEPTCDE",
			expectedGlobalMods: 1,
		},
		{
			name:              "Isotope labeling",
			proforma:          "<15N>PEPTIDE",
			expectedGlobalMods: 1,
		},
		{
			name:              "Multiple global mods",
			proforma:          "<15N><[Carbamidomethyl]@C>PEPTCDE",
			expectedGlobalMods: 2,
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			seq, err := FromProforma(tt.proforma)
			if err != nil {
				t.Fatalf("Failed to parse ProForma '%s': %v", tt.proforma, err)
			}

			globalMods := seq.GetGlobalMods()
			if len(globalMods) != tt.expectedGlobalMods {
				t.Errorf("Expected %d global modifications, got %d", tt.expectedGlobalMods, len(globalMods))
			}
		})
	}
}

func TestChimericSequences(t *testing.T) {
	tests := []struct {
		name               string
		proforma           string
		expectedPeptidoforms int
		expectedFirstSeq   string
		expectedSecondSeq  string
	}{
		{
			name:               "Basic chimeric",
			proforma:           "PEPTIDE/2+ANOTHER/3",
			expectedPeptidoforms: 2,
			expectedFirstSeq:   "PEPTIDE",
			expectedSecondSeq:  "ANOTHER",
		},
		{
			name:               "Complex chimeric with modifications",
			proforma:           "[Acetyl]-PEP[+79.966]TIDE-[Amidated]/2[+Na+]+S[Phospho]EQ/3",
			expectedPeptidoforms: 2,
			expectedFirstSeq:   "PEPTIDE",
			expectedSecondSeq:  "SEQ",
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			seq, err := FromProforma(tt.proforma)
			if err != nil {
				t.Fatalf("Failed to parse ProForma '%s': %v", tt.proforma, err)
			}

			if !seq.IsChimeric() {
				t.Errorf("Expected sequence to be chimeric")
			}

			peptidoforms := seq.GetPeptidoforms()
			if len(peptidoforms) != tt.expectedPeptidoforms {
				t.Errorf("Expected %d peptidoforms, got %d", tt.expectedPeptidoforms, len(peptidoforms))
			}

			if len(peptidoforms) >= 1 {
				firstSeq := peptidoforms[0].ToStrippedString()
				if firstSeq != tt.expectedFirstSeq {
					t.Errorf("Expected first sequence '%s', got '%s'", tt.expectedFirstSeq, firstSeq)
				}
			}

			if len(peptidoforms) >= 2 {
				secondSeq := peptidoforms[1].ToStrippedString()
				if secondSeq != tt.expectedSecondSeq {
					t.Errorf("Expected second sequence '%s', got '%s'", tt.expectedSecondSeq, secondSeq)
				}
			}
		})
	}
}

func TestMultiChainSequences(t *testing.T) {
	proforma := "PEPTIDE//SEQUENCE//THIRD"
	seq, err := FromProforma(proforma)
	if err != nil {
		t.Fatalf("Failed to parse multi-chain ProForma '%s': %v", proforma, err)
	}

	if !seq.IsMultiChain() {
		t.Errorf("Expected sequence to be multi-chain")
	}

	chains := seq.GetChains()
	if len(chains) != 3 {
		t.Errorf("Expected 3 chains, got %d", len(chains))
	}

	expectedChains := []string{"PEPTIDE", "SEQUENCE", "THIRD"}
	for i, chain := range chains {
		stripped := chain.ToStrippedString()
		if stripped != expectedChains[i] {
			t.Errorf("Expected chain %d to be '%s', got '%s'", i, expectedChains[i], stripped)
		}
	}
}

func TestRoundTripConversion(t *testing.T) {
	tests := []string{
		"PEPTIDE",
		"ELVIS[Phospho]K",
		"[Acetyl]-PEPTIDE-[Amidated]",
		"PEPTIDE/2",
		"PEPTIDE/2[+Na+]",
		"<[Carbamidomethyl]@C>PEPTCDE",
		// Note: Some complex cases might have minor formatting differences
	}

	for _, proforma := range tests {
		t.Run(proforma, func(t *testing.T) {
			seq, err := FromProforma(proforma)
			if err != nil {
				t.Fatalf("Failed to parse ProForma '%s': %v", proforma, err)
			}

			regenerated := seq.ToProforma()
			if regenerated == "" {
				t.Errorf("Failed to regenerate ProForma for '%s'", proforma)
			}

			// Parse the regenerated ProForma to ensure it's valid
			seq2, err := FromProforma(regenerated)
			if err != nil {
				t.Errorf("Failed to parse regenerated ProForma '%s': %v", regenerated, err)
			}

			// Compare the stripped sequences
			if seq.ToStrippedString() != seq2.ToStrippedString() {
				t.Errorf("Stripped sequences don't match: '%s' vs '%s'", 
					seq.ToStrippedString(), seq2.ToStrippedString())
			}
		})
	}
}

func TestInfoTags(t *testing.T) {
	proforma := "ELVIS[Phospho|INFO:newly discovered]K"
	seq, err := FromProforma(proforma)
	if err != nil {
		t.Fatalf("Failed to parse ProForma '%s': %v", proforma, err)
	}

	// Check that the modification has info tags
	aa := seq.seq[4] // Position of S
	mods := aa.GetMods()
	if len(mods) == 0 {
		t.Errorf("Expected modification with info tag, but found no modifications")
		return
	}

	mod := mods[0]
	infoTags := mod.GetInfoTags()
	if len(infoTags) == 0 {
		t.Errorf("Expected info tags, but found none")
		return
	}

	expectedTag := "newly discovered"
	if infoTags[0] != expectedTag {
		t.Errorf("Expected info tag '%s', got '%s'", expectedTag, infoTags[0])
	}
}

func TestUtilityFunctions(t *testing.T) {
	seq, err := FromProforma("PEPTIDE")
	if err != nil {
		t.Fatalf("Failed to parse ProForma: %v", err)
	}

	// Test Count function
	count := seq.Count("P", 0, 7)
	if count != 2 {
		t.Errorf("Expected 2 P residues, got %d", count)
	}

	// Test Gaps function
	gapSeq, err := FromProforma("PEP-TIDE")
	if err != nil {
		t.Fatalf("Failed to parse ProForma with gap: %v", err)
	}

	gaps := gapSeq.Gaps()
	expectedGaps := []bool{false, false, false, true, false, false, false, false}
	for i, gap := range gaps {
		if i < len(expectedGaps) && gap != expectedGaps[i] {
			t.Errorf("Expected gap at position %d to be %v, got %v", i, expectedGaps[i], gap)
		}
	}

	// Test GetLength
	if seq.GetLength() != 7 {
		t.Errorf("Expected length 7, got %d", seq.GetLength())
	}

	// Test ToMap
	seqMap := seq.ToMap()
	if seqMap["sequence"] != "PEPTIDE" {
		t.Errorf("Expected sequence 'PEPTIDE' in map, got %v", seqMap["sequence"])
	}
}

func TestComplexProFormaExamples(t *testing.T) {
	// These are examples from the TypeScript package README
	tests := []struct {
		name     string
		proforma string
		shouldParse bool
	}{
		{
			name:     "Joint representation",
			proforma: "ELVIS[U:Phospho|+79.966331]K",
			shouldParse: true,
		},
		{
			name:     "Observed mass",
			proforma: "ELVIS[U:Phospho|Obs:+79.978]K",
			shouldParse: true,
		},
		{
			name:     "Crosslinks",
			proforma: "PEPTK[XL:DSS#XL1|+138.068|INFO:reaction=NHS]IDE",
			shouldParse: true,
		},
		{
			name:     "Gap notation",
			proforma: "RTAAX[+367.0537]WT",
			shouldParse: true,
		},
		{
			name:     "Multiple info tags",
			proforma: "ELVIS[Phospho|INFO:newly discovered|INFO:Created on 2021-06]K",
			shouldParse: true,
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			seq, err := FromProforma(tt.proforma)
			if tt.shouldParse {
				if err != nil {
					t.Errorf("Expected to parse '%s', but got error: %v", tt.proforma, err)
				} else if seq == nil {
					t.Errorf("Expected to parse '%s', but got nil sequence", tt.proforma)
				}
			} else {
				if err == nil && seq != nil {
					t.Errorf("Expected parsing to fail for '%s', but it succeeded", tt.proforma)
				}
			}
		})
	}
}